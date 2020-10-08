## library
suppressPackageStartupMessages(c(
    library(BiocManager),
    options(repos = BiocManager::repositories()),
    library(shiny),
    library(shinyWidgets),
    library(tidyverse),
    library(SummarizedExperiment),
    library(plot3D),
    library(data.table),
    library(rgl),
    library(htmltools),
    library(dplyr),
    library(rgl),
    library(htmltools),
    library(grid),
    library(gridBase),
    library(ggkm),
    library(sqldf),
    library(rstan),
    library(autoimage)
))

## Load data
umap_cancertype <- fread("umap_3d_coors.tsv")
load("kelly.colours.rda")
load("ditto.colours.rda")
load("tumor_samples.Rda")
load("supplemental_gene.Rda")
load("gdc_reactome_path_analysis.Rda")
load("cancerdat.rda")
load("survmodel-ethbkd.stansave")
load("dfclas.rda")

cluster.colours <- c(kelly.colours[c(3:12)], "grey") ## Color scheme
names(cluster.colours) <- c(1:11) ## Assign cluster to color

sample.colours <- umap_cancertype %>% dplyr::select(sample_id, clust_knn) %>% column_to_rownames(var = "sample_id")
 
## Generating selection inputs
cancer_type <- tumor_samples %>% .$project_code %>% unique %>% ## get cancer type info
    .[order(.)] ## Order alphabetically

genename <- colnames(supplemental_gene) ## get gene names
genename <- genename[order(genename)] ## order alphabetically

pathwayname <- rownames(gdc_reactome_path_analysis) ## get pathway names
pathwayname <- pathwayname[order(pathwayname)] ## order alphabetically

## User interface setup
ui <- (fluidPage(
    ## App title
    titlePanel("Interactive Supplemental File"),
    h3("Nguyen et al. 2020"),
    ## Interactive settings
    playwidgetOutput("player"),
    
    mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel(title = "Pan-cancer", h5("Pan and zoom enabled"),
                             selectInput("projection.selected", label = "Projection Type:",
                                         choices = list("Cluster" = "clust_knn",
                                                        "Cancer type" = "project_code",
                                                        "Stage" = "ajcc_pathologic_tumor_stage",
                                                        "Metastasis" = "ajcc_metastasis_pathologic_pm",
                                                        "Log mutation count" = "mut_count"),
                                         selected = "Cluster"),
                            # rglwidgetOutput("pancan",  width = 500, height = 500)), #rglwidgetOutput("pancan",  width = 500, height = 500)
                    fluidRow(
                        column(width = 8, rglwidgetOutput("pancan",  width = 500, height = 500)),
                        column(width = 10, plotOutput("pancanlegend")))),
                    tabPanel(title = "Cancer Type",
                             pickerInput("cancer_type", label = "Cancer type:", ## Cancer type
                                         choices = cancer_type, selected = "BRCA", multiple = TRUE,
                                         options = list(`actions-box` = TRUE)),
                             plotOutput("cancertype")),
                    tabPanel(title = "Gene", 
                             selectInput("genename", label = "Gene:", ## Gene 
                                         choices = genename, selected = "TP53"),
                             pickerInput("cancer_type_gene", label = "Cancer type:", ## Cancer type
                                         choices = cancer_type, selected = cancer_type, multiple = TRUE,
                                         options = list(`actions-box` = TRUE)),
                             plotOutput("gene")),
                    tabPanel(title = "Pathway", 
                             selectInput("pathwayname", label = "Pathway:", ## Pathway
                                         choices = pathwayname, selected = "Cellular response to hypoxia"),
                             pickerInput("cancer_type_pathway", label = "Cancer type:", ## Cancer type
                                         choices = cancer_type, selected = cancer_type, multiple = TRUE,
                                         options = list(`actions-box` = TRUE)),
                             plotOutput("pathway")),
                    tabPanel(title = "Survival",
                             pickerInput("sim.cancertype", label = "Cancer type:", ## Cancer type
                                         choices = tisnames, selected = "BRCA", multiple = TRUE,
                                         options = list(`actions-box` = TRUE)),
                             pickerInput("sim.cluster", label = "Cluster:", ## Cluster
                                         choices = c(as.character(1:10), "All"), selected = "All", multiple = TRUE,
                                         options = list(`actions-box` = TRUE)),
                             pickerInput("sim.sex", label = "Sex:", ## Sex
                                         choices = list("MALE", "FEMALE"), selected = "FEMALE"),
                             pickerInput("sim.ethnic", label = "Ethnicity:", ## Ethnicity
                                         choices = list("European descent" = "WHITE", 
                                                        "African descent"= "BLACK OR AFRICAN AMERICAN"), selected = "European descent"),
                             sliderInput("sim.range", "Age range:", ## Age range
                                         min = 1, max = 100,
                                         value = c(30,70)),
                             sliderInput("sim.population", label = "Population Size:", ## Population size
                                         min = 5, max = 500,
                                         value = 100),
                             checkboxInput("show.sim.k.eq.1", label = "Show Control",
                                           value = FALSE),
                             plotOutput("survival"))
                    )
        )
    )
)


server <- (function(input, output, session) {
    ## Pancancer plot code
    
    output$pancan <- renderRglwidget ({ 
        
        projection.df <- umap_cancertype %>% dplyr::select(plot_x, plot_y, plot_z, input$projection.selected) %>%
            'colnames<-' (c("plot_x", "plot_y", "plot_z", "projection_type"))
        
        if(input$projection.selected == "ajcc_pathologic_tumor_stage") {
            projection.df <- projection.df %>% filter(projection_type != "Not Available")
            color_anno <- c('darkgoldenrod3','dodgerblue3','yellow3','tomato2')
            names(color_anno) <- c("Stage I", "Stage II", "Stage III", "Stage IV") ## Assign stage to color
        }
        if(input$projection.selected == "ajcc_metastasis_pathologic_pm") {
            projection.df <- projection.df %>% filter(!is.na(projection_type))
            color_anno <- c('#a6cee3','#1f78b4')
            names(color_anno) <- c("M0", "M1")
        }
        if(input$projection.selected == "project_code") {
            cancernames <- unique(projection.df$projection_type)
            cancernames <- cancernames[order(cancernames)]
            color_anno <- ditto_colours[2:(length(cancernames)+1)]
            names(color_anno) <- cancernames
        }
        if(input$projection.selected == "clust_knn") {
            clusternames <- unique(projection.df$projection_type)
            clusternames <- clusternames[order(clusternames)]
            color_anno <- kelly.colours[3:(length(clusternames)+2)]
            names(color_anno) <- clusternames
        }
        if(input$projection.selected == "mut_count") {
            projection.df$projection_type <- log10(projection.df$projection_type + 1) ## adding pseudocount for 0's
            max.value <- max(projection.df$projection_type)
            min.value <- min(projection.df$projection_type)
            mycols <- colorRampPalette(c("#fee0d2", "#de2d26"))
        }
        
        
        open3d(useNULL=TRUE)
        if(input$projection.selected == "mut_count") {
            plot3d(projection.df[, c("plot_x", "plot_y", "plot_z"),],
                   type = "s", size = 2.5,
                   xlab = "UMAP 1", ylab = "UMAP 2", zlab = "UMAP 3",
                   aspect = c(8.3, 10.5, 7.6),
                   col = mycols(10))["data"]
        } else{
        plot3d(projection.df[, c("plot_x", "plot_y", "plot_z"),],
               type = "s", size = 2.5,
               xlab = "UMAP 1", ylab = "UMAP 2", zlab = "UMAP 3",
               aspect = c(8.3, 10.5, 7.6),
               col = color_anno[as.character(projection.df$projection_type)])["data"]
        }
        rglwidget()
        
    })
    
    output$pancanlegend <- renderPlot({
        if(input$projection.selected == "clust_knn") {
        cluster.colours <- c("grey", kelly.colours[c(3:12)]) ## Color scheme for cluster
        names(cluster.colours) <- c(0:10) ## assigning cluster to the colors
        MyOrder = matrix(1:10, nrow = 2, ncol = 5, byrow = T)
        legend("topleft", inset = c(0.15), names(cluster.colours[2:11])[MyOrder], fill=cluster.colours[2:11][MyOrder], ncol=5)
        }
        if(input$projection.selected == "ajcc_metastasis_pathologic_pm") {
            color_anno <- c('#a6cee3','#1f78b4')
            names(color_anno) <- c("M0", "M1")
            MyOrder = matrix(1:2, nrow = 1, ncol = 2, byrow = T)
            legend("topleft", inset = c(0.15), legend=names(color_anno[1:2])[MyOrder], fill=color_anno[1:2][MyOrder], ncol=2)
        }
        if(input$projection.selected == "project_code") {
            cancernames <- unique(umap_cancertype$project_code)
            cancernames <- cancernames[order(cancernames)]
            color_anno <- ditto_colours[2:(length(cancernames)+1)]
            names(color_anno) <- cancernames
            MyOrder = matrix(1:24, nrow = 6, ncol = 4, byrow = T)
            legend("topleft", inset = c(0.15), legend = names(color_anno[1:23])[MyOrder[1:23]], fill=color_anno[1:23][MyOrder[1:23]], ncol=4)
        }
        if(input$projection.selected == "ajcc_pathologic_tumor_stage") {
            color_anno <- c('darkgoldenrod3','dodgerblue3','yellow3','tomato2')
            names(color_anno) <- c("Stage I", "Stage II", "Stage III", "Stage IV") ## Assign stage to color
            MyOrder = matrix(1:4, nrow = 1, ncol = 4, byrow = T)
            legend("topleft", inset = c(0.15), legend=names(color_anno[1:4])[MyOrder], fill=color_anno[1:4][MyOrder], ncol=4)
        }
        if(input$projection.selected == "mut_count") {
            mycols <- colorRampPalette(c("#fee0d2", "#de2d26"))
            gl <- grid.layout(nrow=10, ncol=5) 
            vp.1 <- viewport(layout.pos.col=2, layout.pos.row=2)
            pushViewport(viewport(layout=gl))
            pushViewport(vp.1)
            par(new=TRUE, fig=gridFIG())
            par(mar = c(0,0,2,0), mai=c(0.1,0.1,0.1,0.1)) ## setting plot parameters
            legend.scale(c(0, 4.3), col = mycols(10))
            popViewport()
        }
        })
    ## Cancer type plot code
        output$cancertype <- renderPlot({
            # setup layout
            gl <- grid.layout(nrow=2, ncol=2) 
            ## Setup viewports - How figures will be laid out
            vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1)
            vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1) 
            vp.3 <- viewport(layout.pos.col=1, layout.pos.row=2) 
            vp.4 <- viewport(layout.pos.col=2, layout.pos.row=2) 
            # init layout
            pushViewport(viewport(layout=gl))
            
            pushViewport(vp.1)
            par(new=TRUE, fig=gridFIG())
            par(mar = c(0,2,0,0), mai=c(0.1,0.1,0.1,0.1)) ## setting plot parameters
            with(umap_cancertype %>% ## using umap_cancertype data
                     filter(project_code %in% input$cancer_type), ## Only include selected cancer type
                 scatter3D(plot_y, -plot_z, plot_x, ## Plot cancer type using umap 3d coors
                           bg = cluster.colours[as.character(clust_knn)], pch = 21, cex = 0.8, lwd = 0.2,
                           theta = 0, phi = 65, scale = F,
                           xlim = c(-6, 4.5), zlim = c(-4.5, 4), ylim = c(-4, 5),
                           xlab = "UMAP 2", ylab = "UMAP 3", zlab = "UMAP 1",
                           colvar = NULL))
            popViewport()
            
            ## Cluster frequency calculation
            class_freq <- umap_cancertype %>% dplyr::select(clust_knn, project_code) %>% ## select relevant columns
                .[.$project_code %in% input$cancer_type,] %>% ## Only include selected cancer type
                mutate(clust_knn = as.character(.$clust_knn)) ## Convert cluster membership to character
            
            freq_by_class <- list() ## created list to store data in the for loop below
            
            for(n in seq_along(1:max(umap_cancertype$clust_knn))) { ## Seq along clusters 1-10
                df <- data.frame(clust_knn = n,
                                 freq = nrow(class_freq[class_freq$clust_knn %in% n, ]), ## count number of tumor in each cluster
                                 stringsAsFactors = FALSE)
                freq_by_class[[n]] <- df ## store data in list
            }
            freq_by_class <- freq_by_class %>% do.call(rbind, .) %>% ## rbind list
                mutate(clust_knn = factor(x = as.character(.$clust_knn), levels = as.character(c(1:max(umap_cancertype$clust_knn))))) ## turn to factor and relevel 

            
            pushViewport(vp.2)
            par(new=TRUE, fig=gridFIG())
            par(mar = c(2,4,2,0))
            barplot(height = freq_by_class$freq / nrow(class_freq), ## freq / tot_n_cancertype 
                    names.arg = c(1:10), col = cluster.colours,
                    ylab = "Relative Proportion",
                    las = 1) ## always horizontal axis
            popViewport()
            
            if(length(input$cancer_type) == 1) { ## if only one cancer type is selected then do this
                
            spec_cancer <- umap_cancertype %>% filter(project_code == input$cancer_type) %>% ## Filter for selected cancer type
                dplyr::select(plot_x, plot_y, plot_z, Subtype_Selected) ## Select relevant column 
            spec_cancer[is.na(spec_cancer)] <- "NA" ## Convert NA to character "NA"
            
            s_name <- unique(spec_cancer$Subtype_Selected) ## What are the subtype names
            s_name <- s_name[order(s_name,na.last=T)] ## Order the names ABC and have NA last
            
            
            s_color <- ditto_colours[c(2:(length(s_name) + 1))]
            names(s_color) <- s_name
            
            pushViewport(vp.3)
            par(new=TRUE, fig=gridFIG())
            par(mar = c(0,2,0,0), mai=c(0.1,0.1,0.1,0.1)) ## setting plot parameters
            ## color scheme
            with(spec_cancer,
                 scatter3D(plot_y, -plot_z, plot_x, 
                           bg = s_color[as.character(Subtype_Selected)], pch = 21, cex = 0.8, lwd = 0.2,
                           theta = 0, phi = 65, scale = F,
                           xlim = c(-6, 4.5), zlim = c(-4.5, 4), ylim = c(-4, 5),
                           xlab = "UMAP 2", ylab = "UMAP 3", zlab = "UMAP 1",
                           colvar = NULL))
            legend("top", inset = c(-0.01), legend=names(s_color[1:2]), col=s_color[1:2], pch = 20, bty = "n", horiz = T, cex = 1.2)
            legend("top", inset = c(0.05), legend=names(s_color[3:4]), col=s_color[3:4], pch = 20, bty = "n", horiz = T, cex = 1.2)
            legend("top", inset = c(0.11), legend=names(s_color[5:6]), col=s_color[5:6], pch = 20, bty = "n", horiz = T, cex = 1.2)
            legend("top", inset = c(0.17), legend=names(s_color[7:8]), col=s_color[7:8], pch = 20, bty = "n", horiz = T, cex = 1.2)
            popViewport()
            
            # Subtype frequency
            class_freq <- umap_cancertype %>% dplyr::select(clust_knn, project_code, Subtype_Selected) %>%
                .[.$project_code %in% input$cancer_type,] %>%
                dplyr::select(-project_code) %>%
                mutate(clust_knn = as.character(.$clust_knn))
            class_freq[is.na(class_freq)] <- "NA"
            
            class_freq$clust_knn <- factor(class_freq$clust_knn, levels = 1:10)
            counts <- table(class_freq$Subtype_Selected, class_freq$clust_knn)
            
            for(n in 1:ncol(counts)) {
                message(n)
                tot_n_clust <- sum(counts[,n])
                counts[,n] <- round(counts[,n] / tot_n_clust, digits = 3)
            }
            sum(counts[,1])
            
            pushViewport(vp.4)
            par(new=TRUE, fig=gridFIG())
            par(mar = c(2,4,2,0))
            barplot(height = counts,
                    names.arg = c(1:10), col = s_color,
                    ylab = "Relative Proportion",
                    las = 1) ## always horizontal axis
            popViewport()
            }
        })
        
    
    ## Gene plot code
    output$gene <- renderPlot({ 
        par(mfrow = c(1,2), mai = c(0.1, 0.1, 0.1, 0.1)) ## creating plot layout 1 x 2
        par(mar = c(2,2,0,0))
        selected.tumor.gene <- umap_cancertype %>% dplyr::select(sample_id, project_code) %>% ## Selecting relevant columns
            filter(project_code %in% c(input$cancer_type_gene)) ## Filter for selected cancers
            
        gene.df <- supplemental_gene %>% dplyr::select(input$genename) %>% ## filtering for gene
            "colnames<-" ("selectedgene") %>% ## rename column
            rownames_to_column(var = "sample_id") %>% ## move sample id to column
            left_join(umap_cancertype[,c("plot_x", "plot_y", "plot_z", "clust_knn", "sample_id")], ., by = "sample_id") %>% ## left join data by sample_id
            filter(selectedgene == 1) %>% ## only keep sample ids with a mutation in the gene
            filter(sample_id %in% selected.tumor.gene$sample_id) ## filter data with selected cancer type
        with(gene.df, ## only select tumor samples that have a mutation in the selected gene
             scatter3D(plot_y, -plot_z, plot_x,
                       bg = cluster.colours[as.character(gene.df$clust_knn)], pch = 21, cex = 0.8, lwd = 0.2,
                       theta = 0, phi = 65, scale = F,
                       xlim = c(-6, 4.5), zlim = c(-4.5, 4), ylim = c(-4, 5),
                       xlab = "UMAP 2", ylab = "UMAP 3", zlab = "UMAP 1",
                       colvar = NULL))
    })
    ## Pathway plot code
    output$pathway <- renderPlot({ 
        par(mfrow = c(1,2), mai = c(0.1, 0.1, 0.1, 0.1)) ## creating plot layout 1 x 2
        par(mar = c(2,2,0,0))
        selected.tumor.path <- umap_cancertype %>% dplyr::select(sample_id, project_code) %>% ## Selecting relevant columns
            filter(project_code %in% c(input$cancer_type_pathway)) ## Filter for selected cancers
       pathway.df <- gdc_reactome_path_analysis %>% t() %>% as.data.frame() %>% ## transpose data for left_joining later
            dplyr::select(input$pathwayname) %>% ## filtering for pathway
            "colnames<-" ("selectedpath") %>% ## rename column
            rownames_to_column(var = "sample_id") %>% ## move sample id to column
            left_join(umap_cancertype[,c("plot_x", "plot_y", "plot_z", "clust_knn", "sample_id")], ., by = "sample_id") %>% ## left join data by sample_id
            filter(selectedpath == 1) %>%
           filter(sample_id %in% selected.tumor.path$sample_id) ## filter data with selected cancer type
        with(pathway.df, ## only select tumor samples that have a mutation in the selected gene
             scatter3D(plot_y, -plot_z, plot_x, 
                       bg = cluster.colours[as.character(pathway.df$clust_knn)], pch = 21, cex = 0.8, lwd = 0.2,
                       theta = 0, phi = 65, scale = F,
                       xlim = c(-6, 4.5), zlim = c(-4.5, 4), ylim = c(-4, 5),
                       xlab = "UMAP 2", ylab = "UMAP 3", zlab = "UMAP 1",
                       colvar = NULL))
    })
    
    ## Survival plot code
    output$survival <- renderPlot({
        
        ## Survival plot code
        set.seed(20200805) #a date based seed for reproducibility
        ## these are scaled in the stan file, so we need to scale here as well
        r20scale <- 1e-5
        ageratescale <- .05

        r20 <- c(0.365209, 0.367732, 0.093299, 0.098039) * r20scale
        names(r20) <- c("WF", "WM", "BF", "BM")

        agerate <- c(1.197896, 1.291932, 1.649803, 1.734934) * ageratescale
        names(agerate) <- c("WF", "WM", "BF", "BM")

        ## Model ##
        lccdf = function(r20,agerate,age,k,t){
            return(-r20*((365*exp((agerate*k*t)/365+age*agerate-20*agerate))/(agerate*k)
                         -(365*exp(age*agerate-20*agerate))/(agerate*k)));
        }


        simdeath = function(r20,agerate,age,k){
            r = runif(1)
            lr = log(r)
            tdays = uniroot(function(t) return(lccdf(r20,agerate,age,k,t)-lr),c(0,365*400))
            return(tdays$root) ## years since diagnosis at "age"
        }

        # cancerdat <- read_csv("clinical_plus_cluster.csv")
        # #cancerdat <- read_csv("~/TTmanu/code/supplemental_file_shiny/clinical_plus_cluster.csv")
        # cancerdat$type <- as.factor(cancerdat$type) ## the tissue
        # tisnames <- levels(cancerdat$type)
        # save(cancerdat, tisnames, file = "~/TTmanu/code/supplemental_file_shiny/cancerdat.rda")
        #load("cancerdat.rda")
        ## Getting vital status and final time after diagnosis
        endpoints <- sqldf("select case when death_days_to is NULL then last_contact_days_to else death_days_to end as finaltime, case when death_days_to is NULL then 0 else 1 end as finalstatus from cancerdat;")

        cancerdat <- cbind(cancerdat,endpoints)
        ## Are there any NA's in finaltime column?
        sum(is.na(cancerdat$finaltime))
        ## Remove data with NA's in finaltime column
        cancerdat <- cancerdat[!is.na(cancerdat$finaltime),]

        #load("../code/survival/survmodel-ethbkd.stansave") ## loads a variable called "sampseth" which is a stanfit
        #load("survmodel-ethbkd.stansave")
        asampseth <- as.array(sampseth)
        dfsamps = as.data.frame(sampseth)
        
        ## Data for kmeans tissue
        tissue.kvalues <- dfsamps[,235:257] ## 
        names(tissue.kvalues) <- tisnames
        #names(tissue.kvalues) <- gsub(pattern = ".*\\[", replacement = "", x = names(tissue.kvalues)) ## removing brackets from colnames because shiny doesn't like
        #names(tissue.kvalues) <- gsub(pattern = "\\].*", replacement = "", x = names(tissue.kvalues))
        

        ## Reactive options
        sim.cancertype.selection <- as.character(input$sim.cancertype)
        #sim.cancertype.selection <- c("OV")
        sim.cluster.selection <- input$sim.cluster
        sim.sex.selection <- input$sim.sex
        #sim.sex.selection <- "FEMALE"
        sim.ethnic.selection <- input$sim.ethnic
        #sim.ethnic.selection <- "WHITE"
        
        if( any(sim.cancertype.selection %in% c("PRAD")) & sim.sex.selection == "FEMALE") { stop('invalid selection [incompatible sex and cancer type]')}
        if( any(sim.cancertype.selection %in% c("OV", "UCEC")) & sim.sex.selection == "MALE") {stop('invalid selection [incompatible sex and cancer type]')}

        agemodel <- paste0(substring(text = sim.ethnic.selection, first = 1, last = 1),
                           substring(text = sim.sex.selection, first = 1, last = 1))
        #################
        ## Actual data ##
        #################
        ## Which level is the specific cancer type
        cancerdat.filtered <- cancerdat[cancerdat$type == sim.cancertype.selection & ## filter by cancer type
                                        cancerdat$gender == sim.sex.selection,] ## filter by sex
        ages <- cancerdat.filtered$age_at_initial_pathologic_diagnosis
        
        ####################
        ## Simulated data ##
        ####################
        mixedages = runif(input$sim.population,input$sim.range[1],input$sim.range[2])
        simulated.model <- list()
        for(type in sim.cancertype.selection) {
          ## Getting the ktissue mean for specified cancer
          kmean.tis <- tissue.kvalues[, names(tissue.kvalues) == type]
          kmean.tis <- mean(kmean.tis)
          # ## Getting the kclass means for specified cancer
          # kclnames <- grep("k\\[",names(dfsamps), value=TRUE)
          # dfclas <- dfsamps[,kclnames]
          # dfclas <- stack(dfclas)
          # names(dfclas) <- c("value","coefname")
          # ## Separate coefname to only include class
          # dfclas$class <- gsub(pattern = ".*,", replacement = "", x = dfclas$coefname) 
          # dfclas$class <- gsub(pattern = "\\].*", replacement = "", x = dfclas$class)
          # ## Separate coefname to only include tissue
          # dfclas$tiss <- gsub(pattern = ",.*", replacement = "", x = dfclas$coefname) 
          # dfclas$tiss <- gsub(pattern = ".*\\[", replacement = "", x = dfclas$tiss) 
          # 
          # convert.tiss2names <- data.frame(tiss = as.character(1:23),
          #                                  tissnames = tisnames)
          # 
          # dfclas <- dfclas %>% left_join(., convert.tiss2names, by = "tiss") %>% ## leftjoin data
          #     dplyr::select(-tiss) ## Remove tissue number system
          # rm(convert.tiss2names)
          # save(dfclas, file = "~/TTmanu/code/supplemental_file_shiny/dfclas.rda")
          #load("dfclas.rmd")
          #load("~/TTmanu/code/supplemental_file_shiny/dfclas.rmd")
          ## Getting the kclass mean for specified cluster
          for(tumorclass in sim.cluster.selection) {
            if(tumorclass == "All") {
              kcl <- dfclas[which(dfclas$tissnames %in% type),]
              kmean.cl <- mean(kcl$value)
            } else{
              kcl <- dfclas[which(dfclas$tissnames %in% type &
                                    dfclas$class %in% tumorclass),]
              kmean.cl <- mean(kcl$value)
            }
            
            ## Getting kmean value ##
            kmean <- kmean.tis * kmean.cl
            ## Getting age info
            
            ##Name for cancer type and cluster
            listname <- paste0(type, " Cluster ", tumorclass, "; k = ", round(kmean, digits = 1))
            
            simulated.model[[listname]] <- sapply(mixedages,function(a) tryCatch(simdeath(r20[agemodel],agerate[agemodel],a,kmean)/365,error= function(e) {0}))
          }
        }
        
        simulated.model <- as.data.frame(do.call(cbind, simulated.model))
        simulated.model$status <- 1
        
        simulated.model <- simulated.model %>% pivot_longer(-c(status),
                                        values_to = "time",
                                        names_to = "condition")
        
        survyrs.ref <-  sapply(mixedages,function(a) tryCatch(simdeath(r20[agemodel],agerate[agemodel],a,1)/365,error= function(e) {0}))
        simulated.ref <- data.frame("General Population; k = 1" = survyrs.ref,
                          status = rep(1, length(survyrs.ref)))
        names(simulated.ref) <- c("General Population; k = 1", "status")
        simulated.ref <- simulated.ref %>% pivot_longer(-c(status),
                                   values_to = "time",
                                   names_to = "condition")
        
        if(input$show.sim.k.eq.1 == FALSE) {
        max.x <- max(simulated.model$time) + 5 
        } else {
          simulated.model <- rbind(simulated.ref, simulated.model)
          max.x <- max(simulated.model$time) + 5
        }
        
        if(input$sim.ethnic == "WHITE") {
            ethnic.group <- "European descent"
        } else {
            ethnic.group <- "African descent"
        }
        
        ggplot(simulated.model, aes(time = time, color = condition, status = status)) +
          geom_km() +
          theme_classic() +
          theme(legend.position = c(.8, .75)) +
          coord_cartesian(xlim = c(0,max.x)) +
          scale_color_manual(values = ditto_colours) +
          labs(x = "Survival Time (years past diagnosis)",
               y = "Probabilty",
               colour = "",
               title = sprintf("Survival for simulated %s %s patients", ethnic.group, tolower(sim.sex.selection))) +
          theme(axis.text = element_text(size = 12),
                axis.title = element_text(size = 14))
    })
    
    
})   
shinyApp(ui = ui, server = server)



