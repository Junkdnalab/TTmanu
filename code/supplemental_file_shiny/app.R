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
    library(car),
    library(grid),
    library(gridBase)
))

## Load data
umap_cancertype <- fread("umap_3d_coors.tsv")
load("kelly.colours.rda")
load("ditto.colours.rda")
load("tumor_samples.Rda")
load("supplemental_gene.Rda")
load("gdc_reactome_path_analysis.Rda")

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
                             rglwidgetOutput("pancan",  width = 500, height = 500)),
                    tabPanel(title = "Cancer Type",
                             
                             pickerInput("cancer_type", label = "Cancer type:", ## Cancer type
                                         choices = cancer_type, selected = "BRCA", multiple = TRUE,
                                         options = list(`actions-box` = TRUE)),
                             plotOutput("cancertype")),
                    tabPanel(title = "Gene", 
                             selectInput("genename", label = "Gene:", ## Gene 
                                         choices = genename, selected = "TP53"),
                             plotOutput("gene")),
                    tabPanel(title = "Pathway", 
                             selectInput("pathwayname", label = "Pathway:", ## Pathway
                                         choices = pathwayname, selected = "Cellular response to hypoxia"),
                             plotOutput("pathway")))
        )
    )
)


server <- (function(input, output, session) {
    ## Pancancer plot code
    output$pancan <- renderRglwidget ({ 
        open3d(useNULL=TRUE)
        plot3d(umap_cancertype[, c("plot_x", "plot_y", "plot_z"),],
               type = "s", size = 2.5,
               xlab = "UMAP 1", ylab = "UMAP 2", zlab = "UMAP 3",
               aspect = c(8.3, 10.5, 7.6),
               col = cluster.colours[as.character(umap_cancertype$clust_knn)])["data"]
        rglwidget()
        
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
                           bg = s_color, pch = 21, cex = 0.8, lwd = 0.2,
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
        gene.df <- supplemental_gene %>% dplyr::select(input$genename) %>% ## filtering for gene
            "colnames<-" ("selectedgene") %>% ## rename column
            rownames_to_column(var = "sample_id") %>% ## move sample id to column
            left_join(umap_cancertype[,c("plot_x", "plot_y", "plot_z", "clust_knn", "sample_id")], ., by = "sample_id") %>% ## left join data by sample_id
            filter(selectedgene == 1)
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
       pathway.df <- gdc_reactome_path_analysis %>% t() %>% as.data.frame() %>% ## transpose data for left_joining later
            dplyr::select(input$pathwayname) %>% ## filtering for pathway
            "colnames<-" ("selectedpath") %>% ## rename column
            rownames_to_column(var = "sample_id") %>% ## move sample id to column
            left_join(umap_cancertype[,c("plot_x", "plot_y", "plot_z", "clust_knn", "sample_id")], ., by = "sample_id") %>% ## left join data by sample_id
            filter(selectedpath == 1)
        with(pathway.df, ## only select tumor samples that have a mutation in the selected gene
             scatter3D(plot_y, -plot_z, plot_x, 
                       bg = cluster.colours[as.character(pathway.df$clust_knn)], pch = 21, cex = 0.8, lwd = 0.2,
                       theta = 0, phi = 65, scale = F,
                       xlim = c(-6, 4.5), zlim = c(-4.5, 4), ylim = c(-4, 5),
                       xlab = "UMAP 2", ylab = "UMAP 3", zlab = "UMAP 1",
                       colvar = NULL))
    })
    
})   
shinyApp(ui = ui, server = server)



