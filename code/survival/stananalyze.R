library(sqldf)
library(readr)
library(tidyverse)
library(rstan)
library(foreach)
library(bayesplot)
library(gridExtra)
library(openxlsx)
library(ggridges)

cancerdat = read_csv("../../data/survival/clinical_plus_cluster.csv");
treatdays = read_csv("../../data/survival/tcga_clinical_data.csv")
# knnclust = read_tsv("../../data/survival/peter_clusters_knn.tsv")

## mortality data from CDC:
namesmort = c("agerng","pdead","nsurv","ndie","pyrs","pyrsabove","expectedlife")
menmort = read_csv("../../data/survival/lifetables/Table02.csv",skip=2)
## drop the final row, it's just a text info statement
menmort = menmort[1:NROW(menmort)-1,] 
wommort = read_csv("../../data/survival/lifetables/Table03.csv",skip=2)
wommort = wommort[1:NROW(wommort)-1,]

names(menmort) = namesmort
names(wommort) = namesmort
yrs=seq(0,100,by=1);
menmort$yrnum = yrs
wommort$yrnum = yrs

ggplot(menmort) + geom_line(aes(yrnum,nsurv/100000),col="blue") +geom_line(aes(yrnum,nsurv/100000),col="red",data=wommort)


sqldf("select count(*) from cancerdat where age_at_initial_pathologic_diagnosis < 20")
sqldf("select count(*) from cancerdat where age_at_initial_pathologic_diagnosis < 10")

sqldf("select type,count(*) from cancerdat group by type");
sqldf("select type,avg(age_at_initial_pathologic_diagnosis) as meanage from cancerdat group by type");
#sqldf("select type,count( as meanage from cancerdat group by type");



## Create a dataset merging various components



deadpat = sqldf("select patient_id,age_at_initial_pathologic_diagnosis as age, gender,race, death_days_to,type,clust_knn as clust, clust as origclust,`ajcc_pathologic_tumor_stage.y` as stage, last_contact_days_to as lastcont from cancerdat where (death_days_to is not null or lastcont is not null) and age is not null and age > 9")

unique(deadpat$race)
deadpat[deadpat$race == "BLACK OR AFRICAN AMERICAN","race"] = "BLACK"
deadpat[grep("\\[",deadpat$race),"race"] = "UNK"
deadpat[grep("PACIFIC",deadpat$race),"race"] = "PACIFIC"
deadpat[grep("INDIAN",deadpat$race),"race"] = "AMERIND"
deadpat$race = as.factor(deadpat$race)
levels(deadpat$race)

ggplot(deadpat,aes(race))+geom_bar()



## we exclude young children (there are none in this dataset at this
## time) because infant mortality has a separate shape to it in the
## first few years of life



deadpat$patient_id=as.factor(deadpat$patient_id)
deadpat$type=as.factor(deadpat$type)
deadpat$clust=deadpat$clust; # stan uses 1 based indexing, our new knn_clust also does
deadpat$tevent=deadpat$death_days_to
deadpat[is.na(deadpat$tevent),"tevent"] = deadpat[is.na(deadpat$tevent),"lastcont"]
deadpat[is.na(deadpat$death_days_to),"eventtype"] = "bfollow"
deadpat[!is.na(deadpat$death_days_to),"eventtype"] = "adeath"
deadpat$eventtype = as.factor(deadpat$eventtype)
deadpat$gender = as.factor(deadpat$gender)
levels(deadpat$gender)

standata = list(Np=NROW(deadpat),Nt=NROW(levels(deadpat$type)),Nc=NROW(unique(deadpat$clust)),
                patient=as.numeric(deadpat$patient_id),age=deadpat$age,tissue=as.numeric(deadpat$type),
                menccdf=menmort$nsurv/100000,
                womccdf=wommort$nsurv/100000,
                gender=as.numeric(deadpat$gender),
                tclass=as.integer(deadpat$clust),tevent = deadpat$tevent,
                eventtype=as.numeric(deadpat$eventtype))

### Fit the survival model with one rate per tissue*class
smodel = stan_model(file="survmodel.stan")


## optimization for initialization doesn't always work... 
start = optimizing(smodel,data=standata);
start2 = optimizing(smodel,data=standata);
start3 = optimizing(smodel,data=standata);

tryCatch(load("survmodel.stansave"),error=function(e){
    samps = sampling(smodel,data=standata,chains=3,warmup=2000,iter=3000,cores=4,control=list(max_treedepth=9)) #,init=list(start,start2,start3));
    save("samps",file="survmodel.stansave")
})





######################
## Fit a factored model, where there's one ktis per tissue, and one kcls per class


smodelf = stan_model(file="survmodelfactor.stan")
start = optimizing(smodelf,data=standata);
start2 = optimizing(smodelf,data=standata);
start3 = optimizing(smodelf,data=standata);


tryCatch(load("survmodel-factor.stansave"),error=function(e){
    sampf = sampling(smodelf,data=standata,chains=3,warmup=500,iter=1000,cores=4,control=list(max_treedepth=12),init=list(start,start2,start3));
    save("sampf", file="survmodel-factor.stansave")
})
load("survmodel-factor.stansave")

sampfsum = summary(sampf)
sampfsum

asampsf = as.array(sampf)


pdf("agebasedmodel-fact.pdf")

mcmc_intervals(asampsf,regex_pars=c("kcls"))
mcmc_intervals(asampsf,regex_pars=c("ktis"))
mcmc_intervals(asampsf,regex_pars=c("r20","agerate"))

r20f = mean(asampsf[,,"r20[1]"])
r20m = mean(asampsf[,,"r20[2]"])
agrf = mean(asampsf[,,"agerate[1]"])
agrm = mean(asampsf[,,"agerate[2]"])


sprintf("r20f mean = %f",r20f)
sprintf("r20m mean = %f",r20m)
sprintf("agrf mean = %f", agrf)
sprintf("agrm mean = %f", agrm)

ourlccdf = function(t,age,agerate,k,r20){
    r20 = r20 * 1e-5
    agerate = agerate*.05
    return(-r20*((365*exp((agerate*k*t)/365+age*agerate-20*agerate))/(agerate*k)
        -(365*exp(age*agerate-20*agerate))/(agerate*k)));
}

ccdfsimpmkr = function(r20,agert){
    function(age){
        exp(ourlccdf(age*365,0,agert,1,r20))
    }
}

mccdf = ccdfsimpmkr(r20m,agrm)
fccdf = ccdfsimpmkr(r20f,agrf)

ggplot(menmort) + geom_line(aes(yrnum,log(nsurv/100000)),col="blue",data=menmort) +
    geom_line(aes(yrnum,log(nsurv/100000)),col="red",data=wommort)+
    stat_function(fun=function(x) log(mccdf(x)),col="green")+
    stat_function(fun=function(x) log(fccdf(x)),col="orange")+labs(title="Log(1-cdf(age)) for Men and Women (CDC and model fit)",x="Age (yrs)",y="log(1-cdf)")


## class 3 and 13 look to be more dangerous, is that because they're
## mostly found in more dangerous tissues?

for(i in 1:13){
    print(ggplot(deadpat[deadpat$clust == i,],aes(type))+geom_bar() +labs(title=paste("Class =",i))+ theme(axis.text.x=element_text(angle=90)))
}
    
dev.off()
#system("evince agebasedmodel-fact.pdf&")


## pdf("factoredModelr0tc.pdf")
## mcmc_dens(asampsf,regex_pars=c("r0c"))+coord_cartesian(xlim=c(0,4));
## mcmc_dens(asampsf,regex_pars=c("r0t"))+coord_cartesian(xlim=c(0,4));
## dev.off();



###### Fit a model split out by black vs white race. These are the
###### only two categories where we have both data in our cancer
###### database, as well as a CDC life table specific to that race.


namesmort = c("agerng","pdead","nsurv","ndie","pyrs","pyrsabove","expectedlife")
wmenmort = read_csv("../../data/survival/lifetables/Table05.csv",skip=2)
## drop the final row, it's just a text info statement
wmenmort = wmenmort[1:NROW(wmenmort)-1,] 
wwommort = read_csv("../../data/survival/lifetables/Table06.csv",skip=2)
wwommort = wwommort[1:NROW(wwommort)-1,]

bmenmort = read_csv("../../data/survival/lifetables/Table05.csv",skip=2)
## drop the final row, it's just a text info statement
bmenmort = bmenmort[1:NROW(bmenmort)-1,] 
bwommort = read_csv("../../data/survival/lifetables/Table06.csv",skip=2)
bwommort = bwommort[1:NROW(bwommort)-1,]




names(wmenmort) = namesmort
names(bmenmort) = namesmort
names(wwommort) = namesmort
names(bwommort) = namesmort
yrs=seq(0,100,by=1);
wmenmort$yrnum = yrs
bmenmort$yrnum = yrs
wwommort$yrnum = yrs
bwommort$yrnum = yrs

ethndat = deadpat[deadpat$race %in% c("WHITE","BLACK"),]
ethndat$race = as.factor(as.character(ethndat$race))


standataeth = list(Np=NROW(ethndat),Nt=NROW(levels(ethndat$type)),Nc=NROW(unique(ethndat$clust)),
                patient=as.numeric(ethndat$patient_id),age=ethndat$age,tissue=as.numeric(ethndat$type),
                menccdf=menmort$nsurv/100000,
                womccdf=wommort$nsurv/100000,
                gender=as.numeric(ethndat$gender), race = as.numeric(ethndat$race),
                tclass=as.integer(ethndat$clust),tevent = ethndat$tevent,
                eventtype=as.numeric(ethndat$eventtype))


standataeth$menccdfw=wmenmort$nsurv/1e5
standataeth$menccdfb=bmenmort$nsurv/1e5
standataeth$womccdfw=wwommort$nsurv/1e5
standataeth$womccdfb=bwommort$nsurv/1e5

sethmodel = stan_model(file="survmodelethbkd.stan")

#start = optimizing(sethmodel,data=standataeth,algorithm="Newton",verbose=TRUE);
#start2 = optimizing(sethmodel,data=standataeth,algorithm="Newton",verbose=TRUE);
#start3 = optimizing(sethmodel,data=standataeth,algorithm="Newton",verbose=TRUE);

#starts = vb(sethmodel,data=standataeth,importance_resampling=TRUE,adapt_iter=10000,iter=10000)

#sampseth = sampling(sethmodel,data=standataeth,chains=3,warmup=500,iter=1000,cores=4,control=list(max_treedepth=12),init=list(start,start2,start3));


tryCatch(load("survmodel-ethbkd.stansave"),error=function(e){
    sampseth = sampling(sethmodel,data=standataeth,chains=3,warmup=2000,iter=4000,thin=2,cores=4,control=list(max_treedepth=12))

    save("sampseth",file="survmodel-ethbkd.stansave")
})
load("survmodel-ethbkd.stansave")

asampseth = as.array(sampseth)

pdf("ethnicparams.pdf")
mcmc_intervals(asampseth,regex_pars=c("ktis"))

for(i in 1:23){
    print(mcmc_intervals(asampseth,regex_pars=paste("k\\[",i,",",sep="")))
}

print(mcmc_intervals(asampseth,regex_pars="(r20)|(agera)"))


#ccdf


dev.off()


sprintf("r20wf mean = %f",mean(asampseth[,,"r20w[1]"]))
sprintf("r20wm mean = %f",mean(asampseth[,,"r20w[2]"]))
sprintf("r20bf mean = %f",mean(asampseth[,,"r20b[1]"]))
sprintf("r20bm mean = %f",mean(asampseth[,,"r20b[2]"]))

sprintf("agrwf mean = %f", mean(asampseth[,,"ageratew[1]"]))
sprintf("agrwm mean = %f", mean(asampseth[,,"ageratew[2]"]))
sprintf("agrbf mean = %f", mean(asampseth[,,"agerateb[1]"]))
sprintf("agrbm mean = %f", mean(asampseth[,,"agerateb[2]"]))










pdf("tissuerates.pdf",width=32,height=16)
nclass = length(unique(deadpat$clust))
ntype = length(levels(deadpat$type))
plotlist = list()
for (i in 1:ntype){
    df = data.frame(r0t=as.vector(asamps[,,sprintf("r0t[%d]",i)]))
    plotlist = append(plotlist,list(ggplot(df,aes(r0t))+geom_density(fill="red")+coord_cartesian(xlim=c(0,6),ylim=c(0,2))+labs(x=levels(deadpat$type)[i])))
    for(j in 1:nclass){
        df = data.frame(rate=as.vector(asamps[,,sprintf("r0[%d,%d]",i,j)]))
        plotlist = append(plotlist,list(ggplot(df,aes(rate))+geom_density(fill="blue")+labs(x=sprintf("%s type %d",levels(deadpat$type)[i],j))+coord_cartesian(xlim=c(0,4),ylim=c(0,2))))
    }
}
grid.arrange(grobs=plotlist,nrow=nclass+1,ncol=ntype,as.table=FALSE)
dev.off()

### Alternate viz

asamps <- rstan::extract(samps) # `extract()` conflicts: rstan & tidyverse

bsamps <- foreach(i = 1:length(levels(factor(deadpat$clust))), .combine = 'rbind') %do% cbind(
    gather(
        as.data.frame(asamps$k[,,i]),
        key = "cancer.type",
        value = "k"),
    clust = rep(
        paste("clust_", as.numeric(levels(factor(deadpat$clust))[i])-1, sep = ''), 
        dim(asamps$k[,,i])[1]))


for (i in 1:length(levels(deadpat$type))) {
    bsamps$cancer.type[bsamps$cancer.type==paste("V", i, sep = '')] <- as.character(levels(deadpat$type)[i])
}


csamps <- subset(bsamps, bsamps$clust != "clust_0")

tisplot2 <- ggplot(csamps) + 
    geom_density_ridges(
        aes(x=k, y=clust, fill=clust),
        panel_scaling = FALSE) +
    facet_wrap(. ~ cancer.type, nrow = 3) +
    theme_minimal() +
    scale_fill_brewer(palette = "Paired") +
    xlim(0, 3) +
    theme(legend.position = "none") +
    ggtitle("tissue by class")


ggsave(filename = "tissuerates_agebased.pdf", plot = tisplot2, height = 8.5, width = 11)
ggsave(filename = "tissuerates_agebased.png", plot = tisplot2, height = 8.5, width = 11, dpi = 300)


### cancer specific 

dsamps <- as.data.frame(asamps$ktis)
names(dsamps) <- levels(deadpat$type)
dsamps <- gather(dsamps, key="cancer_type", value="ktis")
tisplot <- ggplot(dsamps, aes(y=cancer_type, x=ktis)) + 
    geom_density_ridges() +
    theme_minimal() +
    xlim(0, 22) + 
    xlab(expression("k"["tis"])) +
    ylab("cancer type") +
    ggtitle("Tissue specific effects")

ggsave(filename = "tissue_effects.pdf", plot = tisplot, height = 6, width = 4)
ggsave(filename = "tissue_effects.png", plot = tisplot, height = 6, width = 4, dpi = 300)





#### debugging fits

load("survmodel_effcl30.stansave")

df = as.data.frame(samps)
dfs = stack(df)

comp = sqldf("select a.rowid,patient_id,age,gender,type,clust,tevent,eventtype,b.`values` as tdeath from deadpat a join dfs b on b.ind == 'tdeath[' || a.rowid || ']'");


myplot = ggplot(comp) + facet_wrap(~type) + stat_ecdf(aes(x=tdeath),color="blue") + stat_ecdf(aes(x=tevent),data=subset(comp,eventtype=="adeath"),col="green") + stat_ecdf(aes(x=tevent),data=subset(comp,eventtype=="bfollow"),col="orange")+coord_cartesian(xlim=c(0,365*20))+labs(title="New Model ECDF Comparisons")

ggsave("newmodel.png",myplot,width=10,height=10)



ggplot(comp) + geom_density(aes(x=tdeath,fill=type),alpha=.2)+coord_cartesian(xlim=c(0,365*20))





myplot = ggplot(comp) + facet_wrap(~type) + geom_density(aes(x=tdeath,fill="Predicted"),alpha=.5) + geom_density(aes(x=tevent,fill="Actual"),alpha=.5)+coord_cartesian(xlim=c(0,365*20))



ggsave("/tmp/densities.png",myplot)
system("gthumb /tmp/densities.png")




## ecdfs:

ggplot(comp) + stat_ecdf(aes(x=tdeath)) + stat_ecdf(aes(x=tevent))+coord_cartesian(xlim=c(0,365*20))+facet_wrap(~type)

ggplot(comp) + stat_ecdf(aes(x=tdeath)) + stat_ecdf(aes(x=tevent),data=comp[comp$eventtype=="bfollow",])+coord_cartesian(xlim=c(0,365*20))+facet_wrap(~type)

ggplot(comp) + stat_ecdf(aes(x=tdeath)) + stat_ecdf(aes(x=tevent),data=comp[comp$eventtype=="adeath",])+coord_cartesian(xlim=c(0,365*20))+facet_wrap(~type)









### simon's suggested plots


library(bayesplot)
library(patchwork)

load("tissue_key.rda")

outdata <- tibble(as.data.frame(standata[c("age","tissue","gender","tclass","tevent","eventtype")]))
y <- outdata$tevent



yrep <- (as.data.frame(samps)) %>% select(starts_with("tdeath"))
yrep <- janitor::clean_names(yrep)
yrep <- as.matrix(yrep)
subset = sample(1:NROW(yrep),25)
plotlist_dead <- list()
for(tumortissue in tissue_key$tissue_num) {
    samples <- which(outdata$tissue == tumortissue & outdata$eventtype == 1)
    mysubplt <- ppc_dens_overlay(y[samples], yrep[subset,samples],n_dens=128) + coord_cartesian(xlim=c(0,9000)) + ggtitle(paste(tissue_key %>% filter(tissue_num == tumortissue) %>% pull(tissue), "n", "=", length(samples)))
                                                                               plotlist_dead <- c(plotlist_dead, list(mysubplt))
    }
myplt1 = wrap_plots(plotlist_dead)
    
plotlist_alive <- list()
for(tumortissue in tissue_key$tissue_num) {
  samples <- which(outdata$tissue == tumortissue & outdata$eventtype != 1)
  mysubplt <- ppc_dens_overlay(y[samples], yrep[subset,samples],n_dens=128) + coord_cartesian(xlim=c(0,9000)) + ggtitle(paste(tissue_key %>% filter(tissue_num == tumortissue) %>% pull(tissue),
                                                                           "n", "=", length(samples)))
  plotlist_alive <- c(plotlist_alive, list(mysubplt))
}
myplt2 = wrap_plots(plotlist_alive)

pdf("comparison.pdf",10,10)
print(myplt1)
print(myplt2)
dev.off()

#### KM tests


testdf = data.frame(t=rgamma(1000,3,2/1000),ID=1:1000,eventtype=c(rep("adeath",times=200),rep("bfollow",times=800)))
testdf[testdf$eventtype == "bfollow","t"] = testdf[testdf$eventtype == "bfollow","t"] - rexp(800,1/10)
testdf$isdeath = sapply(testdf$eventtype,FUN=function(x) if(x == "bfollow"){return(0)}else{return(1)})

dp = ggplot(testdf[testdf$eventtype == "adeath",])+geom_km(aes(time=t,status=isdeath))+labs(title="Deaths Only")+coord_cartesian(xlim=c(0,10000))

ap = ggplot(testdf) + geom_km(aes(time=t,status=isdeath))+labs(title="Deaths and Followups")+coord_cartesian(xlim=c(0,10000))

wrap_plots(dp,ap)



### check out the older model:

load("survmodel_knn3525_old.stansave")
ls()
df=as.data.frame(samps)
dfs = stack(df)

comp = sqldf("select a.rowid,patient_id,age,gender,type,clust,tevent,eventtype,b.`values` as tdeath from deadpat a join dfs b on b.ind == 'tdeath[' || a.rowid || ']'");


myplot = ggplot(comp) + facet_wrap(~type) + stat_ecdf(aes(x=tdeath),color="blue") + stat_ecdf(aes(x=tevent),data=subset(comp,eventtype=="adeath"),col="green") + stat_ecdf(aes(x=tevent),data=subset(comp,eventtype=="bfollow"),col="orange")+coord_cartesian(xlim=c(0,365*20))+labs(title="Old Model ECDF Comparisons")

ggsave("oldmodel.png",myplot,width=10,height=10)
