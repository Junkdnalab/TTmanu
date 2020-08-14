library(survival)
library(ggfortify)
library(gridExtra)
library(readr)
library(sqldf)

# run the following if you don't have ggkm
# install.packages("devtools")
# devtools::install_github("sachsmc/ggkm")
library(ggkm)


set.seed(20200805) #a date based seed for reproducibility

## read in estimates from the stan runs, and then generate pseudo K-M curves


## ethnicity based model output for mean r20 values, as printed in the
## Rout file:

## these are scaled in the stan file, so we need to scale here as well
r20scale =1e-5
ageratescale = .05
    

r20wf = 0.365209 * r20scale
r20wm = 0.367732 * r20scale
r20bf = 0.093299 * r20scale
r20bm = 0.098039 * r20scale

ageratewf = 1.197896 * ageratescale
ageratewm = 1.291932 * ageratescale
ageratebf = 1.649803 * ageratescale
ageratebm = 1.734934 * ageratescale


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



## do a survival plot for the age distribution of each cancer using
## ktis from our stan run.  male and female white (black is less well
## estimated due to fewer data points)

load("survmodel-ethbkd.stansave") ## loads a variable called
                                  ## "sampseth" which is a stanfit
## object
class(sampseth)

asampseth = as.array(sampseth)

cancerdat = read_csv("../../data/survival/clinical_plus_cluster.csv");
cancerdat$type = as.factor(cancerdat$type) ## the tissue
endpoints = sqldf("select case when death_days_to is NULL then last_contact_days_to else death_days_to end as finaltime, case when death_days_to is NULL then 0 else 1 end as finalstatus from cancerdat;")

cancerdat = cbind(cancerdat,endpoints);
sum(is.na(cancerdat$finaltime))
cancerdat = cancerdat[!is.na(cancerdat$finaltime),]



plotkmcomp = function(sim,ref,act,actstatus){
    data.frame(
        stime = sim,
        rtime = ref,
        acttime = act,
        status = rep(1, length(ref)),
        actstatus = actstatus
    ) %>%
        pivot_longer(-c(status, actstatus),
            values_to = "time",
            names_to = "condition"
        ) %>%
            mutate(status = case_when(condition != "acttime" ~ as.integer(status),
            condition == "acttime" ~ as.integer(actstatus))) %>%
        ggplot(aes(time = time, color = condition, status = status)) +
        geom_km() +
        theme_classic() +
        theme(legend.position = c(.8, .75))
}
    


plotcancertype = function(cancerdat,i){

    lev = which(levels(cancerdat$type) == i)[[1]]
    
    
    women = cancerdat[cancerdat$type == i & cancerdat$gender=="FEMALE",]
    men = cancerdat[cancerdat$type == i & cancerdat$gender=="MALE",]
    agesw = women$age_at_initial_pathologic_diagnosis
    agesm = men$age_at_initial_pathologic_diagnosis
    kmean = mean(asampseth[,,sprintf("ktis[%d]",lev)])
#    print(sprintf("Cancer %s kmean = %f\n",i,kmean))
#    try(hist(agesw),silent=TRUE)
#    try(hist(agesm),silent=TRUE)

    phis = ggplot(data.frame(age=c(agesw,agesm),
                             sex=as.factor(c(rep("Female",length(agesw)),
                                             rep("Male",length(agesm)))))) +
        geom_histogram(aes(age,color=sex))
    
    survyrswf = sapply(agesw,function(a) tryCatch(simdeath(r20wf,ageratewf,a,kmean)/365,error= function(e) {0}))
    survyrswfref = sapply(agesw,function(a) tryCatch(simdeath(r20wf,ageratewf,a,1)/365,error= function(e) {0}))

    survyrswm = sapply(agesm,function(a) tryCatch(simdeath(r20wm,ageratewm,a,kmean)/365,error= function(e) {0}))
    survyrswmref = sapply(agesm,function(a) tryCatch(simdeath(r20wm,ageratewm,a,1)/365,error= function(e) {0}))
    p1 <- tryCatch(plotkmcomp(survyrswf, survyrswfref, women$finaltime / 365, women$finalstatus) +
        coord_cartesian(xlim = c(0, 80)) +
        labs(
            x = "Survival Time (years past diagnosis)",
            y = "Probability",
            title = sprintf("Survival for %s simulated white female patients", i)
        ) +
        scale_color_manual(
            values = c(
                rtime = "#0072b2",
                stime = "#d55e00",
                acttime = "#009e73"
            ),
            breaks = c("rtime", "stime", "acttime"),
            labels = c(
                "Simulated General Population; k = 1",
                sprintf("Simulated Patient Population; k = %.1f", kmean),
                "Actual Patient Population"
            )
        ),
    error = function(e) {
        print(e)
        ggplot()
    }
    )
    p2 = tryCatch(plotkmcomp(survyrswm,survyrswmref,men$finaltime/365,men$finalstatus)+
                  coord_cartesian(xlim=c(0,80)) +
                  labs(x="Survival Time (years past diagnosis)",y="Probability",
                       title=sprintf("Survival for %s simulated white male patients",i)) +
                  annotate("text",x=80, y=.9,label="Simulated General Population k=1",color="orange",hjust="right") +
                  annotate("text",x=80, y=.7,label=sprintf("Simulated Patient Population k=%.1f",kmean),color="red",hjust="right")+
                  annotate("text",x=80, y=.5,label=sprintf("Actual Patient Population",kmean),color="blue",hjust="right")                  
                 ,
                  error = function(e) {print(e);ggplot()})

    list(p1,p2,phis)
    
}




## pdf("kmexamples.pdf")
## j=0
## for (i in levels(cancerdat$type)){
##     j = j+1
##     plotcancertype(cancerdat,i)
## }
## dev.off()

## work on creating a plot using patchwork


## plot STAD simulated/actual male/female
## plot THCA simulated/actual male/female

