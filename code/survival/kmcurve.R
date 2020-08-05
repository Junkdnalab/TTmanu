library(survival)
library(ggfortify)
library(gridExtra)
library(readr)


install.packages("devtools")
devtools::install_github("sachsmc/ggkm")
library(ggkm)

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

curve(lccdf(r20wf,ageratewf,30,1,x)-log(runif(1)),0,365*100)

simdeath = function(r20,agerate,age,k){
    r = runif(1)
    lr = log(r)
    tdays = uniroot(function(t) return(lccdf(r20,agerate,age,k,t)-lr),c(0,365*400))
    return(tdays$root) ## years since diagnosis at "age"
}

## create example data for a fixed 30 year old without vs with a k=2
## cancer


patientsurvivalsbkd = replicate(100,simdeath(r20wf,ageratewf,30,1))/365
patientsurvivalsk2 = replicate(100,simdeath(r20wf,ageratewf,30,2))/365


fit1 <- survfit(Surv(patientsurvivalsbkd) ~ 1) 
fit2 <- survfit(Surv(patientsurvivalsk2) ~ 1) 
p1 = autoplot(fit1) + coord_cartesian(xlim=c(0,80)) + labs(x="Survival Time (years past 30th birthday)",y="Probability",title="Survival For 30 year old White Females (baseline)")
p2 = autoplot(fit2) + coord_cartesian(xlim=c(0,80)) + labs(x="Survival Time (years past 30th birthday)",y="Probability",title="Survival For 30 year old White Females with k = 2 cancer") 
grid.arrange(p1,p2)


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

plotcancertype = function(cancerdat,i){
    agesw = cancerdat[cancerdat$type == i & cancerdat$gender=="FEMALE","age_at_initial_pathologic_diagnosis"]$age_at_initial_pathologic_diagnosis
    agesm = cancerdat[cancerdat$type == i & cancerdat$gender=="MALE","age_at_initial_pathologic_diagnosis"]$age_at_initial_pathologic_diagnosis
    kmean = mean(asampseth[,,sprintf("ktis[%d]",j)])
    print(sprintf("Cancer %s kmean = %f\n",i,kmean))
    try(hist(agesw),silent=TRUE)
    try(hist(agesm),silent=TRUE)

    survyrswf = sapply(agesw,function(a) tryCatch(simdeath(r20wf,ageratewf,a,kmean)/365,error= function(e) {0}))
    survyrswfref = sapply(agesw,function(a) tryCatch(simdeath(r20wf,ageratewf,a,1)/365,error= function(e) {0}))
    survyrswm = sapply(agesm,function(a) tryCatch(simdeath(r20wm,ageratewm,a,kmean)/365,error= function(e) {0}))
    survyrswmref = sapply(agesm,function(a) tryCatch(simdeath(r20wm,ageratewm,a,1)/365,error= function(e) {0}))
#    kmf = try(survfit(Surv(survyrswf)~1),silent=TRUE)
#    kmfref = try(survfit(Surv(survyrswfref)~1),silent=TRUE)
#    kmm = try(survfit(Surv(survyrswm)~1),silent=TRUE)
#    kmmref = try(survfit(Surv(survyrswmref)~1),silent=TRUE)
    p1 = tryCatch(
        ggplot(data.frame(stime=survyrswf,rtime=survyrswfref,status=rep(1,length(survyrswf)))) + geom_km(aes(time=stime,status=status),color="red") + geom_km(aes(time=rtime,status=status),color="orange") + coord_cartesian(xlim=c(0,80)) +
        labs(x="Survival Time (years past diagnosis)",
             y="Probability",
             title=sprintf("Survival for %s simulated white female patients",i)) +
        annotate("text",x=50, y=.8,label="General Population k=1",color="orange") +
        annotate("text",x=50, y=.7,label=sprintf("Patient Population k=%.1f",kmean),color="red"),
        error= function(e) {ggplot()})
    p2 = tryCatch(
        ggplot(data.frame(stime=survyrswm,rtime=survyrswmref,status=rep(1,length(survyrswm)))) + geom_km(aes(time=stime,status=status),color="red") + geom_km(aes(time=rtime,status=status),color="orange") + coord_cartesian(xlim=c(0,80)) +
        labs(x="Survival Time (years past diagnosis)",y="Probability",
             title=sprintf("Survival for %s simulated white male patients",i)) +
        annotate("text",x=50, y=.8,label="General Population k=1",color="orange") +
        annotate("text",x=50, y=.7,label=sprintf("Patient Population k=%.1f",kmean),color="red"),
        error = function(e) {ggplot()})
    
    tryCatch(print(grid.arrange(p1,p2)),error=function(e) {print(ggplot())})
    
}




pdf("kmexamples.pdf")
j=0
for (i in levels(cancerdat$type)){
    j = j+1
    plotcancertype(cancerdat,i)
}
dev.off()
