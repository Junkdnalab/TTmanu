// A survival based model, using race to select the background risk rate

functions {
  real ourmodel_lpdf(real t, real age, real agerate, real k, real r20 ){
    return( (-(365*r20*exp((agerate*k*t)/365+age*agerate-20*agerate))/(agerate*k))
	    +(agerate*k*t)/365+log(r20)+(365*exp(age*agerate-20*agerate)*r20)/(agerate*k)
	    +age*agerate-20*agerate);
    
  }
  real ourmodel_lccdf(real t, real age, real agerate, real k, real r20){
    return(-r20*((365*exp((agerate*k*t)/365+age*agerate-20*agerate))/(agerate*k)
            -(365*exp(age*agerate-20*agerate))/(agerate*k)));
    
  }
}

// we only really have usable data on white and black races, so we use
// race = 1 for black, and race = 2 for white and split out


data{
  int Np;
  int Nt;
  int Nc;
  int age[Np];
  int race[Np];
  int patient[Np];
  int tissue[Np];
  int tclass[Np];
  real tevent[Np];
  int eventtype[Np];
  real menccdfw[101];
  real womccdfw[101];
  real menccdfb[101];
  real womccdfb[101];
  int gender[Np];
}

parameters{
  real<lower=0.0> r20w[2];
  real<lower=0.0> r20b[2];
  real<lower=0.0> k[Nt,Nc];
  real<lower=0.0> ktis[Nt];
  real<lower=0.0> ageratew[2];
  real<lower=0.0> agerateb[2];
}

model{
  ageratew ~ gamma(2.0,1.0/1.0);
  agerateb ~ gamma(2.0,1.0/1.0);
  r20w ~ gamma(1.5,0.5/1.0);
  r20b ~ gamma(1.5,0.5/1.0);
  for(tis in 1:Nt){
    ktis[tis] ~ gamma(1.2,0.2/1.0);
    for( cl in 1:Nc){
      k[tis,cl] ~ gamma(3.5,2.5/1.0);
    }
  }

  for(ii in 1:6){
    int i;
    i = 65+ii*5;
  // calibrate the age specific risk rate to CDC published life tables, WHITE
    log(womccdfw[i+1]) ~ normal(ourmodel_lccdf(i*365 | 0,ageratew[1]*.05,1,r20w[1]*1e-5),log(1.05));
    log(menccdfw[i+1]) ~ normal(ourmodel_lccdf(i*365 | 0,ageratew[2]*.05,1,r20w[2]*1e-5),log(1.05));
    // calibrate the age specific risk rate to CDC published life tables, BLACK
    log(womccdfb[i+1]) ~ normal(ourmodel_lccdf(i*365 | 0,agerateb[1]*.05,1,r20b[1]*1e-5),log(1.05));
    log(menccdfb[i+1]) ~ normal(ourmodel_lccdf(i*365 | 0,agerateb[2]*.05,1,r20b[2]*1e-5),log(1.05));
  }

  
  
  
  for(i in 1:Np){
    int tis = tissue[i];
    int cl = tclass[i];

    if(race[i] == 2){ // white
      if(eventtype[i] == 1){
	target += ourmodel_lpdf(tevent[i]| age[i],ageratew[gender[i]] * 0.05 ,ktis[tis]*k[tis,cl],r20w[gender[i]] * 1e-5);
      }else{
	target += ourmodel_lccdf(tevent[i]| age[i],ageratew[gender[i]] * 0.05 ,ktis[tis]*k[tis,cl],r20w[gender[i]] * 1e-5);
      }
    }else{ // black
      if(eventtype[i] == 1){
	target += ourmodel_lpdf(tevent[i]| age[i],agerateb[gender[i]] * 0.05 ,ktis[tis]*k[tis,cl],r20b[gender[i]] * 1e-5);
      }else{
	target += ourmodel_lccdf(tevent[i]| age[i],agerateb[gender[i]] * 0.05 ,ktis[tis]*k[tis,cl],r20b[gender[i]] * 1e-5);
      }      
    }
  }
}



