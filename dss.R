
print('hi dss there');
dss<-function(ic50,slope,max,min.conc.tested,max.conc.tested,y=10,DSS.type=2,concn_scale=1e-9){
#rdata should be in in format containing IC50, SLOPE, MAX,MIN.Concentration,MAX.Concentration

  a=as.numeric(unname(max))

  b=as.numeric(unname(slope))
  d=0 # min response
  ic50 = as.numeric(unname(ic50))
  min.conc.tested = as.numeric(unname(min.conc.tested))
  max.conc.tested = as.numeric(unname(max.conc.tested))
  Min.Conc<- log10(min.conc.tested*concn_scale) #
  Max.Conc<- max.conc.tested
  x2<-log10(Max.Conc*concn_scale)  


  if(is.na(ic50)||is.na(b)||is.na(a)||is.na(Min.Conc)||is.na(Max.Conc)){
    dss<-NA
  }
  else if(isTRUE(ic50>=Max.Conc)){
    dss<-0
  }
  else if(isTRUE(b==0)){
    dss<-0
  }
  else{
  if(a>100){ a<-100  }
  if(isTRUE(b<0)){ b<--b  }
    c<-log10(ic50*concn_scale)
    if(a>y){
      if(y!=0){
        x1<-(c - ((log(a-y)-log(y-d))/(b*log(10))))
        if(isTRUE(x1 < Min.Conc)){x1<-Min.Conc}
  		else if(isTRUE(x1 > x2)){x1<-x2}
        }
      else {x1<-Min.Conc}

    # This is a logistic function used in Dotmatics.com
    # y = d+(a-d)/(1+10^(b*(c-x)))
    #inverse function
    # x = c - ((log(a-y)-log(d-y))/(b*log(10)))
      
    int_y=(((((a-d)*log(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1)) - (y*(x2-x1))

       total_area<-(x2-Min.Conc)*(100-y)

    if(DSS.type==1){
    	norm_area<-((int_y/total_area)*100)#DSS1
    }
    if(DSS.type==2){
#	if(a>100){a<-100}
    	norm_area<-((int_y/total_area)*100)/log10(a)#DSS2 #AUC1
        if(isTRUE(norm_area > 50)){ norm_area <- 0}
    }
    if(DSS.type==3){
#	if(a>100){a<-100}
    	norm_area<-((int_y/total_area)*100)*(log10(100)/log10(a))*((x2-x1)/(x2-Min.Conc)) #DSS3 #AUC5
    }
    if(isTRUE(norm_area < 0|norm_area > 100)){
        dss<-0
      }else{
        dss<-round(norm_area,digits=4)}
    } else {dss<-0} 
   } 
return (dss)
}

