# load packages
lapply(c("drc","ggplot2","reshape2","xgboost","parallel", "openxlsx", "dplyr", "NNLM", "modeest", "mlrMBO", "smoof", "kriging"), library, character.only = !0)

# remove outliers
outlier_remove <- compiler::cmpfun(function(xTmp, iqr_ = 1.5){
  qq <- unname(quantile(xTmp, probs=c(.25, .75), na.rm = T))
  outlier_detector <- iqr_ * IQR(xTmp, na.rm = T)
  xTmp < (qq[1] - outlier_detector) | xTmp > (qq[2] + outlier_detector)
})


# fit single agent dose-response curve
CALC_IC50_EC50_DSS = function(xpr_tbl, DSS_typ, readoutCTX = F, drug_name="drug_name")
{
  tryCatch({
    
    mat_tbl <- data.frame(inhibition=as.numeric(xpr_tbl), dose = as.numeric(names(xpr_tbl))); #dose = 10**(1:length(xpr_tbl)));
    mat_tbl$logconc = log10(mat_tbl$dose); mat_tbl$viability = 100 - mat_tbl$inhibition;
    mat_tbl$inhibition2 = mat_tbl$inhibition; mat_tbl$viability2 = mat_tbl$viability;
    mat_tbl <- mat_tbl[order(mat_tbl[,"dose"]),] 
    
    if(any(duplicated(mat_tbl$inhibition))) mat_tbl$inhibition <- seq(from = 0, length.out = length(mat_tbl$inhibition), by = 0.01) + mat_tbl$inhibition; 
    
    
    estimate_param <- tryCatch({drm(inhibition ~ logconc, data = mat_tbl, fct = LL.4(fixed = c(NA, NA, NA,NA),names = c("SLOPE","MIN","MAX","IC50")),logDose=10,control = drmc(errorm = F))}, 
                               warning=function(w){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)},
                               error=function(e){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)})
    coef_estim <- coef(estimate_param); names(coef_estim) <- c("SLOPE","MIN","MAX","IC50"); coef_estim["SLOPE"] <- coef_estim["SLOPE"]*-1 
    coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=T), max(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
    coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=T),coef_estim["IC50"]); coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
    coef_estim["IC50"] <- log10(coef_estim["IC50"]); coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
    coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition<0),max(mat_tbl$logconc,na.rm=T),coef_estim["IC50"]); coef_estim["MIN"] <- 0; coef_estim["MAX"] <- max(mat_tbl$inhibition,na.rm=T)
    min_lower <- ifelse(min(mat_tbl$inhibition,na.rm=T) > 0,min(mat_tbl$inhibition,na.rm=T),0); min_lower <- ifelse(min_lower >= 100,99,min_lower)
    coef_estim["MAX"] <- ifelse(coef_estim["MAX"]>100,100,coef_estim["MAX"]); coef_estim["MAX"] <- ifelse(coef_estim["MAX"]<0,100,coef_estim["MAX"])
    max_lower <- ifelse(max(mat_tbl$inhibition,na.rm=T)>100,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T)); max_lower <- ifelse(max_lower < 0,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T)); 
    max_lower <- ifelse(max_lower < 0,0,max_lower); max_lower <- ifelse(max_lower > 100,100,max_lower); run_avg <- caTools::runmean(mat_tbl$inhibition, 10); max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)]>run_avg[nrow(mat_tbl)]),max(mat_tbl$inhibition[run_avg>run_avg[nrow(mat_tbl)]]),coef_estim["MAX"])
    max_upper <- ifelse(any(mat_tbl$inhibition > max_upper),mean(mat_tbl$inhibition[mat_tbl$inhibition > max_upper])+5,max_upper)
    max_upper <- ifelse(max_upper < 0,coef_estim["MAX"],max_upper); max_upper <- ifelse(max_upper > 100,100,max_upper) #coef_estim["MAX"]
    max_upper <- ifelse(max_lower > max_upper,coef_estim["MAX"],max_upper);mean_inh_last = mean(tail(mat_tbl$inhibition,2),na.rm=T)
    if(mean_inh_last < 60) {
      if(mean_inh_last > 25) coef_estim["IC50"] <- mean(mat_tbl$logconc,na.rm=T) else if(mean_inh_last < 25) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)}
    if(mean(mat_tbl$inhibition[1:3],na.rm=T)<5) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)
    if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + 0.001
    
    #adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
    nls_result_ic50_old <- function(){
      tryCatch({
        nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]),lower=list(SLOPE=0,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),upper=list(SLOPE=4,MIN=0,MAX=100, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      }, error = function(e) {minpack.lm::nlsLM(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,start=list(SLOPE=1, MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]),lower=c(SLOPE=0, MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),upper=c(SLOPE=4, MIN=0,MAX=100, IC50=max(mat_tbl$logconc)))
      })} 
    # IC50 first
    nls_result_ic50 <- nls_result_ic50_old();
    nls_result_ic50_2 <- tryCatch({
      nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]], IC50=median(mat_tbl$logconc)),lower=list(SLOPE=0,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),upper=list(SLOPE=4,MIN=0,MAX=100, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
    },warning = function(w) {nls_result_ic50_old()},error = function(e) {nls_result_ic50_old()})
    
    tryCatch({
      aaa=tryCatch({summary(nls_result_ic50)},error=function(e){summary(nls_result_ic50_2)})
      bbb=tryCatch({summary(nls_result_ic50_2)},error=function(e){summary(nls_result_ic50)})
      
      sumIC50 = list(aaa,bbb)
      ic50std_resid <- round(sqrt(sum((sumIC50[[1]]$residuals)^2)/(length(sumIC50[[1]]$residuals)-1)),1);
      ic50std_resid2 <- round(sqrt(sum((sumIC50[[2]]$residuals)^2)/(length(sumIC50[[2]]$residuals)-1)),1);
      # continue with the best
      switch_ = which.min(c(ic50std_resid, ic50std_resid2)); nls_result_ic50 = list(nls_result_ic50, nls_result_ic50_2)[[switch_]]
    }, error = function(e){})
    
    if(coef(nls_result_ic50)["SLOPE"] <= 0.2){if(mean_inh_last > 60) coef_estim["IC50"] <- min(mat_tbl$logconc,na.rm=T);
    nls_result_ic50 <- nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),lower=list(SLOPE=0.1,MIN=min_lower,MAX=max_lower,IC50=min(mat_tbl$logconc)),upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
    }
    #Calculate the standard error scores
    sumIC50 = summary(nls_result_ic50); ic50std_Error <- sumIC50$coefficients["IC50","Std. Error"];
    ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2)/(length(sumIC50$residuals)-1)),1);
    max_signal <- max(mat_tbl$dose,na.rm=T); min_signal <- min(mat_tbl$dose,na.rm=T)
    
    #prepare final data and convert IC50 back from log scale (inverse)
    coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE","MAX","MIN")]; coef_ic50["IC50"] <- 10^coef_ic50["IC50"]; coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"]<0,max_signal,coef_ic50["IC50"])
    coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<0,max_signal,coef_ic50["IC50"]);coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<10,max_signal,coef_ic50["IC50"])
    coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"]<0,0,coef_ic50["MAX"]);coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition,na.rm=T),min(mat_tbl$inhibition,na.rm=T))>50),min_signal,coef_ic50["IC50"])
    x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=100); yic <- predict(nls_result_ic50, data.frame(logconc=x))
    perInh <- t(matrix(mat_tbl[,"inhibition"],dimnames=list(paste0(rep("D", length(mat_tbl[,"inhibition"])), 1:length(mat_tbl[,"inhibition"])))))
    coef_tec50 = coef_ic50; 
    coef_tec50["IC50"] <- ifelse(coef_tec50["MAX"] > 25, coef_tec50["IC50"], max(mat_tbl$dose,na.rm=T))
    if(readoutCTX){names(coef_tec50) <- c("TC50","SLOPE","MAX","MIN"); ytec <- yic; perViaTox <- perInh;} else{
      names(coef_tec50) <- c("EC50","SLOPE","MAX","MIN"); coef_tec50["SLOPE"] = -1 * coef_tec50["SLOPE"]; # min - 0, max - 77 in ec50 it is max - 100, min - 23
      tmp = coef_tec50["MAX"]; coef_tec50["MAX"] = 100 - coef_tec50["MIN"]; coef_tec50["MIN"] = 100 - tmp; ytec <- 100 - yic;
      perViaTox <- 100 - perInh;
    }
    ############################# 
    #############    DSS
    dss_score <- 100#round(as.numeric(dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=as.integer(DSS_typ))),1);
    coef_ic50 <- c(coef_ic50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,IC50_std_error=ic50std_Error,DSS=dss_score)
    coef_tec50 <- c(coef_tec50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,TEC50_std_error=ic50std_Error)
    IC50_df <- data.frame(DRUG_NAME="drug_name",ANALYSIS_NAME="IC50", t(as.matrix(coef_ic50)), perInh,GRAPH=NA, DSS = as.numeric(dss_score), sDSS = "", SE_of_estimate = as.numeric(ic50std_resid))
    #TEC50_df <- data.frame(DRUG_NAME="drug_name",ANALYSIS_NAME=TEC50,t(as.matrix(coef_tec50)), perViaTox, GRAPH=NA, DSS = as.numeric(dss_score), sDSS = "", SE_of_estimate = as.numeric(ic50std_resid))
    #numeric_cols <- sapply(IC50_df, is.numeric); IC50_df[,numeric_cols] <- round(IC50_df[,numeric_cols],1); numeric_cols <- sapply(TEC50_df, is.numeric); TEC50_df[,numeric_cols] <- round(TEC50_df[,numeric_cols],1)
    
    #  plot if needed
    # icpl <- ggplot2::ggplot(mat_tbl, aes(logconc, inhibition2)) + scale_x_continuous(breaks=mat_tbl$logconc,labels=mat_tbl$dose) +
    #   geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = x, y = yic), aes(x, yic), color="blue", size = 0.8) +
    #   geom_vline(xintercept = log10(coef_ic50["IC50"]), colour="grey", size = 0.8) + ggtitle(paste0(drug_name,"\n"))+
    #   theme_bw() + labs(y = "% inhibition", x = "conc(nM)")  +  ylim(-25, 125) +
    #   geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ic50["IC50"])*0.95, y2=115, text2="IC50"), color="grey", parse=T) +
    #   theme_minimal() +  theme(plot.title = element_text(hjust = 0.5))
    # 
    # print(icpl);
    
    return (list(coef_ic50=coef_ic50,nls_result_ic50=nls_result_ic50));
  })
}



#############################################################################################################
#### Plots

# plot dose-response matrix (args: matrix, title, agent names)

PlotRespMatr <- function(matr, name = "none", d1N = "", d2N = ""){
  
  data.plot <- melt(matr)
  colnames(data.plot) <- c("y","x","Inhibition")
  data.plot$Inhibition <- round(c(matr), 1)
  data.plot$x <- as.factor(data.plot$x)
  data.plot$y <- as.factor(data.plot$y)
  
  axis.x.text <- round(as.numeric(colnames(matr)), 1)
  axis.y.text <- round(as.numeric(rownames(matr)), 1)
  dose.response.p <- ggplot(data.plot, aes_string(x = "x", y = "y")) + geom_tile(aes_string(fill = "Inhibition"), color = "#FCFCFC", size=1) + 
    theme(title = element_text(face = "bold", size = 10)) + 
    geom_text(aes_string(fill = "Inhibition", label = "Inhibition"), size = 3.5) + 
    scale_fill_gradient2(low = "grey", high = "#4682b4", midpoint = 0, name = paste0("inhibition", " (%)"),na.value="white") + 
    scale_x_discrete(labels = axis.x.text) + scale_y_discrete(labels = axis.y.text) + 
    labs(x = d2N, y = d1N) + theme(plot.title = element_text(hjust = 0.5, size = 16)) + guides(fill=F)
  
  dose.response.p <- dose.response.p + theme(axis.text.x = element_text(color = "black", face = "bold", size = 11))
  dose.response.p <- dose.response.p + theme(axis.text.y = element_text(color = "black", face = "bold", size = 11))
  dose.response.p <- dose.response.p + theme(axis.title = element_text(size = 14))
  dose.response.p <- dose.response.p + ggtitle(paste0("\nDose-response matrix (",name,")\n"));
  list(pl = dose.response.p)
}

# plot interaction matrix
PlotIntMatr <- function(matr, d1N, d2N){ 
  
  # single drug deviations
  D1Len = matr[,1]; names(D1Len)[1] = 1e-6; D2Len = matr[1,]; names(D2Len)[1] = 1e-6;
  d1 = tryCatch({ 
    predict(CALC_IC50_EC50_DSS(D1Len, DSS_typ = 2, drug_name = "")$nls_result_ic50) 
  }, error = function(e){D1Len})
  d2 = tryCatch({ 
    predict(CALC_IC50_EC50_DSS(D2Len, DSS_typ = 2, drug_name = "")$nls_result_ic50)
  }, error = function(e){D2Len})
  
  matr[,1] = d1; matr[1,] = d2;
  
  # calculate Bliss approximation
  bliss.mat = matr
  for (k in 2:nrow(matr))
    for (j in 2:ncol(matr))
      bliss.mat[k, j] <- matr[k,1] + matr[1,j] - matr[k,1] *  matr[1,j] / 100
  
  
  syn = matr - bliss.mat
  
  synergy.score = list();
  synergy.score$drug.pairs = data.frame(drug.row = d1N,  drug.col = d2N, concUnit = "nM", blockIDs = 1)
  synergy.score$scores = append(synergy.score$scores, list(syn))
  
  ss  = calcsyn(syn, synergy.score$drug.pairs)
  PlotSynergyShiny(ss, graphnumber = 2, gridsize2 =1)
}

##################################################################################
##### Plot 2D and 3D synergy interaction maps
###################################################################################

PlotSynergyShiny <- compiler::cmpfun(function (data, type = "2D", graphnumber = 1, brushx = NULL, brushy = NULL, gridsize = 1, gridsize2 = 0, 
                                               savee2D = NULL, savee3D = NULL, newscore = NULL, name_3D = NULL, method_ = "Bliss", synScoresMtx = NULL, mostsynarea = 1) 
{
  print("plotinside")
  
  !is.list(data) && {stop("Input data is not a list format!")}
  if (gridsize == -1) {colmap = !0; gridsize = 1} else { colmap = !1 }
  
  summary.score <- data$summary.score; cMat <- data$c
  drug.row <- data$drug.row; drug.col <- data$drug.col
  x.conc <- data$x.conc; y.conc <- data$y.conc
  start.point <- data$start.point; end.point <- data$end.point
  
  if (method_ == "ZIP") {
    if (!is.null(newscore)) {plot.title =  plot.title2 = bquote(~delta ~ " - score: " ~ .(newscore))
    } else { plot.title <- bquote(~delta ~ " - score: " ~ .(summary.score)); plot.title2 <- paste0("delta score: ", summary.score)}
    title3D = paste0(drug.row, " & ", drug.col, " <br>\U03B4 - score: ", summary.score)
  }
  else {
    if (!is.null(newscore)) { plot.title = plot.title2 = paste0(method_, " synergy score: ", newscore)
    } else { plot.title <- paste0(method_, " synergy score: ", summary.score); 
    plot.title2 <- paste0(method_, " synergy score: ",  summary.score);}
    title3D = paste0(drug.row, " & ", drug.col, " <br> ", method_, " synergy score: ", summary.score)
  }
  print(paste0("size- ", gridsize))
  
  if (graphnumber == 2){
    plot2d = melt(cMat);
    myPalette <- colorRampPalette(c("green2", "white", "red1"))(100)
    names(plot2d) <- c("x","y","z")
    gplot2d <- 
      ggplot(plot2d) + aes(x, y, z = z, fill = z)  + geom_raster(interpolate = !0) + 
      geom_contour(color = "white", alpha = 0.5) + 
      scale_fill_gradientn(expression(delta ~ -score), colours = myPalette, limits = c(start.point, end.point), 
                           values = rescale(c(-3, -1, 0, 1, 3))) + 
      scale_x_continuous(drug.col, expand = c(0, 0), 
                         breaks = seq(min(plot2d$x), max(plot2d$x), by = (max(plot2d$x) - min(plot2d$x))/(length(x.conc) - 1)), 
                         labels = round(x.conc, 2)) + 
      scale_y_continuous(drug.row, expand = c(0, 0), 
                         breaks = seq(min(plot2d$y), max(plot2d$y), by = (max(plot2d$y) - min(plot2d$y))/(length(y.conc) - 1)), 
                         labels = round(y.conc, 2)) + 
      theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      theme(axis.text = element_text(size = 10)) + 
      theme(title = element_text(vjust = 12)) 
    
    
    byx = (max(plot2d$x) - min(plot2d$x))/(length(x.conc) - 1);
    byy = (max(plot2d$y) - min(plot2d$y))/(length(y.conc) - 1);
    
    if(gridsize2 != 0)
    {     
      gplot2d <- gplot2d + geom_vline(xintercept = seq(min(plot2d$x), max(plot2d$x), by = byx), linetype = "dotted") + 
        geom_hline(yintercept = seq(min(plot2d$y), max(plot2d$y), by = byy), linetype = "dotted") 
      
      gplot2d + ggtitle(plot.title) + coord_cartesian(xlim = c(brushx[1], brushx[2]), ylim = c(brushy[1], brushy[2]))  + 
        theme(plot.title = element_text(hjust = 0.5))
    }
    
  }
})


###################################################################################
##### Calculate synergy surfaces out of scores
###################################################################################

calcsyn <- compiler::cmpfun(function(scores, drug.pairs)
{
  scores.dose <- t(scores)  
  combscores = scores.dose[-1, -1]; combscores[nrow(combscores),ncol(combscores)] <- 'NA'
  summary.score <- round(mean(as.numeric(combscores), na.rm = !0), 3)
  x.conc <- as.numeric(rownames(scores.dose))
  y.conc <- as.numeric(colnames(scores.dose))
  conc.unit <- drug.pairs$concUnit
  unit.text <- paste0("(", conc.unit, ")")
  drug.row <- paste0(drug.pairs$drug.row, " ", unit.text)
  drug.col <- paste0(drug.pairs$drug.col, " ", unit.text)
  color.range <- round(max(abs(max(as.numeric(combscores), na.rm = !0)), abs(min(as.numeric(combscores), na.rm = !0))) + 5, -1)
  start.point <- -color.range; end.point <- color.range  
  pixels.num = 5 * (length(x.conc) - 1) + 2;
  
  # only for visualization  (max BzCl)
  scores.dose[nrow(scores.dose),ncol(scores.dose)] <- max(scores.dose[nrow(scores.dose)-1,ncol(scores.dose)],
                                                          scores.dose[nrow(scores.dose),ncol(scores.dose)-1],
                                                          scores.dose[nrow(scores.dose)-1,ncol(scores.dose)-1])    
  
  kriged =  tryCatch({
    tmp <- cbind(expand.grid(c(0:(length(x.conc) - 1)), c(0:(length(y.conc) - 1))), c(as.matrix(scores.dose)))        
    kriging(tmp[, 1],tmp[, 2], tmp[, 3], lags = ifelse(dim(scores.dose)[1] < 8, 2 ,3), 
            pixels = pixels.num, model = "spherical")
  },error = function(e){
    appro <- function(x, n) approx(x, n=n)$y
    tryCatch({
      m = apply(t(apply(scores.dose, 1, function(x) appro(x, ncol(scores.dose)*2))), 2, function(x) appro(x, nrow(scores.dose)*2))
      tmp2<- cbind(expand.grid(c(0:(nrow(m)-1)),c(0:(ncol(m)-1))), c(as.matrix(m)))
      kriging(tmp2[, 1], tmp2[, 2], tmp2[, 3], lags = ifelse(dim(m)[1] < 8, 2 ,3), pixels = pixels.num, model = "spherical")
    },error = function(e){ 
      m = apply(t(apply(scores.dose, 1, function(x) appro(x, ncol(scores.dose)*3))), 2, function(x) appro(x, nrow(scores.dose)*3))
      tmp2<- cbind(expand.grid(c(0:(nrow(m)-1)),c(0:(ncol(m)-1))), c(as.matrix(m)))
      kriging(tmp2[, 1], tmp2[, 2], tmp2[, 3], lags = ifelse(dim(m)[1] < 8, 2 ,3), pixels = pixels.num, model = "spherical")
    })
  })
  
  xseq <- round(kriged[["map"]]$x/kriged$pixel)
  yseq <- round(kriged[["map"]]$y/kriged$pixel)
  a <- min(xseq):max(xseq); b <- min(yseq):max(yseq)
  na <- length(a); nb <- length(b)
  res1 <- as.double(rep(0, na * nb))
  res2 <- as.integer(rep(0, na * nb))
  # no more Rcpp here
  #res = krig(xseq = xseq, kriged = kriged[["map"]]$pred, a = a, yseq = yseq, b = b,  
  #           z_len = length(kriged[["map"]]$pred), res1 = res1, na = na, nb = nb, res2 = res2)
  #res1 = res[[1]]; res1[res[[2]] == 0] <- NA
  
  z.len = length(kriged[["map"]]$pred); #
  for(idx1 in 1:na) { #
    for(idx2 in 1:nb) { #
      for(idx3 in 1:z.len) { #
        if(xseq[idx3] == a[idx1] && yseq[idx3] == b[idx2]) { #
          indx_ = idx2+(idx1-1)*nb; #
          res1[indx_] <- kriged[["map"]]$pred[idx3] #
          res2[indx_] <- 1 #
          break #
        } #
      } #
    } #
  } #
  
  res1[res2 == 0] <- NA #
  cMat <- matrix(res1, na, nb, byrow = !0)
  
  # most synergystic region
  max_ = r_ = c_ = -999;
  for(i in 1:(ncol(scores.dose)-2)){
    for(j in 1:(nrow(scores.dose)-2)){
      mean_ = mean(scores.dose[j:(j+2),i:(i+2)], na.rm = !0)
      if(mean_ > max_) {
        max_ = mean_; r_ = j; c_ = i;
      }
    }
  }
  
  return(list(c = cMat, conc.unit = conc.unit, drug.row = drug.row, 
              drug.col = drug.col, start.point = start.point, end.point = end.point, 
              summary.score = summary.score, x.conc = x.conc, y.conc = y.conc, pixels.num = pixels.num, r_ = r_, c_ = c_, max_ = round(max_,3)))
})



# combine multiple plots
multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}


# outliers_ = c();
# unsure_ =c()
# 
# #sub_ = sample(unique(MGdataFull$PairIndex), 50)
# for(ind_Grin in 1:60){ print(ind_Grin)
#   
#   ############## <<--- for mathew ... --->> ##############
#   
#   data_cell_pair <- dplyr::arrange(MGdataFull[MGdataFull$PairIndex==ind_Grin, ], Conc1, Conc2)
#   
#   #data_cell_pair = data_cell_pair[data_cell_pair$Conc1 %in% unique(data_cell_pair$Conc1)[1:4] & data_cell_pair$Conc2 %in% unique(data_cell_pair$Conc2)[1:4], ]
#   
#   data_cell_pair$Response[data_cell_pair$Conc1==max(data_cell_pair$Conc1)&data_cell_pair$Conc2==max(data_cell_pair$Conc2)] =
#    min(c(data_cell_pair$Response[data_cell_pair$Conc1==sort(unique(data_cell_pair$Conc1),decreasing = T)[2]&data_cell_pair$Conc2==max(data_cell_pair$Conc2)],
#          data_cell_pair$Response[data_cell_pair$Conc2==sort(unique(data_cell_pair$Conc2),decreasing = T)[2]&data_cell_pair$Conc1==max(data_cell_pair$Conc1)]))
#   
#   c1 = unique(data_cell_pair$Conc1); c2 = unique(data_cell_pair$Conc2)
#   openxlsx::write.xlsx(data_cell_pair,paste0("Mathews Griner/Real/",ind_Grin,".xlsx"))
#   
#   # hack concenttions (need for rowInd)
#    data_cell_pair$Conc1 = sapply(data_cell_pair$Conc1, function(i) which(i == unique(data_cell_pair$Conc1))-1)
#    data_cell_pair$Conc2 = sapply(data_cell_pair$Conc2, function(i) which(i == unique(data_cell_pair$Conc2))-1)
#   
#  # conc_ic50 = CALC_IC50_EC50_DSS(xpr_tbl = 100-reshape2::acast(data_cell_pair, Conc1~Conc2, value.var = "Response")[-1,][,1])$coef_ic50["IC50"] 
#   # sbol=randtoolbox::sobol(5,dim=2,init=!0,scrambling=0,seed=4711,normal=!1); MatrTr=reshape2::acast(data_cell_pair,Conc1~Conc2,value.var="Response")
#   # rowIndComb=ceiling((nrow(MatrTr)-1)*sbol[,1]);colIndComb=ceiling((ncol(MatrTr)-1)*sbol[,2]) #plot(rowIndComb, colIndComb)
#   # IndComb = sapply(1:5, function(i){
#   #   which(data_cell_pair$Conc1 == rowIndComb[i] & data_cell_pair$Conc2 == colIndComb[i])
#   # })
#   # rowInd = c(which(data_cell_pair$Conc1==0 | data_cell_pair$Conc2 == 0), IndComb)
#  # concRow = which.min(abs(as.numeric(conc_ic50) - c1));
#   
#   #rowInd = which(data_cell_pair$Conc1==0 | data_cell_pair$Conc2 == 0 | data_cell_pair$Conc1==(concRow-1))
#   
#   #rowInd = which(data_cell_pair$Conc1==0 | data_cell_pair$Conc2 == 0 | data_cell_pair$Conc2 == (5-data_cell_pair$Conc1))
#   rowInd = which(data_cell_pair$Conc1==0 | data_cell_pair$Conc2 == 0 | data_cell_pair$Conc1 == floor(max(data_cell_pair$Conc1) / 2 + 1))
#   #rowInd = which(data_cell_pair$Conc1==0 | data_cell_pair$Conc2 == 0 | data_cell_pair$Conc1 %in% c(4,8))
#   
#   #rowInd = which(data_cell_pair$Conc1==0 | data_cell_pair$Conc2 == 0 | data_cell_pair$Conc2 == data_cell_pair$Conc1)
#   rowInd2 <- c(1:nrow(data_cell_pair))[!(1:64 %in% rowInd)];
#   
#   
#   # return back
#   data_cell_pair$Conc1 = sapply(data_cell_pair$Conc1, function(i) c1[which(i == unique(data_cell_pair$Conc1))])
#   data_cell_pair$Conc2 = sapply(data_cell_pair$Conc2, function(i) c2[which(i == unique(data_cell_pair$Conc2))])
#   
#   data_cell_pair$R1 <- data_cell_pair$Response[as.numeric(sapply(seq.int(1, 64, by = 8), function(i) rep(i,8)))]; 
#   data_cell_pair$R2 <- rep(data_cell_pair$Response[1:8], 8);
#   # data_cell_pair[, c("Conc1","Conc2", "R1", "R2")] = as.data.frame(scale(data_cell_pair[, c("Conc1","Conc2", "R1", "R2")], center=!1, scale=colSums(data_cell_pair[, c("Conc1","Conc2", "R1", "R2")])))
#   
#   data_cell_pair_Training <- data_cell_pair[rowInd, c("Conc1","Conc2", "R1", "R2", "Response")]; data_cell_pair_TrainingFull = data_cell_pair_Training;
#   data_cell_pair_Test <- data_cell_pair[rowInd2, c("Conc1","Conc2", "R1", "R2")];   data_cell_pairCp = data_cell_pair;
#   data_cell_pair_Training$Response[data_cell_pair_Training$Response < 0] = 0.0001;
#   PlotRespMatr(100 -reshape2::acast(data_cell_pair_Training, Conc1~Conc2, value.var = "Response"))
#   
# 
#   #library(e1071); library(caret); library(mlr)
#   set.seed(42)
#   data_cell_pair_Test <- data_cell_pair[rowInd2, c("Conc1","Conc2", "R1", "R2")]; 
#   data_cell_pair_Training2 = data_cell_pair_Training[!is.na(data_cell_pair_Training$Response),]
#   influentPoint = NULL # for now
#   
#   #data_cell_pair_Test2 = data_cell_pair_Test; data_cell_pair_Test2$Response=NA; #data_cell_pair_TrainingTmp$Response[bag_] = NA
#   MatrTr = reshape2::acast(data_cell_pair_Training, Conc1~Conc2, value.var = "Response")
#   
#   # check [0,0] conc. 
#   if(MatrTr[1,1] < max(MatrTr[2,1], MatrTr[1,2])) MatrTr[1,1] = 100
#   
#   # calculate Bliss approximation
#   MatrTr = 100 - MatrTr; bliss.mat = MatrTr
#     for (k in 2:nrow(MatrTr))
#     for (j in 2:ncol(MatrTr))
#       bliss.mat[k, j] <- MatrTr[k,1] + MatrTr[1,j] - MatrTr[k,1] *  MatrTr[1,j] / 100
#   
#   # single drug deviations
#   D1Len = MatrTr[,1]; names(D1Len)[1] = 1e-6; D2Len = MatrTr[1,]; names(D2Len)[1] = 1e-6;
#   d1 = tryCatch({ 
#       predict(CALC_IC50_EC50_DSS(D1Len, DSS_typ = 2, drug_name = "")$nls_result_ic50) 
#   }, error = function(e){D1Len})
#   d2 = tryCatch({ 
#     predict(CALC_IC50_EC50_DSS(D2Len, DSS_typ = 2, drug_name = "")$nls_result_ic50)
#   }, error = function(e){D2Len})
#   
#   # deviations
#   devD1 = abs(d1 - (MatrTr[,1])); devD2 = abs(d2 - (MatrTr[1,])); dev_ = abs(bliss.mat - MatrTr); 
#   MatrOutl = matrix(!1, nrow = nrow(MatrTr), ncol = ncol(MatrTr))
#   
#   # check Bliss deviation
#   MatrOutl[-1,-1] = outlier_remove(abs(dev_[-1,-1] - median(dev_[-1,-1], na.rm = T)), iqr_ = 5) & (dev_[-1,-1] > 25)
#   MatrOutl[,1] = as.logical(colSums(MatrOutl, na.rm = T)) & devD1 > 10 | devD1 > 15
#   MatrOutl[1,] = as.logical(rowSums(MatrOutl, na.rm = T)) & devD2 > 10 | devD1 > 15
# 
#   # check if unsure
#   devPrevD1 = abs(sapply(2:length(MatrTr[1,]), function(x) MatrTr[1,x] - MatrTr[1,x-1]))
#   devPrevD2 = abs(sapply(2:length(MatrTr[,1]), function(x) MatrTr[x,1] - MatrTr[x-1,1]))
#   if(any(devPrevD1>75) && any(devPrevD2>75)) unsure_ = c(unsure_, ind_Grin)
#   
#   #
#   MatrOutl[is.na(MatrOutl)] = !1; MatrTr[MatrOutl] = NA
#   if(any(MatrOutl)){ influentPoint = T; warning("Possible outliers were removed"); outliers_ = c(outliers_, ind_Grin) }
#   
#   MatrTr = 100 - MatrTr
#   
#   
#   matr_Out = reshape2::melt( MatrTr ); colnames(matr_Out) <- c("Conc1", "Conc2", "Response"); 
#   matr_Out = dplyr::arrange(matr_Out, Conc1, Conc2)[rowInd,]
#   data_cell_pair_Training2$Response = matr_Out$Response
#   
#   
#   NNMF20 = do.call("cbind", mclapply(1:100, function(i){
#     
#     #bag_ = sample(1:nrow(data_cell_pair_Training), sample(2:round(nrow(data_cell_pair_Training)/3.5),1));
#     data_cell_pair_Test2 = data_cell_pair_Test; data_cell_pair_Test2$Response=NA; #data_cell_pair_TrainingTmp$Response[bag_] = NA
#     MatrTr = reshape2::acast(rbind(data_cell_pair_Training2, data_cell_pair_Test2), Conc1~Conc2, value.var="Response")
#     MatrTr = MatrTr + matrix(runif(1, -0.001, 0.001), nrow = nrow(MatrTr), ncol = ncol(MatrTr))
#     
# 
#     if(length(influentPoint) != 0){
#       if(is.na(MatrTr[1,1])){ MatrTr[1,1] = 100 }; if(is.na(MatrTr[nrow(MatrTr),1])) {MatrTr[nrow(MatrTr),1] = MatrTr[nrow(MatrTr)-1,1]}
#       if(is.na(MatrTr[1, ncol(MatrTr)])) {MatrTr[1, ncol(MatrTr)] = MatrTr[1, ncol(MatrTr)-1]};
#       MatrTr[,1] = zoo::na.approx(MatrTr[,1], rule=2); MatrTr[1,] = zoo::na.approx(MatrTr[1,], rule=2);
#     }
#     
# 
#     #MatrTr = 100 - BaselineCorrectionSD(100 - MatrTr)$corrected.mat
#     nsclc2.nmf <- NNLM::nnmf(as.matrix(MatrTr), sample(2:3,1), verbose = F, beta = sample(seq(.1,1,.01), 3), alpha = sample(seq(.1,1,.01), 3), #rep(.001,3), alpha = rep(.001,3),
#                              max.iter = 500L, loss = "mse", check.k = F);
#     nsclc2.hat.nmf <- with(nsclc2.nmf, W %*% H);
# 
# 
# 
#     matr_Out =reshape2::melt( nsclc2.hat.nmf ); colnames(matr_Out) <- c("Conc1", "Conc2", "Response"); 
#     
#     NNMF20 =  dplyr::arrange(matr_Out, Conc1, Conc2)$Response
#     
#   }, mc.cores = 4))
#  
#   if(sum(colSums(NNMF20 == 0)==0)>1) NNMF20 = NNMF20[,colSums(NNMF20 == 0) == 0]
#   NNMF20 = rowMeans(NNMF20)
#   
# 
#   
#   # average and write results
#   data_cell_pair_Test$Response = NNMF20[rowInd2]; data_cell_pair_Test$Response[data_cell_pair_Test$Response > 100] = 100;
#   dataNNMFout <- dplyr::arrange(gtools::smartbind(data_cell_pair_Test, data_cell_pair_TrainingFull), Conc1, Conc2);
#   data_cell_pair$predNNMF <- dataNNMFout$Response
#   data_cell_pair$S_mseNNMF = sqrt(sum((data_cell_pair_Test$Response - data_cell_pair$Response[rowInd2])^2, na.rm = T) /(length(data_cell_pair$Response[rowInd2])))
#   data_cell_pair$S_mseNNMF   
#   
# 
#    # XGBoostRun <<- list();
#   data_cell_pair_Training2 <- data_cell_pair_Training2[!is.na(data_cell_pair_Training2$Response),];  
#   data_cell_pair_Training2$Response = data_cell_pair_Training2$Response + sapply(1:nrow(data_cell_pair_Training2), function(i) runif(1, -0.0001, 0.0001))
#     
#     obj.fun = makeSingleObjectiveFunction(
#       name = "XGBoost",
#       fn = function(x) {
#         
#         logNtree = x[1]; lambda = x[2]; alpha = x[3]; maxdepth = x[4]; subsample = x[5]; colsample_bytree = x[6]; eta = x[7];  #StackedTrain = c();
#         
#         # repeated CV
#         #for(repCv in 1:3){
#         MAD_ <- mclapply(1:4, function(repCv){
#           
#           MAD_i = 0
#           flds <- caret::createFolds(data_cell_pair_Training2$Response, k = 3, list = T, returnTrain = F);
#           
#           for(k in 1:length(flds)){
#             testData <- data_cell_pair_Training2[flds[[k]], ]; trainData <- data_cell_pair_Training2[-flds[[k]], ]
#             
#             fit = xgboost(data=as.matrix(trainData[,c("R1","R2","Conc1","Conc2")]),label = trainData$Response, verbose = F, 
#                           nrounds=round(2**logNtree), nthread = 1, save_name = paste0("xgboost",repCv,"_",k,".model"),
#                           params=list(objective = "reg:linear", max.depth=maxdepth, eta=eta, lambda = lambda, alpha = alpha, 
#                                       subsample=subsample, colsample_bytree = colsample_bytree))
#             
#             
#             ypred = predict(fit, as.matrix(testData[,c("R1","R2","Conc1","Conc2")])); 
#             MAD_i <- MAD_i + mean(abs(ypred - testData$Response), na.rm = T);
#             
#           #  testData$Response = ypred; StackedTrain = rbind(StackedTrain, testData)
#           }
#           MAD_i
#         }, mc.cores = 4)
#         #XGBoostRun <<- append(XGBoostRun, list(StackedTrain))
#         Reduce(sum, MAD_)
#         
#       },
#       par.set = makeParamSet(
#         makeNumericVectorParam("logNtree", len = 1, lower = 4, upper = 9),
#         makeNumericVectorParam("lambda", len = 1, lower = 0, upper = 3),
#         makeNumericVectorParam("alpha", len = 1, lower = 0, upper = 3),
#         makeIntegerVectorParam("maxdepth", len = 1, lower = 1, upper = 6),
#         makeNumericVectorParam("subsample", len = 1, lower = .4, upper = 1),
#         makeNumericVectorParam("colsample_bytree", len = 1, lower = .4, upper = 1),
#         makeNumericVectorParam("eta", len = 1, lower = .001, upper = .1)
#       ),
#       
#       minimize = !0
#     )
#     
#     des = generateDesign(n = 10, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
#     
#     des$y = apply(des, 1, obj.fun) #as.numeric(mclapply(1:nrow(des), function(i) obj.fun(des[i,]), mc.cores = 6))
#     surr.km = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2", control = list(trace = F))
#     
#     # library(parallelMap)
#     # parallelStartMulticore(cpus = 4, show.info = T)
#     modelsXGBoost = mbo(obj.fun, design = des, learner = surr.km, show.info = !0,
#                         control = setMBOControlInfill(setMBOControlTermination(makeMBOControl(), iters = 30),
#                                                       crit = makeMBOInfillCritEI()))$opt.path$env[["path"]]
#     
# 
#     orderMAD = order(modelsXGBoost$y); models = modelsXGBoost[orderMAD, ]; #run = XGBoostRun[orderMAD]
#     
#     
#     # Fit with 20 models with best parameters
#     pred_ <- do.call("cbind", mclapply(1:4, function(i){
#       fit <- xgboost(as.matrix(data_cell_pair_Training2[,c("R1","R2","Conc1","Conc2")]), label = data_cell_pair_Training2$Response,
#                      verbose = F, nrounds=round(2**models[i,"logNtree"]) , nthread = 1,
#                      params=list(objective = "reg:linear", max.depth=models[i,"maxdepth"], eta=models[i,"eta"], lambda = models[i,"lambda"], 
#                                  alpha = models[i,"alpha"], subsample=models[i,"subsample"], colsample_bytree = models[i,"colsample_bytree"]))
#       predict(fit, as.matrix(data_cell_pair_Test[,c("R1","R2","Conc1","Conc2")])); 
#     }, mc.cores = 4))
#     
#     
#     data_cell_pair_Test$Response = rowMeans(pred_[,1:4]); data_cell_pair_Test$Response[data_cell_pair_Test$Response > 100] = 100;
#     dataXGBoostout <- dplyr::arrange(gtools::smartbind(data_cell_pair_Test, data_cell_pair_TrainingFull), Conc1, Conc2);
#     data_cell_pair$predXGBoost <- dataXGBoostout$Response
#     data_cell_pair$S_mseXGBoost = sqrt(sum((data_cell_pair_Test$Response - data_cell_pair$Response[rowInd2])^2, na.rm = T) /(length(data_cell_pair$Response[rowInd2])))
#     data_cell_pair$S_mseXGBoost
#     
#     #if(length(influentPoint) == 0) final_ = (data_cell_pair$predXGBoost + data_cell_pair$predNNMF)/2 else final_ = data_cell_pair$predNNMF
#     #final_ = data_cell_pair$predXGBoost 
#     final_ =  (data_cell_pair$predXGBoost + data_cell_pair$predNNMF)/2 
#     data_cell_pair$predfinal = final_
#     data_cell_pair$S_msefinal = sqrt(sum((final_[rowInd2] - data_cell_pair$Response[rowInd2])^2, na.rm = T) /(length(data_cell_pair$Response[rowInd2])))
#     data_cell_pair$S_msefinal
#     
#     
#     
#     
#     data_cell_pair_Training2 <- data_cell_pair_Training; data_cell_pair_Training2$Response = data_cell_pair_Training2$Response / 100;
#     data_cell_pair_Training2[is.na(data_cell_pair_Training2$Response),"Response"] = "Nan"
#     
#     simpleFile <- data_cell_pair_Test[, c("Conc1", "Conc2")]; simpleFile$Response = "Nan";
#     simpleFile <- dplyr::arrange(gtools::smartbind(simpleFile, data_cell_pair_Training2[, c("Conc1", "Conc2", "Response")]), Conc1, Conc2);
#     colnames(simpleFile)[grepl("Response", colnames(simpleFile))] <- "Var4";
#     
#     write.table(simpleFile, 'aaa.csv',sep = ",", row.names = !1)
#     charsAAA = gsub('',"",gsub('"', "", readChar('aaa.csv',  file.info('aaa.csv')$size)))
#     cat(charsAAA, file='aaa.csv')
#     
#     system("/home/local/MATLAB/R2017a/bin/matlab -nodesktop -nojvm -nosplash -noFigureWindows -nodisplay -wait -r \"run('./Example.m'); exit\"", wait = !0);
#     
#     datatmp <- dplyr::arrange(data_cell_pair[, c("Conc1","Conc2", "R1", "R2", "Response")], Conc1, Conc2)
#     datatmp$Response = (round(dplyr::arrange(read.csv('aaa-predictions.csv', header = T), Conc1, Conc2)$Var4, 6)) * 100
#     data_cell_pair$predDOSE <- datatmp$Response;
#     
#     # calculate error
#     data_cell_pair$S_mseDOSE = sqrt(sum((datatmp$Response[rowInd2] - data_cell_pair$Response[rowInd2])^2, na.rm = T) /(length(data_cell_pair$Response[rowInd2])))
#     data_cell_pair$S_mseDOSE 
#     
#     
#     
#     # calculate synergy
#     dtSyn = list(dose.response.mats = list(), drug.pairs = data.frame(drug.row="d1",drug.col="d2",concUnit="nM",blockIDs=1,stringsAsFactors=!1));
#     dtSyn$dose.response.mats = c(dtSyn$dose.response.mats, list(data.matrix(100 - reshape2::acast(data_cell_pair, Conc1~Conc2, value.var="Response"))))
#     matr_Out=reshape2::melt(CalculateSynergy(dtSyn,method="Bliss",correction=!1,Emin=NA,Emax=NA)$scores[[1]]); colnames(matr_Out)<-c("Conc1","Conc2","SynScore");
#     data_cell_pair$Bliss=dplyr::arrange(matr_Out,Conc1,Conc2)$SynScore;data_cell_pair$BlissReal=mean(data_cell_pair$Bliss[data_cell_pair$Bliss!=0])
#     matr_Out=reshape2::melt(CalculateSynergy(dtSyn,method="Bliss",correction=!0,Emin=NA,Emax=NA)$scores[[1]]); colnames(matr_Out)<-c("Conc1","Conc2","SynScore");
#     data_cell_pair$BlissCor=dplyr::arrange(matr_Out,Conc1,Conc2)$SynScore;data_cell_pair$BlissCorReal=mean(data_cell_pair$BlissCor[data_cell_pair$BlissCor!=0])
#     # calculate synergy NNMF
#     dtSyn = list(dose.response.mats = list(), drug.pairs = data.frame(drug.row="d1",drug.col="d2",concUnit="nM",blockIDs=1,stringsAsFactors=!1));
#     dtSyn$dose.response.mats = c(dtSyn$dose.response.mats, list(data.matrix(100 - reshape2::acast(data_cell_pair, Conc1~Conc2, value.var="predNNMF"))))
#     matr_Out=reshape2::melt(CalculateSynergy(dtSyn,method="Bliss",correction=!1,Emin=NA,Emax=NA)$scores[[1]]); colnames(matr_Out)<-c("Conc1","Conc2","SynScore");
#     data_cell_pair$BlissNNMF=dplyr::arrange(matr_Out,Conc1,Conc2)$SynScore;data_cell_pair$BlissRealNNMF=mean(data_cell_pair$BlissNNMF[data_cell_pair$BlissNNMF!=0])
#     matr_Out=reshape2::melt(CalculateSynergy(dtSyn,method="Bliss",correction=!0,Emin=NA,Emax=NA)$scores[[1]]); colnames(matr_Out)<-c("Conc1","Conc2","SynScore");
#     data_cell_pair$BlissCorNNMF=dplyr::arrange(matr_Out,Conc1,Conc2)$SynScore;data_cell_pair$BlissCorRealNNMF=mean(data_cell_pair$BlissCorNNMF[data_cell_pair$BlissCorNNMF!=0])
#     # calculate synergy XGBoost
#     dtSyn = list(dose.response.mats = list(), drug.pairs = data.frame(drug.row="d1",drug.col="d2",concUnit="nM",blockIDs=1,stringsAsFactors=!1));
#     dtSyn$dose.response.mats = c(dtSyn$dose.response.mats, list(data.matrix(100 - reshape2::acast(data_cell_pair, Conc1~Conc2, value.var="predXGBoost"))))
#     matr_Out=reshape2::melt(CalculateSynergy(dtSyn,method="Bliss",correction=!1,Emin=NA,Emax=NA)$scores[[1]]); colnames(matr_Out)<-c("Conc1","Conc2","SynScore");
#     data_cell_pair$BlissXGBoost=dplyr::arrange(matr_Out,Conc1,Conc2)$SynScore;data_cell_pair$BlissRealXGBoost=mean(data_cell_pair$BlissXGBoost[data_cell_pair$BlissXGBoost!=0])
#     matr_Out=reshape2::melt(CalculateSynergy(dtSyn,method="Bliss",correction=!0,Emin=NA,Emax=NA)$scores[[1]]); colnames(matr_Out)<-c("Conc1","Conc2","SynScore");
#     data_cell_pair$BlissCorXGBoost=dplyr::arrange(matr_Out,Conc1,Conc2)$SynScore;data_cell_pair$BlissCorRealXGBoost=mean(data_cell_pair$BlissCorXGBoost[data_cell_pair$BlissCorXGBoost!=0])
#     # calculate synergy Averaged
#     dtSyn = list(dose.response.mats = list(), drug.pairs = data.frame(drug.row="d1",drug.col="d2",concUnit="nM",blockIDs=1,stringsAsFactors=!1));
#     dtSyn$dose.response.mats = c(dtSyn$dose.response.mats, list(data.matrix(100 - reshape2::acast(data_cell_pair, Conc1~Conc2, value.var="predfinal"))))
#     matr_Out=reshape2::melt(CalculateSynergy(dtSyn,method="Bliss",correction=!1,Emin=NA,Emax=NA)$scores[[1]]); colnames(matr_Out)<-c("Conc1","Conc2","SynScore");
#     data_cell_pair$Blissfinal=dplyr::arrange(matr_Out,Conc1,Conc2)$SynScore;data_cell_pair$BlissRealfinal=mean(data_cell_pair$Blissfinal[data_cell_pair$Blissfinal!=0])
#     matr_Out=reshape2::melt(CalculateSynergy(dtSyn,method="Bliss",correction=!0,Emin=NA,Emax=NA)$scores[[1]]); colnames(matr_Out)<-c("Conc1","Conc2","SynScore");
#     data_cell_pair$BlissCorfinal=dplyr::arrange(matr_Out,Conc1,Conc2)$SynScore;data_cell_pair$BlissCorRealfinal=mean(data_cell_pair$BlissCorfinal[data_cell_pair$BlissCorfinal!=0])
#     
#     data_cell_pair$average = Biobase::rowMedians(cbind(data_cell_pair$predNNMF, data_cell_pair$predDOSE, data_cell_pair$predXGBoost))
#     data_cell_pair$average = fields::Tps(data_cell_pair[,c("Conc1",  "Conc2")], data_cell_pair$average, lambda= 1e-8)$fitted.values  # smooth a bit
#     
#     # calculate error
#     data_cell_pair$S_mseaverage = sqrt(sum((data_cell_pair$average[rowInd2] - data_cell_pair$Response[rowInd2])^2, na.rm = T) /(length(data_cell_pair$Response[rowInd2])))
#     data_cell_pair$S_mseaverage 
# 
#     
#     pdf(paste0(ind_Grin, ".pdf"), width = 15, height = 8)
# 
#     multiplot(PlotRespMatr(100 -reshape2::acast(data_cell_pair, Conc1~Conc2, value.var = "Response"), "Experiment", data_cell_pair$Drug1[[1]], data_cell_pair$Drug2[[1]])$pl,
#               PlotRespMatr(100 -reshape2::acast(data_cell_pair, Conc1~Conc2, value.var = "average"), "DECREASE", data_cell_pair$Drug1[[1]], data_cell_pair$Drug2[[1]])$pl,
#               PlotRespMatr(100 -reshape2::acast(data_cell_pair, Conc1~Conc2, value.var = "predDOSE"), "Dose model", data_cell_pair$Drug1[[1]],data_cell_pair$Drug2[[1]])$pl,
#               PlotIntMatr(100 -reshape2::acast(data_cell_pair, Conc1~Conc2, value.var = "Response"), data_cell_pair$Drug1[[1]], data_cell_pair$Drug2[[1]]),
#               PlotIntMatr(100 -reshape2::acast(data_cell_pair, Conc1~Conc2, value.var = "average"), data_cell_pair$Drug1[[1]], data_cell_pair$Drug2[[1]]),
#               PlotIntMatr(100 -reshape2::acast(data_cell_pair, Conc1~Conc2, value.var = "predDOSE"), data_cell_pair$Drug1[[1]], data_cell_pair$Drug2[[1]]),cols = 3)
#     
#     
#     dev.off()
#     data_cell_pair466Outl = rbind(data_cell_pair466Outl, data_cell_pair)
# }
# 
# saveRDS(data_cell_pair466Outl, "data_cell_pairDiagDen.RDS")
#     
# 
# 
# 
# 
# 
