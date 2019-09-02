# DSRT_set = "FO4B"; path = "/fs/projects/breeze/db/reports/6701_wrete2_aianevsk"; report.author = "aianevsk"; headerpath = "~/Desktop/projects/DSRT3 (2)/headerQC.txt"; tailerpath = "~/Desktop/projects/DSRT3 (2)/tailerQC.txt"; 
print("HI9")
lapply(c("plotly", "scales", "parallel", "foreach", "gridExtra", "grid", "graphics", "gplots", "ggplot2", "raster", "xtable"), library, character.only = !0)

# #dir.create("./Results", showWarnings=F, recursive=T)
# # headerpath = "/projects/breeze/code/DSRT3/headerQC.txt"; tailerpath = "/projects/breeze/code/DSRT3/tailerQC.txt"
# # 
# # create HTML report, instead of nozzle
# cm <- '<html><meta http-equiv="refresh" content="0; url=./Results/HTMLreport/qc.html" /></html>';
# writeChar(cm, "report.html", nchar(cm, type = "chars"))
# file.copy(from="/projects/breeze/code/DSRT3/HTMLreport/", to="./Results",  recursive = T)
# num_cores = parallel::detectCores()


################################################################################################ 
################# Matrix normalization (3-4 sigma's squeezing)

heatNorm <- compiler::cmpfun(function(datamat){
  mean_ = mean(datamat, na.rm = T); sd_ = sd(datamat, na.rm = T)
  m3sdp = mean_+3*sd_; m3sdn = mean_-3*sd_;
  ind_more = datamat > (m3sdp); ind_less = datamat < (m3sdn)
  ind_more[which(is.na(ind_more))]=FALSE
  ind_less[which(is.na(ind_less))]=FALSE
  datamat[ind_more] = m3sdp+(datamat[ind_more]-m3sdp)/((max(datamat,na.rm=T)-m3sdp)/sd_)
  datamat[ind_less] = m3sdn-(m3sdn-datamat[ind_less])/((m3sdn-min(datamat,na.rm=T))/sd_)
  
  list(datamat = datamat, outliers = ind_more | ind_less)
})

################################################################################################### 
################# outlier removal function.

outlier_remove <- compiler::cmpfun(function(x){
  qq <- unname(quantile(x, probs=c(.25, .75), na.rm = T))
  outlier_detector <- 1.5 * IQR(x, na.rm = T)
  x[x < (qq[1] - outlier_detector) | x > (qq[2] + outlier_detector)] <- NA
  x
})  


################################################################################################### 
################# read any plate of any format. in theory :)

read_plates <- function(file_name = NULL, sheet = 1, ncol_ = 24, nrow_ = 16){
  if (!file.exists(file_name)) stop("File dosn't exist")
  
  # read into array
  ext <- toupper(tools::file_ext(file_name)); gsub_ = " |:|-|;|,|\r|\n|\t";
  if(ext == "XLSX"){
    data_ <- data.frame(lapply(openxlsx::read.xlsx(file_name, colNames = F, rowNames = F, sheet = sheet), as.character), stringsAsFactors=F)
    data_ <- data_[!apply(is.na(data_) | data_ == "", 1, all),]
    data_ <- data_[,!apply(is.na(data_) | data_ == "", 2, all)]
    data_ <- gsub(gsub_,"", as.vector(t((na.omit(as.vector(t(data_)))))))
  } else if (ext %in% c("TXT", "CSV")){
    data_ = readChar(file_name, file.info( file_name )$size , useBytes=T)
    data_ = as.data.frame(strsplit(gsub(gsub_," ", data_), " "), stringsAsFactors = F)
    data_ <- data_[!apply(is.na(data_) | data_ == "", 1, all),]
  }
  
  # find first A letter and "1". STOP if either not found 
  ind_1 = base::match("1", data_); ind_A = base::match("A", data_); if(is.na(ind_1) | is.na(ind_A)) stop("No matrix detected")
  
  ## find date and Bar code. otherwise return 0
  #read_date = ifelse(length(date_Ind <- grep("Date", data_)), data_[date_Ind+1], 0)
  #BAR_code =  ifelse(!is.na(BAR_index <- match("ID1", data_)), data_[BAR_index+1], 0)
  
  # cut matrix from whole data_ and remove everything except raw data
  data_ = data_[ind_A:(base::match(LETTERS[nrow_], data_) + ncol_)]
  data_ = stringr::str_replace_all(data_, "[^[:digit:]]", "")
  data_ = as.numeric(data_[!(data_ == "")])
  
  # return matrix
  matrix(data_, nrow = nrow_, ncol = ncol_, byrow = T)
}

################################################################################################ 
################# Generate RDA 
################################################################################################ 


createeset <- function(file_names,annoframe,infoframe)
{
  infoframe$read_date=""
  
  #extract all the data from the plates
  screen_table = do.call(rbind, lapply(1:length(file_names), function(i) {
    
    
    tryCatch({
    
    file_infoname = gsub("./Data/","", as.character(file_names))[i]; print(file_infoname);  
    file_name = paste0("./www/Data/",file_infoname)
    file_info <- infoframe[match(file_infoname,as.character(infoframe$File_name)),]
    
    if(toupper(file_info$Reader) == "PARADIGM"){
      rawdatamat = read_plates(file_name, sheet = 2)
      read_date=toupper(format(as.Date(Sys.Date(),"%d-%m-%Y"),format="%d-%b-%y"))
      barcode=paste0(read_date, gsub("[[:punct:]]|\\s","",Sys.time()))	
    }
    else if (toupper(file_info$Reader) == "ENSIGHT")
    {
      #if(toupper(tools::file_ext(file_name)) == "CSV")
      #{
      rawdatamat = read_plates(file_name, sheet = 1)
      read_date=toupper(format(as.Date(Sys.Date(),"%d-%m-%Y"),format="%d-%b-%y"))
      barcode=paste0(read_date, gsub("[[:punct:]]|\\s","",Sys.time()))	
      #system(paste0('iconv -t ascii//TRANSLIT ', file_name, " > ",file_name,"NEW.csv"), wait = T)
      #system(paste0("mv ", file_name,"NEW.csv ",file_name), wait = T)
      
      #tbl_dataframe = read.csv(file_name,header=!1, blank.lines.skip=!0, fill=!0,skip=0,nrows=50,
      #                         sep = ",",colClasses=NA,stringsAsFactors=!1)
      #read_date=toupper(format(as.Date(substr(tbl_dataframe[grep("Screen Started:",tbl_dataframe$V1),3],1,10), origin = "1899-12-30"),format="%d-%b-%y"))
      #barcode=as.character(read_date); rowstart=match("A",tbl_dataframe$V1)[1]
      #rawdatamat=data.matrix(tbl_dataframe[rowstart:(rowstart + 15),2:25])
      #} 
    }
    else if (toupper(file_info$Reader) == "CYTATION5")
    {
      if(toupper(tools::file_ext(file_name)) == "XLSX")
      {
        tbl_dataframe = as.data.frame(openxlsx::read.xlsx(file_name, rows = 1:65, cols = 1:30, colNames = !1), stringsAsFactors = !1)
        rowstart = match("A",tbl_dataframe[,2])
        rawdatamat = data.matrix(tbl_dataframe[rowstart:(rowstart + 15),3:26])
        datec=grep("Date",tbl_dataframe); dater=grep("Date",tbl_dataframe[,datec]) #Date: 12.12.2014
        read_date=toupper(format(as.Date(as.numeric(as.character(tbl_dataframe[dater,datec+1])), origin = "1899-12-30"),format="%d-%b-%y"))
        barcodec=grep("Barcode",tbl_dataframe); barcode=unlist(tbl_dataframe[grep("Barcode",tbl_dataframe[,barcodec]),barcodec + 1])
        if(gtools::invalid(barcode)) barcode=read_date
      } else {	
        tbl_dataframe= read.csv(file_name,header=F, strip.white = T,blank.lines.skip=F,sep = "\t",colClasses=NA,stringsAsFactors=F)
        rawdatamat= data.matrix(read.csv(file_name,header=F, strip.white = T,blank.lines.skip=F,
                                         skip=grep("Results",tbl_dataframe[,1]) + 1, sep= "\t", colClasses=NA,stringsAsFactors=F)[1:16,2:25])
        datec=grep("Date",tbl_dataframe); dater=grep("Date",tbl_dataframe[,datec]) 
        read_date=toupper(format(as.Date(gsub("[.]","-",unlist(tbl_dataframe[dater,datec+1])),"%d-%m-%Y"),format="%d-%b-%y"))
        barcode=read_date
      }
    }
    else if (toupper(file_info$Reader) %in% "PHERASTAR")
    {
      if(toupper(tools::file_ext(file_name)) == "CSV")
      {
        rawdatamat = read_plates(file_name, sheet = 1)
        read_date=toupper(format(as.Date(Sys.Date(),"%d-%m-%Y"),format="%d-%b-%y"))
        barcode=paste0(read_date, gsub("[[:punct:]]|\\s","",Sys.time()))	
        #   tbl_dataframe= read.csv(file_name,header=F, blank.lines.skip=T,skip=0,nrows=16,sep = ",",colClasses=NA,stringsAsFactors=F)					
        #   
        #   if(length(grep("Time",tbl_dataframe$V3))!=0) #converted
        #   {	
        #     rawdatamat = read.csv(file_name,header=F, blank.lines.skip=F,nrows=16,
        #                           skip=grep("Time \\[",tbl_dataframe$V1) + 2,sep = ",",colClasses=NA,stringsAsFactors=F)
        #     rawdatamat[1] = NULL					
        #     datec=grep("Date:",tbl_dataframe) #Date: 12.12.2014
        #     read_date = toupper(format(as.Date(gsub("[.]","-",gsub(" ","", tbl_dataframe[grep("Date:",tbl_dataframe[,datec]),datec + 1], fixed = T)),"%d-%m-%Y"),format="%d-%b-%y"))
        #     barcodec=grep("ID1:",tbl_dataframe);barcode = gsub(" ","", tbl_dataframe[grep("ID1:",tbl_dataframe[,barcodec]),barcodec + 1], fixed = T) 
        #   }
        #   else if(length(grep("Time",tbl_dataframe$V1)) != 0) {  # with USA time
        #     rawdatamat = read_plates(file_name);
        #     datec=grep("Date:",tbl_dataframe); 
        #     read_date = toupper(format(as.Date(gsub("Date:|\\s","",strsplit(tbl_dataframe[grep("Date:",tbl_dataframe[,datec]),datec], "Time")[[1]][[1]]), format = "%m/%d/%Y"),format="%d-%b-%y"))
        #     barcodec=grep("ID1:",tbl_dataframe);barcode = gsub("ID1: ","", tbl_dataframe[grep("ID1:",tbl_dataframe[,barcodec]),barcodec], fixed = T) 
        #     
        #     
        #   }	
        #   else { 	
        #     rawdatamat=read.csv(file_name,header=F, blank.lines.skip=F,nrows=16,
        #                         skip=grep("Time \\[",tbl_dataframe$V1) + 2,sep = ",", colClasses=NA, stringsAsFactors=F)
        #     rawdatamat$V1=unlist(gsub("(\\w+\\:\\s+)","",rawdatamat$V1,perl=T))
        #     datec=grep("Date: ",tbl_dataframe) #Date: 12.12.2014
        #     read_date=toupper(format(as.Date(gsub("[.]","-",unlist(strsplit(tbl_dataframe[grep("Date: ",tbl_dataframe[,datec]),datec]," "))[2]),"%d-%m-%Y"),format="%d-%b-%y"))
        #     barcodec=grep("ID1: ",tbl_dataframe); barcode=unlist(strsplit(tbl_dataframe[grep("ID1: ",tbl_dataframe[,barcodec]),barcodec], ": "))[2]
        #   }
        #   if(gtools::invalid(barcode)) barcode=read_date
        #   rawdatamat=data.matrix(rawdatamat)
      }
      else if(toupper(tools::file_ext(file_name)) %in% c("XLS","XLSX"))
      {
        # rawdatamat = read_plates(file_name, sheet = 1)
        # read_date=toupper(format(as.Date(Sys.Date(),"%d-%m-%Y"),format="%d-%b-%y"))
        # barcode=paste0(read_date, gsub("[[:punct:]]|\\s","",Sys.time()))	
        if(toupper(tools::file_ext(file_name)) == "XLS") stop(".XLS is not supported in DSRTv2. Please use .XLSX");
        tbl_dataframe = as.data.frame(openxlsx::read.xlsx(file_name, rows = 1:65, cols = 1:25, colNames = !1), stringsAsFactors = !1)

        if(length(grep("Raw Data",tbl_dataframe[,2]))!=0) # Original
        {
          rowstart=grep("Raw Data",tbl_dataframe[,2]) + 2
          rawdatamat = data.matrix(tbl_dataframe[rowstart:(rowstart + 15),2:25])
          datec=grep("Date: ",tbl_dataframe)
          read_date = toupper(format(as.Date(gsub("[.]","-",unlist(strsplit(tbl_dataframe[grep("Date: ",tbl_dataframe[,datec]),datec], ": "))[2]),"%d-%m-%Y"),format="%d-%b-%y"))
          barcodec=grep("ID1: ",tbl_dataframe); barcode=unlist(strsplit(tbl_dataframe[grep("ID1: ",tbl_dataframe[,barcodec]),barcodec], ": "))[2]
        }
        else
        {
          rowstart=grep("Time \\[",tbl_dataframe[,1]) + 2
          rawdatamat = data.matrix(tbl_dataframe[rowstart:(rowstart + 15),2:25])
          datec=grep("Date:",tbl_dataframe) #Date: 12.12.2014
          read_date= toupper(format(as.Date(gsub("[.]","-", as.Date(as.numeric(as.character(tbl_dataframe[grep("Date:",tbl_dataframe[,datec]),datec + 1])), origin = "1899-12-30")),"%Y-%m-%d"),format="%d-%b-%y"))
          barcodec=grep("ID1:",tbl_dataframe); barcode= gsub(" ","", tbl_dataframe[grep("ID1:",tbl_dataframe[,barcodec]),barcodec + 1], fixed = T)
        }
        if(gtools::invalid(barcode)) barcode=read_date
        rawdatamat=data.matrix(rawdatamat)
      }
    }
    ########################################
    ####After getting raw data matrix.
    ####Code is same for all file types####
    ########################################
    
    dimnames(rawdatamat) <- list(LETTERS[1:16], 1:24); data_tbl <- reshape2::melt(as.matrix(rawdatamat))
    colnames(data_tbl) <- c("Row","Column","rawIntensity")
    
    cols_ <- list("Barcode", "read_date", "screen_id", "readout", "DWell"); 
    rows_ <- list(barcode, read_date, file_info$screen_id, file_info$Readout, paste0(data_tbl$Row, data_tbl$Column))
    invisible(lapply(1:length(cols_), function(i) data_tbl[, cols_[[i]]] <<- rows_[[i]]))
    
    data_tbl <- merge(data_tbl,annoframe[annoframe$Plate==file_info$Plate,],by="DWell", all.x = T)
    
    
    
    # positive and negative controls without outliers
    pos_ctrl <- outlier_remove(data_tbl$rawIntensity[data_tbl$Content %in% "pos"])
    neg_ctrl <- outlier_remove(data_tbl$rawIntensity[data_tbl$Content %in% "neg"])
    
    #Calculate percent inhibition and activation
    avg_low <- mean(pos_ctrl,na.rm=T); avg_high <- mean(neg_ctrl,na.rm=T)
    data_tbl$inhibition_percent <- ((avg_high-data_tbl$rawIntensity)/(avg_high-avg_low))*100
    
    data_tbl
    }, error = function(e) {
      write(paste0("ERROR!!! Something wrong with a plate: ",file_names[i],"or plate is missing"), 'analysis-output.txt'); write(7,"progress_bar.txt");
      stop();
    })
    
    }))
  
  
  
  # browser()
  # if (is.null(screen_table$inhibition_percent))  screen_table$inhibition_percent=data_tbl$inhibition_percent
  # 
  screen_table$Content[screen_table$Column==24 & screen_table$ProductName=="BzCl"] <- "pos"
  return(screen_table)
}


################################################################################################ 
################# QC statistics 

stderr <- compiler::cmpfun(function(x){sqrt(var(x,na.rm=T)/length(na.omit(x)))})
lowsd <- compiler::cmpfun(function(x){return(mean(x)-stderr(x))})
highsd <- compiler::cmpfun(function(x){return(mean(x)+stderr(x))})
pop.sd <- compiler::cmpfun(function(x)(sqrt(var(x)*(length(x)-1)/length(x)))); 
ssmd <- compiler::cmpfun(function(x,y)round((mean(x)-mean(y))/sqrt(var(x)+var(y))));
zfactor <- compiler::cmpfun(function(x,y)round((1-(3 * (pop.sd(x)+pop.sd(y)))/(abs(mean(x)-mean(y)))),2));
robustzfactor <- compiler::cmpfun(function(x,y)round((1-(3 * (mad(x)+mad(y)))/(abs(median(x)-median(y)))),2));
mycol=c("sample"="gray50","cells"="darkgreen","cellsTR"="green1","internal1"="darkorchid1","internal2"="darkorchid2","neg1"="red1","neg"="red1","DMSO"="red1","neg2"="red4","neg3"="firebrick1","neg4"="firebrick2","mneg1"="tomato1","mneg2"="tomato3","pos1"="blue1","pos"="blue1","BzCl"="blue1","bzt"="magenta4","pos2"="slateblue1","mpos1"="steelblue1","mpos2"="steelblue2","blanks"="cyan")
myshape=c("sample"=1,"cells"=16,"cellsTR"=16,"internal1"=16,"internal2"=16,"neg1"=16,"neg"=16,"DMSO"=16,"neg2"=16,"neg3"=16,"neg4"=16,"mneg1"=16,"mneg2"=16,"pos1"=16,"pos"=16,"BzCl"=16,"bzt"=16,"pos2"=16,"mpos1"=16,"mpos2"=16,"blanks"=16)


QC_statistics <- function()
{
  # browser()
  
  data_ <- unique(screen_table[c("screen_id", "Plate")]);
  
  qcstat = do.call(rbind, lapply(1:dim((data_))[1], function(i){
    plate_table=screen_table[screen_table$screen_id == data_[i,1] & screen_table$Plate == data_[i,2], ]
    # browser()
    plate_table=plate_table[with(plate_table, order(Column)), ]
    raw_datamat <- matrix(plate_table$rawIntensity,nrow=16,ncol=24,byrow=F)
    
    productUp = toupper(plate_table$ProductName);
    pos = na.omit(plate_table$rawIntensity[productUp %in% c("BZCL","POS") & plate_table$DCol!=24]); neg = na.omit(plate_table$rawIntensity[productUp %in% c("DMSO","NEG")]);
    cells_1=plate_table$rawIntensity[plate_table$ProductName=="cells" & plate_table$Column==1]
    cells_24=plate_table$rawIntensity[plate_table$ProductName=="cells" & plate_table$Column==24]
    
    outer_negcontrols=plate_table$rawIntensity[productUp %in% c("CELLS","DMSO","NEG","DMSO") & (plate_table$Column %in% c(1,24) | plate_table$DRow %in% LETTERS[c(1,16)])]
    inner_negcontrols=plate_table$rawIntensity[productUp %in% c("CELLS","DMSO","NEG","DMSO") & (plate_table$Column %in% c(2:23) & plate_table$DRow %in% LETTERS[2:15])]
    Out_To_In_Controls=round(median(outer_negcontrols)/median(inner_negcontrols),2)
    
    if(max(plate_table$readout) %in% c("CTX","SyTxG")){
      plate_ssmd <- ssmd(pos,neg);Signal_Vs_BG=round(mean(pos)/mean(neg),1)
    }else{
      plate_ssmd <- ssmd(neg,pos);Signal_Vs_BG=round(mean(neg)/mean(pos),1)
    }
    
    data.frame(Plate = paste0(data_[i,1],"_",data_[i,2]), Z_Prime = zfactor(neg,pos), Robust_Z_Prime = robustzfactor(neg,pos),SSMD = plate_ssmd, Signal_Vs_BG = Signal_Vs_BG, Out_To_In_Controls=Out_To_In_Controls,Mean_Neg = round(mean(neg),1), SD_Neg = round(pop.sd(neg),1), CV_Neg = round(raster::cv(neg),1), Mean_Pos = round(mean(pos),1), SD_Pos = round(pop.sd(pos),1), CV_Pos = round(raster::cv(pos),1), Mean_Cells_Col1 = round(mean(cells_1),1), SD_Cells_Col1 = round(pop.sd(cells_1),1), CV_Cells_Col1 = round(raster::cv(cells_1),1), Mean_Cells_Col24 = round(mean(cells_24),1), SD_Cells_Col24 = round(pop.sd(cells_24),1), CV_Cells_Col24 = round(raster::cv(cells_24),1))
  }))
  
  save(qcstat,file ="./www/RDA/QC.rda")
  openxlsx::write.xlsx(qcstat,file.path(getwd(), "www","Results", "QC_All_Plates.xlsx"), sheetName="Data",row.names=F)
  as.data.frame(qcstat)
}



################################################################################################ 
################# PDF REPORT
################################################################################################ 

PDF_QC <- function(){
  
  pdf("./www/Results/All_QC.pdf",width=24, height=12,family="Helvetica-Narrow")
  ########################################################################################################  
  ################### STAT. TABLE 
  
  par(mfrow=c(1,1),mar=c(0,0,2,0))
  maxrow = 28; npages = ceiling(nrow(qcstat)/maxrow);
  for (i in 1:npages) 
  {
    if(i<npages) qcs=qcstat[seq(1+((i-1)*maxrow), i*maxrow),] else qcs=qcstat[seq(1+((i-1)*maxrow), nrow(qcstat)),]
    
    tg <- gridExtra::tableGrob(qcs, theme = gridExtra::ttheme_default(padding = grid::unit(c(4, 10), "mm")))
    
    colcs=qcs
    colcs[c("Signal_Vs_BG","Mean_Neg", "SD_Neg", "Mean_Pos", "SD_Pos", "Mean_Cells_Col1", "SD_Cells_Col1", "Mean_Cells_Col24", "SD_Cells_Col24")] <- "blue"
    colcs$Plate<-"black"; colcs$CV_Cells_Col24 <- ifelse(colcs$CV_Cells_Col24 > 10 ,"red","blue")
    colcs$Z_Prime <- ifelse(colcs$Z_Prime < 0.5,"red","blue"); colcs$Robust_Z_Prime <- ifelse(colcs$Robust_Z_Prime < 0.5,"red","blue")
    colcs$Out_To_In_Controls <- ifelse((colcs$Out_To_In_Controls < 0.9 | colcs$Out_To_In_Controls > 1.1) ,"red","blue")
    colcs$SSMD <- ifelse(colcs$SSMD < 6,"red","blue"); colcs$CV_Neg <- ifelse(colcs$CV_Neg > 10 ,"red","blue"); 
    colcs$CV_Pos <- ifelse(colcs$CV_Pos > 10 ,"red","blue"); colcs$CV_Cells_Col1 <- ifelse(colcs$CV_Cells_Col1 > 10 ,"red","blue")
    
    gplots::textplot(qcs,col.data=as.matrix(colcs),show.rownames = F,show.colnames = T, col.colnames="red", valign="top",halign="left")
    title("QC Summary:",cex.main=2,col.main="red",adj=0)
  }
  
  ########################################################################################################  
  ################### Bar plots 
  
  multiplot <- compiler::cmpfun(function(plotlist=NULL, file, cols=1, layout=NULL) {
    numPlots = length(plotlist)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) 
      # Make the panel # ncol: Number of columns of plots # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
    
    # Set up the page
    grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = T))
      print(plotlist[[i]], vp = grid::viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  })
  
  par(mfrow=c(1,1),oma=par("mar"), mar=c(2,3,6,2))
  plot_list=lapply(unique(screen_table$screen_id),function(screen)
  {
    set_table=screen_table[screen_table$screen_id==screen,]
    ggplot2::ggplot(data=set_table[set_table$Content %in% c("pos","neg","cells","NA","empty"),], ggplot2::aes(x=Content, y=rawIntensity,fill=factor(Plate))) + ggplot2::stat_summary(fun.y=mean, geom="bar", position= ggplot2::position_dodge(), colour='white')+ggplot2::stat_summary(fun.y=mean, fun.ymin=lowsd, fun.ymax=highsd, geom="errorbar", position=ggplot2::position_dodge(.9),color = 'black', size=.5, width=0.2) + ggplot2::ggtitle(paste("Raw Intensity Controls",screen,sep="::")) + ggplot2::theme_bw() + ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),panel.grid.major = ggplot2::element_blank())
  })
  multiplot(plotlist = plot_list, cols=4)
  
  
  ########################################################################################################  
  ################### Heatmaps 
  
  heatmapfun <- compiler::cmpfun(function(datamat,heatmaptitle, normalize_ = !0){  
    
    dimnames(datamat) <- list(LETTERS[1:nrow(datamat)], 1:ncol(datamat))
    datamat <- datamat[nrow(datamat):1,]
    
    if(normalize_) datamat <- heatNorm(datamat)$datamat
    
    image(x = 1:dim(datamat)[2], y = 1:dim(datamat)[1], z = t(datamat), axes = F, xlab = '', ylab = '', 
          main=heatmaptitle, col=colorRampPalette(c("blue", "white", "red"))(100))
    axis(3,1:dim(datamat)[2], labels=colnames(datamat), las = 3, line = -0.5, tick = 0,cex.axis =1); 
    axis(2,1:dim(datamat)[1], labels=rownames(datamat), las = 2, line = -0.5, tick = 0,cex.axis =.8)  
  })
  
  nr=nc=3; totno=nr*nc
  
  par(mfrow=c(3,3),oma=par("mar"), mar=c(2,3,6,2))
  
  for(screen in as.character(unique(screen_table$screen_id)))
  {
    print(screen)
    set_table=screen_table[screen_table$screen_id==screen,]
    no=0;
    for(plate in unique(set_table$Plate))
    {
      no=no+1; print(plate)
      plate_table=set_table[set_table$Plate==plate,]
      plate_table=plate_table[with(plate_table, order(Column)), ]
      raw_datamat <- matrix(plate_table$rawIntensity,nrow=16,ncol=24,byrow=F)
      plate_id=paste(plate_table$screen_id[1],plate,sep="_")
      heatmapfun(raw_datamat,paste0("Raw data : ",plate_id), normalize_ = T)
    }
    
    while(no < totno){
      no=no+1; frame()
    }
  }
  
  ########################################################################################################  
  ################### Scatter plots 
  
  for(screen in as.character(unique(screen_table$screen_id)))
  {	
    set_table=screen_table[screen_table$screen_id==screen,]
    set_table=set_table[with(set_table, order(Plate,Row,Column)),]
    #p <- ggplot2::ggplot(data=set_table, ggplot2::aes(x=DRow, y=rawIntensity,ymin =0)) + ggplot2::facet_wrap(~Plate,nrow=3,scales="free") + ggplot2::geom_point(position = ggplot2::position_jitter(h=1,w=0.3), ggplot2::aes(colour=Content,shape=Content),show_guide = T) +  ggplot2::scale_shape_manual(values=myshape) +  ggplot2::scale_colour_manual(values=mycol) + ggplot2::scale_x_continuous(breaks=1:length(unique(set_table$DRow)),labels=LETTERS[1:length(unique(set_table$DRow))])+ ggplot2::ggtitle(paste("Raw Intensity: Row Wise  ",screen,sep="::"))  + ggplot2::theme_bw()
    p <- ggplot2::ggplot(data=set_table, ggplot2::aes(x=DRow, y=rawIntensity,ymin =0)) + ggplot2::facet_wrap(~Plate,nrow=3,scales="free") + ggplot2::geom_point(position = ggplot2::position_jitter(h=1,w=0.3), ggplot2::aes(colour=Content,shape=Content),show_guide = T, size = 7) +  ggplot2::scale_shape_manual(values=myshape) +  ggplot2::scale_colour_manual(values=mycol) + ggplot2::ggtitle(paste("Raw Intensity: Row Wise  ",screen,sep="::"))  + ggplot2::theme_bw()
    
    set_table=set_table[with(set_table, order(Plate,Column,Row)),]
    p2 <-ggplot2::ggplot(data=set_table, ggplot2::aes(x=factor(Column), y=rawIntensity,ymin =0)) + ggplot2::facet_wrap(~Plate,nrow=3,scales="free") + ggplot2::geom_point(position = ggplot2::position_jitter(h=1,w=0.3), ggplot2::aes(colour=Content,shape=Content),show_guide = T, size = 7) + ggplot2::scale_shape_manual(values=myshape) +  ggplot2::scale_colour_manual(values=mycol) + ggplot2::ggtitle(paste("Raw Intensity : Column Wise",screen,sep="::")) + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major=ggplot2::element_line(color="gray50"))
    
    p3 <-ggplot2::ggplot(data=set_table[set_table$Content %in% c("pos","neg","cells"),], ggplot2::aes(x=Column, y=rawIntensity,ymin =0)) + ggplot2::facet_wrap(~Plate,nrow=3,scales="free") + ggplot2::geom_point(position = ggplot2::position_jitter(h=1,w=0.3), ggplot2::aes(colour=Content,shape=Content),show_guide = T) + ggplot2::geom_smooth(method=lm,fullrange=F, ggplot2::aes(colour=Content)) + ggplot2::scale_x_continuous(breaks=1:length(unique(set_table$Column))) +  ggplot2::scale_shape_manual(values=myshape)  + ggplot2::scale_colour_manual(values=mycol) + ggplot2::ggtitle(paste("Raw Intensity Controls : Column Wise",screen,sep="::")) + ggplot2::theme_bw()+ ggplot2::theme(panel.grid.major=ggplot2::element_line(color="gray50"))
    
    # browser()
    
    #   print(c(list(p), list(p2)))
    
    print(c(list(p), list(p2), list(p3)))
  }
  dev.off()
}


################################################################################################ 
################# DB PUSH
###############################################################################################

# dbPush <- function(){
#   
#   # connect to DB
#   dbConnection <- function(){
#     library(ROracle); drv <- dbDriver("Oracle"); Sys.setenv(TZ = "GMT"); Sys.setenv(ORA_SDTZ = "GMT")
#     host <- "192.168.0.219"; port <- 1521; sid <- "rora"
#     connect.string <- paste0("(DESCRIPTION=", "(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))", "(CONNECT_DATA=(SID=", sid, ")))")
#     dbConnect(drv, username = "breezeadmin", password = "fdsiQmR4u", dbname = connect.string)
#   }
#   
#   tryCatch({
#     
#     con <<- dbConnection();
#     
#     # get correct DSRT set name for DB
#     DSRTsets = c("FO3D_Plate_16_Onwards", "FO3D_Plate_1_15", "FO3C_Plate_16_30", "FO2Baq", "FIMM_01B", "FIMM_01E")
#     dbDSRTsets = c("FO3D_16--", "FO3D_1-15", "FO3C_16-30", "FO2AB1-aq", "Fimm-01B", "Fimm-01E")
#     db_DSRT_set = dbDSRTsets[which(DSRT_set == DSRTsets)]; if(gtools::invalid(db_DSRT_set)) db_DSRT_set = DSRT_set
#     
#     # prepare infoframe for pushing
#     infoframe$sample_info = infoframe$screen_sample_id <- "Record Already Exists"
#     rownames(infoframe)= gsub("\\|",".",rownames(infoframe))
#     
#     # extract job id and name from path
#     db_job_id=gsub("(\\D+)","",path);
#     a=unlist(strsplit(path,"/")); db_job_name=sub("^(.*)[_].*", "\\1", a[length(a)])
#     
#     # get ids for pipeline, job status, db_status
#     db_pipline_id <- dbGetQuery(con, paste0("select PIPELINE_ID from BREEZEADMIN.CV_PIPELINE where DESCRIPTION = '","DSRT Rpipeline","'"))[[1]]
#     db_job_status_id <- dbGetQuery(con, paste0("select JOB_STATUS_ID from BREEZEADMIN.CV_JOB_STATUS where DESCRIPTION = '","SUCCESSFUL","'"))[[1]]
#     db_status_id <- dbGetQuery(con, paste0("select STATUS_ID from BREEZEADMIN.CV_STATUS_ID where DESCRIPTION = '","VALID","'"))[[1]]
#     
#     df <- data.frame(PK_JOB_ID=NA,NAME=db_job_name,DETAILS=db_job_name,PROJECT_ID=1,PIPELINE_ID=db_pipline_id,CLIENT_TYPE_ID=1,FIMM_USERNAME=report.author,STARTED_TIME=format(Sys.time(), "%d-%b-%y %I.%M.%OS6 %p"),
#                      ENDED_TIME=format(Sys.time(), "%d-%b-%y %I.%M.%OS6 %p"),CONTACT_EMAIL= paste0(report.author,"@mars.genome.helsinki.fi"), JOB_STATUS_ID=db_job_status_id,SGE_ID=1,PARAMETERS_XML="",STATUS_ID=db_status_id, USER_STAMP="spotdar",TIME_STAMP="",VERSION=NA)
#     # fill JOB Table
#     dbWriteTable(con,"JOB", df, row.names = F, overwrite = F,append = T, ora.number = T,schema = "BREEZEADMIN")
#     
#     # get CV tables
#     db_CV_CONTENT <-dbGetQuery(con, "select * from CV_CONTENT"); db_CV_PRODUCT_NAME <- dbGetQuery(con, "select * from CV_PRODUCT_NAME");
#     db_CV_DRUG_NOTATION_ID <- dbGetQuery(con, "select * from CV_DRUG_NOTATION_ID"); db_CV_CONCENTRATION_UNIT_ID <- dbGetQuery(con, "select * from CV_CONCENTRATION_UNIT_ID");
#     
#     screens <- unique(screen_table$screen_id)
#     # screens to delete
#     screen_del<- dbGetQuery(con, 
#                             paste0("SELECT PK_SCREEN_SAMPLE_ID, IDENTIFIER, SCREEN_ID FROM SCREEN_SAMPLE WHERE DESCRIPTION IN ('",paste(screens,collapse="','"),"')"))
#     
#     #delete from raw result, curve fit and Screen_sample
#     dbSendQuery(con, paste0("DELETE FROM RAW_RESULT_SCREEN WHERE SCREEN_SAMPLE_ID IN ('",paste(unique(screen_del$PK_SCREEN_SAMPLE_ID),collapse="','"),"')"))
#     dbSendQuery(con, paste0("DELETE FROM CURVE_FITTING WHERE SCREEN_ID IN ('",paste(unique(screen_del$SCREEN_ID),collapse="','"),"')"))
#     dbSendQuery(con, paste0("DELETE FROM SCREEN_SAMPLE WHERE PK_SCREEN_SAMPLE_ID IN ('",paste(unique(screen_del$PK_SCREEN_SAMPLE_ID),collapse="','"),"')"))
#     dbCommit(con);  print("DELETED")
#     
#     #Check if screens exist in the database
#     screen_tbl <- dbGetQuery(con,"select PK_SCREEN_ID, IDENTIFIER from SCREEN")
#     already_in_db_screens <- screen_tbl$PK_SCREEN_ID[screen_tbl$IDENTIFIER %in% screens]
#     
#     # create empty tables which will pe pushed to DB after subsequent filling
#     screen_smp_t = c(); screen_t = c(); raw_result_t = c();
#     
#     # write screen table
#     for (screen in screens)
#     {
#       #get screen table
#       set_table=screen_table[screen_table$screen_id==screen,]
#       if(gtools::invalid(max(set_table$read_date))) scr_date <- toupper(format(as.Date(Sys.Date(),"%d-%m-%Y"),format="%d-%b-%y")) else scr_date <- max(set_table$read_date)
#       
#       if(!(ifelse(screen %in% screen_tbl$IDENTIFIER,screen_tbl$PK_SCREEN_ID[screen_tbl$IDENTIFIER %in% screen],screen) %in% already_in_db_screens)){
#         db_readout <-dbGetQuery(con, paste0("select SCREEN_READOUT_ID from BREEZEADMIN.CV_SCREEN_READOUT where DESCRIPTION = '", unique(set_table$readout),"'"))[[1]]
#         screen_t = rbind(screen_t, data.frame(PK_SCREEN_ID=NA,IDENTIFIER=screen,DESCRIPTION=screen,DOTM_SOURCE_ID=NA,DOTM_PROTOCOL_ID=NA,SOURCE_IDENTIFIER="RPipeline",ENTITY_ID=249,BREEZE_ALIAS=screen,SAMPLE_TYPE_ID=7,DISEASE_SUB_TYPE_ID=5,MEDIA_TYPE_ID=4,HISTOLOGY_ID=12,DISEASE_STATE_ID=5,
#                                               EXPERIMENT_TYPE_ID=1,PLATE_COUNT=length(unique(set_table$Plate)),DISEASE_GRADE_ID=11,DISEASE_STAGE_ID=11,STATUS_ID=1,USER_STAMP=report.author,TIME_STAMP="",VERSION=NA,SCREEN_READOUT_ID=db_readout,SCREEN_CREATION_DATE=scr_date,ACCESS_POLICY=0,SCREEN_OWNER=report.author))
#       }
#     }
#     if(!gtools::invalid(screen_t)) dbWriteTable(con,"SCREEN", screen_t, row.names = F, overwrite = F,append = T, ora.number = T, schema = "BREEZEADMIN")
#     
#     # write screen sample table
#     for (screen in screens)
#     {
#       #get screen table
#       set_table=screen_table[screen_table$screen_id==screen,]
#       db_SCREEN_ID=dbGetQuery(con, paste0("select PK_SCREEN_ID from BREEZEADMIN.SCREEN where IDENTIFIER = '",screen,"'"))$PK_SCREEN_ID
#       
#       for (sample in paste0(screen,'_',unique(set_table$Plate)))
#       {
#         db_ORDER_NO=sub('.*_\\s*', '', sample); db_IDENTIFIER=paste0(screen,"_D",db_ORDER_NO)
#         screen_smp_t=rbind(screen_smp_t, data.frame(PK_SCREEN_SAMPLE_ID=NA,IDENTIFIER=db_IDENTIFIER,DESCRIPTION=screen,SCREEN_ID=db_SCREEN_ID,ORDER_NO=db_ORDER_NO,STATUS_ID=1,USER_STAMP=report.author,TIME_STAMP="",VERSION=NA))
#       }
#     }
#     dbWriteTable(con, "SCREEN_SAMPLE", screen_smp_t, row.names = F, overwrite = F, append = T, ora.number = T, schema = "BREEZEADMIN")
#     
#     print("Go to raw results...")
#     # write raw result table
#     tryCatch({
#       
#       for (screen in screens)
#       {
#         #get screen table
#         set_table=screen_table[screen_table$screen_id==screen,]
#         
#         #for (sample in paste0(screen,'_',unique(set_table$Plate)))
#         for (db_ORDER_NO in unique(set_table$Plate))
#         {
#           #print(paste0("screen: ",screen,";  sample: ", sample))
#           # to insert into raw result table
#           plate_data=set_table[set_table$Plate==db_ORDER_NO,]
#           sample=paste0(screen,'_',db_ORDER_NO)
#           
#           # add ids for rawresults, job, screen_sample, content, Drug notation, products' ids and names ...
#           #db_ORDER_NO=sub('.*_\\s*', '', sample); 
#           db_IDENTIFIER=paste0(screen,"_D",db_ORDER_NO)
#           db_CV_CONTENT$DESCRIPTION[db_CV_CONTENT$DESCRIPTION=="SAMPLES"]="SAMPLE" 
#           plate_data$JOB_ID <-max(dbGetQuery(con, paste0("select PK_JOB_ID from BREEZEADMIN.JOB where NAME = '",db_job_name,"'")))[[1]]
#           plate_data$SCREEN_SAMPLE_ID <-dbGetQuery(con, paste0("select PK_SCREEN_SAMPLE_ID from BREEZEADMIN.SCREEN_SAMPLE where IDENTIFIER = '",db_IDENTIFIER,"'"))[[1]]
#           plate_data$CONTENT_ID=db_CV_CONTENT$CONTENT_ID[match(toupper(plate_data$Content),db_CV_CONTENT$DESCRIPTION)]
#           plate_data$DRUG_NOTATION_ID=db_CV_DRUG_NOTATION_ID$DRUG_NOTATION_ID[match(plate_data$ProductId,db_CV_DRUG_NOTATION_ID$DESCRIPTION)]
#           plate_data$FIMM_PRODUCT_ID=db_CV_PRODUCT_NAME$FIMM_PRODUCT_ID[match(plate_data$DRUG_NOTATION_ID,db_CV_PRODUCT_NAME$DRUG_NOTATION_ID)]
#           plate_data$PRODUCT_NAME_ID=db_CV_PRODUCT_NAME$PK_PRODUCT_NAME_ID[match(plate_data$DRUG_NOTATION_ID,db_CV_PRODUCT_NAME$DRUG_NOTATION_ID)]
#           
#           Product_Upper = toupper(plate_data$Content)
#           plate_data[Product_Upper=="BLANKS", c("ProductName", "ProductId")] <- "BLANK"; plate_data[Product_Upper %in% c("CELLS","EMPTY"), c("ProductName", "ProductId")] <- "EMPTY"
#           plate_data[Product_Upper=="POS", c("ProductName", "ProductId")] <- "LOW"; plate_data[Product_Upper=="NEG", c("ProductName", "ProductId")] <- "HIGH"
#           plate_data$PLATE_IDENTIFIER= paste0(unique(plate_data$Barcode),"_D", db_ORDER_NO)	
#           plate_data$CONCENTRATION_UNIT="nM";plate_data$PK_RAW_RESULT_ID="";plate_data$STATUS_ID=1;plate_data$USER_STAMP=report.author;plate_data$TIME_STAMP="";plate_data$VERSION=""
#           
#           # replace column names to match DB
#           colnames(plate_data)[sapply(c("DWell","Concentration","DRow","rawIntensity","inhibition_percent"), function(i) which(i == colnames(plate_data)))] <- 
#             c("WELL","CONCENTRATION_AMOUNT","RROW","RAW_INTENSITY_VALUE","PERCENT_INHIBITION")
#           plate_data$RROW = as.integer(unclass(factor(plate_data$RROW))); plate_data$PERCENT_INHIBITION = round(plate_data$PERCENT_INHIBITION, 2); 
#           plate_data$PERCENT_INHIBITION[plate_data$PERCENT_INHIBITION>200] = 200; plate_data$PERCENT_INHIBITION[plate_data$PERCENT_INHIBITION<=-200] = -200; 
#           
#           db_DTM_UNIQUE_layout_plates <- dbGetQuery(con, paste0("select DTM_UNIQUE_PLATE_ID,DESCRIPTION  from BREEZEADMIN.CV_TM_STORAGE_PLATE where BREEZE_ALIAS = '",db_DSRT_set,"' AND STATUS_ID = 1 AND UPPER(description) like '%",db_ORDER_NO,"%'"))
#           db_DTM_UNIQUE_PLATE_ID <- db_DTM_UNIQUE_layout_plates$DTM_UNIQUE_PLATE_ID[grep(paste0("N",db_ORDER_NO),db_DTM_UNIQUE_layout_plates$DESCRIPTION)]		
#           plate_data$DTM_UNIQUE_PLATE_ID=rep(as.character(db_DTM_UNIQUE_PLATE_ID),nrow(plate_data))
#           plate_data$CONCENTRATION_UNIT_ID=db_CV_CONCENTRATION_UNIT_ID$CONCENTRATION_UNIT_ID[match(plate_data$CONCENTRATION_UNIT,db_CV_CONCENTRATION_UNIT_ID$TEXT)]
#           
#           # get final DB RAW_RESULT_SCREEN table
#           plate_data=plate_data[,c("PK_RAW_RESULT_ID","JOB_ID","SCREEN_SAMPLE_ID","FIMM_PRODUCT_ID","PRODUCT_NAME_ID","CONTENT_ID","DTM_UNIQUE_PLATE_ID","RROW","WELL","RAW_INTENSITY_VALUE","PERCENT_INHIBITION","PLATE_IDENTIFIER","CONCENTRATION_AMOUNT","CONCENTRATION_UNIT_ID","STATUS_ID","USER_STAMP","TIME_STAMP","VERSION")]
#           
#           # replace missing
#           filt <- is.na(plate_data$FIMM_PRODUCT_ID)
#           plate_data$FIMM_PRODUCT_ID[filt] <- 20035; plate_data$PRODUCT_NAME_ID[filt] <- 20387; 
#           #plate_data$CONTENT_ID[filt] <- 4
#           
#           raw_result_t = rbind(raw_result_t, plate_data);
#           
#           sample_id_small=paste0(screen,".",db_ORDER_NO)
#           infoframe$sample_info[match(sample_id_small,rownames(infoframe))]=sample
#           infoframe$screen_sample_id[match(sample_id_small,rownames(infoframe))]=plate_data$SCREEN_SAMPLE_ID[[1]]
#         }
#       }
#     }, 
#     error = function(e) {message(e); print("I am failed at filling screen table")}
#     )
#     dbWriteTable(con, "RAW_RESULT_SCREEN", raw_result_t, row.names = F, overwrite = F, append = T, ora.number = T, schema = "BREEZEADMIN")
#     openxlsx::write.xlsx(infoframe, file = "./Results/Screen_Database.xlsx")
#   }, 
#   error = function(e) {message(e); print("I am failed at QC pushing")},
#   finally = {
#     tryCatch({ dbDisconnect(con); },
#              error = function(e) { print("Can't connect to DB")}); 
#     print("QC push end")
#   }
#   )
# }


# 
# if(registered == !0){
# dbwriteinfoframe=infoframe;dbwriteinfoframe$DSRT_Set = DSRT_set;dbwriteinfoframe$ID=NULL;colnames(dbwriteinfoframe)[which(names(dbwriteinfoframe) == "screen_id")] <- "screen_name";
# mydb = dbConnect(MySQL(),user='breezeuser',password='Breezeshinyatlas123!!!',dbname='breezedata', host='127.0.0.1')
# dbwriteinfoframe$USER_ID=DBI::dbGetQuery(mydb, paste0('SELECT USER_ID FROM  user_and_access WHERE USER_EMAIL = "',name,'"'))$USER_ID[[1]]
# for (i in 1:length(unique(dbwriteinfoframe$screen_name))) {
# dbSendQuery(mydb, paste0(c("DELETE FROM screen_table WHERE screen_name = '", unique(dbwriteinfoframe$screen_name)[i], "' AND USER_ID = '", dbwriteinfoframe$USER_ID[1], "';"), collapse = ""))
# }
# dbWriteTable(conn = mydb, "screen_table", dbwriteinfoframe,append = T,overwrite=F,row.names = FALSE)
# dbDisconnect(mydb)}


HTML_QC_report <- function(){
  # combine pushing QC to DB with generating PDF.
  #toRun <- list(pdf_ = PDF_QC, db_ = dbPush)
  #mclapply(toRun, function(x) x(), mc.cores = 2)
  if(exists("DBQCCF")) dbPush(); 
  PDF_QC();
  
  ########################################################################################################  
  ################### Saving excel files and html table(which is used in qc.html, instead of noozle)
  
  screen_table=screen_table[with(screen_table, order(screen_id,Plate,Column,Row)),]
  if(!is.null(screen_table$Barcode)){
    st=screen_table[,c("DWell","DRow","DCol","Plate","ProductId","ProductName","Concentration","Content","screen_id","Barcode","read_date","readout","rawIntensity","inhibition_percent")]
    st$Content <- as.character(st$Content); st$Content[st$ProductName=="dmso"] <- "DMSO"
    st$Content[st$ProductName=="BzCl"] <- "BzCl"; st$Content[st$ProductName=="cells"] <- "cells"
    openxlsx::write.xlsx(st,paste0("./www/Results/Raw_Data",".xlsx"),sheetName="Raw_Data")
  }
  
  
  # create table and convert it to HTML
  htable <- print(xtable::xtable(qcstat[complete.cases(qcstat$Plate), ], caption = "QC Summary:", digits = c(0, 0, 2, 2, 0, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)), include.rownames = FALSE, 
                  caption.placement = "top", type = "html", print.results = F, html.table.attributes = getOption("xtable.html.table.attributes","border=0"))
  
  # change table alignment
  htable <- gsub("right", "center", htable)
  h2table <- sub("<table", "<table id=\"stattab\" data-row-style=\"rowStyle\" data-toggle=\"table\" data-show-refresh=\"true\" data-show-toggle=\"true\" data-show-columns=\"true\" data-search=\"true\" data-select-item-name=\"toolbar1\" data-pagination=\"true\" data-sort-name=\"name\" data-sort-order=\"desc\" ", htable)
  h2table <- gsub("\\s", "", h2table); h2table <- gsub("align", " align", h2table)
  h2table <- gsub("id", " id", h2table); h2table <- gsub("QCSummary", "", h2table)
  h2table <- gsub('<tr><th>', '<thead><th>', h2table); h2table <- gsub('</th></tr>', '</th></thead>', h2table);
  
  
  ########################################################################################################  
  ################### HTML report. JS functions and css are in headerQC.txt and tailerQC.txt.
  
  headerpath="./headerQC.txt"
  tailerpath="./tailerQC.txt"
  
  
  
  reportName = "./www/Results/All_QC.html"
  
  ################################################################################################ 
  ########## HTML table
  
  htable <- sub("<table", "<table id=\"stattab\"", htable); htable <- gsub("\\s", "", htable)
  htable <- gsub('align="center"', '', htable); htable <- gsub("id", " id", htable)
  htable <- gsub("QCSummary", "QC Summary", htable)
  # if neccesary#htable <- gsub('<tr><th>', '<thead><th>', htable); htable <- gsub('</th></tr>', '</th></thead>', htable);
  
  # load prepared and minified css and js
  header <- base::readChar(headerpath, file.info(headerpath)$size)
  tailer <- base::readChar(tailerpath, file.info(tailerpath)$size)
  
  # prepare body
  body <- ""
  screens = unique(screen_table$screen_id)
  
  ################################################################################################ 
  ########  Search in HTML file (search for drugs)
  
  # search information (plate + well)
  Searchinfo = paste0(screen_table$DWell, ", ", screen_table$screen_id, "_", screen_table$Plate)
  # create search table based on this info
  searchTable = list(ProductName = screen_table$ProductName, Searchinfo = Searchinfo)
  searchTable = split(searchTable$Searchinfo, searchTable$ProductName)
  
  # create html for dropdown used to search.  dropDownSearchDrug - html for dropdown(including a list of drugs),# sT - is a hover message showing the (plate + well) information found.
  dropDownSearchDrug = "<select id='dsearch' class='selectpicker' title='Choose Product...' data-style='btn-info' data-live-search='true'>"
  sT = "var sT=["
  names_ = names(sapply(searchTable, names))
  
  for (k in 1:length(names_)) {
    sT = paste0(sT, "[\"<p>", paste0(searchTable[[names_[k]]], collapse = "</p><p>"), "</p>\"],")
    dropDownSearchDrug = paste0(dropDownSearchDrug, "<option value='", k, "'>", names_[k], "</option>")
  }
  sT = substr(sT, 1, nchar(sT) - 1); sT = paste0(sT, "]")
  dropDownSearchDrug = paste0(dropDownSearchDrug, "</select>")
  
  ################################################################################################ 
  ######### Bar graphs
  print("bar")
  
  bd <- do.call(paste0, lapply(1:length(screens), function(i){
    screen = screens[i]
    set_table = screen_table[screen_table$screen_id == screen, ]
    
    # create bar charts ggplot
    g <- ggplot(data = set_table[set_table$Content %in% c("pos", "neg", "cells", "NA", "empty"), ], aes(x = Content, y = rawIntensity, fill = factor(Plate))) + 
      stat_summary(fun.y = mean, geom = "bar", position = position_dodge(), colour = "white") + stat_summary(fun.y = mean, fun.ymin = lowsd, fun.ymax = highsd, 
                                                                                                             geom = "errorbar", position = position_dodge(0.9), color = "black", size = 0.5, width = 0.2) + ggtitle(paste0("<b>", "Raw Intensity Controls::", screen, 
                                                                                                                                                                                                                           "</b>")) + theme_bw() + theme(panel.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + guides(fill = guide_legend(title = "Plate", 
                                                                                                                                                                                                                                                                                                                                                                                                  title.position = "top")) + labs(list(x = "", y = ""))
    save.image("philtestr.RDATA")
    # convert to widget and adjust some parameters
    LyBar <- plotly::as.widget(plotly::ggplotly(g) %>% plotly::config(showLink = F) %>% layout(titlefont = list(size = 14), xaxis = list(title = "content", titlefont = list(size = 12)), 
                                                                                               yaxis = list(title = "raw intensity", titlefont = list(size = 12)), legend = list(tracegroupgap = 3), margin = list(l = 100, b = 35, t = 67)))
    # create JSON
    JsonBar <- htmlwidgets:::toJSON(LyBar)[1]
    # rewrite JSON end
    obj <- paste0(substr(JsonBar, 1, regexpr("\"showLink\":false", JsonBar)[1]), "showLink\":false}},\"evals\":[],\"jsHooks\":[]}}</script>")
    # optimize JSON by rounding each double to 2 dec.  !!!!!!!!!!!!!! do not comment for high performance of HTML!
    obj <- gsubfn::gsubfn("[[:digit:]]+\\.+[[:digit:]]+[[:digit:]]", ~format(round(as.numeric(x), 2), nsmall = 2), obj)

    
    # make container name unique by extracting current proc.time to make it unique (used to put json inside of it)
    widContN <- sub("\\.", "", proc.time()[3])
    widContN <- paste0(gsubfn::gsubfn("[[:digit:]]", ~format(letters[as.numeric(x) + 1]), widContN), letters[i])
    
    # create container which is regulated throug .popka2 class and put JSON inside
    widContainer <- paste0("<div id=\"", widContN, "\" class=\"popka\">", "</div>")
    barGraph <- paste0("<script type=\"application/json\" data-for=\"", widContN, "\">", obj, "</script>")
    renderJs <- paste0("<script> var Data = document.querySelector(\"script[data-for='", widContN, "']\"); var x = JSON.parse(Data.textContent || Data.text).x;Plotly.plot('", 
                       widContN, "', x.data, x.layout, x.config);</script>")
    
    paste(widContainer, barGraph, renderJs, sep = "\n")
  }));
  
  body<- paste0(body, bd)
  
  # substitute some 'imperfections' in obtained JSON
  body <- gsub("factor\\(Plate\\)", "Plate", body)
  body <- gsub("rawIntensity", "Raw intensity", body)
  body <- gsub("\"text\":\"Plate\",\"x\":1.02,\"y\":1", "\"text\":\"Plate\",\"x\":1.02,\"y\":0.95", body)
  
  # create space before next type of plots goes
  body <- paste0(body, "<br><br>")
  
  ################################################################################################ Heat maps
  print("heat")
  for (screen in screens) {
    set_table = screen_table[screen_table$screen_id == screen, ]
    plates = unique(set_table$Plate)
    
    bd <- do.call(paste0, mclapply(1:length(plates), function(i){
      
      plate = plates[i]
      # get data
      plate_table = set_table[set_table$Plate == plate, ]
      plate_table = plate_table[with(plate_table, order(Column)), ]
      plate_id = paste(plate_table$screen_id[1], plate, sep = "_")
      
      # data to display on hover
      raw_datamatdisp <- matrix(paste0("Product name: ", plate_table$ProductName, "<br>Product id: ", plate_table$ProductId, "<br>Raw intensity: ", plate_table$rawIntensity, 
                                       "<br>Concentration: ", plate_table$Concentration), nrow = 16, ncol = 24, byrow = F)
      
      # matrix dispayed in heat map.
      raw_datamat <- matrix(plate_table$rawIntensity, nrow = 16, ncol = 24, byrow = F)
      raw_datamat_content <- matrix(plate_table$Content, nrow = 16, ncol = 24, byrow = F)
      
      # reverse them (needed for correct plotly build)
      raw_datamat_corrected <- apply(t(raw_datamat), 1, rev)
      raw_datamat_cont_corrected <- apply(t(raw_datamat_content), 1, rev)
      pr_name_corrected <- apply(t(raw_datamatdisp), 1, rev)
      data_norm <- heatNorm(raw_datamat_corrected); raw_datamat_corrected <- data_norm$datamat; 
      
      #exclude pos and neg from outliers and reshape
      noPosandNeg = !matrix(grepl(("POS|NEG"),toupper(raw_datamat_cont_corrected)), nrow = 16, ncol = 24, byrow = F)
      outl_ = data_norm$outliers  & noPosandNeg; outl_ =  which(outl_, arr.ind = T)
      # df_ = reshape2::melt(raw_datamat_corrected); colnames(df_) = c("row","col","dt"); 
      
      # create heat map
      plotheat <- subplot(  plotly:: plot_ly(z = raw_datamat_corrected,type = "heatmap",
                                             # colorscale = list(list(0,'rgb(0,0,255)'), list(0.5,'rgb(255,255,255)'), list(1,'rgb(255,0,0)')),
                                             colors = c("blue", "white", "red"),
                                             colorbar = list(outlinewidth = 0, title = "raw <br>intensity", len = 0.9, thickness = 15,
                                                             tickcolor = "#fff", ticks = "inside"), hoverinfo = "y+x+text", text = pr_name_corrected) %>%
                              # add_trace(y=outl_[,1], x=outl_[,2],  data = "possible outlier", showlegend=F) %>%
                              plotly::config(showLink = F) %>%
                              layout(title = paste0("Raw data: ", plate_id, "<br>"), margin = list(t = 100, b = 50, l = 20, r = 10),
                                     xaxis = list(title = "", tickmode = "array", side = "top", ticklen = 1, linecolor = "#ffffff", linewidth = 4,
                                                  tickvals = seq.int(0, ncol(raw_datamat_corrected), by = 1), tickfont = list(family = "serif", size = 10),
                                                  ticks = "outside", ticktext = as.character(seq.int(1, ncol(raw_datamat_corrected), by = 1))),
                                     yaxis = list(title = "", tickmode = "array", ticklen = 1, linecolor = "#ffffff", linewidth = 4,
                                                  tickvals = seq.int(0, nrow(raw_datamat_corrected), by = 1), tickfont = list(family = "serif", size = 10),
                                                  ticks = "outside", ticktext = rev(LETTERS[1:nrow(raw_datamat_corrected)]))))
      #htmltools::save_html(as.widget(plotheat), file = "test2.html")
      
      # create JSON
      JsonHeat <- htmlwidgets:::toJSON(plotly::as_widget(plotly_build(plotheat)))[1]
      
      # modify JSON end
      obj <- gsubfn::gsubfn("[[:digit:]]+\\.+[[:digit:]]+[[:digit:]]", ~format(round(as.numeric(x), 2), nsmall = 2), JsonHeat)
      # obj <- paste0(substr(JsonHeat, 1, regexpr("\"showLink\":false", JsonHeat)[1]), "showLink\":false}},\"evals\":[],\"jsHooks\":[]}")
      
      # optimize JSON by rounding each double to 2 dec.  
      # obj <- gsubfn::gsubfn("[[:digit:]]+\\.+[[:digit:]]+[[:digit:]]", ~format(round(as.numeric(x), 2), nsmall = 2), obj)
      
      # make container name unique by extracting current proc.time to make it unique (used to put json inside of it)
      widContN <- sub("\\.", "", proc.time()[3])
      widContN <- paste0(gsubfn::gsubfn("[[:digit:]]", ~format(letters[as.numeric(x) + 1]), widContN), letters[i])
      widContN <- paste0(widContN,screen,plate);
      
      # (regulate size throug .popka2) wrapContN <- paste0('w', widContN)
      widContainerWide <- paste0("<div id=\"plate", as.character(plate), "\" class=\"plate popka2\">")
      widContainer <- paste0("<div id=\"", widContN, "\" >", "</div>")
      
      heatMap_ <- paste0("<script type=\"application/json\" data-for=\"", widContN, "\">", obj, "</script>")
      
      renderJs <- paste0("<script> var Data = document.querySelector(\"script[data-for='", widContN, "']\"); var x = JSON.parse(Data.textContent || Data.text).x;Plotly.plot('", 
                         widContN, "', x.data, x.layout, x.config);</script>")
      
      paste(widContainerWide, widContainer, heatMap_, renderJs, "</div>", sep = "\n")
      
    }));
    body<- paste0(body, bd, "<br>")
    
  }
  
  ################################################################################################ First scatters
  
  for (screen in screens) {
    set_table = screen_table[screen_table$screen_id == screen, ]
    set_table = set_table[with(set_table, order(Plate, Row, Column)), ]
    
    plates = unique(set_table$Plate)
    
    bd <- do.call(paste0, mclapply(1:length(plates), function(i){
      
      plate = plates[i]
      
      plate_table = set_table[set_table$Plate == plate, ]
      
      scatpl <- ggplot2::ggplot(data=plate_table, ggplot2::aes(x=DRow, y=rawIntensity,ymin =0)) + ggplot2::facet_wrap(~Plate,nrow=3,scales="free") + ggplot2::geom_point(position = ggplot2::position_jitter(h=1,w=0.3), ggplot2::aes(colour=Content,shape=Content),show_guide = T) +  ggplot2::scale_shape_manual(values=myshape) +  ggplot2::scale_colour_manual(values=mycol) + ggplot2::ggtitle(paste("Raw Intensity: Row Wise  ",screen,sep="::"))  + ggplot2::theme_bw() + labs(x = " ", y = " ") + theme(plot.title = element_text(size = 10), legend.title = element_text(size = 10))
      plb <- plotly::plotly_build(scatpl)
      
      plb$x$data[[4]]$marker$opacity = 0.5
      # 
      for (i in 1:length(plb$x$data)) {
        #    plb$x$data[[i]]$name = substr(plb$x$data[[i]]$name, 2, regexpr(",", plb$x$data[[i]]$name)[1] - 1)  #(cells,cells) to cells
        pltCell = plate_table[plate_table$Content == plb$x$data[[i]]$name, ]
        
        j = 1:length(plb$x$data[[i]]$text)
        plb$x$data[[i]]$text[j] = paste0("(", pltCell$Column[j], ", ", pltCell$Row[j], ")<br>", "Product name: ", pltCell$ProductName[j], "<br>", "Product id: ", 
                                         pltCell$ProductId[j], "<br>", "Raw Intensity: ", pltCell$rawIntensity[j], "<br>", "Concentration: ", pltCell$Concentration[j])
      }
      
      plb$x$layout$xaxis$title = "plate row"; plb$x$layout$yaxis$title = "raw intensity"; plb$x$layout$xaxis$titlefont$size = 12; plb$x$layout$yaxis$titlefont$size = 12; plb$x$layout$annotations[[1]]$x = 1.07; plb$x$layout$annotations[[1]]$y = 0.955
      # plb$layout$annotations[[1]]$text = 'Content:'
      
      JsonHeat <- htmlwidgets:::toJSON(plotly::as_widget(plb))[1]
      #obj <- paste0(substr(JsonHeat, 1, regexpr("\"elementId\":null,\"preRenderHook\":null,", JsonHeat)[1]), "elementId\":null,\"preRenderHook\":null,\"evals\":[],\"jsHooks\":[]}")
      obj <- gsubfn::gsubfn("[[:digit:]]+\\.+[[:digit:]]+[[:digit:]]", ~format(round(as.numeric(x), 2), nsmall = 2), JsonHeat)
      # obj <- JsonHeat
      
      # make container name unique by translating current proc.time into string.
      widContN <- sub("\\.", "", proc.time()[3])
      widContN <- paste0(gsubfn::gsubfn("[[:digit:]]", ~format(letters[as.numeric(x) + 1]), widContN), letters[i])
      widContN <- paste0(widContN,screen,plate);
      
      # widContainer <- paste0('<div id='', widContN, '' class='plotly html-widget'></div>')
      
      
      # (regulate size throug .popka2) wrapContN <- paste0('w', widContN)
      widContainerWide <- paste0("<div id=\"plate", as.character(plate), "\" class=\"plate popka2\">")
      widContainer <- paste0("<div id=\"", widContN, "\">", "</div>")
      
      scatF <- paste0("<script type=\"application/json\" data-for=\"", widContN, "\">", obj, "</script>")
      
      renderJs <- paste0("<script> var Data = document.querySelector(\"script[data-for='", widContN, "']\"); var x = JSON.parse(Data.textContent || Data.text).x;Plotly.plot('", 
                         widContN, "', x.data, x.layout, x.config);</script>")
      
      paste(widContainerWide, widContainer, scatF, renderJs, "</div>", sep = "\n")
    }));
    
    body<- paste0(body, bd, "<br>")
  }
  
  
  ### THE COMMENTS FOR FOLLOWING SECTIONS ARE OMMITED (THEY ARE THE SAME AS FOR THE 'FIRST SCATTERS' SECTION)
  
  print("secscat")
  ################################################################################################ Second scatters
  for (screen in screens) {
    set_table = screen_table[screen_table$screen_id == screen, ]
    set_table = set_table[with(set_table, order(Plate, Row, Column)), ]
    
    plates = unique(set_table$Plate)
    
    bd <- do.call(paste0, mclapply(1:length(plates), function(i){
      
      plate = plates[i]
      
      plate_table = set_table[set_table$Plate == plate, ]
      
      scatpl <- ggplot(data = plate_table, aes(x = factor(Column), y = rawIntensity, ymin = 0)) + geom_point(position = position_jitter(h = 1, w = 0.3), 
                                                                                                             aes(colour = Content, shape = Content), show.legend = TRUE, size = 1) + scale_shape_manual(values = myshape) + scale_colour_manual(values = mycol) + 
        theme(panel.grid.major = element_line(color = "gray50")) + ggtitle(paste("Raw Intensity : Column Wise", screen, sep = "::")) + theme_bw() + labs(x = " ", 
                                                                                                                                                         y = " ") + theme(plot.title = element_text(size = 10), legend.title = element_text(size = 10))
      
      
      plb <- plotly::plotly_build(scatpl)
      
      plb$x$data[[4]]$marker$opacity = 0.5
      
      for (i in 1:length(plb$x$data)) {
        #    plb$x$data[[i]]$name = substr(plb$x$data[[i]]$name, 2, regexpr(",", plb$x$data[[i]]$name)[1] - 1)  #(cells,cells) to cells
        pltCell = plate_table[plate_table$Content == plb$x$data[[i]]$name, ]
        
        j = 1:length(plb$x$data[[i]]$text)
        plb$x$data[[i]]$text[j] = paste0("(", pltCell$Column[j], ", ", pltCell$Row[j], ")<br>", "Product name: ", pltCell$ProductName[j], "<br>", "Product id: ",
                                         pltCell$ProductId[j], "<br>", "Raw Intensity: ", pltCell$rawIntensity[j], "<br>", "Concentration: ", pltCell$Concentration[j])
      }
      
      plb$x$layout$xaxis$title = "plate column";plb$x$layout$yaxis$title = "raw intensity";plb$x$layout$xaxis$titlefont$size = 12; plb$x$layout$yaxis$titlefont$size = 12; plb$x$layout$annotations[[1]]$x = 1.07; plb$x$layout$annotations[[1]]$y = 0.955
      # plb$layout$annotations[[1]]$text = 'Content:'
      
      JsonHeat <- htmlwidgets:::toJSON(plotly::as_widget(plb))[1]
      # obj <- paste0(substr(JsonHeat, 1, regexpr("\"elementId\":null,\"preRenderHook\":null,", JsonHeat)[1]), "elementId\":null,\"preRenderHook\":null,\"evals\":[],\"jsHooks\":[]}")
      obj <- gsubfn::gsubfn("[[:digit:]]+\\.+[[:digit:]]+[[:digit:]]", ~format(round(as.numeric(x), 2), nsmall = 2), JsonHeat)
      # obj <- JsonHeat
      
      # make container name unique by translating current proc.time into string.
      widContN <- sub("\\.", "", proc.time()[3])
      widContN <- paste0(gsubfn::gsubfn("[[:digit:]]", ~format(letters[as.numeric(x) + 1]), widContN), letters[i])
      widContN <- paste0(widContN,screen,plate);
      
      # widContainer <- paste0('<div id='', widContN, '' class='plotly html-widget'></div>')
      
      
      # (regulate size throug .popka2) wrapContN <- paste0('w', widContN)
      widContainerWide <- paste0("<div id=\"plate", as.character(plate), "\" class=\"plate popka2\">")
      widContainer <- paste0("<div id=\"", widContN, "\">", "</div>")
      
      scatF <- paste0("<script type=\"application/json\" data-for=\"", widContN, "\">", obj, "</script>")
      
      renderJs <- paste0("<script> var Data = document.querySelector(\"script[data-for='", widContN, "']\"); var x = JSON.parse(Data.textContent || Data.text).x;Plotly.plot('", 
                         widContN, "', x.data, x.layout, x.config);</script>")
      
      paste(widContainerWide, widContainer, scatF, renderJs, "</div>", sep = "\n")
    }));
    
    body<- paste0(body, bd, "<br>")
  }
  
  ################################################################################################ Third scatters
  for (screen in screens) {
    set_table = screen_table[screen_table$screen_id == screen, ]
    set_table = set_table[with(set_table, order(Plate, Row, Column)), ]
    
    plates = unique(set_table$Plate)
    
    bd <- do.call(paste0, mclapply(1:length(plates), function(i){
      
      plate = plates[i]
      
      plate_table = set_table[set_table$Plate == plate, ]
      
      scatpl <- ggplot2::ggplot(data = plate_table[plate_table$Content %in% c("pos", "neg", "cells"), ], ggplot2::aes(x = Column, y = rawIntensity, ymin = 0)) + ggplot2::geom_point(position = ggplot2::position_jitter(h = 1, 
                                                                                                                                                                                                                         w = 0.3), ggplot2::aes(colour = Content, shape = Content), show.legend = TRUE, size = 1) + ggplot2::geom_smooth(method = lm, fullrange = FALSE, ggplot2::aes(colour = Content)) + 
        ggplot2::scale_x_continuous(breaks = 1:length(unique(set_table$Column))) + ggplot2::scale_shape_manual(values = myshape) + ggplot2::scale_colour_manual(values = mycol) + 
        ggplot2::ggtitle(paste("Raw Intensity Controls : Column Wise", screen, sep = "::")) + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray50")) + 
        ggplot2::labs(x = " ", y = " ") + ggplot2::theme(plot.title = ggplot2::element_text(size = 10), legend.title = ggplot2::element_text(size = 10))
      
      
      plb <- plotly::plotly_build(scatpl)
      # 
      for (i in 1:3) {
        #   plb$x$data[[i]]$name = substr(plb$x$data[[i]]$name, 2, regexpr(",", plb$x$data[[i]]$name)[1] - 1)  #(cells,cells) to cells
        pltCell = plate_table[plate_table$Content == plb$x$data[[i]]$name, ]
        
        j = 1:length(plb$x$data[[i]]$text)
        
        plb$x$data[[i]]$text[j] = paste0("(", pltCell$Column[j], ", ", pltCell$Row[j], ")<br>", "Product name: ", pltCell$ProductName[j], "<br>", "Raw Intensity: ",
                                         pltCell$rawIntensity[j], "<br>")
      }
      
      for (i in 4:length(plb$x$data)) {
        plb$x$data[[i]]$hoverinfo = "y"
        plb$x$data[[i]]$text = ""
      }
      
      plb$x$layout$xaxis$title = "plate column"; plb$x$layout$yaxis$title = "raw intensity"; plb$x$layout$xaxis$titlefont$size = 12; plb$x$layout$yaxis$titlefont$size = 12; plb$x$layout$annotations[[1]]$x = 1.07; plb$x$layout$annotations[[1]]$y = 0.955
      # plb$layout$annotations[[1]]$text = 'Content:'
      
      JsonHeat <- htmlwidgets:::toJSON(plotly::as_widget(plb))[1]
      # obj <- paste0(substr(JsonHeat, 1, regexpr("\"elementId\":null,\"preRenderHook\":null,", JsonHeat)[1]), "elementId\":null,\"preRenderHook\":null,\"evals\":[],\"jsHooks\":[]}")
      obj <- gsubfn::gsubfn("[[:digit:]]+\\.+[[:digit:]]+[[:digit:]]", ~format(round(as.numeric(x), 2), nsmall = 2), JsonHeat)
      # obj <- JsonHeat
      
      # make container name unique by translating current proc.time into string.
      widContN <- sub("\\.", "", proc.time()[3])
      widContN <- paste0(gsubfn::gsubfn("[[:digit:]]", ~format(letters[as.numeric(x) + 1]), widContN), letters[i])
      widContN <- paste0(widContN,screen,plate);
      # widContainer <- paste0('<div id='', widContN, '' class='plotly html-widget'></div>')
      
      
      # (regulate size throug .popka2) wrapContN <- paste0('w', widContN)
      widContainerWide <- paste0("<div id=\"plate", as.character(plate), "\" class=\"plate popka2\">")
      widContainer <- paste0("<div id=\"", widContN, "\">", "</div>")
      
      scatF <- paste0("<script type=\"application/json\" data-for=\"", widContN, "\">", obj, "</script>")
      
      renderJs <- paste0("<script> var Data = document.querySelector(\"script[data-for='", widContN, "']\"); var x = JSON.parse(Data.textContent || Data.text).x;Plotly.plot('", 
                         widContN, "', x.data, x.layout, x.config);</script>")
      
      paste(widContainerWide, widContainer, scatF, renderJs, "</div>", sep = "\n")
    }));
    body<- paste0(body, bd, "<br>")
  }
  
  # create html cat(paste(header, body, tailer, sep = ''), file=reportName, sep='') dropDown = paste0(dropDown, '</select>')
  header <- sub("Product search:</h5>", paste0("Product search:</h5>", dropDownSearchDrug), header)
  tailer <- sub("var sT=", sT, tailer)
  
  
  reportData <- paste0(header, htable, "</div><br><br>", body, tailer)
  reportData = paste0(gsub(",(\"source\":).*", ',"source":"A","config":{"modeBarButtonsToRemove":["sendDataToCloud"]},"base_url":"https://plot.ly"},"width":null,"height":null,"sizingPolicy":{"defaultWidth":null,"defaultHeight":null,"padding":5,"viewer":{"defaultWidth":null,"defaultHeight":null,"padding":null,"fill":true,"suppress":false,"paneHeight":null},"browser":{"defaultWidth":null,"defaultHeight":null,"padding":null,"fill":true},"knitr":{"defaultWidth":null,"defaultHeight":null,"figure":true}},"dependencies":null,"elementId":null,"preRenderHook":null,"evals":[],"jsHooks":[]}</script>', strsplit(reportData, "}</script>")[[1]]), collapse = "")
  reportData <- gsub('}\"}],\"cloud\":false,\"showLink\":false}},\"evals\":\\[],\"jsHooks\":\\[]}</script>','"}],"cloud":false},"source":"A","config":{"modeBarButtonsToRemove":["sendDataToCloud"]},"base_url":"https://plot.ly"},"width":null,"height":null,"sizingPolicy":{"defaultWidth":null,"defaultHeight":null,"padding":5,"viewer":{"defaultWidth":null,"defaultHeight":null,"padding":null,"fill":true,"suppress":false,"paneHeight":null},"browser":{"defaultWidth":null,"defaultHeight":null,"padding":null,"fill":true},"knitr":{"defaultWidth":null,"defaultHeight":null,"figure":true}},"dependencies":null,"elementId":null,"preRenderHook":null,"evals":[],"jsHooks":[]}</script>',reportData)

  writeChar(reportData, reportName, nchar(reportData, type = "chars"))
  
  return(h2table)
}


