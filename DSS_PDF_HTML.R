lapply(c("plotly", "scales", "parallel", "foreach", "gridExtra", "grid", "graphics", "gplots", "ggplot2", "raster", "xtable","plyr","openxlsx"), library, character.only = !0)
graphics.off()

###########################################################################################
############ Distance functions (for clustering)

na.dist <- function(t.dist,...) {
  t.dist <- as.matrix(t.dist); t.limit <- 1.1*max(t.dist,na.rm=T); if (min(t.dist,na.rm=T)==-99999){t.dist[which(t.dist==-99999)]<-t.limit}else{t.dist[is.na(t.dist)] <- t.limit}
  return(as.dist(t.dist))}

mydist <- function(x,dmethod,usecor){
  if(dmethod %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) return (na.dist(dist(x,dmethod)));if(dmethod %in% c("euclid","cosangle","abscosangle")) 	return(na.dist(as.dist(as.matrix(distancematrix(x,d=dmethod)))));
  if(dmethod %in% c("pearson","abspearson", "kendall", "abskendall", "spearman", "absspearman")){if(dmethod=='pearson')	return(na.dist(as.dist(1-cor(t(x), method=dmethod, use=usecor))));if(dmethod=='abspearson') return(na.dist(as.dist(1-abs(cor(t(x), method='pearson', use=usecor))))); if(dmethod=='kendall') return(na.dist(as.dist(1-cor(t(x), method=dmethod, use=usecor))));if(dmethod=='abskendall') return(na.dist(as.dist(1-abs(cor(t(x), method='kendall', use=usecor)))));if(dmethod=='spearman') return(na.dist(as.dist(1-cor(t(x), method=dmethod, use=usecor))));if(dmethod=='absspearman') return(na.dist(as.dist(1-abs(cor(t(x), method='spearman', use=usecor)))));}}


###########################################################################################
############ files and annotations
#directory_name <- getwd()

# read combi_table
load("final_tbl_conc_drugs.rds"); Original_DSS2 <- final_tbl_conc_drugs; rm(final_tbl_conc_drugs); DSS_cols = grep("DSS", colnames(Original_DSS2)); 
screens_=Original_DSS2[, grep("Experiment_id", colnames(Original_DSS2))]; 
if(typeof(screens_) != "character") screens_ = unname(apply(screens_, 2, function(c) unique(na.omit(c)))) else screens_ = unique(na.omit(screens_));

Name = paste0('Waterfall_plots_of_50 Active Drugs (Original DSS',DSS_typ,"_",Sys.Date(),')');
PDFname = paste0(Name,".pdf"); HTMLname = paste0(Name,".html"); HEATname = paste0('Heatmap (Original DSS',DSS_typ,"_",Sys.Date(),').html');HEATname2 = paste0('Heatmap (Selective sDSS',DSS_typ,"_",Sys.Date(),').html');Treename= paste0("OriginalDSStree",".html");HEATauc = paste0('Heatmap (AUC DSS',DSS_typ,"_",Sys.Date(),').html');

# Adding of Mechanism/Targets and  Class.explained 
if(layout_ == "fimm"){
  fo5a=openxlsx::read.xlsx("./datasets/Drug_annotations_table/FO5A annotations.xlsx");colnames(fo5a)[1]="ID";
  merged_and_ordered <- join(Original_DSS2[,c("ID","Mechanism/Targets")], fo5a[,c("ID","Mechanism/Targets")],by="ID",match="first")
  merged_and_ordered2 <- join(Original_DSS2[,c("ID","Class.explained")], fo5a[,c("ID","Class.explained")],by="ID",match="first")
  row.names(merged_and_ordered)=row.names(Original_DSS2); row.names(merged_and_ordered2)=row.names(Original_DSS2)
  merged_and_ordered[,2]=NULL; merged_and_ordered2[,2]=NULL
  Original_DSS2$`Mechanism/Targets`=merged_and_ordered$`Mechanism/Targets`
  Original_DSS2$`Class.explained`=merged_and_ordered2$`Class.explained`
  }else{
  fo5a=openxlsx::read.xlsx("./datasets/Drug_annotations_table/FO5A annotations.xlsx")
  merged_and_ordered <- join(Original_DSS2[,c("DRUG_NAME","Mechanism/Targets")], fo5a[,c("DRUG_NAME","Mechanism/Targets")],by="DRUG_NAME",match="first")
  merged_and_ordered2 <- join(Original_DSS2[,c("DRUG_NAME","Class.explained")], fo5a[,c("DRUG_NAME","Class.explained")],by="DRUG_NAME",match="first")
  row.names(merged_and_ordered)=row.names(Original_DSS2); row.names(merged_and_ordered2)=row.names(Original_DSS2)
  merged_and_ordered[,2]=NULL; merged_and_ordered2[,2]=NULL
  Original_DSS2$`Mechanism/Targets`=merged_and_ordered$`Mechanism/Targets`
  Original_DSS2$`Class.explained`=merged_and_ordered2$`Class.explained`}


if(sDSS){
  denrow=openxlsx::read.xlsx(paste0("./www/Results/DSS/","Selective_DSS",as.character(DSS_typ),"_",Sys.Date(),".xlsx"))
  Original_sDSS2 = Original_DSS2; Original_sDSS2[,DSS_cols] = denrow[,-1:-2]
  };


source("d3heatmap.R")
###########################################################################################
############ Bar plots

if(barplot_){
  
  ##########################
  #### PDF plots
  
  pdf(paste0("./www/Results/DSS/",PDFname), width = 12, height = 10,family="Helvetica-Narrow")
  
  for (i in 1:length(DSS_cols)){
    
    # top 50 DSS for visualization (in decr. order)
    toPlot = Original_DSS2[,c(grep("DRUG_NAME", colnames(Original_DSS2))[[1]],DSS_cols[i])]; toPlot$DSS = as.numeric(toPlot$DSS);
    toPlot <- toPlot[complete.cases(toPlot),]; toPlot = toPlot[order(toPlot[,2], decreasing = T)[1:50],]; 
    x = names(toPlot)[1]; y = names(toPlot)[2]; title = screens_[i];
    
    p <- ggplot(toPlot, aes_string(x = names(toPlot)[1], y = names(toPlot)[2])) + geom_bar(stat="identity", fill = "red", width=.7) + 
      scale_x_discrete(limits = toPlot[[x]]) + ggtitle(title[[1]]) + theme(axis.ticks.length=unit(0.1, "cm"),  axis.ticks.x=element_blank(),
                                                                           axis.text.x = element_text(angle=90, margin=unit(c(-0.5,0,0,0), "cm")), panel.background = element_blank(),
                                                                           axis.text.y = element_text(margin=unit(c(0,0,0,0), "cm")), plot.title = element_text(size = 14, face="bold", hjust = 0.5)) + xlab("") + ylab("DSS")
    
    print(p)
  }
  dev.off()
}

warning("HI")
##########################
#### HTML plots

if(HTMLreport)
{
  
  Htmlplots = paste0(sapply(1:length(DSS_cols), function(i){
    print(i)
    # top 50 DSS for visualization (in decr. order)
    toPlot = Original_DSS2[,c(1:4, c(5+11*(i-1)):c((4+11*i)))]; #toPlot = Original_DSS2[,c(1, -3+i*10 ,2,3,  10*i-5, 10*i-4, 10*i+3)];
    TEC50 = ifelse(any(grepl("TC50", colnames(toPlot))), "TC50", "EC50"); toPlot$DSS = as.numeric(toPlot$DSS)
    # toPlot = Original_DSS2[,c("ID", "DRUG_NAME", ,3,  10*i-5, 10*i-4, 10*i+3)]; toPlot$DSS = as.numeric(toPlot$DSS)
    toPlot = toPlot[order(toPlot[,"DSS"], decreasing = T)[1:50],]; #toPlot <- toPlot[complete.cases(toPlot),]
    
    p <- plot_ly(x = toPlot$DRUG_NAME, y = toPlot$DSS, type = "bar", hoverinfo = "x+text", color = I("red"),
                 text = paste0("Name: <b>",toPlot$DRUG_NAME,"</b><br>DSS: <b>",toPlot$DSS,"</b><br>AUC: <b>",toPlot$AUC,"</b><br>IC50/EC50: <b>",toPlot[[TEC50]],"</b><br>ID: <b>",
                               toPlot$ID, '</b><br>Mechanism/Targets: <b>',toPlot$`Mechanism/Targets`,'</b><br>Class explained: <b>', toPlot$Class.explained, '</b>')) %>% 
      layout(title = as.character(screens_[i]), xaxis = list(title = "", categoryarray = y, categoryorder = "array", tickfont = list(size = 10.5), 
                                                             ticks = "outside",tickangle = -45), yaxis = list(title = as.character("DSS"), tickfont = list(size = 10.5), ticks = "outside"),
             margin = list(b = 100, l = 40), height = 650)
    p$x$data$hoverinfo = "x+text"; p$x$data$inherit = !1; p$x$data$marker = NULL; p$x$data$marker$color = "#FC8D62"; p$x$data$name = "red";
    p$x$data$text = p$x$attrs[[1]]$text; p$x$data$type="bar"; p$x$data$x=p$x$attrs[[1]]$x; p$x$data$y=p$x$attrs[[1]]$y; p$x$layout$height = 650;
    p$x$layout$margin = NULL; p$x$layout$margin$b= 100; p$x$layout$margin$l= 40; p$x$layout$margin$t= 25; p$x$layout$margin$r= 10;
    p$x$layout$title=NULL; p$x$layout$title=p$x$layoutAttrs[[1]]$title;
    
    obj <- htmlwidgets:::toJSON(as.widget(p))[1]
    # obj <- gsubfn::gsubfn("[[:digit:]]+\\.+[[:digit:]]+[[:digit:]]", ~format(round(as.numeric(x), 2), nsmall = 2), obj)
    
    GraphsEChtml = do.call(paste0, lapply(1:nrow(toPlot), function(j){
      paste0('<div class = "FigSt" id="', gsub("\\s|[[:punct:]]", "", toPlot$DRUG_NAME[j]),'">', toPlot$GRAPH[j], '</div>')}))
    
    # make container name unique by extracting current proc.time (used to put json inside of it)
    widContN <- paste0(gsubfn("[[:digit:]]", ~format(letters[as.numeric(x) + 1]), sub("\\.", "", proc.time()[3])), letters[i])
    
    # (regulate size throug .screen) wrapContN <- paste0('w', widContN)
    widContainerWide <- paste0("<div id=\"screen", as.character(i), "\" class=\"screen\">")
    widContainer <- paste0("<div id=\"", widContN, "\" >", "</div>")
    
    barchart <- paste0("<script type=\"application/json\" data-for=\"", widContN, "\">", obj, "</script>")
    renderJs <- paste0("<script> var Data = document.querySelector(\"script[data-for='", widContN, "']\"); var x = JSON.parse(Data.textContent || Data.text).x; x.data = [x.data]; x.config = undefined; x.base_url = undefined;x.height = undefined;x.source = undefined;x.url = undefined;x.width = undefined; Plotly.plot('", widContN, "', x.data, x.layout, x.config);", "var myPlot = document.getElementById('",widContN,"');myPlot.on('plotly_hover',function(data){$(\"#screen",as.character(i)," > #\"+data.points[0].x.replace(/[^a-zA-Z0-9]/g,'')).fadeTo('fast', 1);}); myPlot.on('plotly_unhover', function(data){$(\"#screen",as.character(i)," > #\"+data.points[0].x.replace(/[^a-zA-Z0-9]/g,'')).fadeTo('fast', 0);});" ,"</script>")
    
    paste0(widContainerWide, widContainer, GraphsEChtml, barchart, renderJs, "</div>")
  }))
  
  header <- base::readChar(headerwater, file.info(headerwater)$size)
  Htmlplots = paste0(header, paste0(Htmlplots, collapse=""),'</body></html>')
  writeChar(Htmlplots, paste0("./www/Results/DSS/",HTMLname), nchar(Htmlplots, type = "chars"))
}

warning("HI2")
warning(length(DSS_cols))
warning("HI3")
###########################################################################################
############ Heat maps DSS
#browser();browser();browser();browser();

if(heatmap_)
{
  
  i_ = 1:length(DSS_cols) # number of screens
  
  if(length(i_) > 1) # we need at least two screens for clustering.
  {
    # remove strings with all NAs
    Original_DSS2 <- Original_DSS2[apply(Original_DSS2[,grep("Experiment_id", colnames(Original_DSS2))], 1, function(x) !all(is.na(x))),]
    
    Experiment_id_cols = grep("Experiment_id", colnames(Original_DSS2));    
    GRAPH_cols = grep("GRAPH", colnames(Original_DSS2));
    NAME_col = grep("DRUG_NAME", colnames(Original_DSS2))[[1]];
    TEC50_cols = grep("EC50$|TC50$", colnames(Original_DSS2));
    Mechanism_cols = grep("Mechanism/Targets$", colnames(Original_DSS2));
    ID_cols = grep("ID$", colnames(Original_DSS2));
    Class_cols = grep("Class.explained$", colnames(Original_DSS2));
    Std_error= grep("TEC50_std_error$", colnames(Original_DSS2));           
    AUC_Value= grep("AUC$", colnames(Original_DSS2));            
    
    for(i in i_) Original_DSS2[,GRAPH_cols[i]] = gsub("\"","\'", Original_DSS2[,GRAPH_cols[i]])
    
    hoverdata =  as.data.frame(sapply(i_, function(i)
      
      as.data.frame(apply(Original_DSS2[,c(Experiment_id_cols[i],NAME_col,DSS_cols[i],AUC_Value[i],TEC50_cols[i],ID_cols,Std_error[i],Mechanism_cols,Class_cols,GRAPH_cols[i])], 1, function(x) paste0('Experiment id: <b>',x[1],'</b><br>Name: <b>',x[2], '</b><br>DSS: <b>', 
                                                                                                                                                                                                       x[3], '</b><br>AUC: <b>',x[4], '</b><br>IC50/EC50: <b>',x[5], '</b><br>ID: <b>',x[6], '</b><br>Standard Error: <b>',x[7], '</b><br>Mechanism/Targets: <b>',x[8],
                                                                                                                                                                                                       '</b><br>Class explained: <b>',x[9], '</b><br>', x[10])), stringsAsFactors = F)), stringsAsFactors = F)
    
    custom_pallete <- c("#FFFFFF","#FF4B4B","#CC0000")
    #colrs = colorRampPalette(RColorBrewer::brewer.pal(9,"Reds"))(2000)[1:1700]
    #colrs =  colorRampPalette(c("white","red"))
    colrs = colorRampPalette(custom_pallete)(2000)[1:1700]
    toPlot <- apply(Original_DSS2[,DSS_cols[i_]],2,as.numeric)
    
    
    # height should also be changed in header and tailer .txt.
    if(dmethod_ == "none")  Rowv = Colv = !1 else Rowv = Colv = !0
    
    p <- d3heatmap(toPlot, Rowv = Rowv,Colv = Colv, scale = "none", dendrogram = "both",color = colrs,
                   distfun = function(x) mydist(x, dmethod= ifelse(exists("dmethod_"), dmethod_, "euclidean"), usecor = "pairwise.complete.obs"),
                   hclustfun = function(x) hclust(x, method="ward.D2"),  width = 1280, yaxis_font_size = 10,
                   labRow = Original_DSS2[,NAME_col], labCol = as.character(screens_),
                   k_col = length(i_), cellnote = hoverdata)
    
    
    # the drug names order has changed due to clustering. match the order of mechanismofactions with new order of drugs.
    data_ = Original_DSS2[,c("DRUG_NAME","Mechanism/Targets")];
    data_ = data_[match(as.data.frame(p$x$matrix$rows, stringsAsFactors = F)[[1]], data_$DRUG_NAME),]
    
    #recursively add spaces to repeated strings, while we get rid off all duplicates 
    while(any(duplicated(data_[,2])))
      data_[min(which(duplicated(data_[,2]) == !0)),2] = paste0(data_[min(which(duplicated(data_[,2]) == !0)),2], " ")
    
    #add mechanism of action to heatmap. (i.e. rows2 in JSON.) rows - drugs, cols - screens, rows2 - mechanisms.
    p$x$matrix$rows2 = data_[,2];
    p$x$options$yaxis_width = 400; p$x$options$xaxis_height = 200;
    
    #add colors for colorscale. (replace 10 with desired number of colors/breaks in color legend)
    colbr = 15;
    p$x$matrix$colrs = colrs[seq.int(1, length(colrs), round(length(colrs)/colbr))]
    
    #add labels for min... middle... and max color.
    p$x$matrix$colrslabs = round(seq.int(0, max(toPlot, na.rm = T),length.out = round(colbr)))
    p$x$matrix$colrsper = sapply(1:colbr, function(i) {
      if(i!=colbr) sum(toPlot > p$x$matrix$colrslabs[i] & toPlot < p$x$matrix$colrslabs[i+1],  na.rm = T)
      else sum(toPlot > p$x$matrix$colrslabs[i],  na.rm = T)})
    p$x$matrix$colrsper = round(p$x$matrix$colrsper/(nrow(toPlot)*ncol(toPlot)),3)
    
    header <- base::readChar(headerpath, file.info(headerpath)$size)
    # header <- gsub(":6700,",":3000,",header)
    tailer <- base::readChar(tailerpath, file.info(tailerpath)$size)
    # tailer <- gsub(":1700,",":900,",header)
    # 
    obj <- rjson::toJSON(p)
    obj <- paste0(substr(obj, 1, regexpr("\"width\":1280,", obj)[1]), "evals\":[],\"jsHooks\":[]}")
    obj <- paste0(header,obj,tailer)
    
    #optimize JSON by rounding each double to 2 dec. 
    obj <- gsubfn::gsubfn("[[:digit:]]+\\.+[[:digit:]]+[[:digit:]]", ~format(round(as.numeric(x), 3), nsmall = 3), obj)
    obj <- stringr::str_trim(obj)
    
    
    
    obj <- gsub("sizeforpdfprint",paste0(nrow(Original_DSS2)*4.1,"mm"),obj)
    
    # obj <- gsub("6700",nrow(Original_DSS2)*20,obj)
    # if(nrow(Original_DSS2)<40){obj <- gsub("6700",nrow(Original_DSS2)*50,obj)}else{obj <- gsub("6700",nrow(Original_DSS2)*20,obj)}
    if(nrow(Original_DSS2)<16 && nrow(Original_DSS2)>11){obj <- gsub("6700",nrow(Original_DSS2)*61,obj)} else if(nrow(Original_DSS2)<9){obj <- gsub("6700",nrow(Original_DSS2)*190,obj)} else if(nrow(Original_DSS2)>=9 && nrow(Original_DSS2)<=11){obj <- gsub("6700",nrow(Original_DSS2)*70,obj)} 
    else if(nrow(Original_DSS2)>=16 && nrow(Original_DSS2)<25){obj <- gsub("6700",nrow(Original_DSS2)*45,obj)} else if(nrow(Original_DSS2)>=25 && nrow(Original_DSS2)<39){obj <- gsub("6700",nrow(Original_DSS2)*32,obj)} else{obj <- gsub("6700",nrow(Original_DSS2)*20,obj)}
    # obj <- gsub(":6700,",":2700,",obj)
    writeChar(obj, paste0("./www/Results/DSS/",HEATname), nchar(obj, type = "chars") + 50)
    
    
    
    
    
    
    
    
    
    
    toPlot <- apply(Original_DSS2[,AUC_cols[i_]],2,as.numeric)
    
    p <- d3heatmap(toPlot, Rowv = Rowv,Colv = Colv, scale = "none", dendrogram = "both",color = colrs,
                   distfun = function(x) mydist(x, dmethod= ifelse(exists("dmethod_"), dmethod_, "euclidean"), usecor = "pairwise.complete.obs"),
                   hclustfun = function(x) hclust(x, method="ward.D2"),  width = 1280, yaxis_font_size = 10,
                   labRow = Original_DSS2[,NAME_col], labCol = as.character(screens_),
                   k_col = length(i_), cellnote = hoverdata)
    
    
    # the drug names order has changed due to clustering. match the order of mechanismofactions with new order of drugs.
    data_ = Original_DSS2[,c("DRUG_NAME","Mechanism/Targets")];
    data_ = data_[match(as.data.frame(p$x$matrix$rows, stringsAsFactors = F)[[1]], data_$DRUG_NAME),]
    
    #recursively add spaces to repeated strings, while we get rid off all duplicates 
    while(any(duplicated(data_[,2])))
      data_[min(which(duplicated(data_[,2]) == !0)),2] = paste0(data_[min(which(duplicated(data_[,2]) == !0)),2], " ")
    
    #add mechanism of action to heatmap. (i.e. rows2 in JSON.) rows - drugs, cols - screens, rows2 - mechanisms.
    p$x$matrix$rows2 = data_[,2];
    p$x$options$yaxis_width = 400; p$x$options$xaxis_height = 200;
    
    #add colors for colorscale. (replace 10 with desired number of colors/breaks in color legend)
    colbr = 15;
    p$x$matrix$colrs = colrs[seq.int(1, length(colrs), round(length(colrs)/colbr))]
    
    #add labels for min... middle... and max color.
    p$x$matrix$colrslabs = round(seq.int(0, max(toPlot, na.rm = T),length.out = round(colbr)))
    p$x$matrix$colrsper = sapply(1:colbr, function(i) {
      if(i!=colbr) sum(toPlot > p$x$matrix$colrslabs[i] & toPlot < p$x$matrix$colrslabs[i+1],  na.rm = T)
      else sum(toPlot > p$x$matrix$colrslabs[i],  na.rm = T)})
    p$x$matrix$colrsper = round(p$x$matrix$colrsper/(nrow(toPlot)*ncol(toPlot)),3)
    
    header <- base::readChar(headerpath2, file.info(headerpath2)$size)
    # header <- gsub(":6700,",":3000,",header)
    tailer <- base::readChar(tailerpath, file.info(tailerpath)$size)
    # tailer <- gsub(":1700,",":900,",header)
    # 
    obj <- rjson::toJSON(p)
    obj <- paste0(substr(obj, 1, regexpr("\"width\":1280,", obj)[1]), "evals\":[],\"jsHooks\":[]}")
    obj <- paste0(header,obj,tailer)
    
    #optimize JSON by rounding each double to 2 dec. 
    obj <- gsubfn::gsubfn("[[:digit:]]+\\.+[[:digit:]]+[[:digit:]]", ~format(round(as.numeric(x), 3), nsmall = 3), obj)
    obj <- stringr::str_trim(obj)
    
    
    
    obj <- gsub("sizeforpdfprint",paste0(nrow(Original_DSS2)*4.1,"mm"),obj)
    
    # obj <- gsub("6700",nrow(Original_DSS2)*20,obj)
    # if(nrow(Original_DSS2)<40){obj <- gsub("6700",nrow(Original_DSS2)*50,obj)}else{obj <- gsub("6700",nrow(Original_DSS2)*20,obj)}
    if(nrow(Original_DSS2)<16 && nrow(Original_DSS2)>11){obj <- gsub("6700",nrow(Original_DSS2)*61,obj)} else if(nrow(Original_DSS2)<9){obj <- gsub("6700",nrow(Original_DSS2)*190,obj)} else if(nrow(Original_DSS2)>=9 && nrow(Original_DSS2)<=11){obj <- gsub("6700",nrow(Original_DSS2)*70,obj)} 
    else if(nrow(Original_DSS2)>=16 && nrow(Original_DSS2)<25){obj <- gsub("6700",nrow(Original_DSS2)*45,obj)} else if(nrow(Original_DSS2)>=25 && nrow(Original_DSS2)<39){obj <- gsub("6700",nrow(Original_DSS2)*32,obj)} else{obj <- gsub("6700",nrow(Original_DSS2)*20,obj)}
    # obj <- gsub(":6700,",":2700,",obj)
    writeChar(obj, paste0("./www/Results/DSS/",HEATauc), nchar(obj, type = "chars") + 50)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # TREE
    # exclude rows with 0
    toPlot <- apply(Original_DSS2[,DSS_cols[i_]],2,as.numeric); ZEROrows = as.numeric(which(rowSums(toPlot, na.rm = T) == 0)) 
    if(length(ZEROrows) > 0){Original_DSS2 = Original_DSS2[-ZEROrows,]; toPlot = toPlot[-ZEROrows,];}  
    Original_DSS2$Class.explained[is.na(Original_DSS2$Class.explained)] = "Unknown"
    
    p2 <- d3heatmap::d3heatmap(toPlot, Rowv = !0, scale = "none", dendrogram = "both",color = colrs, k_row = 4,
                               distfun = function(x) mydist(x, dmethod= ifelse(exists("dmethod_"), dmethod_, "euclidean"), usecor = "pairwise.complete.obs"),
                               hclustfun = function(x) hclust(x, method="ward.D2"),  width = 1280, yaxis_font_size = 10,
                               labRow = paste0(Original_DSS2[,NAME_col], '||', Original_DSS2$Class.explained, '||',as.numeric(rowMeans(Original_DSS2[,DSS_cols], na.rm = !0))), labCol = as.character(screens_),
                               k_col = length(i_), cellnote = hoverdata)
    
    objTree = rjson::toJSON(p2$x$rows); drugClasses = unique(Original_DSS2$Class.explained); 
    treepath = "tree_.html"; tree_ <- base::readChar(treepath, file.info(treepath)$size)
    tree_ <- sub("replaceJSONhere",objTree,tree_); tree_ <- sub("replaceClasseshere",paste0('"', paste0(drugClasses, collapse = '","'), '"'),tree_);
    writeChar(tree_, paste0("./www/Results/DSS/",Treename), nchar(tree_, type = "chars") + 50)
    
    
    
    
    if(sDSS){
      # remove strings with all NAs
      
      Original_sDSS2 <- Original_sDSS2[apply(Original_sDSS2[,grep("Experiment_id", colnames(Original_sDSS2))], 1, function(x) !all(is.na(x))),]
      
      hoverdata =  do.call("cbind", lapply(i_, function(i)
        sapply(1:nrow(Original_sDSS2), function(x){
          paste0('Experiment id: <b>',Original_sDSS2[x,Experiment_id_cols[i]],'</b><br>Name: <b>',Original_sDSS2[x,NAME_col], '</b><br>sDSS: <b>',Original_sDSS2[x,DSS_cols[i]],'</b><br>AUC: <b>',Original_sDSS2[x,AUC_Value[i]], '</b><br>IC50/EC50: <b>',Original_sDSS2[x,TEC50_cols[i]], 
                 '</b><br>ID: <b>',Original_sDSS2[x,ID_cols], '</b><br>Mechanism/Targets: <b>',Original_sDSS2[x,Mechanism_cols], '</b><br>Class explained: <b>',Original_sDSS2[x,Class_cols],
                 '</b><br>DSS: <b>',Original_DSS2[Original_DSS2$DRUG_NAME == Original_sDSS2[x,]$DRUG_NAME,DSS_cols[i]])})
      ))
      
      
        
      toPlot <- apply(Original_sDSS2[,DSS_cols[i_]],2,as.numeric); toPlotCp = toPlot;
      colrs = colorRampPalette(c("#0000CC","#FFFFFF","#CC0000"))(2000)
      
      
      toPlot[which(!is.na(toPlot) & (toPlot < 0))] = scales::rescale(toPlot[which(!is.na(toPlot) & (toPlot < 0))], c(-1, 0))
      toPlot[which(!is.na(toPlot) & (toPlot > 0))] = scales::rescale(toPlot[which(!is.na(toPlot) & (toPlot > 0))], c(0, 1))
      
      if(dmethod_ == "none")  Rowv = Colv = !1 else Rowv = Colv = !0
      
      p <- d3heatmap(toPlot, Rowv = Rowv,Colv = Colv, scale = "none", dendrogram = "both",color = colrs,
                     distfun = function(x) mydist(x, dmethod= ifelse(exists("dmethod_"), dmethod_, "euclidean"), usecor = "pairwise.complete.obs"),
                     hclustfun = function(x) hclust(x, method="ward.D2"),  width = 1280, yaxis_font_size = 10,
                     labRow = Original_sDSS2[,NAME_col], labCol = as.character(screens_),
                     k_col = length(i_), cellnote = hoverdata)
      
      
      # the drug names order has changed due to clustering. match the order of mechanismofactions with new order of drugs.
      data_ = Original_sDSS2[,c("DRUG_NAME","Mechanism/Targets")];
      data_ = data_[match(as.data.frame(p$x$matrix$rows, stringsAsFactors = F)[[1]], data_$DRUG_NAME),]
      
      #recursively add spaces to repeated strings, while we get rid off all duplicates 
      while(any(duplicated(data_[,2])))
        data_[min(which(duplicated(data_[,2]) == !0)),2] = paste0(data_[min(which(duplicated(data_[,2]) == !0)),2], " ")
      
      #add mechanism of action to heatmap. (i.e. rows2 in JSON.) rows - drugs, cols - screens, rows2 - mechanisms.
      p$x$matrix$rows2 = data_[,2];
      p$x$options$yaxis_width = 400; p$x$options$xaxis_height = 200;
      
      #add colors for colorscale. (replace 10 with desired number of colors/breaks in color legend)
      colbr = 7;
      # p$x$matrix$colrs = colrs[seq.int(1, length(colrs), round(length(colrs)/colbr))]
      sumsum2=sum(!is.na(toPlot) & (toPlot < 0))
      sumsum3=sum(!is.na(toPlot) & (toPlot > 0))
      answer2ans=(sumsum2*colbr)/(sumsum2+sumsum3)-1
      answer3ans=(sumsum3*colbr)/(sumsum2+sumsum3)+1
      
      
      
      
      
      colors1=colrs[1:1000]
      colors2=colrs[1001:2000]
      
      devided2=(length(colors1))
      # sup=round(devided2/answer2ans)
      # sup2=round(devided2/answer3ans)
      # sup3=sum(sup,sup2)
      
      seq11=seq.int(1, (length(colors1)), round(devided2/answer2ans))
      seq21=seq.int(1, (length(colors1)), round(devided2/answer3ans))
      
      # seq41 = c(seq11, seq31)
      
      p$x$matrix$colrs1 = colors1[seq11]
      p$x$matrix$colrs2 = colors2[seq21]
      
      p$x$matrix$colrs = c( p$x$matrix$colrs1, p$x$matrix$colrs2)
      
      # p$x$matrix$colrs = colrs[seq.int(1, (length(colrs))/2, round(devided2/answer2ans))]
      # p$x$matrix$colrs = colrs[seq.int(1, (length(colrs)), sup3)]
      #add labels for min... middle... and max color.
      toPlot = toPlotCp;
      
      
      p$x$matrix$colrslabs = round(seq.int(min(toPlot, na.rm = T), max(toPlot, na.rm = T),length.out = round(colbr)))
      p$x$matrix$colrsper = sapply(1:colbr, function(i) {
        if(i!=colbr) sum(toPlot > p$x$matrix$colrslabs[i] & toPlot < p$x$matrix$colrslabs[i+1],  na.rm = T)
        else sum(toPlot > p$x$matrix$colrslabs[i],  na.rm = T)})
      p$x$matrix$colrsper = round(p$x$matrix$colrsper/(nrow(toPlot)*ncol(toPlot)),3)
      
      header <- base::readChar(headerpath, file.info(headerpath)$size)
      tailer <- base::readChar(tailerpath, file.info(tailerpath)$size)
      
      obj <- rjson::toJSON(p)
      obj <- paste0(substr(obj, 1, regexpr("\"width\":1280,", obj)[1]), "evals\":[],\"jsHooks\":[]}")
      obj <- paste0(header,obj,tailer)
      
      #optimize JSON by rounding each double to 2 dec. 
      obj <- gsubfn::gsubfn("[[:digit:]]+\\.+[[:digit:]]+[[:digit:]]", ~format(round(as.numeric(x), 3), nsmall = 3), obj)
      obj <- stringr::str_trim(obj)
      obj <- gsub("sizeforpdfprint",paste0(nrow(Original_sDSS2)*4.1,"mm"),obj)
      if(nrow(Original_sDSS2)<16 && nrow(Original_sDSS2)>11){obj <- gsub("6700",nrow(Original_sDSS2)*61,obj)} else if(nrow(Original_sDSS2)<9){obj <- gsub("6700",nrow(Original_sDSS2)*190,obj)} else if(nrow(Original_sDSS2)>=9 && nrow(Original_sDSS2)<=11){obj <- gsub("6700",nrow(Original_sDSS2)*70,obj)} 
      else if(nrow(Original_sDSS2)>=16 && nrow(Original_sDSS2)<25){obj <- gsub("6700",nrow(Original_sDSS2)*45,obj)} else if(nrow(Original_sDSS2)>=25 && nrow(Original_sDSS2)<39){obj <- gsub("6700",nrow(Original_sDSS2)*32,obj)} else{obj <- gsub("6700",nrow(Original_sDSS2)*20,obj)}
      
      
      
      writeChar(obj, paste0("./www/Results/DSS/",HEATname2), nchar(obj, type = "chars") + 50)
      
    }
    
  }
}




####
# JSON for final report, index.html
if(sDSS){listToJSON <- c(numberofscreens = numberofscreens, doneDSS = T,sDSS=T, PDFname = PDFname, HTMLname = HTMLname, HEATname = HEATname,HEATname2 = HEATname2,Treename=Treename, XLSXname = list.files(path = "./www/Results/DSS", pattern="Original_DSS2_[0-9]{4}-[0-9]{2}-[0-9]{2}[.]xlsx$"), XLSXname2 = list.files(path = "./www/Results/DSS", pattern = "Selective_DSS2_[0-9]{4}-[0-9]{2}-[0-9]{2}[.]xlsx$")[[1]], AUCname = list.files(path = "./www/Results/DSS", pattern = "Original_AUC_[0-9]{4}-[0-9]{2}-[0-9]{2}[.]xlsx$")[[1]])}
if(!sDSS){listToJSON <- c(numberofscreens = numberofscreens, doneDSS = T,sDSS=F, PDFname = PDFname, HTMLname = HTMLname, HEATname = HEATname,Treename=Treename, XLSXname = list.files(path = "./www/Results/DSS", pattern="Original_DSS2_[0-9]{4}-[0-9]{2}-[0-9]{2}[.]xlsx$")[[1]] , AUCname = list.files(path = "./www/Results/DSS", pattern = "Original_AUC_[0-9]{4}-[0-9]{2}-[0-9]{2}[.]xlsx$")[[1]],HEATauc = HEATauc)}

jsonL <- rjson::toJSON(listToJSON)
#setwd(directory_name);
write(jsonL, "./www/Results/HTMLreport/DSS.json")

rm(final_tbl_conc_drugs); rm(screen_table); rm(drug_annot_tbl_full)


