
############################################ NanoStatsG ############################################


#' @title Plots statistics 
#' @description NanoStatsG plots statistics parsing the metadata table returned by NanoTableG
#' @param NanoGList Object of class list returned by NanoPrepareG
#' @param NanoGTable Metadata table returned by NanoTableG
#' @param DataOut Where results will be saved. Use the same directory specified as "DataOut" for NanoTableG
#' @param KeepGGObj Store data.frames for ggplot plots in folder. Useful for personalize plot colors. Default to FALSE.
#' @return Plots: \cr 
#' - Yield.pdf (accumulation of reads and bps); \cr 
#' - RBLQ.pdf (# reads, # bps, length and quality overview every 30 minutes of experiment); \cr 
#' - LvQ.pdf (length and quality compared jointly); \cr 
#' - PFGC.pdf (passed and failed reads, GC content if previously computed); \cr 
#' - Activity.pdf (channels activity (# bps). Inactive channels are grey-colored) \cr
#' Tables: \cr
#' - metadata.fltrd.txt (metadata table for high-quality passed sequence)
#' - ShortSummary.txt (table with major statistics for the experiment)
#' @examples
#' DataOut <- "/path/to/DataOut"
#' # Need a list previously generated with NanoPreparGM()
#' # Need a table previously generated with NanoTableG()
#' # If List from NanoPrepareG() and Table from NanoTableG() exist: 
#' # Do not save ggplot2 tables:
#' NanoStatsG(List,Table, DataOut=DataOut)
#' # Save ggplot2 tables:
#' NanoStatsG(List,Table, DataOut=DataOut,KeepGGObj=TRUE)

## should be in description file
library(reshape2)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)

NanoStatsG <- function(NanoGList, NanoGTable, DataOut, KeepGGObj = FALSE) {
  
  Label <- NanoGList[[3]]    
  
  Dirs <- file.path(DataOut, Label)
  Directory <- file.path(DataOut, Label)  #here is the metadata.txt
  ## seems to appear twice - Dirs and Directory are identical
  
  TableInDirectory <- list.files(Directory,pattern="metadata.txt")
    if(length(TableInDirectory) == 0) {
    stop("Use the same directory specified for NanoTableG function")
  }
  
  CumulativeInDirectory <- list.files(Directory,pattern=("Yield.pdf"))
    if(length(CumulativeInDirectory) != 0) {
    stop("Can't use a directory that already contains previous results")  
  }
  ## functions ##
  
  #recreate the same plot of NanoStatsM, but requires different parsing of metadata table
  Increment <- function(x)
  {
    return(x+8)
  }
  
  Increment2 <- function(x)
  {
    return(x+32)
  }
  
  Increment3<-function(x)
  {
    return(x+128)
  }
  
  rotate<-function(x) 
  { 
    t(apply(x, 1, rev))
  }
  
  define_region <- function(row,col)
  {
    viewport(layout.pos.row = row, layout.pos.col = col)
  }

  ################
  options(scipen=9999)
  #setwd(Directory)

          
  Really_Pass_File <- which(as.numeric(NanoGTable[,6]) >= 7)
  List.Files.HDF5_Fail_Length <- length(which(as.numeric(NanoGTable[,6]) < 7))
  List.Files.HDF5_Pass.length <- nrow(NanoGTable)-List.Files.HDF5_Fail_Length
  
  NanoTable2 <- NanoGTable[Really_Pass_File,]
  
  write.table(NanoTable2, file.path(Directory, 'metadata.fltrd.txt'), col.names=T, row.names=F, quote=F, sep="\t")
   
  Table_HDF5_Def <- NanoTable2[,1:6]   
  Time_2 <- as.numeric(Table_HDF5_Def[,4])
  Run_Duration <- round(as.numeric(difftime(as.POSIXct(max(Time_2),origin="1/1/1970"), as.POSIXct(min(Time_2), origin="1/1/1970"), units="hours")))
  Relative_Time <- scales::rescale(Time_2, to=c(0,Run_Duration))  
  #Table_HDF5_Def<-cbind(Table_HDF5,Time_Rescaled)
  #colnames(Table_HDF5_Def)<-c("Flowcell ID","Read Id","Channel Number","Relative Time","Length of Read","Quality","Relative Experimental Time")
  
  
  #Relative_Time<-as.numeric(Table_HDF5_Def[,7])
  Relative_Time_Per_Hours<-seq(from=min(round(Relative_Time)), to=max(round(Relative_Time)), by=0.5)
  Template_Length<-as.numeric(Table_HDF5_Def[,5])
  Quality_Score<-as.numeric(Table_HDF5_Def[,6])
    
  message("Analyzing...")
  
  lenRTPH <- length(Relative_Time_Per_Hours)
  Reads_Per_Hour <- numeric(lenRTPH)
  Base_Pairs_Per_Hour <- numeric(lenRTPH)
  Max_Length_Per_Hour <- numeric(lenRTPH)
  Mean_Length_Per_Hour <- numeric(lenRTPH)
  Min_Length_Per_Hour <- numeric(lenRTPH)
  Min_Quality_Score_Per_Hour <- numeric(lenRTPH)
  Mean_Quality_Score_Per_Hour <- numeric(lenRTPH)
  Max_Quality_Score_Per_Hour <- numeric(lenRTPH)
  
  
  for (ii in seq_len(lenRTPH)) {
    
    if (ii < lenRTPH) {
      Index_Hours<-which(Relative_Time >= Relative_Time_Per_Hours[ii] & Relative_Time < Relative_Time_Per_Hours[ii+1])
      if (length(Index_Hours) == 0) {
        next
      }
      else
      {
        Reads_Per_Hour[ii]<-length(Index_Hours)
        Base_Pairs_Per_Hour[ii]<-sum(Template_Length[Index_Hours])
        Mean_Length_Per_Hour[ii]<-mean(Template_Length[Index_Hours])
        Max_Length_Per_Hour[ii]<-max(Template_Length[Index_Hours])
        Min_Length_Per_Hour[ii]<-min(Template_Length[Index_Hours])
        Mean_Quality_Score_Per_Hour[ii]<-mean(Quality_Score[Index_Hours])
        Min_Quality_Score_Per_Hour[ii]<-min(Quality_Score[Index_Hours])
        Max_Quality_Score_Per_Hour[ii]<-max(Quality_Score[Index_Hours])
      }
    }
    else {
      
      Index_Hours<-which(Relative_Time == Relative_Time_Per_Hours[ii])
      
      if (length(Index_Hours) == 0) {
        next
      }
      else {
        Reads_Per_Hour[ii]<-length(Index_Hours)
        Base_Pairs_Per_Hour[ii]<-sum(Template_Length[Index_Hours])
        Mean_Length_Per_Hour[ii]<-Mean_Length_Per_Hour[ii-1]
        Max_Length_Per_Hour[ii]<-Max_Length_Per_Hour[ii-1]
        Min_Length_Per_Hour[ii]<-Min_Length_Per_Hour[ii-1]
        Mean_Quality_Score_Per_Hour[ii]<-Mean_Quality_Score_Per_Hour[ii-1]
        Min_Quality_Score_Per_Hour[ii]<-Min_Quality_Score_Per_Hour[ii-1]
        Max_Quality_Score_Per_Hour[ii]<-Max_Quality_Score_Per_Hour[ii-1]
      }
    }
  }

  
  
  Cumulative_Reads<-cumsum(Reads_Per_Hour)
  Cumulative_Basepairs<-cumsum(Base_Pairs_Per_Hour)

  if (is.na(Table_HDF5_Def[1,3])) { #skip muxes reordering

    Channel_Vector<-as.numeric(Table_HDF5_Def[,2])
    #Mux_Vector<-as.numeric(Table_HDF5_Def[,3])
    Channels_Number<-c(1:512)
    
    
    Base_Pairs_Per_Channel <- NULL
    
    Table_HDF5_Re <- vector(mode = "list", length(Channels_Number))
    
    for (iii in seq_len(length(Channels_Number))) {
      Ind_Chann<-which(Channel_Vector == Channels_Number[iii])
      #Mux_Associated<-sort(Mux_Vector[Ind_Chann], index.return=TRUE)$ix
      if (length(Ind_Chann) == 0) {
        next
      }

      if (length(Ind_Chann) == 1) {

        Table_HDF5_Re[[iii]]<-Table_HDF5_Def[Ind_Chann,]

      }

      else {
        Table_HDF5_Re[[iii]]<-Table_HDF5_Def[Ind_Chann,]
        #Table_HDF5_Reordered<-rbind(Table_HDF5_Reordered,Table_HDF5_Re) slow for very large tables
      }
      Base_Pairs_Per_Channel[iii]<-sum(Template_Length[Ind_Chann])
    }

    Table_HDF5_Reordered<-do.call(rbind,Table_HDF5_Re) # a lot faster

  }

  else { #muxes are present, reorder

    Channel_Vector<-as.numeric(Table_HDF5_Def[,2])
    Mux_Vector<-as.numeric(Table_HDF5_Def[,3])
    Channels_Number<-c(1:512)
    
    
    Base_Pairs_Per_Channel <- NULL
    
    Table_HDF5_Re <- vector(mode = "list", length(Channels_Number))
    
    for (iii in seq_len(length(Channels_Number))) {
      Ind_Chann<-which(Channel_Vector == Channels_Number[iii])
      Mux_Associated<-sort(Mux_Vector[Ind_Chann], index.return=TRUE)$ix
      if (length(Ind_Chann) == 0) {
        next
      }

      if (length(Ind_Chann) == 1) {

        Table_HDF5_Re[[iii]]<-Table_HDF5_Def[Ind_Chann,]

      }

      else {

        Table_HDF5_Re[[iii]]<-Table_HDF5_Def[Ind_Chann,][Mux_Associated,]
        #Table_HDF5_Reordered<-rbind(Table_HDF5_Reordered,Table_HDF5_Re) slow for very large tables
      }
      Base_Pairs_Per_Channel[iii]<-sum(Template_Length[Ind_Chann])
    }

    Table_HDF5_Reordered<-do.call(rbind,Table_HDF5_Re) # a lot faster

  }
  
  #rownames(Table_HDF5_Reordered)<-c()
  
  
  #PLOT CUMULATIVE READS/BP


  message("Plotting...")

  x <- Relative_Time_Per_Hours
  y0.1 <- Cumulative_Reads
  data0.1 <- data.frame('x' = x,'y' = y0.1)
  data0.1$group <- "reads yield"
  
  
  Cumulative_Reads_Plot<-ggplot(data0.1, aes(x=x, y=y, col=group)) +
    geom_line(size=.5, colour = "dodgerblue4") + 
    scale_x_continuous(name="time(hrs)", breaks=(seq(0,Run_Duration,2)))+
    scale_y_continuous(name="# reads")+
    #geom_ribbon(data=subset(data0.1,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
    #scale_color_manual(name='', values=c("reads yield" = "dodgerblue4"))+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_text(size=11),axis.title.y = element_text(size=11))
    #theme(legend.position="bottom")
    #ggtitle("Cumulative Reads")
 


  y0.2 <- Cumulative_Basepairs
  data0.2 <- data.frame('x' = x,'y' = y0.2)
  data0.2$group <- "bps yield"
  
  Cumulative_Base_Pairs_Plot <- ggplot(data0.2, aes(x=x, y=y, col=group)) +
    geom_line(size=.5, colour = "firebrick4") + 
    scale_x_continuous(name="time(hrs)", breaks=(seq(0,Run_Duration,2)))+
    scale_y_continuous(name="# bps")+
    #geom_ribbon(data=subset(data0.2,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
    #scale_color_manual(name='', values=c("bps yield" = "firebrick4"))+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_text(size=11),axis.title.y = element_text(size=11))
    #theme(legend.position="bottom")
    #ggtitle("Cumulative Base Pairs")
  
  Cumulative_Plot <- grid.arrange(Cumulative_Reads_Plot, Cumulative_Base_Pairs_Plot, nrow=2, ncol=1, top = textGrob("Reads and base pairs yield",gp=gpar(fontsize=20)))
  
  ggsave(file.path(Directory, "Yield.pdf"), device="pdf", Cumulative_Plot, height=10,width=15)


  #PLOT PER-HOUR READS/BPs/QUALITY/LENGTH

  y1 <- Reads_Per_Hour
  data1 <- data.frame('x' = x, 'y' = y1)
  data1$group <- "# reads / 30 mins"
  
  Reads_Per_Hour_Plot<-ggplot(data1, aes(x=x, y=y, col=group)) +
    geom_line(size=.5, colour = "dodgerblue4") + 
    scale_x_continuous(name="time (hrs)", breaks=(seq(0,Run_Duration,2)))+
    scale_y_continuous(name="# reads")+
    #geom_ribbon(data=subset(data1,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
    #scale_color_manual(name='', values=c("# reads / 30 mins" = "dodgerblue4"))+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_text(size=11),axis.title.y = element_text(size=11))
    #theme(legend.position="bottom")
    #ggtitle("Reads")
  
  
  y2 <- Base_Pairs_Per_Hour
  data2 <- data.frame('x' = x, 'y' = y2)
  data2$group <- "# bps / 30 mins"
  
  
  Base_Pairs_Per_Hour_Plot <- ggplot(data2, aes(x=x, y=y, col=group)) +
    geom_line(size=.5, colour = "firebrick4") + 
    scale_x_continuous(name="time (hrs)", breaks=(seq(0,Run_Duration,2)))+
    scale_y_continuous(name="# bps")+
    #geom_ribbon(data=subset(data2,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
    #scale_color_manual(name='', values=c("# bps / 30 mins" = "firebrick4"))+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_text(size=11),axis.title.y = element_text(size=11))
    #theme(legend.position="bottom")
    #ggtitle("Base Pairs")
  
  
  y3.0<-log10(Mean_Length_Per_Hour)
  data3.0.0<-data.frame('x'=x,'y'=Mean_Length_Per_Hour)
  data3.0<-data.frame('x'=x,'y'=y3.0)
  data3.0$group<-"avg length"
  
  y3.1<-log10(Max_Length_Per_Hour)
  data3.1<-data.frame('x'=x,'y'=y3.1)
  data3.1$group<-"max length"
  
  y3.2<-log10(Min_Length_Per_Hour)
  data3.2<-data.frame('x'=x,'y'=y3.2)
  data3.2$group<-"min length"
  
  data3<-rbind(data3.0, data3.1, data3.2)
  
  Length_Per_Hour_Plot<-ggplot(data3, aes(x=x, y=y, col=group)) +
    geom_line(size=.5) + 
    scale_x_continuous(name="time (hrs)", breaks=(seq(0,Run_Duration,2)))+
    scale_y_continuous(name="length (bps)", breaks=c(2,3,4,5,6), labels=c("100","1000","10000","100000","1000000"))+
    #geom_ribbon(data=subset(data3,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0, alpha=.7) +
    scale_color_manual(name='', values=c("avg length" = "forestgreen", "min length" = "green", "max length" = "darkolivegreen"))+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_text(size=11),axis.title.y = element_text(size=11))+
    theme(legend.position="bottom")
    #ggtitle("Length")
  
  y4.0<-Mean_Quality_Score_Per_Hour
  data4.0<-data.frame('x'=x,'y'=y4.0)
  data4.0$group<-"avg quality"
  
  y4.1<-Min_Quality_Score_Per_Hour
  data4.1<-data.frame('x'=x,'y'=y4.1)
  data4.1$group<-"min quality"
  
  y4.2<-Max_Quality_Score_Per_Hour
  data4.2<-data.frame('x'=x,'y'=y4.2)
  data4.2$group<-"max quality"
  
  data4<-rbind(data4.0,data4.1,data4.2)
  
  Quality_Score_Per_Hour_Plot<-ggplot(data4, aes(x=x, y=y, col=group)) +
    geom_line(size=.5) + 
    scale_x_continuous(name="time (hrs)", breaks=(seq(0,Run_Duration,2)))+
    scale_y_continuous(name="quality (phred)")+
    #geom_ribbon(data=subset(data4,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,alpha=.7) +
    scale_color_manual(name='', values=c("avg quality" = "darkorange", "min quality" = "gold", "max quality" = "orangered"))+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_text(size=11),axis.title.y = element_text(size=11))+
    theme(legend.position="bottom")
    #ggtitle("Quality")
  
  
  Others_Plot<-grid.arrange(Reads_Per_Hour_Plot,Base_Pairs_Per_Hour_Plot,Length_Per_Hour_Plot,Quality_Score_Per_Hour_Plot, nrow=2, ncol=2, top = textGrob("Reads number, base pairs number, reads length and quality calculated every 30 minutes", gp=gpar(fontsize=20)))
  
  ggsave(file.path(Directory, "RBLQ.pdf"), device="pdf", Others_Plot, height=10, width=15)
  
  #PASS/FAIL
  
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold"),
      axis.text.x=element_blank()
    )
  
  Data_Pass_Fail_Percentage <- data.frame(
    group = c("% passed", "% failed"),
    value = c((List.Files.HDF5_Pass.length/(List.Files.HDF5_Pass.length+List.Files.HDF5_Fail_Length)), ((List.Files.HDF5_Fail_Length)/(List.Files.HDF5_Pass.length+List.Files.HDF5_Fail_Length)))
  )
  
  
  Data_Pass_Fail_Percentage_Plot<-ggplot(Data_Pass_Fail_Percentage, aes(x="", y=value, fill=group))+
    geom_bar(width = 1, size = 1, color = "white", stat = "identity")+
    geom_text(aes(label = scales::percent(value)), position = position_stack(vjust = 0.5))+
    coord_polar("y", start=0) +
    scale_fill_brewer(palette="Dark2") +
    blank_theme+
    theme(legend.title=element_blank(),legend.position="bottom")
    #ggtitle("% Passed and Failed/Skipped ")
  
  
  Data_Pass_Fail_Tot <- data.frame(
    group = c("# passed", "# failed"),
    value = c(List.Files.HDF5_Pass.length, List.Files.HDF5_Fail_Length)
  )
  
  Data_Pass_Fail_Tot_Plot<-ggplot(Data_Pass_Fail_Tot, aes(x="", y=value, fill=group))+
    geom_bar(width = 1, size = 1, color = "white", stat = "identity")+
    geom_text(aes(label = value), position = position_stack(vjust = 0.5))+
    coord_polar("y", start=0) +
    scale_fill_brewer(palette="Accent") +
    blank_theme+
    theme(legend.title=element_blank(), legend.position="bottom")
    #ggtitle("# Passed and Failed/Skipped")
  
  #WILL BE SAVED WITH GC CONTENT
  
  
  ###################### LENGTH_VS_QUALITY ######################
    

    
  limit <- 500000 #limit number of point to plot. Too slow otherwise

  if (length(Template_Length) <= limit) {

    Tot <- data.frame(cbind(Template_Length,Quality_Score))
    colnames(Tot) <- c("Template_Length","Quality_Score")

  }

  else { 

    message("Too many points for LvsQ scatterplot: rescaling, but mantaining proportions ...")
    tmp <- cbind(Template_Length,Quality_Score)
    Sample <- data.frame(tmp[sample(nrow(tmp), limit),])
    indmin_l <- which.min(Template_Length) 
    indmax_l <- which.max(Template_Length)
    indmin_q <- which.min(Quality_Score)
    indmax_q <- which.max(Quality_Score)
    Bound <- tmp[c(indmin_l, indmax_l, indmin_q, indmax_q),]
    Tot<-data.frame(rbind(Sample,Bound))
    colnames(Tot) <- c("Template_Length","Quality_Score")  

  }

  ScatterTheme <- list(labs(x="length (bps)",y="quality (phred)"),theme_bw(), theme(legend.position=c(1,0),legend.justification=c(1,0), legend.background=element_blank(),legend.direction="horizontal", legend.title=element_text(face="bold.italic")))
    
    
  hist_top_mean_length<-ggplot(data.frame(Template_Length), aes(x=Template_Length))+theme_bw()+ geom_histogram(aes(y = ..count../1000),col="darkolivegreen", fill="forestgreen",boundary = min(Template_Length), bins=30)+labs(x="",y=expression("Count"["(10^3)"]))+scale_x_continuous(limits=c(min(Template_Length),max(Template_Length)))
  hist_right_mean_quality<-ggplot(data.frame(Quality_Score), aes(x=Quality_Score))+ theme_bw()+ geom_histogram(aes(y = ..count../1000),col="orangered", fill="darkorange",boundary = min(Quality_Score), bins=30)+ labs(x="",y=expression("Count"["(10^3)"]))+coord_flip()+scale_x_continuous(limits=c(min(Quality_Score),max(Quality_Score)))
  empty <- ggplot()+geom_point(aes(1,1), colour="white")+
    theme(axis.ticks=element_blank(), 
          panel.background=element_blank(), 
          axis.text.x=element_blank(), axis.text.y=element_blank(),           
          axis.title.x=element_blank(), axis.title.y=element_blank())
    
  scatter <- ggplot(Tot,aes(x=Tot$Template_Length, y=Tot$Quality_Score))+geom_point(col="grey27", size=.09, alpha=.4)+scale_x_continuous(limits=c(min(Tot$Template_Length),max(Tot$Template_Length)))+scale_y_continuous(limits=c(min(Tot$Quality_Score),max(Tot$Quality_Score)))+stat_density2d(aes(col=..level.., alpha=..level..)) + scale_color_continuous(low="darkblue",high="darkred") +geom_smooth(method=lm,linetype=2, size=.5,col="black",se=F) + guides(alpha="none",col=guide_legend(title="Density"))+ ScatterTheme
          
    
  Length_VS_Quality_Plot<-grid.arrange(hist_top_mean_length, empty, scatter, hist_right_mean_quality, ncol=2, nrow=2, widths=c(4,1), heights=c(1, 4))
    
    
  ggsave(file.path(Directory, "LvsQ.pdf"), device="pdf", Length_VS_Quality_Plot, height=10,width=15)
  
  
  
  
  if (!is.na(Table_HDF5_Def[1,3])) { #if muxes are present

    Mux_Numbers<-c(1:4)
    
    Chan<-as.numeric(Table_HDF5_Reordered[,2])
    Mu<-as.numeric(Table_HDF5_Reordered[,3])
    Le<-as.numeric(Table_HDF5_Reordered[,5])
    
    List_Of_Mux <- vector(mode = "list", length(Channels_Number))
    
    for (iii in seq_len(length(Channels_Number))) {
      Ind_Chann<-which(Chan == Channels_Number[iii])
      Mux_Associated_Number<-sort(Mu[Ind_Chann])
      Table_Mux <- NULL
      for (lll in 1:length(Mux_Numbers)) {
        Ind_Mux<-which(Mux_Associated_Number == Mux_Numbers[lll])
        Chan_Mux<-Chan[Ind_Chann][Ind_Mux]
        if (length(Chan_Mux) == 0) {
          Mux<-NA
          Lenght_Per_Mux<-NA
          Table_Mu<-cbind(NA, NA, NA)
          Table_Mux<-rbind(Table_Mux,Table_Mu)
        }
        else {
          Mux<-Mu[Ind_Chann][Ind_Mux]
          Lenght_Per_Mux<-sum(Le[Ind_Chann][Ind_Mux])
          Table_Mu<-cbind(unique(Chan_Mux), unique(Mux),Lenght_Per_Mux)
          Table_Mux<-rbind(Table_Mux,Table_Mu)
        }
      }
      List_Of_Mux[[iii]]<-Table_Mux
    }
    
    Table_Mux_Def<-do.call(rbind,List_Of_Mux)

  }

  ### else skip muxes rapresentation. Table_Mux_Def does not exist in the environment
  
  
  #PLOT CORRELATION MATRIXES (CHANNEL (AND MUXES))
  
  m1<-matrix(Base_Pairs_Per_Channel[1:32], ncol=8, nrow=4, byrow=TRUE)
  m2<-matrix(Base_Pairs_Per_Channel[449:480], ncol=8, nrow=4, byrow=TRUE)
  m3<-matrix(Base_Pairs_Per_Channel[385:416], ncol=8, nrow=4, byrow=TRUE)
  m4<-matrix(Base_Pairs_Per_Channel[321:352], ncol=8, nrow=4, byrow=TRUE)
  m5<-matrix(Base_Pairs_Per_Channel[257:288], ncol=8, nrow=4, byrow=TRUE)
  m6<-matrix(Base_Pairs_Per_Channel[193:224], ncol=8, nrow=4, byrow=TRUE)
  m7<-matrix(Base_Pairs_Per_Channel[129:160], ncol=8, nrow=4, byrow=TRUE)
  m8<-matrix(Base_Pairs_Per_Channel[65:96], ncol=8, nrow=4, byrow=TRUE)
  mdef3<-rbind(m1,m2,m3,m4,m5,m6,m7,m8)
  m9<-rotate(matrix(Base_Pairs_Per_Channel[33:64], ncol=8, nrow=4, byrow=TRUE))
  m10<-rotate(matrix(Base_Pairs_Per_Channel[481:512], ncol=8, nrow=4, byrow=TRUE))
  m11<-rotate(matrix(Base_Pairs_Per_Channel[417:448], ncol=8, nrow=4, byrow=TRUE))
  m12<-rotate(matrix(Base_Pairs_Per_Channel[353:384], ncol=8, nrow=4, byrow=TRUE))
  m13<-rotate(matrix(Base_Pairs_Per_Channel[289:320], ncol=8, nrow=4, byrow=TRUE))
  m14<-rotate(matrix(Base_Pairs_Per_Channel[225:256], ncol=8, nrow=4, byrow=TRUE))
  m15<-rotate(matrix(Base_Pairs_Per_Channel[161:192], ncol=8, nrow=4, byrow=TRUE))
  m16<-rotate(matrix(Base_Pairs_Per_Channel[97:128], ncol=8, nrow=4, byrow=TRUE))
  mdef4<-rbind(m9,m10,m11,m12,m13,m14,m15,m16)
  Matrixbpchannel<-cbind(mdef3,mdef4)
  


  if (exists("Table_Mux_Def", envir=environment())) {


    BasePairs_Per_Mux<-as.numeric(Table_Mux_Def[,3])
    
    First_Eight_Disposition<-c(3,4,1,2,6,5,8,7)
    Second_Eight_Disposition<-Increment(First_Eight_Disposition)
    Third_Eight_Disposition<-Increment(Second_Eight_Disposition)
    Fouth_Eight_Disposition<-Increment(Third_Eight_Disposition)
    First_Line<-c(First_Eight_Disposition,Second_Eight_Disposition,Third_Eight_Disposition,Fouth_Eight_Disposition)
    Second_Line<-Increment2(First_Line)
    Third_Line<-Increment2(Second_Line)
    Fourth_Line<-Increment2(Third_Line)
    First_Block<-c(First_Line,Second_Line,Third_Line,Fourth_Line)
    Second_Block<-Increment3(First_Block)
    Third_Block<-Increment3(Second_Block)
    Fourth_Block<-Increment3(Third_Block)
    Fifth_Block<-Increment3(Fourth_Block)
    Sixth_Block<-Increment3(Fifth_Block)
    Seventh_Block<-Increment3(Sixth_Block)
    Eight_Block<-Increment3(Seventh_Block)
    Ninth_Block<-Increment3(Eight_Block)
    Tenth_Block<-Increment3(Ninth_Block)
    Eleventh_Block<-Increment3(Tenth_Block)
    Twelfth_Block<-Increment3(Eleventh_Block)
    Thirtheenth_Block<-Increment3(Twelfth_Block)
    Fourtheenth_Block<-Increment3(Thirtheenth_Block)
    Fiftheenth_Block<-Increment3(Fourtheenth_Block)
    Sixtheenth_Block<-Increment3(Fiftheenth_Block)
    
    M1<-matrix(BasePairs_Per_Mux[First_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M2<-matrix(BasePairs_Per_Mux[Fiftheenth_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M3<-matrix(BasePairs_Per_Mux[Thirtheenth_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M4<-matrix(BasePairs_Per_Mux[Eleventh_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M5<-matrix(BasePairs_Per_Mux[Ninth_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M6<-matrix(BasePairs_Per_Mux[Seventh_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M7<-matrix(BasePairs_Per_Mux[Fifth_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M8<-matrix(BasePairs_Per_Mux[Third_Block], ncol= 32, nrow= 4, byrow=TRUE)
    Mdef3<-rbind(M1,M2,M3,M4,M5,M6,M7,M8)
    M9<-rotate(matrix(BasePairs_Per_Mux[Second_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M10<-rotate(matrix(BasePairs_Per_Mux[Sixtheenth_Block],ncol= 32, nrow= 4, byrow=TRUE))
    M11<-rotate(matrix(BasePairs_Per_Mux[Fourtheenth_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M12<-rotate(matrix(BasePairs_Per_Mux[Twelfth_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M13<-rotate(matrix(BasePairs_Per_Mux[Tenth_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M14<-rotate(matrix(BasePairs_Per_Mux[Eight_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M15<-rotate(matrix(BasePairs_Per_Mux[Sixth_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M16<-rotate(matrix(BasePairs_Per_Mux[Fourth_Block], ncol= 32, nrow= 4, byrow=TRUE))
    Mdef4<-rbind(M9,M10,M11,M12,M13,M14,M15,M16)
    MatrixMuxActivity<-cbind(Mdef3,Mdef4)

  }


  #PLOTTING "FALSE" CORRELATION MATRIXES
  
  Palette <- colorRampPalette(brewer.pal(9, "Reds"))
  
  adjMatrixbpchannel<-melt(rotate(t(Matrixbpchannel)))
  Plot_Channel_Activity<-ggplot(data=adjMatrixbpchannel, aes(x=Var1, y=Var2)) +
    geom_tile(aes(fill=value), color="white", size=2)+
    scale_fill_gradientn(colours= Palette(4), na.value="grey70",limits=c(min(adjMatrixbpchannel[,3], na.rm=TRUE), max(adjMatrixbpchannel[,3], na.rm=TRUE))) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="bottom",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),legend.text=element_text(size=10),legend.text.align = 0)+
    guides(fill = guide_colorbar(barwidth= 61,barheight=.5,title="channels activity (# bps)", title.position="top",title.hjust=0.5))
    #ggtitle("Channels Activity")
  
  #ggsave("Activity.pdf", device="pdf",Plot_Channel_Activity, height=10, width=18)



  if (exists("MatrixMuxActivity", envir=environment())) {


    adjMatrixMuxActivity<-melt(rotate(t(MatrixMuxActivity)))
    Plot_Mux_Activity<-ggplot(data=adjMatrixMuxActivity, aes(x=Var1, y=Var2)) +
      geom_tile(aes(fill=value), color="white", size=2)+
      scale_fill_gradientn(colours=Palette(4), na.value="grey70",limits=c(min(adjMatrixMuxActivity[,3], na.rm=TRUE), max(adjMatrixMuxActivity[,3], na.rm=TRUE))) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="bottom",
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank(),legend.text=element_text(size=10))+
      guides(fill = guide_colorbar(barwidth= 61, barheight=.5,title="muxes activity (# bps)", title.position="top",title.hjust=0.5))
      #ggtitle("Muxes Activity")
    
    
    Plot_Tot<-grid.arrange(Plot_Channel_Activity,Plot_Mux_Activity,nrow=2, ncol=1, widths=15, heights=c(12,12))

    ggsave(file.path(Directory, "Activity.pdf"), device="pdf",Plot_Tot, height=10, width=18)


  }


  else {

    ggsave("Activity.pdf", device="pdf",Plot_Channel_Activity, height=10, width=18)

  }
  
  
  ###PLOT GC CONTENT#######
  
  
  
  if (NanoTable2[1,7] == "GC_Content") {
    
    pdf('PFGC.pdf', height=10, width=15, onefile=TRUE)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1, 2)))
    print(Data_Pass_Fail_Percentage_Plot, vp=define_region(1, 1))
    print(Data_Pass_Fail_Tot_Plot, vp = define_region(1, 2))
    dev.off()

  }
  
  else {
    
    GC_Content_To_Plot<-as.numeric(NanoGTable[,7][!is.na(NanoGTable[,7])])
    Hist_GC_Content<-ggplot(data.frame(GC_Content_To_Plot), aes(x=GC_Content_To_Plot))+theme_bw()+ geom_histogram(col="midnightblue", fill="cyan4", bins=30)+labs(x="",y="Count")+ggtitle("GC Content")+theme(plot.title = element_text(hjust = 0.5,size=16))
    pdf('PFGC.pdf', height=10, width=15, onefile=TRUE)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2)))
    print(Data_Pass_Fail_Percentage_Plot, vp=define_region(1, 1))
    print(Data_Pass_Fail_Tot_Plot, vp = define_region(1, 2))
    print(Hist_GC_Content, vp = define_region(2, 1:2))
    dev.off()
  }
  
  ###Get major stats###

  Longest <- max(Template_Length)
  Shortest <- min(Template_Length)
  MeanDim <- mean(Template_Length)
  MedianDim <- median(Template_Length)
  HQ <- max(Quality_Score)
  LQ <- min(Quality_Score)
  MeanQ <- mean(Quality_Score)
  MedianQ <- median(Quality_Score)
  Passed_num<-List.Files.HDF5_Pass.length
  Failed_num<-List.Files.HDF5_Fail_Length
  Passed_rat<-List.Files.HDF5_Pass.length/(List.Files.HDF5_Pass.length+List.Files.HDF5_Fail_Length)
  Failed_rat<-List.Files.HDF5_Fail_Length/(List.Files.HDF5_Pass.length+List.Files.HDF5_Fail_Length)


  dfval <- rbind(Label, Longest,Shortest,MeanDim,MedianDim,HQ, LQ, MeanQ, MedianQ, Cumulative_Basepairs[length(Cumulative_Basepairs)],Passed_num, Failed_num, Passed_rat, Failed_rat)
  rownames(dfval) <- c('Sample','Longest sequence (bps)', 'Shortest sequence (bps)', 'Mean length (bps)', 'Median length (bps)', 'Highest quality sequence (phred)', 'Lowest quality sequence (phred)', 'Mean quality (phred)', 'Median quality (phred)', 'Passed reads throughput (bps)', '# passed reads', '# failed/skipped reads', 'Ratio passed', 'Ratio failed/skipped')


  write.table(dfval,file.path(Directory, 'ShortSummary.txt'),col.names=F, sep="\t", row.names=T, quote=F)

  #save data for NanoCompare

  Comparison_Directory<-file.path(Directory, 'DataForComparison')
  dir.create(Comparison_Directory, showWarnings = FALSE)

  write.table(data1, file.path(Comparison_Directory, "Reads.txt"),col.names=T, sep="\t")
  write.table(data2, file.path(Comparison_Directory, "Bases.txt"),col.names=T, sep="\t")
  write.table(data3.0.0, file.path(Comparison_Directory, "Length.txt"),col.names=T, sep="\t")
  write.table(data4.0, file.path(Comparison_Directory, "Quality.txt"),col.names=T, sep="\t")
  
  if (KeepGGObj == TRUE) {

      GG_Directory<-file.path(Directory, 'GGTables')
      dir.create(GG_Directory, showWarnings = FALSE)
      write.table(data0.1,file.path(GG_Directory, 'Yield_reads.txt'),col.names=T, sep="\t", row.names=F, quote=F)
      write.table(data0.2,file.path(GG_Directory, 'Yield_bps.txt'),col.names=T, sep="\t", row.names=F, quote=F)
      write.table(data1,file.path(GG_Directory, 'RBLQ_reads.txt'),col.names=T, sep="\t", row.names=F, quote=F)
      write.table(data2,file.path(GG_Directory, 'RBLQ_bps.txt'),col.names=T, sep="\t", row.names=F, quote=F)
      write.table(data3,file.path(GG_Directory, 'RBLQ_length.txt'),col.names=T, sep="\t", row.names=F, quote=F)
      write.table(data4,file.path(GG_Directory, 'RBLQ_quality.txt'),col.names=T, sep="\t", row.names=F, quote=F)
      write.table(Data_Pass_Fail_Percentage,file.path(GG_Directory, 'PFGC_percentage.txt'),col.names=T, sep="\t", row.names=F, quote=F)
      write.table(Data_Pass_Fail_Tot,file.path(GG_Directory, 'PFGC_number.txt'),col.names=T, sep="\t", row.names=F, quote=F)
      write.table(Tot,file.path(GG_Directory, 'LvsQ_scatter.txt'),col.names=T, sep="\t", row.names=F, quote=F)
      write.table(adjMatrixbpchannel,file.path(GG_Directory, 'Activity_channels.txt'),col.names=T, sep="\t", row.names=F, quote=F)

      if (NanoTable2[1,7] != "GC_Content") {

        write.table(data.frame(GC_Content_To_Plot),file.path(GG_Directory, 'PFGC_GCC.txt'),col.names=T, sep="\t", row.names=F, quote=F)

      }


      if (exists("adjMatrixMuxActivity", envir=environment())) {

        write.table(adjMatrixbpchannel,file.path(GG_Directory, 'Activity_muxes.txt'),col.names=T, sep="\t", row.names=F, quote=F)


      }
  }

  message("Done")    
  
}

