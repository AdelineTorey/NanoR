############################################ NanoPrepareG ############################################



#' @title Prepares sequencing summary and .fastq files
#' @description NanoPrepareG prepares MinION and GridION X5 sequencing summary and .fastq files for other functions from this package. Name of the function is inherited from NanoR previous version as MinION and GridION X5 had different default-output formats and only GridION X5 outputted sequencing summary and .fastq files. From MinION release 18.12 and GridION 18.12.1 outputs are the same.
#' @param DataSummary Path to sequencing summary file/s folder
#' @param DataFastq Path to passed .fastq folder. 
#' @param Cores Number of cores to accelerate sequencing summary files reading when dealing with multiple ones. Default to 1
#' @param Label Label to identify the experiment. A folder with this name will be created in DataOut directory.
#' @details NanoPreareG can find desired input files recursively, so be careful to specify path to folder contaning a unique data type. Old releases of GridION X5 stored all the .fastq files (passed and failed) in the same directory. In this case, this directory can be given as "DataFastq" input.
#' @return Object of class list
#' @examples
#' #do not run
#' DataSummary<-'path/to/sequencing_summary'
#' DataFastq<-'path/to/fastq_pass'
#' Label<-'Exp'
#' #new behaviour
#' List<-NanoPrepareG(DataSummary, DataFastq, Label=Label)
#' #old behaviour (same folder for DataSummary and DataFastq)
#' Data<-'path/to/sequencing_summary'<-'path/to/fastq_pass'
#' List<-NanoPrepareG(Data, Data, Cores=5, Label=Label)


NanoPrepareG <- function(DataSummary, DataFastq, Label) {

  FastqFiles <- list.files(DataFastq, pattern = ".fastq", full.names = TRUE, recursive = TRUE)
  SummariesFiles <- list.files(DataSummary, full.names = TRUE, pattern = "sequencing_summary", recursive = TRUE)
  label <- as.character(Label)         
  
  ## print message with number of files in folder
  message(length(FastqFiles), " passed .fastq files in folder")
  message(length(SummariesFiles), " sequencing summary files in folder")
  message("Reading and organizing sequencing summary ...")
  

    Read_Table_Summary <- function(File){

    Table <- read.table(File, header = FALSE, sep="\t", skip=1)

    RelativeTimeToAdd <- (as.numeric(Table[,11])+as.numeric(Table[,13]))## calculate a relative time that will be rescaled
    Read_Id <- as.character(Table[,3])
    Channel <- as.numeric(Table[,5])
    Mux <- as.numeric(Table[,6])    
    Length <- as.numeric(Table[,14])
    Qscore <- as.numeric(Table[,15])
    Table <- cbind(Read_Id, Channel, Mux, RelativeTimeToAdd, Length, Qscore)
    return(Table)
  }

  SummaryTable <- Read_Table_Summary(SummariesFiles) #if only one table, assume is the new format
  colnames(SummaryTable) <- c("Read Id", "Channel Number", "Mux Number", "Relative Time", "Length of Read", "Quality")
  List <- list(FastqFiles, SummaryTable, label)

  message("Done")
  return(List)
}