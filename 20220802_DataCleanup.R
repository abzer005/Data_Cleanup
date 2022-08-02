install.packages("dplyr")
install.packages("ggplot2")
install.packages("IRdisplay")

library(dplyr)
library(ggplot2)
library(IRdisplay)

# setting the current directory as the working directory
Directory <- normalizePath(readline("Enter the path of the folder with input files: "),"/",mustWork=FALSE)
setwd(Directory)

dirs <- dir(path=paste(getwd(), sep=""), full.names=TRUE, recursive=TRUE)
folders <- unique(dirname(dirs))
files <- list.files(folders, full.names=TRUE)
files_1 <- basename((files))
files_2 <- dirname((files))
# Creating a Result folder
dir.create(path=paste(files_2[[1]], "_Results", sep=""), showWarnings = TRUE)
fName <-paste(files_2[[1]], "_Results", sep="")

print(files_1)

#--------------------------------------------------------------------------------
  #Gets the extension of each file. Ex:csv
  pattern <- c()
  for (i in files_1){
    sep_file <- substr(i, nchar(i)-2,nchar(i))
    pattern <- rbind(pattern,sep_file)
  }
  
input <- as.double(unlist(strsplit(readline("Specify the file index of gapfilled feature-file, metadata separated by commas:"), split=",")))
ft <- read.csv(files_1[input[1]],sep = ifelse(pattern[input[1]]!="csv","\t",","), header=TRUE,check.names = FALSE) # By applying 'row.names = 1', the 1st column 'ID' becomes the row names
md <-read.csv(files_1[input[2]], sep = ifelse(pattern[input[2]]!="csv","\t",","), header=TRUE,check.names = FALSE)

  #If you have non gapfilled feature file:
  if(readline("Do you have non gap-filled feature table? Y/N:")=="Y"){
    x <- as.double(readline("Enter the ID number of non-gap-filled feature file:"))
    nft<- read.csv(files_1[x],sep=ifelse(pattern[x]!="csv","\t",","), header = TRUE,check.names = FALSE)
}
#--------------------------------------------------------------------------------
  ##Removing Peak area extensions
colnames(ft) <- gsub(' Peak area','',colnames(ft))
md$filename<- gsub(' Peak area','',md$filename)

#Removing if any NA columns present in the md file
md <- md[,colSums(is.na(md))<nrow(md)]

#Changing the row names of the files
rownames(md) <- md$filename
md <- md[,-1]
rownames(ft) <- paste(ft$'row ID',round(ft$'row m/z',digits = 3),round(ft$'row retention time',digits = 3), sep = '_')

#Picking only the files with column names containing 'mzML'
ft <- ft[,grep('mzML',colnames(ft))]

# if nft file exists, we perform all the above for nft as well
if(exists("nft")==T){
  colnames(nft) <- gsub(' Peak area','',colnames(nft))
  rownames(nft) <- paste(nft$'row ID',round(nft$'row m/z',digits = 3),round(nft$'row retention time',digits = 3), sep = '_')
  nft <- nft[,grep('mzML',colnames(nft))]
}

#--------------------------------------------------------------------------------------------------------------------------------------------------
new_ft<- ft[,order(colnames(ft))] #ordering the ft by its column names
new_md <-md[order(rownames(md)),] #ordering the md by its row names

#lists the colnames(ft) which are not present in md
unmatched_ft <- colnames(new_ft)[which(is.na(match(colnames(new_ft),rownames(new_md))))] 
cat("These", length(unmatched_ft),"columns of feature table are not present in metadata:")
if((length(unmatched_ft) %% 2) ==0)
{matrix(data=unmatched_ft,nrow=length(unmatched_ft)/2,ncol=2)}else
{matrix(data=unmatched_ft,nrow=(length(unmatched_ft)+1)/2,ncol=2)}

flush.console()
Sys.sleep(0.2)

#lists the rownames of md which are not present in ft
unmatched_md <- rownames(new_md)[which(is.na(match(rownames(new_md),colnames(new_ft))))] 
cat("These", length(unmatched_md),"rows of metadata are not present in feature table:")
if((length(unmatched_md) %% 2) ==0)
{matrix(data=unmatched_md,nrow=length(unmatched_md)/2,ncol=2)}else
{matrix(data=unmatched_md,nrow=(length(unmatched_md)+1)/2,ncol=2)}

#Removing those unmatching columns and rows:
if(length(unmatched_ft)!=0){new_ft <- subset(ft, select = -c(which(is.na(match(colnames(ft),rownames(md))))) )}
if(length(unmatched_md)!=0){new_md <- md[-c(which(is.na(match(rownames(md),colnames(ft))))),]}

#checking the dimensions of our new ft and md:
cat("The number of rows and columns in our original ft is:",dim(ft),"\n")
cat("The number of rows and columns in our new ft is:",dim(new_ft),"\n")
cat("The number of rows and columns in our new md is:",dim(new_md))

new_ft<- new_ft[,order(colnames(new_ft))] #ordering the ft by its column names
new_md <-new_md[order(rownames(new_md)),] #ordering the md by its row names

#checking if they are the same
if(identical(colnames(new_ft),rownames(new_md))==T)
{print("The column names of ft and rownames of md are the same")}else{print("The column names of ft and rownames of md are not the same")}

#-Subsetting the metdata by user-defined conditions:----------------------------------------------------------------------------------------------------
  #Storing the ft and md in different names
  Meta_Filter <- md
  input_data <- ft
  
  # Visualising the different levels in each attribute of the metadata
  lev <- c()
  for(i in 1:ncol(Meta_Filter)){
    x <- levels(as.factor(Meta_Filter[,i]))
    if(is.double(Meta_Filter[,i])==T){x=round(as.double(x),2)}
    x <-toString(x)
    lev <- rbind(lev,x)
  }
  Info_mat <- data.frame(INDEX=c(1:ncol(Meta_Filter)),ATTRIBUTES=colnames(Meta_Filter),LEVELS=lev,row.names=NULL)
  Info_mat
  
  #Splitting batch-wise:
  Batch_info <- ifelse(readline("Do you want to split your data batch-wise? Y/N: ")=="Y",as.double(readline("Please enter the index of the attributes having batch information:")),"No")
  if(is.numeric(Batch_info)==TRUE){
    Levels_Batch <- levels(as.factor(Meta_Filter[,Batch_info]))
    IRdisplay::display(data.frame(INDEX=c(1:length(Levels_Batch)),LEVELS=Levels_Batch))
    
    flush.console()  
    Sys.sleep(0.2)
    
    Cdtn <-as.double(unlist(strsplit(readline("Enter the index number of condition(s) you want to KEEP (separated by commas):"), split=',')))
    
    #Getting all the rows of metadata that satisfies each element of the condition and storing it as an element in Split_data list
    Split_data <-list()
    for (j in 1:length(Cdtn)){
      Split_data[[j]] <- Meta_Filter[(Meta_Filter[,Batch_info]==Levels_Batch[Cdtn[j]]),]
    }
    
    Batch_data <-do.call(rbind, Split_data) # unlists the Split data and combines them by row
    flush.console()  
    Sys.sleep(0.2)
    
    IRdisplay::display(head(Batch_data)) #Visualising the Batch_data
    dim(Batch_data)
  }
 
  #-----------------------------------
  #If batch data exists, it will take it as "data", else take Meta_filter as "data"
  if(exists("Batch_data")==T){data <-Batch_data}else{data <-Meta_Filter}
  
  Info_mat
  Condition <- as.double(unlist(readline("Enter the index number of the attribute to split sample and blank:")))
  
  Levels_Cdtn <- levels(as.factor(data[,Condition[1]]))
  print(matrix(Levels_Cdtn,length(Levels_Cdtn)))
  
  #Among the shown levels of an attribute, select the ones to keep
  Blk_id <- as.double(unlist(readline("Enter the index number of your BLANK:")))
  paste0('You chosen blank is:',Levels_Cdtn[Blk_id])
  
  #Splitting the data into control and samples based on the metadata
  md_Blank <- data[(data[,Condition] == Levels_Cdtn[Blk_id]),]
  Blank <- input_data[,which(colnames(input_data)%in%rownames(md_Blank))] 
  md_Samples <- data[(data[,Condition] != Levels_Cdtn[Blk_id]),]
  Samples <- input_data[,which(colnames(input_data)%in%rownames(md_Samples))] 

#-Function:Frequency plot-------------------------------------------------------------------------------
  FrequencyPlot <- function(x1,x2){
    
    bins <- c(-1,0,(1 * 10^(seq(0,10,1)))) #creating bins from -1 to 10^10 using sequence function seq()
    scores_x1 <- cut(as.matrix(x1),bins,labels = c('0','1','10','1E2','1E3','1E4','1E5','1E6','1E7','1E8','1E9','1E10')) #cut function cuts the give table into its appropriate bins
    
    #transform function convert the tables into a column format: easy for visualization 
    Table_x1<-transform(table(scores_x1)) #contains 2 columns: "scores_x1", "Freq"
    
    #Repeating the same steps for x2
    scores_x2 <- cut(as.matrix(x2),bins,labels = c('0','1','10','1E2','1E3','1E4','1E5','1E6','1E7','1E8','1E9','1E10'))
    Table_x2<-transform(table(scores_x2))
    
    #Getting the names of x1 and x2
    arg1 <- deparse(substitute(x1))
    arg2 <-deparse(substitute(x2))
    
    #Creating a data frame for plotting
    data_plot <- as.data.frame(c(Table_x1$Freq,Table_x2$Freq)) #Concatenating the frequency info of both tables rowwise
    colnames(data_plot) <- "Freq" #naming the 1st column as 'Freq'
    data_plot$Condition <- c(rep(arg1,12),rep(arg2,12)) #adding a 2nd column 'Condition', which just repeats the name of x1 and x2 accordingly
    data_plot$Range_bins <- rep(Table_x1$scores_x1,2) #Adding 3rd column 'Range Bins'
    data_plot$Log_Freq <- log(data_plot$Freq+1) #Log scaling the frequency values
    
    ## GGPLOT2
    BarPlot <- ggplot(data_plot, aes(Range_bins, Log_Freq, fill = Condition)) + 
      geom_bar(stat="identity", position = "dodge", width=0.4) + 
      scale_fill_brewer(palette = "Set1") +
      ggtitle(label="Frequency plot") +
      xlab("Range") + ylab("(Log)Frequency") + labs(fill = "Data Type") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
      theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
      theme(plot.title = element_text(hjust = 0.5)) # centering the plot title
    
    print(BarPlot)
  }  
  
#---Blank Removal:-----------------------------------------------------------------------------
  if(readline('Do you want to perform Blank Removal- Y/N:')=='Y'){
    
    #When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
    Cutoff <- as.numeric(readline('Enter Cutoff value between 0.1 & 1:')) # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3
    
    #Getting mean for every feature in blank and Samples
    Avg_blank <- rowMeans(Blank, na.rm= FALSE, dims = 1) # set na.rm = FALSE to check if there are NA values. When set as TRUE, NA values are changed to 0
    Avg_samples <- rowMeans(Samples, na.rm= FALSE, dims = 1)
    
    #Getting the ratio of blank vs Sample
    Ratio_blank_Sample <- (Avg_blank+1)/(Avg_samples+1)
    
    # Creating a bin with 1s when the ratio>Cutoff, else put 0s
    Bg_bin <- ifelse(Ratio_blank_Sample > Cutoff, 1, 0 )
    Blank_removal <- cbind(Samples,Bg_bin)
    
    # Checking if there are any NA values present. Having NA values in the 4 variables will affect the final dataset to be created
    temp_NA_Count <-cbind(Avg_blank ,Avg_samples,Ratio_blank_Sample,Bg_bin)
    
    print('No of NA values in the following columns:')
    print(colSums(is.na(temp_NA_Count)))
    
    #Calculating the number of background features and features present
    print(paste("No.of Background or noise features:",sum(Bg_bin ==1,na.rm = TRUE)))
    print(paste("No.of features after excluding noise:",(nrow(Samples) - sum(Bg_bin ==1,na.rm = TRUE)))) 
    
    Blank_removal <- Blank_removal %>% filter(Bg_bin == 0) # Taking only the feature signals
    Blank_removal <- as.matrix(Blank_removal[,-ncol(Blank_removal)]) # removing the last column Bg_bin 
  }
#--------------------------------------------------------------------------------
  GapFilled <-Blank_removal
  if(exists("nft")==T){NotGapFilled <-nft}
  
  if(readline('Do you want to perform Imputation with minimum value of NonGapFilled table? - Y/N:')=='Y'){
    plot<- FrequencyPlot(GapFilled,NotGapFilled)
    
    Arg1 = plot$data$Condition[1]
    Arg2 = plot$data$Condition[13]
    
    # getting the 2nd minimum value of non-gap filled data. (The first minimum value in the data table is usually zero)
    RawLOD <- round(min(NotGapFilled[NotGapFilled!=min(NotGapFilled)]))
    
    print(paste0("The minimum value greater than 0 for ",Arg1,":", round(min(GapFilled[GapFilled!=min(GapFilled)]))))
    print(paste0("The minimum value greater than 0 for ",Arg2,":", RawLOD))
    
    Imputed <- GapFilled
    Imputed[Imputed<RawLOD] <- RawLOD # Replacing values<RawLOD with RawLOD
    #head(Imputed) 
  }
  
  write.csv(Imputed, file.path(fName,paste0('Imputed_QuantTable_with_MinValue_',RawLOD,'_Ecoli.csv')),row.names =TRUE) 

#Imputation with Cutoff LOD:-------------------------------------------------------------------------------
  FrequencyPlot(GapFilled,GapFilled)
  if(readline('Do you want to perform Imputation with a Cutoff LOD? - Y/N:')=='Y'){
    Cutoff_LOD <-as.numeric(readline("Enter your Cutoff LOD as seen in the plot:"))  #Enter the LOD value as seen in the frequency plot
    Imputed <- GapFilled
    Imputed[Imputed <Cutoff_LOD] <- Cutoff_LOD
    dim(Imputed)
  }
  
  Imputed<-Imputed[rowMeans(Imputed)!= RawLOD,]  
  write.csv(Imputed,file.path(fName, paste0('Imputed_QuantTable_filled_with_',Cutoff_LOD,'_CutOff_Used_',Cutoff,'_Bsub','.csv')),row.names =TRUE)

#Normalisation:--------------------------------------------------------------------------------
  if (readline("Do you want to perform Normalization: Y/N:") == 'Y'){
    #Getting column-wise sums of the input-data
    sample_sum <- colSums(Imputed, na.rm= TRUE, dims = 1)
    
    #Dividing each element of a particular column with its column sum
    Normalized_data <- c()
    for (i in 1:ncol(Imputed)){
      x <- Imputed[,i] / sample_sum[i]
      Normalized_data <- cbind(Normalized_data, x)
    }
    colnames(Normalized_data) <- names(sample_sum)
    
    dim(Normalized_data)
  } else return(head(Imputed))
  
  
  print(paste('No.of NA values in Normalized data:',sum(is.na(Normalized_data)== TRUE)))
  
  write.csv(Normalized_data, file.path(fName,'Normalised_Quant_table_Ecoli.csv'),row.names =TRUE) 
#END---------------------------------------------------------------------------------------------------