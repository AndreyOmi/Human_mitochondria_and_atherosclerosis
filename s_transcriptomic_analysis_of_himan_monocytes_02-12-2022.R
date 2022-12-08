### Автор: Андрей Омельченко
### Дата:  31-08-2021
### Название: Скрипт для транскриптомного анализа моноцитов человека.
###


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
 #BiocManager::install("biomaRt")
 #BiocManager::install("DESeq2")

#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)

require(DOSE) ####### для отрисовки Dotplot

library(organism, character.only = TRUE)

library(stringi)
library(stringr)
library(reshape)
library(nortest)
library(ggplot2)
library(openxlsx)
library(tidyr) ##### From long to wide data format 
library(biomaRt)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(DOSE) ####### для отрисовки Dotplot
library(pathview)


rm(list = ls())

PathToResults <- "../Results"
PathToSources <- "../Sources"
PathToReferences <- "../References"
PathToDerivedData <- "../DerivedData"


PathToWorkingDirectory <- "C:\\Users\\Андрей\\OneDrive\\Documents\\Work\\НИИМЧ\\Транскриптомика атеросклероза\\Scripts"
#PathToWorkingDirectory <- "C:\\Users\\306-workstation\\Documents\\Andrey_Omelchenko\\RIHM_306_workstation\\Транскриптомика атеросклероза\\Scripts"

setwd(PathToWorkingDirectory)

PathToSource_R_Functions <- "C:\\Users\\Андрей\\OneDrive\\Documents\\Work\\R_Projects"

source(file.path(PathToSource_R_Functions, "f_R_functions_for_Science_and_Develoment_AO_13-01-2022.R"))

################# ПЕРВИЧНАЯ РАБОТА С ДАННЫМИ ###########################

##### открытие файла Excel с нормированными данными

FileName <- file.path(PathToSources, "all.RPKM.counts.new.xlsx")
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$all.2.RPKM.counts

InputData <- ListOfFileSheets$all.2.RPKM.counts

##### Подготовка данных для анализа

Transcripts_Table <- data.frame()

names(InputData)

##### Первая линия экспериментов

First_Sub_Table <- InputData[,c("gNam", "gLen", "71_7_Control", "72_7_Native_LDL", "73_7_Atherogenic_LDL", "74_7_Desialylated_LDL", "75_7_Acc_LDL", "76_7_Latex", "77_Ox_LDL")]
names(First_Sub_Table)
str(First_Sub_Table)

melt_First_Sub_Table <- melt(First_Sub_Table, id=c("gNam", "gLen"))
melt_First_Sub_Table$variable <- as.character(melt_First_Sub_Table$variable)
melt_First_Sub_Table$"Experiment_Line" <- 1

names(melt_First_Sub_Table)
str(melt_First_Sub_Table)

names(melt_First_Sub_Table) <- c("Transcript_name", "Transcript_length", "Impact_type", "RPKM_value", "Experiment_Line")

Transcripts_Table <- rbind.data.frame(Transcripts_Table, melt_First_Sub_Table)

##### Вторая линия экспериментов

Second_Sub_Table <- InputData[,c("gNam", "gLen", "81_8_Control", "82_8_Native_LDL", "83_8_Atherogenic_LDL", "84_8_Desialyated_LDL", "85_8_Acc_LDL", "86_Latex", "87_8_Ox_LDL")]
names(Second_Sub_Table)
str(Second_Sub_Table)

melt_Second_Sub_Table <- melt(Second_Sub_Table, id=c("gNam", "gLen"))
melt_Second_Sub_Table$variable <- as.character(melt_Second_Sub_Table$variable)
melt_Second_Sub_Table$"Experiment_Line" <- 2

names(melt_Second_Sub_Table)
str(melt_Second_Sub_Table)

names(melt_Second_Sub_Table) <- c("Transcript_name", "Transcript_length", "Impact_type", "RPKM_value", "Experiment_Line")

Transcripts_Table <- rbind.data.frame(Transcripts_Table, melt_Second_Sub_Table)

##### Третья линия экспериментов

Third_Sub_Table <- InputData[,c("gNam", "gLen", "91_9_Control", "92_9_Native_LDL", "93_9_Atherogenic_LDL", "96_9_Larex", "95_9_AccLDL", "94_9_DesialylatedLDL", "97_9_OxLDL")]
names(Third_Sub_Table)
str(Third_Sub_Table)

melt_Third_Sub_Table <- melt(Third_Sub_Table, id=c("gNam", "gLen"))
melt_Third_Sub_Table$variable <- as.character(melt_Third_Sub_Table$variable)
melt_Third_Sub_Table$"Experiment_Line" <- 3

names(melt_Third_Sub_Table)
str(melt_Third_Sub_Table)

names(melt_Third_Sub_Table) <- c("Transcript_name", "Transcript_length", "Impact_type", "RPKM_value", "Experiment_Line")

Transcripts_Table <- rbind.data.frame(Transcripts_Table, melt_Third_Sub_Table)

str(Transcripts_Table)

######### Удаление префиксов из вида модификации

unique(Transcripts_Table$Impact_type)

unique(str_remove(Transcripts_Table$Impact_type, "[0-9]{2}_[[0-9]_]*"))

Transcripts_Table$Impact_type <- str_remove(Transcripts_Table$Impact_type, "[0-9]{2}_[[0-9]_]*")

######## Исправление опечаток 

Transcripts_Table$Impact_type[which(Transcripts_Table$Impact_type == "Larex")] <- "Latex"
Transcripts_Table$Impact_type[which(Transcripts_Table$Impact_type == "Desialyated_LDL")] <- "Desialylated_LDL"
Transcripts_Table$Impact_type[which(Transcripts_Table$Impact_type == "DesialylatedLDL")] <- "Desialylated_LDL"
Transcripts_Table$Impact_type[which(Transcripts_Table$Impact_type == "AccLDL")] <- "Acc_LDL"
Transcripts_Table$Impact_type[which(Transcripts_Table$Impact_type == "OxLDL")] <- "Ox_LDL"

str(Transcripts_Table)
    
unique(Transcripts_Table$Impact_type)
unique(Transcripts_Table$Experiment_Line)

############ Сохранение исходных данных в файл .Rdata

Transcripts_Table_FileName <- file.path(PathToSources, paste("all.RPKM.counts.new.Rdata", sep = ''))
saveRDS(Transcripts_Table, file = Transcripts_Table_FileName)

list.files(PathToSources)

##### открытие файла Excel с сырыми данными

FileName <- file.path(PathToSources, "all.counts.new.xlsx")
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$all.2.counts

InputData <- InputData_Excel_Sheet

##### Подготовка данных для анализа

Transcripts_Table <- data.frame()

names(InputData)

##### Первая линия экспериментов

First_Sub_Table <- InputData[,c("gNam", "71_7_Control", "72_7_Native_LDL", "73_7_Atherogenic_LDL", "74_7_Desialylated_LDL", "75_7_Acc_LDL", "76_7_Latex", "77_Ox_LDL")]
names(First_Sub_Table)
str(First_Sub_Table)

melt_First_Sub_Table <- melt(First_Sub_Table, id=c("gNam"))
melt_First_Sub_Table$variable <- as.character(melt_First_Sub_Table$variable)
melt_First_Sub_Table$"Experiment_Line" <- 1

names(melt_First_Sub_Table)
str(melt_First_Sub_Table)

names(melt_First_Sub_Table) <- c("Transcript_name", "Impact_type", "Reads_count", "Experiment_Line")

Transcripts_Table <- rbind.data.frame(Transcripts_Table, melt_First_Sub_Table)

##### Вторая линия экспериментов

Second_Sub_Table <- InputData[,c("gNam", "81_8_Control", "82_8_Native_LDL", "83_8_Atherogenic_LDL", "84_8_Desialyated_LDL", "85_8_Acc_LDL", "86_Latex", "87_8_Ox_LDL")]
names(Second_Sub_Table)
str(Second_Sub_Table)

melt_Second_Sub_Table <- melt(Second_Sub_Table, id=c("gNam"))
melt_Second_Sub_Table$variable <- as.character(melt_Second_Sub_Table$variable)
melt_Second_Sub_Table$"Experiment_Line" <- 2

names(melt_Second_Sub_Table)
str(melt_Second_Sub_Table)

names(melt_Second_Sub_Table) <- c("Transcript_name", "Impact_type", "Reads_count", "Experiment_Line")

Transcripts_Table <- rbind.data.frame(Transcripts_Table, melt_Second_Sub_Table)

##### Третья линия экспериментов

Third_Sub_Table <- InputData[,c("gNam", "91_9_Control", "92_9_Native_LDL", "93_9_Atherogenic_LDL", "96_9_Larex", "95_9_AccLDL", "94_9_DesialylatedLDL", "97_9_OxLDL")]
names(Third_Sub_Table)
str(Third_Sub_Table)

melt_Third_Sub_Table <- melt(Third_Sub_Table, id=c("gNam"))
melt_Third_Sub_Table$variable <- as.character(melt_Third_Sub_Table$variable)
melt_Third_Sub_Table$"Experiment_Line" <- 3

names(melt_Third_Sub_Table)
str(melt_Third_Sub_Table)

names(melt_Third_Sub_Table) <- c("Transcript_name", "Impact_type", "Reads_count", "Experiment_Line")

Transcripts_Table <- rbind.data.frame(Transcripts_Table, melt_Third_Sub_Table)

str(Transcripts_Table)

######### Удаление префиксов из вида модификации

unique(Transcripts_Table$Impact_type)

unique(str_remove(Transcripts_Table$Impact_type, "[0-9]{2}_[[0-9]_]*"))

Transcripts_Table$Impact_type <- str_remove(Transcripts_Table$Impact_type, "[0-9]{2}_[[0-9]_]*")

######## Исправление опечаток 

Transcripts_Table$Impact_type[which(Transcripts_Table$Impact_type == "Larex")] <- "Latex"
Transcripts_Table$Impact_type[which(Transcripts_Table$Impact_type == "Desialyated_LDL")] <- "Desialylated_LDL"
Transcripts_Table$Impact_type[which(Transcripts_Table$Impact_type == "DesialylatedLDL")] <- "Desialylated_LDL"
Transcripts_Table$Impact_type[which(Transcripts_Table$Impact_type == "AccLDL")] <- "Acc_LDL"
Transcripts_Table$Impact_type[which(Transcripts_Table$Impact_type == "OxLDL")] <- "Ox_LDL"

str(Transcripts_Table)

unique(Transcripts_Table$Impact_type)
unique(Transcripts_Table$Experiment_Line)

############ Сохранение исходных данных в файл .Rdata

Transcripts_Table_FileName <- file.path(PathToSources, paste("all.counts.new.Rdata", sep = ''))
saveRDS(Transcripts_Table, file = Transcripts_Table_FileName)

list.files(PathToSources)

######## Открытие файла .txt

InputData <- read.table(file.path(PathToSources, paste("Nikita_all_counts.txt", sep = '')), header = TRUE)

head(InputData)

str(InputData)

InputData$geneNams <- as.character(InputData$geneNams)

names(InputData) <- c("Transcript_name", "Control", "Native_LDL", "Ox_LDL", "Acc_LDL", "HDL")

##### модификация вида таблицы по накоплению холестерина

InputData_long <- gather(data = InputData, experiment_type, reads_count, c("Control", "Native_LDL", "Ox_LDL", "Acc_LDL", "HDL"))

InputData_long$storage_type <- ifelse(InputData_long$experiment_type %in% c("Control", "Native_LDL", "HDL"), "no_storage", "storage")

head(InputData_long)

unique(InputData_long$storage_type)

unique(InputData_long$experiment_type)

###### Преобразование таблицы в широкую форму

InputData_wide <- spread(data = InputData_long, key = storage_type, value = reads_count)

############ Сохранение исходных данных в файл .Rdata

Transcripts_Table_FileName <- file.path(PathToSources, paste("Nikita_all_counts.Rdata", sep = ''))
saveRDS(InputData, file = Transcripts_Table_FileName)

list.files(PathToSources)

################# КОНЕЦ ПЕРВИЧНОЙ РАБОТЫ С ДАННЫМИ ###########################


################ Открытие файла с подготовленными исходными нормализованными данными

InputData_FileName <- file.path(PathToSources, "all.RPKM.counts.new.Rdata")
InputData <- readRDS(file = InputData_FileName)

names(InputData)
str(InputData)

unique(InputData$Impact_type)
unique(InputData$Experiment_Line)
unique(InputData$Impact_type)

################ Открытие файла с подготовленными исходными сырыми данными

InputData_FileName <- file.path(PathToSources, "all.counts.new.Rdata")
InputData <- readRDS(file = InputData_FileName)

names(InputData)
str(InputData)

unique(InputData$Impact_type)
unique(InputData$Experiment_Line)

head(InputData,  n = 100 )

################ Открытие файла с подготовленными данными первого эксперимента по секвенированию РНК

InputData_FileName <- file.path(PathToSources, "Nikita_all_counts.Rdata")
InputData <- readRDS(file = InputData_FileName)

names(InputData)
str(InputData)

head(InputData)

############ Некоторая проверка данных

InputData[InputData$Transcript_name == "ENSG00000067715",]
Max_RPKM_DataSet[Max_RPKM_DataSet$Transcript_name == "ENSG00000067715",]

InputData[InputData$RPKM_value[InputData$Impact_type == "Control"] >= 1000,]

############ Расчеты для выбора транскриптов

summary(InputData$RPKM_value[InputData$Impact_type == "Control"])
quantile(InputData$RPKM_value[InputData$Impact_type == "Control"])

unique(InputData$Transcript_name[InputData$Impact_type == "Control" & InputData$RPKM_value >= 100])

length(unique(InputData$Transcript_name[InputData$Impact_type == "Control" & InputData$RPKM_value >= 100]))

Max_RPKM_DataSet <- subset.data.frame(InputData, InputData$RPKM_value[InputData$Impact_type == "Control"] >= 100)

head(Max_RPKM_DataSet)
unique(Max_RPKM_DataSet$Impact_type)

quantile(Max_RPKM_DataSet$RPKM_value[Max_RPKM_DataSet$Impact_type == "Control"])

RPKM_DataSet <- data.frame()


for (Line_counter in unique(InputData$Experiment_Line))
{
  
message("For experiment № ", Line_counter, " process is beginning...")
  
Unique_Transcript_name <- unique(InputData$Transcript_name[InputData$Impact_type == "Control" & InputData$RPKM_value >= 100 & InputData$Experiment_Line == Line_counter])

Max_RPKM_DataSet <- subset.data.frame(InputData, (InputData$Transcript_name %in% Unique_Transcript_name & InputData$Experiment_Line == Line_counter))

Max_RPKM_DataSet$"Max_Diff" <- NA
Max_RPKM_DataSet$"Max_Diff<=5%" <- NA

Process_counter <- 0

for (Trancript_counter in unique(Max_RPKM_DataSet$Transcript_name))
  {
    temp_TranscriptData <- subset.data.frame(Max_RPKM_DataSet, Max_RPKM_DataSet$Transcript_name == Trancript_counter)
    #temp_TranscriptData$SD <- sd(temp_TranscriptData$RPKM_value[temp_TranscriptData$Impact_type %in% c("Control", "Native_LDL", "Atherogenic_LDL")])
    #temp_TranscriptData$SD <- sd(temp_TranscriptData$RPKM_value)
    if("Control" %in%  unique(temp_TranscriptData$Impact_type))
     {
      Mean_Control_Value <- mean(temp_TranscriptData$RPKM_value[temp_TranscriptData$Impact_type == "Control"], na.rm = TRUE)
      max_diff_value <- max(abs(Mean_Control_Value - temp_TranscriptData$RPKM_value[temp_TranscriptData$Impact_type %in% c("Native_LDL", "Atherogenic_LDL")]))
      temp_TranscriptData$Max_Diff <- max_diff_value
      temp_TranscriptData$`Max_Diff<=5%` <- ifelse( (Mean_Control_Value*0.05) >= max_diff_value, TRUE, FALSE)
    
      Max_RPKM_DataSet$Max_Diff[Max_RPKM_DataSet$Transcript_name == Trancript_counter]  <- temp_TranscriptData$Max_Diff 
      Max_RPKM_DataSet$`Max_Diff<=5%`[Max_RPKM_DataSet$Transcript_name == Trancript_counter]  <- temp_TranscriptData$`Max_Diff<=5%`
    
      Process_counter <- Process_counter + 1
    
      message(Process_counter, " transcript from ", length(unique(Max_RPKM_DataSet$Transcript_name)), " is done!")
    } else
         {
           message(Process_counter, " transcript from ", length(unique(Max_RPKM_DataSet$Transcript_name)), " has't control value !!!!!!") 
         }
}

RPKM_DataSet <- rbind.data.frame(RPKM_DataSet, Max_RPKM_DataSet)

}

summary(RPKM_DataSet$Max_Diff)



RPKM_DataSet[RPKM_DataSet$RPKM_value >= 100,]
RPKM_DataSet[RPKM_DataSet$`Max_Diff<=5%` == TRUE,]

Result_DataSet <- RPKM_DataSet[RPKM_DataSet$`Max_Diff<=5%` == TRUE,]

unique(Result_DataSet$Transcript_name)

length(unique(Result_DataSet$Transcript_name))

########## Просмотр результатов в исходном файле Excel

InputData_Excel_Sheet[InputData_Excel_Sheet$gNam %in% unique(Result_DataSet$Transcript_name),]

########## тестовая выборка тысячи значений

seq()

x <- runif(1000, 8, 10)
x
sd(x)

Test_Table <- data.frame()

for (mean_counter in seq(from = 5, to = 2000, by = 10))
  {
   values <- runif(n = 1000, min = mean_counter, max = (mean_counter*1.1))
   temp_df <- data.frame("Mean" = mean(values),
                         "SD" = sd(values),
                         "Mean.SD" = mean(values)/sd(values)
                        )
   Test_Table <- rbind.data.frame(temp_df, Test_Table)
}

Test_Table

plot(x = Test_Table$Mean, y = Test_Table$Mean.SD, data = Test_Table)



x2 <- sample(x = c(9.0:10.0), size = 1000, replace = TRUE)
x2
sd(x2)

x3 <- seq(from = 1, to = 10, by = 0.1)
x3
sd(x3)

######### конец тестовой вставки

######## Запрос к базе данных Ensembl по ensembl_transcript_id

listEnsembl()

listDatasets(useMart("ensembl"))

Ensembl_DataSet <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

listFilters(mart = Ensembl_DataSet)

listAttributes(mart = Ensembl_DataSet)

######## Собственно запрос

affyids <- unique(Result_DataSet$Transcript_name)

affyids <- "ENSG00000067715"

Ensembl_Query_Results <- getBM(attributes = c( "ensembl_gene_id",
                                               "hgnc_symbol",
                                               "hgnc_trans_name",
                                               "transcript_biotype",
                                               "description",
                                               "external_synonym",
                                               "mim_gene_description",
                                               "phenotype_description"
                                              ),
                              filters = "ensembl_gene_id",
                              values = affyids, 
                              mart = Ensembl_DataSet,
                              useCache = FALSE
                              )

Ensembl_Query_Results

unique(Ensembl_Query_Results$hgnc_symbol)
length(unique(Ensembl_Query_Results$hgnc_symbol))

######### Добавление идентификатора для кодирующих генов упоминающихся первый раз

names(Ensembl_Query_Results)
unique(Ensembl_Query_Results$transcript_biotype)

Ensembl_Query_Results$"select" <- NA

Ensembl_Query_Results$select[!(duplicated(Ensembl_Query_Results$hgnc_symbol)) & Ensembl_Query_Results$transcript_biotype == "protein_coding"] <- 1
Protein_coding_IDs <- Ensembl_Query_Results$ensembl_gene_id[!(duplicated(Ensembl_Query_Results$hgnc_symbol)) & Ensembl_Query_Results$transcript_biotype == "protein_coding"]
Protein_coding_IDs

Ensembl_Query_Results$select[!(duplicated(Ensembl_Query_Results$hgnc_symbol)) & Ensembl_Query_Results$transcript_biotype == "nonsense_mediated_decay"] <- 2
nonsense_mediated_decay_IDs <- Ensembl_Query_Results$ensembl_gene_id[!(duplicated(Ensembl_Query_Results$hgnc_symbol)) & Ensembl_Query_Results$transcript_biotype == "nonsense_mediated_decay"]
nonsense_mediated_decay_IDs

Ensembl_Query_Results$select[!(duplicated(Ensembl_Query_Results$hgnc_symbol)) & Ensembl_Query_Results$transcript_biotype == "retained_intron"] <- 3
retained_intron_IDs <- Ensembl_Query_Results$ensembl_gene_id[!(duplicated(Ensembl_Query_Results$hgnc_symbol)) & Ensembl_Query_Results$transcript_biotype == "retained_intron"]
retained_intron_IDs

Ensembl_Query_Results$select[!(duplicated(Ensembl_Query_Results$hgnc_symbol)) & Ensembl_Query_Results$transcript_biotype == "processed_transcript"] <- 4
processed_transcript_IDs <- Ensembl_Query_Results$ensembl_gene_id[!(duplicated(Ensembl_Query_Results$hgnc_symbol)) & Ensembl_Query_Results$transcript_biotype == "processed_transcript"]
processed_transcript_IDs

######### Сохранение результатов в Excel file

Input_Table <- InputData_Excel_Sheet[InputData_Excel_Sheet$gNam %in% unique(Result_DataSet$Transcript_name),]
Result_Table <- Ensembl_Query_Results

str(Input_Table)
str(Result_Table)

Results_List <- list(Input_Table, Result_Table)
names(Results_List) <- c("Из исходной таблицы", "Из базы данных Ensembl")
ResultFileName <- file.path(PathToResults, paste("Кандидатные гены из транскриптома моноцитов_07-09-2021.xlsx", sep = ''))
save_to_excel_file_by_openxlsx(ResultFileName, Results_List, save_row_names = TRUE)


###################### ТРАНСКРИПЦИОННЫЙ АНАЛИЗ #########################

####### Анализ дифференциальной экспрессии

head(InputData)

#### Добавление названия гена в таблицу исходных данных

######## Запрос к базе данных Ensembl по ensembl_transcript_id

listEnsembl()

listDatasets(useMart("ensembl"))

Ensembl_DataSet <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

listFilters(mart = Ensembl_DataSet)

listAttributes(mart = Ensembl_DataSet)

######## Собственно запрос

affyids <- unique(InputData$Transcript_name)

Ensembl_Query_Results <- getBM(attributes = c( "ensembl_gene_id",
                                               "hgnc_symbol"
                                               #"hgnc_trans_name",
                                               #"transcript_biotype",
                                               #"description",
                                               #"external_synonym",
                                               #"mim_gene_description",
                                               #"phenotype_description"
),
filters = "ensembl_gene_id",
values = affyids, 
mart = Ensembl_DataSet,
useCache = FALSE
)

Ensembl_Query_Results

unique(Ensembl_Query_Results$hgnc_symbol)
length(unique(Ensembl_Query_Results$hgnc_symbol))

unique(InputData$Transcript_name)
length(unique(InputData$Transcript_name))

# ##### Собственно добавление названия генов
# 
# InputData$"hgnc_symbol" <- NA
# 
# index <- match(InputData$Transcript_name, Ensembl_Query_Results$ensembl_gene_id)
# InputData[!is.na(index), "hgnc_symbol"] <- Ensembl_Query_Results[index[!is.na(index)], "hgnc_symbol"]
# 
# ######## Неизвестные транскрипты
# 
# InputData[is.na(InputData$hgnc_symbol),]

############ Создание матрицы с численностями ридов

InputData$Impact_type_Experiment_Line <- paste(InputData$Impact_type, InputData$Experiment_Line, sep = "_")

str(InputData)

InputData_spread <- spread(data = InputData[, !(names(InputData) %in% c("Impact_type", "Experiment_Line"))], key = Impact_type_Experiment_Line, value = Reads_count)

row.names(InputData_spread) <- InputData_spread$Transcript_name

InputData_Matrix <- as.matrix(x = InputData_spread[, !(names(InputData_spread) %in% "Transcript_name")])

Metadata_Matrix <- matrix(, nrow = length(colnames(InputData_Matrix)), ncol = 1)

row.names(Metadata_Matrix) <- colnames(InputData_Matrix)

colnames(Metadata_Matrix) <- "Experiment_type"

Metadata_Matrix[, "Experiment_type"] <- str_remove(string = row.names(Metadata_Matrix), pattern = "_[1-3]")

dds <- DESeqDataSetFromMatrix(countData = InputData_Matrix, colData = Metadata_Matrix, design = ~ Experiment_type)

str(dds)

DE_Table <- DESeq(dds)

plotDispEsts(DE_Table)

res <- results(DE_Table)
head(res)

plotMA(res)

DE_results_dataframe <- as.data.frame(res)

ggplot(data = DE_results_dataframe, aes(x = log2FoldChange, y = pvalue)) +
geom_point() +
theme_AO

########### Добавление маркеров для повышенной и пониженной экспрессии

# add a column of NAs
DE_results_dataframe$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
DE_results_dataframe$diffexpressed[DE_results_dataframe$log2FoldChange > 0.6 & DE_results_dataframe$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
DE_results_dataframe$diffexpressed[DE_results_dataframe$log2FoldChange < -0.6 & DE_results_dataframe$pvalue < 0.05] <- "DOWN"

########## Добавление названий генов в таблицу результатов

######## Запрос к базе данных Ensembl по ensembl_transcript_id

listEnsembl()

listDatasets(useMart("ensembl"))

Ensembl_DataSet <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

listFilters(mart = Ensembl_DataSet)

listAttributes(mart = Ensembl_DataSet)

######## Собственно запрос

affyids <- unique(rownames(DE_results_dataframe))

Ensembl_Query_Results <- getBM(attributes = c( "ensembl_gene_id",
                                               "hgnc_symbol"
                                               #"hgnc_trans_name",
                                               #"transcript_biotype",
                                               #"description",
                                               #"external_synonym",
                                               #"mim_gene_description",
                                               #"phenotype_description"
),
filters = "ensembl_gene_id",
values = affyids, 
mart = Ensembl_DataSet,
useCache = FALSE
)

Ensembl_Query_Results

unique(Ensembl_Query_Results$hgnc_symbol)
length(unique(Ensembl_Query_Results$hgnc_symbol))

unique(unique(rownames(DE_results_dataframe)))
length(unique(unique(rownames(DE_results_dataframe))))

##### Собственно добавление названия генов

DE_results_dataframe$Transcript_name <- rownames(DE_results_dataframe)
DE_results_dataframe$hgnc_symbol <- NA

index <- match(DE_results_dataframe$Transcript_name, Ensembl_Query_Results$ensembl_gene_id)
DE_results_dataframe[!is.na(index), "hgnc_symbol"] <- Ensembl_Query_Results[index[!is.na(index)], "hgnc_symbol"]

######## Неизвестные транскрипты

DE_results_dataframe[is.na(DE_results_dataframe$hgnc_symbol),]


########## Отрисовка графика дифференцтальной экспрессии генов

Plot_File_Name <- file.path(PathToResults, paste("Differential gene expression analysis", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

Figure <- ggplot(data = DE_results_dataframe, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = hgnc_symbol)) +
  geom_point(size = 4) +
  #geom_text(label = DE_results_dataframe$hgnc_symbol[DE_results_dataframe$diffexpressed == "UP" | DE_results_dataframe$diffexpressed == "DOWN"], nudge_x = 0.25, nudge_y = 0.4, check_overlap = T) +
  geom_text(label = DE_results_dataframe$hgnc_symbol, nudge_x = 0.25, nudge_y = 0.3, check_overlap = T) +
  #geom_label(data = DE_results_dataframe[(DE_results_dataframe$diffexpressed == "UP" | DE_results_dataframe$diffexpressed == "DOWN"),], aes(label = hgnc_symbol)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme_AO

print(Figure)

dev.off()

############ Выбор из таблицы целевых генов

DE_Genes_Set_Selected <- DE_results_dataframe[DE_results_dataframe$hgnc_symbol %in% c(
  "ANXA1",
  "F2RL1",
  "IL15",
  "PERK",
  "TSPYL2",
  "DUSP1",
  "IL7",
  "IL7R",
  "IL8",
  "TIGIT"
),]

############ Выбор генов с дифференциальной экспрессией

DE_Genes_Set <- DE_results_dataframe[(DE_results_dataframe$diffexpressed == "UP" | DE_results_dataframe$diffexpressed == "DOWN"),]

########### Добавление описания генов с дифференциальной экспрессией

affyids <- unique(DE_Genes_Set$Transcript_name)

Ensembl_Query_Results <- getBM(attributes = c( "ensembl_gene_id",
                                               "hgnc_symbol",
                                               "hgnc_trans_name",
                                               "transcript_biotype",
                                               "description",
                                               "external_synonym",
                                               "mim_gene_description",
                                               "phenotype_description"
),
filters = "ensembl_gene_id",
values = affyids, 
mart = Ensembl_DataSet,
useCache = FALSE
)

Ensembl_Query_Results

names(Ensembl_Query_Results)[which(names(Ensembl_Query_Results) == "ensembl_gene_id")] <- names(DE_Genes_Set)[which(names(DE_Genes_Set) == "Transcript_name")]

DE_Genes_Set_Desc <- merge(DE_Genes_Set, Ensembl_Query_Results, by = "Transcript_name")

DE_Genes_Set_Desc$hgnc_symbol.y <- NULL

#### Сохранение результатов  в Excel файл

str(DE_results_dataframe)
str(DE_Genes_Set)
str(DE_Genes_Set_Selected)

rownames(DE_results_dataframe) <- NULL
rownames(DE_Genes_Set) <- NULL
rownames(DE_Genes_Set_Selected) <- NULL

List_Of_Result_Table_All <- List_Of_Result_Table[1]
names(List_Of_Result_Table_All) <- c("Все гены")
Table_ResultFileName <- file.path(PathToResults, paste("Результаты анализа дифференциальной экспрессии всех генов", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table_All)

CSV_ResultFileName <- file.path(PathToResults, paste("Результаты анализа дифференциальной экспрессии всех генов", Sys.Date(), ".csv", sep = "_"))
write.csv(List_Of_Result_Table[[1]], file = CSV_ResultFileName, row.names = FALSE)

CSV_ResultFileName <- file.path(PathToResults, paste("Результаты анализа дифференциальной экспрессии отличающихся генов", Sys.Date(), ".csv", sep = "_"))
write.csv(List_Of_Result_Table[[2]], file = CSV_ResultFileName, row.names = FALSE)

CSV_ResultFileName <- file.path(PathToResults, paste("Результаты анализа дифференциальной экспрессии 10 генов Василия", Sys.Date(), ".csv", sep = "_"))
write.csv(List_Of_Result_Table[[3]], file = CSV_ResultFileName, row.names = FALSE)

####### Сохраниение листа в файл RData

list.files(PathToSources)

Table_FileName <- file.path(PathToSources, paste("Результаты анализа дифференциальной экспрессии генов", Sys.Date(), ".Rdata", sep = "_"))
saveRDS(List_Of_Result_Table, file = Table_FileName)

###### Открытие файла с таблицами дифференциальной экспрессии генов

InputData_FileName <- file.path(PathToSources, paste("Результаты анализа дифференциальной экспрессии генов_2022-01-23_.Rdata", sep = "_"))
List_Of_Result_Table <- readRDS(file = InputData_FileName)

str(List_Of_Result_Table)

names(List_Of_Result_Table)

List_Of_Result_Table$`Все гены`

############ Работа с данными первого опыта

row.names(InputData) <- InputData$Transcript_name

InputData_Matrix <- as.matrix(x = InputData[, !(names(InputData) %in% "Transcript_name")])

Metadata_Matrix <- matrix(, nrow = length(colnames(InputData_Matrix)), ncol = 1)

row.names(Metadata_Matrix) <- colnames(InputData_Matrix)

colnames(Metadata_Matrix) <- "Experiment_type"

Metadata_Matrix[which(row.names(Metadata_Matrix) %in% c("Control", "Native_LDL", "HDL")), "Experiment_type"] <- "no_accumulation"
Metadata_Matrix[which(row.names(Metadata_Matrix) %in% c("Ox_LDL", "Acc_LDL")), "Experiment_type"] <- "accumulation"

dds <- DESeqDataSetFromMatrix(countData = InputData_Matrix, colData = Metadata_Matrix, design = ~ Experiment_type)

str(dds)

DE_Table <- DESeq(dds)

plotDispEsts(DE_Table)

res <- results(DE_Table)
head(res)

plotMA(res)

DE_results_dataframe <- as.data.frame(res)

ggplot(data = DE_results_dataframe, aes(x = log2FoldChange, y = pvalue)) +
  geom_point() +
  theme_AO

########### Добавление маркеров для повышенной и пониженной экспрессии

# add a column of NAs
DE_results_dataframe$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
DE_results_dataframe$diffexpressed[DE_results_dataframe$log2FoldChange > 0.6 & DE_results_dataframe$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
DE_results_dataframe$diffexpressed[DE_results_dataframe$log2FoldChange < -0.6 & DE_results_dataframe$pvalue < 0.05] <- "DOWN"

########## Добавление названий генов в таблицу результатов

######## Запрос к базе данных Ensembl по ensembl_transcript_id

listEnsembl()

listDatasets(useMart("ensembl"))

Ensembl_DataSet <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

listFilters(mart = Ensembl_DataSet)

listAttributes(mart = Ensembl_DataSet)

######## Собственно запрос

affyids <- unique(rownames(DE_results_dataframe))

Ensembl_Query_Results <- getBM(attributes = c( "ensembl_gene_id",
                                               "hgnc_symbol"
                                               #"hgnc_trans_name",
                                               #"transcript_biotype",
                                               #"description",
                                               #"external_synonym",
                                               #"mim_gene_description",
                                               #"phenotype_description"
),
filters = "ensembl_gene_id",
values = affyids, 
mart = Ensembl_DataSet,
useCache = FALSE
)

Ensembl_Query_Results

unique(Ensembl_Query_Results$hgnc_symbol)
length(unique(Ensembl_Query_Results$hgnc_symbol))

unique(unique(rownames(DE_results_dataframe)))
length(unique(unique(rownames(DE_results_dataframe))))

##### Собственно добавление названия генов

DE_results_dataframe$Transcript_name <- rownames(DE_results_dataframe)
DE_results_dataframe$hgnc_symbol <- NA

index <- match(DE_results_dataframe$Transcript_name, Ensembl_Query_Results$ensembl_gene_id)
DE_results_dataframe[!is.na(index), "hgnc_symbol"] <- Ensembl_Query_Results[index[!is.na(index)], "hgnc_symbol"]

######## Неизвестные транскрипты

DE_results_dataframe[is.na(DE_results_dataframe$hgnc_symbol),]

#### Добавление номера транскрипта к генам с неизвестным названием гена

DE_results_dataframe$hgnc_symbol[which(is.na(DE_results_dataframe$hgnc_symbol))] <- DE_results_dataframe$Transcript_name[which(is.na(DE_results_dataframe$hgnc_symbol))]


########## Отрисовка графика дифференцтальной экспрессии генов

Plot_File_Name <- file.path(PathToResults, paste("Differential gene expression analysis from first experience", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

Figure <- ggplot(data = DE_results_dataframe, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = hgnc_symbol)) +
  geom_point(size = 4) +
  #geom_text(label = DE_results_dataframe$hgnc_symbol[DE_results_dataframe$diffexpressed == "UP" | DE_results_dataframe$diffexpressed == "DOWN"], nudge_x = 0.25, nudge_y = 0.4, check_overlap = T) +
  geom_text(label = DE_results_dataframe$hgnc_symbol, nudge_x = 0.25, nudge_y = 0.3, check_overlap = T) +
  #geom_label(data = DE_results_dataframe[(DE_results_dataframe$diffexpressed == "UP" | DE_results_dataframe$diffexpressed == "DOWN"),], aes(label = hgnc_symbol)) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme_AO

print(Figure)

dev.off()

############ Выбор из таблицы целевых генов

DE_Genes_Set_Selected <- DE_results_dataframe[DE_results_dataframe$hgnc_symbol %in% c(
  "ANXA1",
  "F2RL1",
  "IL15",
  "PERK",
  "TSPYL2",
  "DUSP1",
  "IL7",
  "IL7R",
  "IL8",
  "TIGIT"
),]

############ Выбор генов с дифференциальной экспрессией

DE_Genes_Set <- DE_results_dataframe[(DE_results_dataframe$diffexpressed == "UP" | DE_results_dataframe$diffexpressed == "DOWN"),]

########### Добавление описания генов с дифференциальной экспрессией

affyids <- unique(DE_Genes_Set$Transcript_name)

Ensembl_Query_Results <- getBM(attributes = c( "ensembl_gene_id",
                                               "hgnc_symbol",
                                               "hgnc_trans_name",
                                               "transcript_biotype",
                                               "description",
                                               "external_synonym",
                                               "mim_gene_description",
                                               "phenotype_description"
),
filters = "ensembl_gene_id",
values = affyids, 
mart = Ensembl_DataSet,
useCache = FALSE
)

Ensembl_Query_Results

names(Ensembl_Query_Results)[which(names(Ensembl_Query_Results) == "ensembl_gene_id")] <- names(DE_Genes_Set)[which(names(DE_Genes_Set) == "Transcript_name")]

DE_Genes_Set_Desc <- merge(DE_Genes_Set, Ensembl_Query_Results, by = "Transcript_name")

DE_Genes_Set_Desc$hgnc_symbol.y <- NULL

#### Сохранение результатов  в Excel файл

str(DE_results_dataframe)
str(DE_Genes_Set)
str(DE_Genes_Set_Selected)

rownames(DE_results_dataframe) <- NULL
rownames(DE_Genes_Set) <- NULL
rownames(DE_Genes_Set_Selected) <- NULL

List_Of_Result_Table <- list(DE_results_dataframe, DE_Genes_Set, DE_Genes_Set_Selected)
names(List_Of_Result_Table) <- c("Все гены", "Отличающаяся экспрессия", "Гены от Василия")
Table_ResultFileName <- file.path(PathToResults, paste("Результаты анализа дифференциальной экспрессии всех генов", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)


##### открытие файла Excel с DEGs

FileName <- file.path(PathToResults, paste("Результаты анализа дифференциальной экспрессии всех генов_2022-01-30_.xlsx"))
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$`Отличающаяся экспрессия`

DEGs_Genes_Table <- InputData_Excel_Sheet

####### открытие файла с генами, ассоциированными с атеросклерозом

Atherosclerosis_Genes_Table <- read.table(file = file.path(PathToSources, paste("Top_100_genes_atherosclerosis_14-02-2022.txt")), sep = "\n")
names(Atherosclerosis_Genes_Table) <- "Atherosclerosis_Genes"

ER_stress_Genes_Table <- read.table(file = file.path(PathToSources, paste("ER_stress_Genes_17-02-2022.txt")), sep = "\n")
names(ER_stress_Genes_Table) <- "ER_stress_Genes"

Vasily_interested_Genes_Table <- read.table(file = file.path(PathToSources, paste("Vasily_interested_genes_24-09-2021.txt")), sep = "\n")
names(Vasily_interested_Genes_Table) <- "Vasily_interested_genes"

###### Сколько DEGs относятся к атеросклерозу

DEGs_Genes_Table$hgnc_symbol[DEGs_Genes_Table$hgnc_symbol %in% Atherosclerosis_Genes_Table$Atherosclerosis_Genes]

DEGs_Genes_Table[DEGs_Genes_Table$hgnc_symbol %in% Atherosclerosis_Genes_Table$Atherosclerosis_Genes,]

###### Сколько DEGs относятся к стрессу ЭР

DEGs_Genes_Table$hgnc_symbol[DEGs_Genes_Table$hgnc_symbol %in% ER_stress_Genes_Table$ER_stress_Genes]

DEGs_Genes_Table[DEGs_Genes_Table$hgnc_symbol %in% ER_stress_Genes_Table$ER_stress_Genes, ]

######## Открытие файлов с таблицами tsv с генами DEGs

list.files(PathToDerivedData)

tsv_files <- list.files(path = file.path(PathToDerivedData, "Atherosclerosis_Genes"))

tsv_files

Atherosclerosis_Network_Table <- data.frame()

for (files_counter in seq_along(tsv_files))
  {
   temp_file_content <- read.csv(file = file.path(PathToDerivedData, "Atherosclerosis_Genes", tsv_files[files_counter]), sep = "\t")
   protein_name <- str_remove(string = str_extract(string =  tsv_files[files_counter], pattern = "^.+_s"), pattern = "_s")
   temp_file_content$"node_of_interes" <- protein_name
   Atherosclerosis_Network_Table <- rbind.data.frame(Atherosclerosis_Network_Table, temp_file_content)
  }

names(Atherosclerosis_Network_Table)[1] <- "node1"

######## Открытие файлов с таблицами tsv с генами ER stress

tsv_files <- list.files(path = file.path(PathToDerivedData, "ER_stress_Genes"))

tsv_files

ER_stress_Network_Table <- data.frame()

for (files_counter in seq_along(tsv_files))
{
  temp_file_content <- read.csv(file = file.path(PathToDerivedData, "ER_stress_Genes", tsv_files[files_counter]), sep = "\t")
  protein_name <- str_remove(string = str_extract(string =  tsv_files[files_counter], pattern = "^.+_s"), pattern = "_s")
  temp_file_content$"node_of_interes" <- protein_name
  ER_stress_Network_Table <- rbind.data.frame(ER_stress_Network_Table, temp_file_content)
}

names(ER_stress_Network_Table)[1] <- "node1"

####### Поиск генов интересов в таблице с атеросклерозом

str(Atherosclerosis_Network_Table)

Atherosclerosis_Network_Table$node1[Atherosclerosis_Network_Table$node1 %in% as.character(Vasily_interested_Genes_Table$Vasily_interested_genes)]
Atherosclerosis_Network_Table$node2[Atherosclerosis_Network_Table$node2 %in% as.character(Vasily_interested_Genes_Table$Vasily_interested_genes)]

####### Поиск генов интересов в таблице с ER stress

str(ER_stress_Network_Table)

ER_stress_Network_Table$node1[ER_stress_Network_Table$node1 %in% as.character(Vasily_interested_Genes_Table$Vasily_interested_genes)]
ER_stress_Network_Table$node2[ER_stress_Network_Table$node2 %in% as.character(Vasily_interested_Genes_Table$Vasily_interested_genes)]

############################   Gene Set Enrichment Analysis ##############################

# reading in data from deseq2

list.files(PathToResults)
CSV_ResultFileName <- file.path(PathToResults, "Результаты анализа дифференциальной экспрессии всех генов_2022-01-23_.csv")
df <- read.csv(CSV_ResultFileName, header = TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$Transcript_name

# omit any NA values 
gene_list <- na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list <- sort(gene_list, decreasing = TRUE)

######## Выбор параметров для аннотации исследуемого генома

keytypes(org.Hs.eg.db)

gse <- gseGO(geneList = gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 1000, 
             minGSSize = 3, 
             maxGSSize = 100, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

gse

###### Отрисовка точечного графика

Plot_File_Name <- file.path(PathToResults, paste("Dotplot of GSEA", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

dev.off()

####### Encrichment Map

Plot_File_Name <- file.path(PathToResults, paste("Encrichment Map of GSEA", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

emapplot(gse, showCategory = 10)

dev.off()


##### Category Netplot categorySize can be either 'pvalue' or 'geneNum'
##### The cnetplot depicts the linkages of genes and biological concepts

Plot_File_Name <- file.path(PathToResults, paste("Category Netplot of GSEA", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

cnetplot(gse, categorySize = "pvalue", foldChange = gene_list, showCategory = 5)

dev.off()

####### Ridgeplot

Plot_File_Name <- file.path(PathToResults, paste("Ridgeplot of GSEA", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

ridgeplot(gse) + labs(x = "enrichment distribution")

dev.off()

####### GSEA Plot
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId

Plot_File_Name <- file.path(PathToResults, paste("GSEA Plot", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

dev.off()

################ PubMed trend of enriched terms

Plot_File_Name <- file.path(PathToResults, paste("PubMed trend of enriched terms", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

terms <- gse$Description[1:5]
pmcplot(terms, 2012:2022, proportion=FALSE)

dev.off()



################# KEGG Gene Set Enrichment Analysis #################

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids <- bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids <- ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 <- df[df$Transcript_name %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y <- dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)


########## Create gseKEGG object ############

kegg_organism <- "hsa"

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 1000,
               minGSSize    = 3,
               maxGSSize    = 100,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

###### Отрисовка точечного графика

Plot_File_Name <- file.path(PathToResults, paste("Dotplot of KEGG", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

dev.off()

####### Encrichment Map

Plot_File_Name <- file.path(PathToResults, paste("Encrichment Map of KEGG", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

emapplot(kk2)

dev.off()

##### Category Netplot categorySize can be either 'pvalue' or 'geneNum'
##### The cnetplot depicts the linkages of genes and biological concepts

Plot_File_Name <- file.path(PathToResults, paste("Category Netplot of KEGG", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list)

dev.off()

####### Ridgeplot

Plot_File_Name <- file.path(PathToResults, paste("Ridgeplot of KEGG", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

ridgeplot(kk2) + labs(x = "enrichment distribution")

dev.off()

####### GSEA Plot
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId

Plot_File_Name <- file.path(PathToResults, paste("GSEA Plot by KEGG", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

gseaplot(x = kk2, by = "all", title = gse$Description[1], geneSetID = 1)

dev.off()

############# Pathview ########

gseKEGG_Table <- gseKEGG(geneList = kegg_gene_list)

gseKEGG_Table

# Produce the native KEGG plot (PNG)
hsa <- pathview(gene.data = kegg_gene_list, pathway.id = "hsa00010", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
hsa <- pathview(gene.data = kegg_gene_list, pathway.id="hsa00010", species = kegg_organism, kegg.native = FALSE)


Plot_File_Name <- file.path(PathToResults, paste("Pathview by KEGG", Sys.Date(), ".png", sep = "_"))

knitr::include_graphics(file.path(PathToWorkingDirectory, "hsa00010.pathview.png"))




