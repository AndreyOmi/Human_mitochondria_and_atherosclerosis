### Автор: Андрей Омельченко
### Дата:  07-09-2021
### Название: Скрипт для исследований провосполительного ответа и толерантности цибридов.
###


library(stringi)
library(stringr)
library(reshape2)
library(nortest)
library(ggplot2)
library(car) ###### для теста на выбросы
library(hablar)   ### изменение типов для полей таблиц
library(magrittr) ### конвеерные операции
library(pgirmess) ### Spatial Analysis and Data Mining for Field Ecologists
library(PMCMR) #### Calculate Pairwise Multiple Comparisons of Mean Rank Sums
library(tidyr) #### конвертация между широкими и длинными форматами таблиц
library(Rmisc) #### for CI calculating
library(ggpubr) ### для отрисовки уровней значимости на графиков сравнения
library(jsonlite) ### чтение и запись в JSON файл
library(psych) ### for PCA and EFA and corr.test()

rm(list = ls())

PathToResults <- "../Results"
PathToSources <- "../Sources"
PathToReferences <- "../References"
PathToDerivedData <- "../DerivedData"


PathToWorkingDirectory <- "C:\\Users\\Андрей\\OneDrive\\Documents\\Work\\НИИМЧ\\Провоспаление и толерантность цибридов\\Scripts"

setwd(PathToWorkingDirectory)

PathToSource_R_Functions <- "C:\\Users\\Андрей\\OneDrive\\Documents\\Work\\R_Projects"

source(file.path(PathToSource_R_Functions, "f_R_functions_for_Science_and_Develoment_AO_13-01-2022.R"))


################# ПЕРВИЧНАЯ РАБОТА С ДАННЫМИ ###########################

##### открытие файла Excel с данными о провоспалении от 19-11-2021

FileName <- file.path(PathToSources, "Рейтинг воспаления - цибридные линии 19.11.2021.xlsx")
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$`Ответ по цитокинам_for_R`

unique(InputData_Excel_Sheet$Номер.эксперимента)

##### Подготовка данных для анализа

names(InputData_Excel_Sheet)
names(InputData_Excel_Sheet) <- c("number",
                                  "date",
                                  "type",
                                  "cell_line",
                                  "1st_LPS_hit",
                                  "2st_LPS_hit",            
                                  "IL-1-beta",
                                  "TNF",
                                  "MCP1",
                                  "IL-6",
                                  "IL-8",
                                  "IL-1-beta_tolerance",
                                  "TNF_tolerance",
                                  "MCP1_tolerance",
                                  "IL-6_tolerance",
                                  "IL-8_tolerance"
)


str(InputData_Excel_Sheet)

unique(InputData_Excel_Sheet$`1st_LPS_hit`)
unique(InputData_Excel_Sheet$`2st_LPS_hit`)

InputData_Excel_Sheet$`2st_LPS_hit`[is.na(InputData_Excel_Sheet$`2st_LPS_hit`)] <- 0

InputData_Excel_Sheet$"stimulation" <- NA

for (rows_counter in 1:nrow(InputData_Excel_Sheet))
{
  if(InputData_Excel_Sheet$`1st_LPS_hit`[rows_counter] == 0 & InputData_Excel_Sheet$`2st_LPS_hit`[rows_counter] == 0)
  {
    InputData_Excel_Sheet$stimulation[rows_counter] <- "no"
  } else if(InputData_Excel_Sheet$`1st_LPS_hit`[rows_counter] == 1 & InputData_Excel_Sheet$`2st_LPS_hit`[rows_counter] == 0)
  {
    InputData_Excel_Sheet$stimulation[rows_counter] <- "first"
  } else if(InputData_Excel_Sheet$`1st_LPS_hit`[rows_counter] == 0 & InputData_Excel_Sheet$`2st_LPS_hit`[rows_counter] == 1)
  {
    InputData_Excel_Sheet$stimulation[rows_counter] <- "second"
  } else if(InputData_Excel_Sheet$`1st_LPS_hit`[rows_counter] == 1 & InputData_Excel_Sheet$`2st_LPS_hit`[rows_counter] == 1)
  {
    InputData_Excel_Sheet$stimulation[rows_counter] <- "both"
  }
}

unique(InputData_Excel_Sheet$stimulation)

Proinflam_Subset <- InputData_Excel_Sheet[, c("number", "type", "cell_line", "stimulation", "IL-1-beta", "TNF", "MCP1", "IL-6", "IL-8")]
melt_Proinflam_Subset <- melt(data = Proinflam_Subset, id = c("number", "type", "cell_line", "stimulation"))

names(melt_Proinflam_Subset)[c(5,6)] <- c("mediator", "concentration")

str(melt_Proinflam_Subset)

melt_Proinflam_Subset$mediator <- as.character(melt_Proinflam_Subset$mediator)

####### Сохранение таблицы данных с провоспалением в файл Rdata

Save_FileName <- file.path(PathToSources, paste("Провоспаление цибридов_19-11-2021.Rdata", sep = ''))
saveRDS(melt_Proinflam_Subset, file = Save_FileName)

ordered(list.files(PathToSources))

##### открытие файла Excel с данными о провоспалении от 14-10-2021

FileName <- file.path(PathToSources, "Провоспаление и толерантность цибридов_14-10-2021.xlsx")
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$Лист1

unique(InputData_Excel_Sheet$Номер.эксперимента)

##### Подготовка данных для анализа

names(InputData_Excel_Sheet)
names(InputData_Excel_Sheet) <- c("number",
                                  "date",
                                  "type",
                                  "cell_line",
                                  "1st_LPS_hit",
                                  "2st_LPS_hit",            
                                  "IL-1-beta",
                                  "TNF",
                                  "MCP1",
                                  "IL-6",
                                  "IL-8"
                                  )


str(InputData_Excel_Sheet)

unique(InputData_Excel_Sheet$`1st_LPS_hit`)
unique(InputData_Excel_Sheet$`2st_LPS_hit`)

InputData_Excel_Sheet$`2st_LPS_hit`[is.na(InputData_Excel_Sheet$`2st_LPS_hit`)] <- 0

InputData_Excel_Sheet$"stimulation" <- NA

for (rows_counter in 1:nrow(InputData_Excel_Sheet))
 {
  if(InputData_Excel_Sheet$`1st_LPS_hit`[rows_counter] == 0 & InputData_Excel_Sheet$`2st_LPS_hit`[rows_counter] == 0)
        {
         InputData_Excel_Sheet$stimulation[rows_counter] <- "no"
        } else if(InputData_Excel_Sheet$`1st_LPS_hit`[rows_counter] == 1 & InputData_Excel_Sheet$`2st_LPS_hit`[rows_counter] == 0)
                {
                 InputData_Excel_Sheet$stimulation[rows_counter] <- "first"
                } else if(InputData_Excel_Sheet$`1st_LPS_hit`[rows_counter] == 0 & InputData_Excel_Sheet$`2st_LPS_hit`[rows_counter] == 1)
                        {
                         InputData_Excel_Sheet$stimulation[rows_counter] <- "second"
                        } else if(InputData_Excel_Sheet$`1st_LPS_hit`[rows_counter] == 1 & InputData_Excel_Sheet$`2st_LPS_hit`[rows_counter] == 1)
                                {
                                 InputData_Excel_Sheet$stimulation[rows_counter] <- "both"
                                }
 }

unique(InputData_Excel_Sheet$stimulation)

Proinflam_Subset <- InputData_Excel_Sheet[, c("number", "type", "cell_line", "stimulation", "IL-1-beta", "TNF", "MCP1", "IL-6", "IL-8")]
melt_Proinflam_Subset <- melt(data = Proinflam_Subset, id = c("number", "type", "cell_line", "stimulation"))

names(melt_Proinflam_Subset)[c(5,6)] <- c("mediator", "concentration")

str(melt_Proinflam_Subset)

melt_Proinflam_Subset$mediator <- as.character(melt_Proinflam_Subset$mediator)



####### Открытие файла с митофагией в клеточных линиях

FileName <- file.path(PathToSources, "Митофагия в цибридах_AO.xlsx")
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$Лист1

names(InputData_Excel_Sheet)
str(InputData_Excel_Sheet)

####### Сохранение таблицы данных с митофагией в клеточных линиях в файл Rdata

Save_FileName <- file.path(PathToSources, paste("Митофагия в цибридах_01-12-2021.Rdata", sep = ''))
saveRDS(InputData_Excel_Sheet, file = Save_FileName)

####### Сохранение таблицы данных с провоспалительными медиаторами в файл Rdata

Save_FileName <- file.path(PathToSources, paste("Провоспаление цибридов_14-10-2021.Rdata", sep = ''))
saveRDS(melt_Proinflam_Subset, file = Save_FileName)

list.files(PathToSources)

##### открытие файла Excel с данными рейтинговой оценки провоспаления

FileName <- file.path(PathToSources, "Рейтинг воспаления - цибридные линии 19.11.2021.xlsx")
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$`Рейтинг воспаления_for_R_2`

str(InputData_Excel_Sheet)

unique(InputData_Excel_Sheet$Линия)
unique(Input_Inflammatory_Data$cell_line)

setdiff(unique(InputData_Excel_Sheet$Линия), unique(Input_Inflammatory_Data$cell_line))
setdiff(unique(Input_Inflammatory_Data$cell_line), unique(InputData_Excel_Sheet$Линия))

Ranking_Data <- InputData_Excel_Sheet

Ranking_Data <- aggregate(InputData_Excel_Sheet[,3:5], by = list(InputData_Excel_Sheet$Линия), FUN = mean_without_NA)
Ranking_Data


names(Ranking_Data) <- c("cell_line", "Inflammation_Rank_Sum", "Inflammation_Rank_First", "Inflammation_Rank_Second")

####### Сохранение таблицы данных с рейтингами в файл Rdata

Save_FileName <- file.path(PathToSources, paste("Рейтинг воспаления обобщенные данные_19-11-2021.Rdata", sep = ''))
saveRDS(Ranking_Data, file = Save_FileName)

ordered(list.files(PathToSources))


########################### ОТКРЫТИЕ ФАЙЛА ДАННЫХ О ПРОВОСПАЛЕНИИ ###########################

InputData_FileName <- file.path(PathToSources, "Провоспаление цибридов_19-11-2021.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
Input_Inflammatory_Data <- ListOfInputData_FileSheets

names(Input_Inflammatory_Data)
str(Input_Inflammatory_Data)

xtabs(~ stimulation + type, data =  Input_Inflammatory_Data)

Input_Inflammatory_Data$type[which(Input_Inflammatory_Data$type == "К")] <- "K"

Input_Inflammatory_Data$type <- str_replace_all(string = Input_Inflammatory_Data$type, pattern = "-", replacement = "_")

##################### удаление эксперимента №1 от 20.03.2020

Input_Inflammatory_Data <- Input_Inflammatory_Data[Input_Inflammatory_Data$number != 1,]

####### Открытие файла данных с различающимися мутациями

InputData_FileName <- file.path(PathToSources, "Таблица различающихся митохондриальных мутаций_27-09-2021.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
Diff_Mutations_Table <- ListOfInputData_FileSheets

names(Diff_Mutations_Table)
str(Diff_Mutations_Table)

Diff_Mutations_Table %<>% convert(
  chr(
      Original_Line_Name,
      Derivative_Line_Name,
      Names_Diff_mutations_Plus,
      Names_Diff_mutations_Minus
     )
)

Diff_Mutations_Table <- as.data.frame(Diff_Mutations_Table)

Diff_Mutations_Table

###### Открытие файла с данными по клеточным линиям, их мутациям и гетероплазмии

list.files(path = PathToSources)

InputData_FileName <- file.path(PathToSources, "Длинный_формат_таблицы_исходных_данных_18-05-2021.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
Input_Data <- ListOfInputData_FileSheets

names(Input_Data)
str(Input_Data)

Input_Data %<>% convert(
  chr(
    ID,
    Line,
    Mutation_name,
    Agent
  )
)

str(Input_Data)
class(Input_Data)

Input_Data_Lines_Mutations_Heteroplasmy <- as.data.frame(Input_Data)
Input_Data_Lines_Mutations_Heteroplasmy_base_cons <- subset.data.frame(x = Input_Data_Lines_Mutations_Heteroplasmy, Input_Data_Lines_Mutations_Heteroplasmy$Agent == "base_cons")
Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line <- str_remove_all(Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line, "Native ")
Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line <- str_remove_all(Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line, "-")
Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line[which(Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line == "HSMAM1")] <- "MAM1"
Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line[which(Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line == "HSMAM2")] <- "MAM2"
Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line[which(Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line == "HSMAM3")] <- "MAM3"
Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line[which(Input_Data_Lines_Mutations_Heteroplasmy_base_cons$Line == "THP1")] <- "THP1"

Input_Data <- as.data.frame(Input_Data)

###### Открытие файла с данными с митофагией в клеточных линиях

list.files(path = PathToSources)

InputData_FileName <- file.path(PathToSources, "Митофагия в цибридах_01-12-2021.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
Input_Data_Mitophagia <- ListOfInputData_FileSheets

names(Input_Data_Mitophagia)
str(Input_Data_Mitophagia) 

###### коррекция названий полей таблицы (если нужно)

names(Input_Data_Mitophagia)[which(names(Input_Data_Mitophagia) == "Cell.line")] <- "cell_line"

###### Выбор мутаций без воздействия и с гетероплазмией больше 5

Input_Data[is.na(Input_Data$Heteroplasia_value), "Heteroplasia_value"] <- 0

Input_Data_Native <- subset(Input_Data, !(Input_Data$Agent %in% c("natLDL", "athLDL", "base_v", "oligo_v", "FCCP_v")) & Input_Data$Agent == "base_cons" & Input_Data$Heteroplasia_value > 5)

unique(Input_Data_Native$Mutation_name)
length(unique(Input_Data_Native$Mutation_name))

unique(Input_Data_Native$Agent)
length(unique(Input_Data_Native$Agent))

unique(Input_Data_Native$Line) 
length(unique(Input_Data_Native$Line))

Input_Data_Native

Input_Data_Lines_Mutations_Heteroplasmy <- Input_Data_Native

###### Открытие файла с рейтингом провоспаления

InputData_FileName <- file.path(PathToSources, "Рейтинг воспаления обобщенные данные_19-11-2021.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
Ranking_Data <- ListOfInputData_FileSheets

names(Ranking_Data)

###### коррекция названий полей таблицы (если нужно)

names(Ranking_Data)[which(names(Ranking_Data) == "Линия")] <- "Line"
names(Ranking_Data)[which(names(Ranking_Data) == "Общий")] <- "Inflammation_Rank_Sum"
names(Ranking_Data)[which(names(Ranking_Data) == "База.Активация")] <- "Первая_стимуляция"
names(Ranking_Data)[which(names(Ranking_Data) == "Толерантность")] <- "Вторая_стимуляция"

str(Ranking_Data)

Ranking_Data

######################## ОТКРЫТИЕ ФАЙЛА С ДАННЫМИ С МИТОФАГИЕЙ И С ПРОВОСПАЛЕНИЕМ ######################################

####### открытие файла Excel

list.files(path = PathToSources)

FileName <- file.path(PathToSources, "Митофагия в цибридах_индивидуальные клетки_АО.xlsx")
ListOfFileSheets <- get_excel_file_content_by_xlconnect(FileName)
Input_Data <- ListOfFileSheets$`Исходные данные`

###### предварительная обработка данных

#### 1 - non-responder
#### 2 - tolerance
#### 3 - non-tolerance

Input_Data$"Inflamantion_type" <- NA

for(rows_counter in 1:nrow(Input_Data))
  {
   switch(Input_Data$Группа.по.типу.воспалительного.ответа..1...нонресп..2.толер..3.нетолер[rows_counter],
          '1'  = Input_Data$Inflamantion_type[rows_counter] <- "non-responder",
          '2' = Input_Data$Inflamantion_type[rows_counter] <- "tolerance",
          '3' = Input_Data$Inflamantion_type[rows_counter] <- "non-tolerance"
         )

  }

############ КОНЕЦ ОТКРЫТИЯ И ПРЕОБРАЗОВАНИЯ ФАЙЛА С ДАННЫМИ С МИТОФАГИЕЙ И С ПРОВОСПАЛЕНИЕМ ##############################

###### Для работы с первой стимуляцией

names(Ranking_Data)[which(names(Ranking_Data) == "Первая_стимуляция")] <- "Inflammation_Rank_First"

###### Для работы со второй стимуляцией

names(Ranking_Data)[which(names(Ranking_Data) == "Вторая_стимуляция")] <- "Inflammation_Rank_Second"

###### Визуализация исходных данных данных для всех клеточных линий вместе

Plot_File_Name <- file.path(PathToResults, "Влияние LPS на накопление провоспалительных медиаторов_28-09-2021.jpg")
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

ggplot(Input_Inflammatory_Data, aes(x = stimulation, y = concentration, fill = mediator)) +
  geom_boxplot(outlier.size = 3, position = position_dodge(width = 0.9)) +
  stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE,) +
  scale_fill_brewer(palette = "Set2", name = "Провоспалительный\n медиатор") +
  #scale_shape_manual(values=c(1:14), name = "Клеточная линия") +
  scale_x_discrete(limits=c("no", "first", "second", "both"), labels = c("Контроль", "только\nпервая\n стимуляция", "только\nвторая\n стимуляция", "две стимуляции")) +
  scale_y_log10(breaks=10^(0:3)) +
  xlab(label = "Стимуляция") +
  ylab(label = "Концентрация медиатора") +
  geom_point(aes(shape = "выброс"),  alpha = 0) +
  geom_point(aes(shape = "медиана"),  alpha = 0) +
  geom_point(aes(shape = "среднее арифметическое"),  alpha = 0) +
  guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
  ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
  theme_AO +
  theme(legend.key.size = unit(1, "cm")) 

dev.off()

Input_Inflammatory_Data$stimulation <- factor(Input_Inflammatory_Data$stimulation, levels = c(c("no", "first", "second", "both")))

Plot_File_Name <- file.path(PathToResults, "Влияние LPS на накопление провоспалительных медиаторов_28-09-2021_2.jpg")
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

ggplot(Input_Inflammatory_Data, aes(x = mediator, y = concentration, fill = stimulation)) +
  geom_boxplot(outlier.size = 3, position = position_dodge(width = 0.9)) +
  stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE,) +
  scale_fill_brewer(palette = "Set3", name = "Стимуляция", limits = c("no", "first", "second", "both"), labels = c("Контроль", "только первая стимуляция", "только вторая стимуляция", "две стимуляции")) +
  #scale_shape_manual(values=c(1:14), name = "Клеточная линия") +
  #scale_x_discrete(limits=c("no", "first", "second", "both"), labels = c("Контроль", "только\nпервая стимуляция", "только\nвторая стимуляция", "две стимуляции")) +
  scale_y_log10(breaks=10^(0:3)) +
  xlab(label = "Провоспалительный медиатор") +
  ylab(label = "Концентрация медиатора") +
  geom_point(aes(shape = "выброс"),  alpha = 0) +
  geom_point(aes(shape = "медиана"),  alpha = 0) +
  geom_point(aes(shape = "среднее арифметическое"),  alpha = 0) +
  guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
  ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
  theme_AO +
  theme(legend.key.size = unit(1, "cm"))


dev.off()

unique(Input_Inflammatory_Data$type)

Input_Inflammatory_Data$type <- factor(Input_Inflammatory_Data$type, levels = c("K", "K_K", "LPS", "K_LPS", "LPS_K", "LPS_LPS"))

Plot_File_Name <- file.path(PathToResults, "Влияние LPS на накопление провоспалительных медиаторов_28-09-2021_3.jpg")
jpeg(filename = Plot_File_Name, width = 1500, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

ggplot(Input_Inflammatory_Data, aes(x = type, y = concentration, fill = mediator)) +
  geom_boxplot(outlier.size = 3, position = position_dodge(width = 0.9)) +
  stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE,) +
  scale_fill_brewer(palette = "Set2", name = "Провоспалительный\n медиатор") +
  #scale_shape_manual(values=c(1:14), name = "Клеточная линия") +
  scale_x_discrete(limits=c("K", "K_K", "LPS", "K_LPS", "LPS_K", "LPS_LPS"), labels = c("Контроль", "Контроль-Контроль", "LPS", "Контроль-LPS", "LPS-Контроль", "LPS-LPS")) +
  scale_y_log10(breaks=10^(0:3)) +
  xlab(label = "Тип воздействия") +
  ylab(label = "Концентрация медиатора") +
  geom_point(aes(shape = "выброс"),  alpha = 0) +
  geom_point(aes(shape = "медиана"),  alpha = 0) +
  geom_point(aes(shape = "среднее арифметическое"),  alpha = 0) +
  guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
  ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
  theme_AO +
  theme(legend.key.size = unit(1, "cm")) 

dev.off()


Plot_File_Name <- file.path(PathToResults, "Влияние LPS на накопление провоспалительных медиаторов_28-09-2021_4.jpg")
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

ggplot(Input_Inflammatory_Data, aes(x = mediator, y = concentration, fill = type)) +
  geom_boxplot(outlier.size = 3, position = position_dodge(width = 0.9)) +
  stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE,) +
  scale_fill_brewer(palette = "Set3", name = "Тип воздействия", limits=c("K", "K_K", "LPS", "K_LPS", "LPS_K", "LPS_LPS"), labels = c("Контроль", "Контроль-Контроль", "LPS", "Контроль-LPS", "LPS-Контроль", "LPS-LPS")) +
  scale_y_log10(breaks=10^(0:3)) +
  xlab(label = "Провоспалительный медиатор") +
  ylab(label = "Концентрация медиатора") +
  geom_point(aes(shape = "выброс"),  alpha = 0) +
  geom_point(aes(shape = "медиана"),  alpha = 0) +
  geom_point(aes(shape = "среднее арифметическое"),  alpha = 0) +
  guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
  ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
  theme_AO +
  theme(legend.key.size = unit(1, "cm"))


dev.off()

########### Проверка на аутлайеры

Outlier_Table <- data.frame()

Loop_Flag <- TRUE

Loop_counter <- 0

rows_number <- nrow(Input_Inflammatory_Data)

Input_Inflammatory_Data$"ID" <- 1:nrow(Input_Inflammatory_Data)

Old_Input_Inflammatory_Data <- Input_Inflammatory_Data

Old_Input_Inflammatory_Data$"Otlier_Flag" <- NA

while (Loop_Flag)
    {
     
     fit <- lm(concentration ~ mediator + stimulation, data = Input_Inflammatory_Data)
     
     temp_Outlier_Results <- outlierTest(fit)

     temp_Outlier_Table <- data.frame(
                                      "Проверяемый_показатель" = "concentration",
                                      "Значение" = Input_Inflammatory_Data[as.numeric(names(temp_Outlier_Results$rstudent)), "concentration"],
                                      "Номер строки в таблице данных" = as.numeric(names(temp_Outlier_Results$rstudent)),
                                      "Номер строки в исходной таблице данных" = Input_Inflammatory_Data$ID[as.numeric(names(temp_Outlier_Results$rstudent))],
                                      "нескорректированное p-значение" = temp_Outlier_Results$p,
                                      "p-значение Бонферрони" = temp_Outlier_Results$bonf.p,
                                      "Достоверность" = temp_Outlier_Results$signif,
                                      "Критическая точка" = temp_Outlier_Results$cutoff
                                    )
     
     if(!is.na(unique(temp_Outlier_Table$Значение)[1]))
      {
       Outlier_Table <- rbind.data.frame(Outlier_Table, temp_Outlier_Table)
       temp_IDs <- Input_Inflammatory_Data$ID[temp_Outlier_Table$Номер.строки.в.таблице.данных]
       Input_Inflammatory_Data <- Input_Inflammatory_Data[-c(temp_Outlier_Table$Номер.строки.в.таблице.данных),]
      
       Old_Input_Inflammatory_Data$Otlier_Flag[which(Old_Input_Inflammatory_Data$ID %in% temp_IDs)] <- 1
       Loop_counter <- Loop_counter + 1
       message("Итерация № ", Loop_counter, " : удалено ", nrow(temp_Outlier_Table), " строк. Осталось ", nrow(Input_Inflammatory_Data), " строк из ", rows_number)
      } else
           {
            message("Итерация № ", Loop_counter, " : выбросов нет.")
            Loop_Flag <- FALSE
           }
}

Outlier_Table$Значение %in% Old_Input_Inflammatory_Data$concentration[Old_Input_Inflammatory_Data$Otlier_Flag == 1]
Old_Input_Inflammatory_Data$concentration[Old_Input_Inflammatory_Data$Otlier_Flag == 1] %in% Outlier_Table$Значение

Old_Input_Inflammatory_Data$ID <- NULL

#### Сохранение результатов  в Excel файл 

List_Of_Result_Table <- list(Outlier_Table, Old_Input_Inflammatory_Data)
names(List_Of_Result_Table) <- c("Таблица выбросов", "Исходная таблица данных")
Table_ResultFileName <- file.path(PathToResults, "Наличие аутлайеров в данных о провоспалительных медиаторах_и_исходная_таблица_14-09-2021.xlsx")
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)


#### Удаление аутлайеров по тесту Бонферонни

nrow(Input_Inflammatory_Data)
nrow(temp_Outlier_Table)

Input_Inflammatory_Data_wo_outliers <- Input_Inflammatory_Data[-c(temp_Outlier_Table$Номер.строки.в.таблице.данных),]

nrow(Input_Inflammatory_Data_wo_outliers)

####### Визуализация исходных данных данных без выбросов

Plot_File_Name <- file.path(PathToResults, ".jpg")
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

ggplot(Input_Inflammatory_Data_wo_outliers, aes(x = stimulation, y = concentration, colour = mediator)) +
  geom_boxplot(outlier.size = 3) +
  scale_colour_brewer(palette = "Set1", name = "Провоспалительный\n медиатор") +
  #scale_shape_manual(values=c(1:14), name = "Клеточная линия") +
  scale_x_discrete(limits=c("no", "first", "second", "both")) +
  xlab(label = "Стимуляция") +
  ylab(label = "Концентрация медиатора") +
  theme_AO

dev.off()

#### Удаление аутлайеров по графику "Ящик с усами" (не подходит, т.к. теряются медиаторы)

Input_Inflammatory_Data_wo_outliers <- Input_Inflammatory_Data[Input_Inflammatory_Data$concentration >= 5500,]


######################### Визуализация данных по каждому медиатору  и каждой клеточной линии в отдельности

unique(Input_Inflammatory_Data$mediator)

for (cell_line_counter in seq_along(unique(Input_Inflammatory_Data$cell_line)))
  {
   for (mediator_counter in seq_along(unique(Input_Inflammatory_Data$mediator)))
     {
      temp_mediator_line_Data <- subset.data.frame(Input_Inflammatory_Data, Input_Inflammatory_Data$mediator == unique(Input_Inflammatory_Data$mediator)[mediator_counter] & Input_Inflammatory_Data$cell_line == unique(Input_Inflammatory_Data$cell_line)[cell_line_counter])
   
      File_Name <- paste("Провоспаление_", unique(Input_Inflammatory_Data$mediator)[mediator_counter], "_", unique(Input_Inflammatory_Data$cell_line)[cell_line_counter],  ".jpg", sep = "")
      Plot_File_Name <- file.path(PathToResults, File_Name)
      jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")
   
      Figure <- ggplot(temp_mediator_line_Data, aes(x = stimulation, y = concentration)) +
                geom_boxplot(outlier.size = 3) +
                stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 1), show.legend = FALSE,) +
                scale_x_discrete(limits = c("no", "first", "second", "both"), labels = c("Контроль", "Tолько 1-ая\n стимуляция", "Только 2-ая\n стимуляция", "Обе\n стимуляции")) +
                xlab(label = "Стимуляция") +
                ylab(label = "Концентрация медиатора") +
                ggtitle(paste("Провоспаление_", unique(Input_Inflammatory_Data$mediator)[mediator_counter], "_", unique(Input_Inflammatory_Data$cell_line)[cell_line_counter], sep = "")) +
                theme_AO
      
     print(Figure)
      
     dev.off() 
   }
 }


############ Тест Манна-Уитни для контрольного уровня и первой стимуляции по медиаторам и клеточным линиям

names(Input_Inflammatory_Data)
unique(Input_Inflammatory_Data$type)

K_LPS_Data <- subset.data.frame(Input_Inflammatory_Data, Input_Inflammatory_Data$type %in% c("K", "LPS"))

names(K_LPS_Data)

xtabs(~ type + stimulation, data = K_LPS_Data)

wilcox.test(K_LPS_Data$concentration[which(K_LPS_Data$type == "K")],
            K_LPS_Data$concentration[which(K_LPS_Data$type == "LPS")]
            )

Diff_Inflam_Table <- data.frame()

for (mediator_counter in unique(K_LPS_Data$mediator))
  {
   for (cell_line_counter in unique(K_LPS_Data$cell_line))
     {
      temp_Inflam_Data <- subset.data.frame(K_LPS_Data, K_LPS_Data$mediator == mediator_counter & K_LPS_Data$cell_line == cell_line_counter)
      
      if(length(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K"]) > 1 & length(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS"] > 1))
       {
        
        
        
        temp_wilcoxon_results <- wilcox.test(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K"],
                                           temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS"],
                                           alternative = "less"
                                           )
      
      #   temp_t_results <- t.test(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K"],
      #                          temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS"]
      #                          )
      #   
      #  Norm_Distrb_Control <- rnorm(n = 100, mean = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K"], na.rm = TRUE) , sd = sd(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K"], na.rm = TRUE))
      #  Norm_Distrb_Control <- Norm_Distrb_Control[Norm_Distrb_Control > 0]
      # 
      # #Norm_test_results_Control <- lillie.test(Norm_Distrb_Control)  
      #   
      # Norm_Distrb_LPS <- rnorm(n = 100, mean = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS"], na.rm = TRUE) , sd = sd(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS"], na.rm = TRUE))
      # Norm_Distrb_LPS <- Norm_Distrb_LPS[Norm_Distrb_LPS > 0]
      # 
      # #Norm_test_results_LPS <- lillie.test(Norm_Distrb_LPS)
      # 
      # Norm_t_results <- t.test(Norm_Distrb_Control, Norm_Distrb_LPS)
      # 
      # Norm_one_way_t_results <- t.test(Norm_Distrb_Control, Norm_Distrb_LPS, alternative = "less")
      
      temp_Diff_Inflam_Table <- data.frame("Mediator" = mediator_counter,
                                           "Cell_line" = cell_line_counter,
                                           "mean_control" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K"], na.rm = TRUE),
                                           "mean_LPS" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS"], na.rm = TRUE),
                                           "Diff_means" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS"], na.rm = TRUE) - mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K"], na.rm = TRUE),
                                           #"t_test_p-value" = temp_t_results$p.value,
                                           "mann-whitney_test_p-value" = temp_wilcoxon_results$p.value,
                                           #"adjusted_t_test_p-value" = Norm_t_results$p.value,
                                           #"adjusted_one_way_t_test_p-value" = Norm_one_way_t_results$p.value,
                                           "result" = ifelse(temp_wilcoxon_results$p.value < 0.05, "Различия есть", "Различий нет")
                                          ) 
        
      } else
           {
             temp_Diff_Inflam_Table <- data.frame("Mediator" = mediator_counter,
                                                  "Cell_line" = cell_line_counter,
                                                  "mean_control" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K"], na.rm = TRUE),
                                                  "mean_LPS" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS"], na.rm = TRUE),
                                                  "Diff_means" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS"], na.rm = TRUE) - mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K"], na.rm = TRUE),
                                                  #"t_test_p-value" = "недостаточно данных",
                                                  "mann-whitney_test_p-value" = "недостаточно данных",
                                                  #"adjusted_t_test_p-value" = "недостаточно данных",
                                                  #"adjusted_one_way_t_test_p-value" = "недостаточно данных",
                                                  "result" = "недостаточно данных"
             )  
           }
      
      
      
      Diff_Inflam_Table <- rbind.data.frame(Diff_Inflam_Table, temp_Diff_Inflam_Table)
        
      # ggplot(temp_Inflam_Data, aes(x = type, y = concentration)) +
      #   geom_boxplot(outlier.size = 3) +
      #   stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 1), show.legend = FALSE,) +
      #   #scale_x_discrete(limits = c("no", "first", "second", "both"), labels = c("Контроль", "Tолько 1-ая\n стимуляция", "Только 2-ая\n стимуляция", "Обе\n стимуляции")) +
      #   xlab(label = "Тип эксперимента") +
      #   ylab(label = "Концентрация медиатора") +
      #   ggtitle(paste("Провоспаление_", mediator_counter, "_", cell_line_counter, sep = "")) +
      #   theme_AO
      
     }
  }

######## Вариация расчетов, когда значимыми считаются все различия

Diff_Inflam_Table_Sig <- Diff_Inflam_Table

######## Изучение данных с значимой разницей

Diff_Inflam_Table_Sig <- subset.data.frame(Diff_Inflam_Table, Diff_Inflam_Table$result == "Различия есть")
Diff_Inflam_Table_Sig


############### Оценка влияния индивидуальных мутаций на скорость поглощения кисролода и на провоспалительный ответ

####### Проверка наименования линий в таблице о провоспалении
Diff_Inflam_Table_Sig$Cell_line <- as.character(Diff_Inflam_Table_Sig$Cell_line)
unique(Diff_Inflam_Table_Sig$Cell_line)

####### Проверка наименования линий в таблице о дыхании

unique(Input_Data_Native$Line)

####### Унификация наименования мутаций

Diff_Inflam_Table_Sig$Cell_line <- str_remove_all(Diff_Inflam_Table_Sig$Cell_line, "-")

Input_Data_Native$Line <- str_remove_all(Input_Data_Native$Line, "-")
Input_Data_Native$Line <- str_remove_all(Input_Data_Native$Line, "Native ")

Input_Inflammatory_Data$cell_line <- str_remove_all(Input_Inflammatory_Data$cell_line, "-")
Input_Inflammatory_Data$cell_line <- str_remove_all(Input_Inflammatory_Data$cell_line, "Native ")

##### Проверка соответствия названий линий в таблицах

setdiff(unique(Diff_Inflam_Table_Sig$Cell_line), Input_Data_Native$Line)
setdiff(unique(Input_Inflammatory_Data$cell_line), Input_Data_Native$Line)

####### Поиск уникальных мутаций для клеточных линий

Diff_Inflam_Table_Sig$"mutations" <- NA 
Diff_Inflam_Table_Sig$"unique_line_mutations" <- NA
Diff_Inflam_Table_Sig$"unique_group_mutations" <- NA

unique(Input_Data_Native$Line)

for(cell_line_counter in seq_along(unique(Diff_Inflam_Table_Sig$Cell_line)))
  {
   if(unique(Diff_Inflam_Table_Sig$Cell_line)[cell_line_counter] %in% Input_Data_Native$Line)
    { 
     line_mutations <- Input_Data_Native$Mutation_name[Input_Data_Native$Line == unique(Diff_Inflam_Table_Sig$Cell_line)[cell_line_counter] & Input_Data_Native$Heteroplasia_value > 5.0]
     
     another_all_mutations <- unique(Input_Data_Native$Mutation_name[Input_Data_Native$Line != unique(Diff_Inflam_Table_Sig$Cell_line)[cell_line_counter] & Input_Data_Native$Heteroplasia_value > 5.0])
     another_group_mutations <- unique(Input_Data_Native$Mutation_name[!(Input_Data_Native$Line %in% unique(Diff_Inflam_Table_Sig$Cell_line)) & Input_Data_Native$Heteroplasia_value > 5.0])
     unique_line_Mutations <- setdiff(line_mutations, another_all_mutations)
     unique_group_Mutations <- setdiff(line_mutations, another_group_mutations)
     Diff_Inflam_Table_Sig$mutations[which(Diff_Inflam_Table_Sig$Cell_line == unique(Diff_Inflam_Table_Sig$Cell_line)[cell_line_counter])] <- paste(line_mutations, collapse = " + ")
     Diff_Inflam_Table_Sig$unique_line_mutations[which(Diff_Inflam_Table_Sig$Cell_line == unique(Diff_Inflam_Table_Sig$Cell_line)[cell_line_counter])] <- ifelse(length(unique_line_Mutations) >= 1, unique_line_Mutations, "уникальных мутаций в линии нет")
     Diff_Inflam_Table_Sig$unique_group_mutations[which(Diff_Inflam_Table_Sig$Cell_line == unique(Diff_Inflam_Table_Sig$Cell_line)[cell_line_counter])] <- ifelse(length(unique_group_Mutations) >= 1, unique_group_Mutations, "уникальных мутаций в группе линий нет")
    } else
         {
           Diff_Inflam_Table_Sig$mutations[which(Diff_Inflam_Table_Sig$Cell_line == unique(Diff_Inflam_Table_Sig$Cell_line)[cell_line_counter])] <- "нет данных"
           Diff_Inflam_Table_Sig$unique_line_mutations[which(Diff_Inflam_Table_Sig$Cell_line == unique(Diff_Inflam_Table_Sig$Cell_line)[cell_line_counter])] <- "нет данных"
           Diff_Inflam_Table_Sig$unique_group_mutations[which(Diff_Inflam_Table_Sig$Cell_line == unique(Diff_Inflam_Table_Sig$Cell_line)[cell_line_counter])] <- "нет данных" 
         }
     }

Diff_Inflam_Table_Sig



####### Сохранение длинного формата таблицы результатов в Еxcel файл

List_Of_Result_Table <- list(Diff_Inflam_Table,
                             Diff_Inflam_Table_Sig
                             )
names(List_Of_Result_Table) <- c("Сравнения всех линий",
                                 "Выявление мутаций"
                                 )
Table_ResultFileName <- file.path(PathToResults, "Сравнительный_анализ_Контроль-LPS_14-10-2021.xlsx")
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

####### Есть ли вообще уникальные мутации в клеточных линиях

Unique_Mutations_in_Lines <- data.frame()

for(line_counter in seq_along(unique(Input_Data_Native$Line)))
 {
  Line_Data <- subset.data.frame(Input_Data_Native, Input_Data_Native$Line == unique(Input_Data_Native$Line)[line_counter] & Input_Data_Native$Heteroplasia_value > 5)
  Another_Line_Data <- subset.data.frame(Input_Data_Native, Input_Data_Native$Line != unique(Input_Data_Native$Line)[line_counter] & Input_Data_Native$Heteroplasia_value > 5)
  
  line_mutations <- unique(Line_Data$Mutation_name)[order(unique(Line_Data$Mutation_name))]
  another_all_mutations <- unique(Another_Line_Data$Mutation_name)[order(unique(Another_Line_Data$Mutation_name))]
  
  Unique_Mutations <- setdiff(line_mutations, another_all_mutations)
  
  temp_Unique_Mutations_in_Lines <- data.frame("Line" = unique(Input_Data_Native$Line)[line_counter],
                                               "Unique_Mutations" = ifelse(length(Unique_Mutations >= 1), paste(Unique_Mutations, collapse = " + "), "уникальных мутаций нет")
                                              )
  Unique_Mutations_in_Lines <- rbind.data.frame(Unique_Mutations_in_Lines, temp_Unique_Mutations_in_Lines)
 }

Unique_Mutations_in_Lines

####### Проверка наименования линий в таблице о провоспалении

unique(Input_Inflammatory_Data$cell_line)

####### Проверка наименования линий в таблице о различных мутациях

unique(Diff_Mutations_Table$Original_Line_Name)
unique(Diff_Mutations_Table$Derivative_Line_Name)

####### Унификация наименования мутаций

Input_Inflammatory_Data$cell_line <- str_remove_all(Input_Inflammatory_Data$cell_line, "-")
Diff_Mutations_Table$Original_Line_Name <- str_remove_all(Diff_Mutations_Table$Original_Line_Name, "-")
Diff_Mutations_Table$Derivative_Line_Name <- str_remove_all(Diff_Mutations_Table$Derivative_Line_Name, "-")
Diff_Mutations_Table$Original_Line_Name <- str_remove_all(Diff_Mutations_Table$Original_Line_Name, "Native ")
Diff_Mutations_Table$Derivative_Line_Name <- str_remove_all(Diff_Mutations_Table$Derivative_Line_Name, "Native ")


##### Проверка соответствия названий линий в таблицах

setdiff(unique(Input_Inflammatory_Data$cell_line), union(Diff_Mutations_Table$Original_Line_Name, Diff_Mutations_Table$Derivative_Line_Name))
Fail_Lines_Name <- setdiff(union(Diff_Mutations_Table$Original_Line_Name, Diff_Mutations_Table$Derivative_Line_Name), unique(Input_Inflammatory_Data$cell_line))
Fail_Lines_Name

##### Исключение линий, которых нет в провоспалительной таблице из таблицы различий

Fail_Lines_Numbers <- which(Diff_Mutations_Table$Original_Line_Name %in% Fail_Lines_Name | Diff_Mutations_Table$Derivative_Line_Name %in% Fail_Lines_Name)
Diff_Mutations_Table[Fail_Lines_Numbers,]
Diff_Mutations_Table[-Fail_Lines_Numbers,]

New_Diff_Mutations_Table <- Diff_Mutations_Table[-Fail_Lines_Numbers,]
New_Diff_Mutations_Table

####### Основные расчеты

###### Определение количества уникальных мутаций

names(New_Diff_Mutations_Table)[which(names(New_Diff_Mutations_Table) == "Names_Diff_mutations_Plus")] <- "Names_Diff_mutations"

Mutations_list <- unique(unlist(strsplit(New_Diff_Mutations_Table$Names_Diff_mutations, split = " + ", fixed = TRUE, )))

New_Diff_Mutations_Table$"Original_Line_Control" <- NA
New_Diff_Mutations_Table$"Derivative_Line_Control" <- NA
New_Diff_Mutations_Table$"Diff_Control" <- NA
New_Diff_Mutations_Table$"Original_Line_LPS" <- NA
New_Diff_Mutations_Table$"Derivative_Line_LPS" <- NA
New_Diff_Mutations_Table$"Diff_LPS" <- NA

Mutations_Table <- data.frame()

Line_Mediator_List <- list()


for (mediator_counter in seq_along(unique(Input_Inflammatory_Data$mediator))) 
  {
   New_Diff_Mutations_Table <- Diff_Mutations_Table[-Fail_Lines_Numbers,]
   
   names(New_Diff_Mutations_Table)[which(names(New_Diff_Mutations_Table) == "Names_Diff_mutations_Plus")] <- "Names_Diff_mutations"
   
   Mutations_list <- unique(unlist(strsplit(New_Diff_Mutations_Table$Names_Diff_mutations, split = " + ", fixed = TRUE, )))
   
   New_Diff_Mutations_Table$"Original_Line_Control" <- NA
   New_Diff_Mutations_Table$"Derivative_Line_Control" <- NA
   New_Diff_Mutations_Table$"Diff_Control" <- NA
   New_Diff_Mutations_Table$"Original_Line_LPS" <- NA
   New_Diff_Mutations_Table$"Derivative_Line_LPS" <- NA
   New_Diff_Mutations_Table$"Diff_LPS" <- NA
  
   for (rows_counter_1 in 1:nrow(New_Diff_Mutations_Table))
     {
      Original_Line_SubSet <- subset.data.frame(Input_Inflammatory_Data, Input_Inflammatory_Data$cell_line == New_Diff_Mutations_Table$Original_Line_Name[rows_counter_1] & Input_Inflammatory_Data$mediator == unique(Input_Inflammatory_Data$mediator)[mediator_counter])
      Derivative_Line_SubSet <- subset.data.frame(Input_Inflammatory_Data, Input_Inflammatory_Data$cell_line == New_Diff_Mutations_Table$Derivative_Line_Name[rows_counter_1] & Input_Inflammatory_Data$mediator == unique(Input_Inflammatory_Data$mediator)[mediator_counter])

      New_Diff_Mutations_Table$Original_Line_Control[rows_counter_1] <- mean(Original_Line_SubSet$concentration[Original_Line_SubSet$type == "K"], na.rm = TRUE)
      New_Diff_Mutations_Table$Derivative_Line_Control[rows_counter_1] <- mean(Derivative_Line_SubSet$concentration[Derivative_Line_SubSet$type == "K"], na.rm = TRUE)
      New_Diff_Mutations_Table$Diff_Control[rows_counter_1] <- New_Diff_Mutations_Table$Derivative_Line_Control[rows_counter_1] - New_Diff_Mutations_Table$Original_Line_Control[rows_counter_1]
  
      New_Diff_Mutations_Table$Original_Line_LPS[rows_counter_1] <- mean(Original_Line_SubSet$concentration[Original_Line_SubSet$type == "LPS"], na.rm = TRUE)
      New_Diff_Mutations_Table$Derivative_Line_LPS[rows_counter_1] <- mean(Derivative_Line_SubSet$concentration[Derivative_Line_SubSet$type == "LPS"], na.rm = TRUE)
      New_Diff_Mutations_Table$Diff_LPS[rows_counter_1] <- New_Diff_Mutations_Table$Derivative_Line_LPS[rows_counter_1] - New_Diff_Mutations_Table$Original_Line_LPS[rows_counter_1]
   }
  
   New_Diff_Mutations_Table$"Check" <- NA
  
   while(NA %in% New_Diff_Mutations_Table$Check) 
      {
       for(rows_counter in 1:nrow(New_Diff_Mutations_Table))
         {
          temp_Diff_Mutations <- unlist(strsplit(New_Diff_Mutations_Table$Names_Diff_mutations[rows_counter], split = " + ", fixed = TRUE, ))
          if(length(temp_Diff_Mutations) == 0)
           {
            New_Diff_Mutations_Table$Check[rows_counter] <- 0
           } else if(length(temp_Diff_Mutations) == 1)
                   {
                    temp_Mutations_Table <- data.frame(
                                                       Mutation_name = temp_Diff_Mutations,
                                                       Original_Line_Name = New_Diff_Mutations_Table$Original_Line_Name[rows_counter],
                                                       Derivative_Line_Name = New_Diff_Mutations_Table$Derivative_Line_Name[rows_counter],
                                                       Number_Diff_mutations = New_Diff_Mutations_Table$Number_Diff_mutations[rows_counter],
                                                       Mediator_name = unique(Input_Inflammatory_Data$mediator)[mediator_counter], 
                                                       Diff_value_Control = round(New_Diff_Mutations_Table$Diff_Control[rows_counter], digits = 1),
                                                       Relative_Diff_value_Control = round(New_Diff_Mutations_Table$Diff_Control[rows_counter] / New_Diff_Mutations_Table$Original_Line_Control[rows_counter], digits = 3),
                                                       Diff_value_LPS = round(New_Diff_Mutations_Table$Diff_LPS[rows_counter], digits = 1),
                                                       Relative_Diff_value_LPS = round(New_Diff_Mutations_Table$Diff_LPS[rows_counter] / New_Diff_Mutations_Table$Original_Line_LPS[rows_counter], digits = 3),
                                                       stringsAsFactors = FALSE
                                                      )
                    Mutations_Table <- rbind.data.frame(Mutations_Table, temp_Mutations_Table)
                    New_Diff_Mutations_Table$Check[rows_counter] <- 1
                    rm(temp_Mutations_Table)
                   } else if(length(temp_Diff_Mutations) > 1)
                           {
                            ###### Расчет средних значений для каждой из уже известной мутации
                            
                            temp_temp_Mutation_Table <- data.frame()       
                     
                            for (Mutation_names_counter in unique(Mutations_Table$Mutation_name))
                              {
                               temp_values_Table <- data.frame(Mutation_name = as.character(Mutation_names_counter),
                                                               Average_value_Control = round(mean(Mutations_Table$Diff_value_Control[Mutations_Table$Mutation_name == Mutation_names_counter]), digits = 1),
                                                               Average_value_LPS = round(mean(Mutations_Table$Diff_value_LPS[Mutations_Table$Mutation_name == Mutation_names_counter]), digits = 1),
                                                               stringsAsFactors = FALSE
                                                                      )
                               temp_temp_Mutation_Table <- rbind.data.frame(temp_temp_Mutation_Table, temp_values_Table)
                              }
                            
                            ###### Исключение известных мутаций из клеточных линий
                            
                            Know_Mutations <- temp_temp_Mutation_Table$Mutation_name
                            
                            Unknow_Mutations <- setdiff(temp_Diff_Mutations, Know_Mutations)
                            
                            if(length(Unknow_Mutations) == 1)
                            {
                              Impact_Know_Mutations_value_Control <- sum(temp_temp_Mutation_Table$Average_value_Control[which(temp_temp_Mutation_Table$Mutation_name %in% Know_Mutations)])
                              Impact_Know_Mutations_value_LPS <- sum(temp_temp_Mutation_Table$Average_value_LPS[which(temp_temp_Mutation_Table$Mutation_name %in% Know_Mutations)])
                              
                              Diff_value_Control <- New_Diff_Mutations_Table$Diff_Control[rows_counter] - Impact_Know_Mutations_value_Control
                              Diff_value_LPS <- New_Diff_Mutations_Table$Diff_LPS[rows_counter] - Impact_Know_Mutations_value_LPS
                              
                              Relative_Diff_value_Control <- ifelse(is.null(New_Diff_Mutations_Table$Original_Line_Control[rows_counter]), NA, Diff_value_Control / New_Diff_Mutations_Table$Original_Line_Control[rows_counter])
                              Relative_Diff_value_LPS <- ifelse(is.null(New_Diff_Mutations_Table$Original_Line_LPS[rows_counter]), NA, Diff_value_LPS / New_Diff_Mutations_Table$Original_Line_LPS[rows_counter])
                              
                              
                              temp_Mutations_Table <- data.frame(
                                                                 Mutation_name = Unknow_Mutations,
                                                                 Original_Line_Name = New_Diff_Mutations_Table$Original_Line_Name[rows_counter],
                                                                 Derivative_Line_Name = New_Diff_Mutations_Table$Derivative_Line_Name[rows_counter],
                                                                 Number_Diff_mutations = New_Diff_Mutations_Table$Number_Diff_mutations[rows_counter],
                                                                 Mediator_name = unique(Input_Inflammatory_Data$mediator)[mediator_counter], 
                                                                 Diff_value_Control = Diff_value_Control,
                                                                 Relative_Diff_value_Control = round(Relative_Diff_value_Control, digits = 3),
                                                                 Diff_value_LPS = round(Diff_value_LPS, digits = 1),
                                                                 Relative_Diff_value_LPS = round(Relative_Diff_value_LPS, digits = 3),
                                                                 stringsAsFactors = FALSE
                                                                 )
                              Mutations_Table <- rbind.data.frame(Mutations_Table, temp_Mutations_Table)
                              New_Diff_Mutations_Table$Check[rows_counter] <- 1
                              rm(temp_temp_Mutation_Table, temp_Mutations_Table)
                            } else
                                 {
                                  New_Diff_Mutations_Table$Check[rows_counter] <- 0  
                                 }
                            
                           }
         }
  }
  
   names(New_Diff_Mutations_Table) <- c("Исходная клеточная линия",
        "Производная клеточная линия",
        "Число отличающихся мутаций",
        "Число появившихся мутаций",
        "Число элиминированных мутаций",
        "Отличающиеся (появившиеся) мутации",
        "Элиминированные мутации",
        "Контроль в исходной линии",
        "Контроль в производной линии",
        "Различия в контроле",
        "LPS в исходной линии",
        "LPS в производной линии",
        "Различия при 1-ой стимуляции",
        "Роль мутаций установлена"
        )
   
   Line_Mediator_List[[unique(Input_Inflammatory_Data$mediator)[mediator_counter]]] <- New_Diff_Mutations_Table
   
  }

names(Mutations_Table) <- c("Мутация",
        "Исходная клеточная линия",
        "Производная клеточная линия",
        "Число отличающихся мутаций",
        "Провоспалительный медиатор",
        "Абсолютные различия в контроле",
        "Относительные различия в контроле",
        "Абсолютные различия при 1-ой стимуляции",
        "Относительные различия при 1-ой стимуляции"
        )
   
   Line_Mediator_List[["Вклад отдельных мутаций"]] <- Mutations_Table

   ########## Сохранение таблиц сходств и различий, а также таблицы влияния мутаций  в Excel файл
   
   Results_List <- Line_Mediator_List
   ResultFileName <- file.path(PathToResults, paste("Таблицы вклада митохондриальных_мутаций в провоспалительный эффект_14-10-2021.xlsx", sep = ''))
   save_to_excel_file_by_XLConnect(ResultFileName, Results_List)
   
   
####### Выявление мутаций ассоциированных с повышенным провоспалительным ответом

###### Сравнительный анализ различий в значимых линиях и медиаторах  для первой стимуляции
   
Kruskal_results_Table_LPS <- data.frame()
Kruskal_posthoc_Table_LPS <- data.frame()
   
Set_Mediators <- as.character(unique(Diff_Inflam_Table_Sig$Mediator))

for(mediator_counter in seq_along(Set_Mediators))
  {
   temp_Lines <- Diff_Inflam_Table_Sig$Cell_line[which(Diff_Inflam_Table_Sig$Mediator == Set_Mediators[mediator_counter])]
   
   if(length(temp_Lines) > 1)
    {
   
     temp_SubSet_Inflam_Data <- subset.data.frame(Input_Inflammatory_Data, Input_Inflammatory_Data$cell_line %in% temp_Lines & Input_Inflammatory_Data$mediator == Set_Mediators[mediator_counter] & Input_Inflammatory_Data$type == "LPS")
   
     kruskal_results <- kruskal.test(concentration ~ cell_line, data = temp_SubSet_Inflam_Data)
   
     temp_Kruskal_results_Table <- data.frame(
      "Медиатор" = Set_Mediators[mediator_counter],
      "Название теста" = kruskal_results$method,
      "Название набора данных" = kruskal_results$data.name,
      "Значение статистики" = kruskal_results$statistic,
      "р-значение" = kruskal_results$p.value,
      "Результат" = ifelse(kruskal_results$p.value > 0.05, "Различий нет", "Различия есть")
       )
   Kruskal_results_Table_LPS <- rbind.data.frame(temp_Kruskal_results_Table, Kruskal_results_Table_LPS)
   
   if(temp_Kruskal_results_Table$р.значение < 0.05)
   {
     kruskal_posthoc_results <- kruskalmc(concentration ~ cell_line, data = temp_SubSet_Inflam_Data, probs = 0.05)
     
     names(kruskal_posthoc_results$dif.com) <- c("Уровень различий", "Критический уровень различий", "Достоверность различий")
     
     temp_kruskal_posthoc_results_Table <- data.frame(
       "Медиатор" = Set_Mediators[mediator_counter],
       "Пары сравнения" = row.names(kruskal_posthoc_results$dif.com),
       kruskal_posthoc_results$dif.com[,1:3],
       "p-значение" = kruskal_posthoc_results$signif.level
     )
     Kruskal_posthoc_Table_LPS <- rbind.data.frame(temp_kruskal_posthoc_results_Table, Kruskal_posthoc_Table_LPS)
   }
  }

 }  

Kruskal_results_Table_LPS
Kruskal_posthoc_Table_LPS

###### Добавление мутаций в различающихся линиях

Kruskal_posthoc_Table_LPS$"Отличающиеся мутации" <- NA

for (rows_counter in 1:nrow(Kruskal_posthoc_Table_LPS))
  {
   One_Line <- str_remove(string = str_extract(string = Kruskal_posthoc_Table_LPS$Пары.сравнения[rows_counter], pattern = ".+-"), pattern = "-")
   Another_Line <- str_remove(string = str_extract(string = Kruskal_posthoc_Table_LPS$Пары.сравнения[rows_counter], pattern = "-.+"), pattern = "-")
   
   One_Line_Mutations <- Input_Data_Native$Mutation_name[Input_Data_Native$Line == One_Line]
   Another_Mutations <- Input_Data_Native$Mutation_name[Input_Data_Native$Line == Another_Line]
   
   Diff_Lines_Mutations <- union(setdiff(One_Line_Mutations, Another_Mutations), setdiff(Another_Mutations, One_Line_Mutations))
   
   Kruskal_posthoc_Table_LPS$`Отличающиеся мутации`[rows_counter] <- paste(Diff_Lines_Mutations, collapse = " + ")
   
  }

unique(Kruskal_posthoc_Table_LPS$`Отличающиеся мутации`)

Mutations_list <- unique(unlist(strsplit(Kruskal_posthoc_Table_LPS$`Отличающиеся мутации`, split = " + ", fixed = TRUE, )))

Mutations_list

length(Mutations_list)
  
###### Сравнительный анализ различий в значимых линиях и медиаторах  для контроля

Kruskal_results_Table_K <- data.frame()
Kruskal_posthoc_Table_K <- data.frame()

Set_Mediators <- as.character(unique(Diff_Inflam_Table_Sig$Mediator))

for(mediator_counter in seq_along(Set_Mediators))
{
  temp_Lines <- Diff_Inflam_Table_Sig$Cell_line[which(Diff_Inflam_Table_Sig$Mediator == Set_Mediators[mediator_counter])]
  
  if(length(temp_Lines) > 1)
  {
    
    temp_SubSet_Inflam_Data <- subset.data.frame(Input_Inflammatory_Data, Input_Inflammatory_Data$cell_line %in% temp_Lines & Input_Inflammatory_Data$mediator == Set_Mediators[mediator_counter] & Input_Inflammatory_Data$type == "K")
    
    kruskal_results <- kruskal.test(concentration ~ cell_line, data = temp_SubSet_Inflam_Data)
    
    temp_Kruskal_results_Table <- data.frame(
      "Медиатор" = Set_Mediators[mediator_counter],
      "Название теста" = kruskal_results$method,
      "Название набора данных" = kruskal_results$data.name,
      "Значение статистики" = kruskal_results$statistic,
      "р-значение" = kruskal_results$p.value,
      "Результат" = ifelse(kruskal_results$p.value > 0.05, "Различий нет", "Различия есть")
    )
    Kruskal_results_Table_K <- rbind.data.frame(temp_Kruskal_results_Table, Kruskal_results_Table_K)
    
    if(temp_Kruskal_results_Table$р.значение < 0.05)
    {
      kruskal_posthoc_results <- kruskalmc(concentration ~ cell_line, data = temp_SubSet_Inflam_Data, probs = 0.05)
      
      names(kruskal_posthoc_results$dif.com) <- c("Уровень различий", "Критический уровень различий", "Достоверность различий")
      
      temp_kruskal_posthoc_results_Table <- data.frame(
        "Медиатор" = Set_Mediators[mediator_counter],
        "Пары сравнения" = row.names(kruskal_posthoc_results$dif.com),
        kruskal_posthoc_results$dif.com[,1:3],
        "p-значение" = kruskal_posthoc_results$signif.level
      )
      Kruskal_posthoc_Table_K <- rbind.data.frame(temp_kruskal_posthoc_results_Table, Kruskal_posthoc_Table_K)
    }
  }
  
}   

Kruskal_results_Table_K
Kruskal_posthoc_Table_K

####### Сохранение длинного формата таблицы результатов в Еxcel файл

List_Of_Result_Table <- list(Kruskal_results_Table_LPS, Kruskal_posthoc_Table_LPS, Kruskal_results_Table_K)
names(List_Of_Result_Table) <- c("медиаторы по линиям при LPS", "отличающиеся медиаторы при LPS", "медиаторы по линиям в контроле")
Table_ResultFileName <- file.path(PathToResults, "Исследование_различий_клеточных_линий_первой_стимуляции_14-10-2021.xlsx")
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

########################### Исcледование толерантности ####################


names(Input_Inflammatory_Data)
unique(Input_Inflammatory_Data$type)

head(Input_Inflammatory_Data)

for(rows_counter in 1:nrow(Input_Inflammatory_Data))
  {

   temp_SubSet <- subset.data.frame(Input_Inflammatory_Data, Input_Inflammatory_Data$cell_line == Input_Inflammatory_Data[rows_counter, "cell_line"] &
                        Input_Inflammatory_Data$mediator == Input_Inflammatory_Data[rows_counter, "mediator"] &
                        !is.na(Input_Inflammatory_Data$concentration)
                     )
   melt
  }


####### исследование толерантности на основании классификации от Никиты Никифорова от 29-11-2021



LPS_LPS_Data <- subset.data.frame(Input_Inflammatory_Data, Input_Inflammatory_Data$type %in% c("K_LPS", "LPS_LPS"))

names(LPS_LPS_Data)

xtabs(~ type + stimulation, data = LPS_LPS_Data)


Diff_Inflam_Table_LPS_LPS <- data.frame()

for (mediator_counter in unique(LPS_LPS_Data$mediator))
{
  for (cell_line_counter in unique(LPS_LPS_Data$cell_line))
  {
    temp_Inflam_Data <- subset.data.frame(LPS_LPS_Data, LPS_LPS_Data$mediator == mediator_counter & LPS_LPS_Data$cell_line == cell_line_counter)
    
    if(length(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K_LPS"]) > 1 & length(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS_LPS"] > 1))
    {
      temp_wilcoxon_results_hyper <- wilcox.test(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K_LPS"],
                                                 temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS_LPS"],
                                                 alternative = "less"
                                                 )
      
      temp_wilcoxon_results_responce <- wilcox.test(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K_LPS"],
                                                 temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS_LPS"],
                                                 alternative = "two.sided"
                                                 )
      
      temp_wilcoxon_results_tolerance <- wilcox.test(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K_LPS"],
                                                    temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS_LPS"],
                                                    alternative = "greater"
      )
      
      temp_t_results <- t.test(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K_LPS"],
                               temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS_LPS"]
      )
      
      Norm_Distrb_Control <- rnorm(n = 100, mean = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K_LPS"], na.rm = TRUE) , sd = sd(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K_LPS"], na.rm = TRUE))
      Norm_Distrb_Control <- Norm_Distrb_Control[Norm_Distrb_Control >= 0]
      
      #Norm_test_results_Control <- lillie.test(Norm_Distrb_Control)  
      
      Norm_Distrb_LPS <- rnorm(n = 100, mean = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS_LPS"], na.rm = TRUE) , sd = sd(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS_LPS"], na.rm = TRUE))
      Norm_Distrb_LPS <- Norm_Distrb_LPS[Norm_Distrb_LPS >= 0]
      
      #Norm_test_results_LPS <- lillie.test(Norm_Distrb_LPS)
      
      Norm_t_results <- t.test(Norm_Distrb_Control, Norm_Distrb_LPS)
      
      Norm_one_way_t_results <- t.test(Norm_Distrb_Control, Norm_Distrb_LPS, alternative = "less")
      
      temp_Diff_Inflam_Table_LPS_LPS <- data.frame("Mediator" = mediator_counter,
                                           "Cell_line" = cell_line_counter,
                                           "mean_K_LPS" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K_LPS"], na.rm = TRUE),
                                           "mean_LPS_LPS" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS_LPS"], na.rm = TRUE),
                                           "Diff_means" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS_LPS"], na.rm = TRUE) - mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K_LPS"], na.rm = TRUE),
                                           "t_test_p-value" = temp_t_results$p.value,
                                           "hyper_mann-whitney_test_p-value" = temp_wilcoxon_results_hyper$p.value,
                                           "responce_mann-whitney_test_p-value" = temp_wilcoxon_results_responce$p.value,
                                           "tolerance_mann-whitney_test_p-value" = temp_wilcoxon_results_tolerance$p.value,
                                           "adjusted_t_test_p-value" = Norm_t_results$p.value,
                                           "adjusted_one_way_t_test_p-value" = Norm_one_way_t_results$p.value,
                                           "hyper_result" = ifelse(temp_wilcoxon_results_hyper$p.value < 0.05, "hyperresponder", "нет различий"),
                                           "responce_result" = ifelse(temp_wilcoxon_results_responce$p.value < 0.05, "есть различия", "non-responder"),
                                           "tolerance_result" = ifelse(temp_wilcoxon_results_tolerance$p.value < 0.05, "tolerance", "нет различий")
      ) 
      
    } else
    {
      temp_Diff_Inflam_Table_LPS_LPS <- data.frame("Mediator" = mediator_counter,
                                           "Cell_line" = cell_line_counter,
                                           "mean_K_LPS" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K_LPS"], na.rm = TRUE),
                                           "mean_LPS_LPS" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS_LPS"], na.rm = TRUE),
                                           "Diff_means" = mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "LPS_LPS"], na.rm = TRUE) - mean(temp_Inflam_Data$concentration[temp_Inflam_Data$type == "K_LPS"], na.rm = TRUE),
                                           "t_test_p-value" = "недостаточно данных",
                                           "hyper_mann-whitney_test_p-value" = "недостаточно данных",
                                           "responce_mann-whitney_test_p-value" = "недостаточно данных",
                                           "tolerance_mann-whitney_test_p-value" = "недостаточно данных",
                                           "adjusted_t_test_p-value" = "недостаточно данных",
                                           "adjusted_one_way_t_test_p-value" = "недостаточно данных",
                                           "hyper_result" = "недостаточно данных",
                                           "responce_result" = "недостаточно данных",
                                           "tolerance_result" = "недостаточно данных"
      )  
    }
    
    
    
    Diff_Inflam_Table_LPS_LPS <- rbind.data.frame(Diff_Inflam_Table_LPS_LPS, temp_Diff_Inflam_Table_LPS_LPS)
    
    
    
    # ggplot(temp_Inflam_Data, aes(x = type, y = concentration)) +
    #   geom_boxplot(outlier.size = 3) +
    #   stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 1), show.legend = FALSE,) +
    #   #scale_x_discrete(limits = c("no", "first", "second", "both"), labels = c("Контроль", "Tолько 1-ая\n стимуляция", "Только 2-ая\n стимуляция", "Обе\n стимуляции")) +
    #   xlab(label = "Тип эксперимента") +
    #   ylab(label = "Концентрация медиатора") +
    #   ggtitle(paste("Провоспаление_", mediator_counter, "_", cell_line_counter, sep = "")) +
    #   theme_AO
    
  }
}

Diff_Inflam_Table_LPS_LPS

######## Изучение данных с значимой разницей

Diff_Inflam_Table_Sig_LPS_LPS <- subset.data.frame(Diff_Inflam_Table_LPS_LPS,
                                                   Diff_Inflam_Table_LPS_LPS$hyper_result == "hyperresponder" |
                                                   Diff_Inflam_Table_LPS_LPS$responce_result == "non-responder" |
                                                   Diff_Inflam_Table_LPS_LPS$tolerance_result == "tolerance"
                                                   )
Diff_Inflam_Table_Sig_LPS_LPS

############### Оценка влияния индивидуальных мутаций на скорость поглощения кисролода и на провоспалительный ответ

####### Проверка наименования линий в таблице о провоспалении
Diff_Inflam_Table_Sig_LPS_LPS$Cell_line <- as.character(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)
unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)

####### Проверка наименования линий в таблице о дыхании

unique(Input_Data_Native$Line)

##### Проверка соответствия названий линий в таблицах

setdiff(unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line), Input_Data_Native$Line)
setdiff(unique(Input_Inflammatory_Data$cell_line), Input_Data_Native$Line)

####### Поиск уникальных мутаций для клеточных линий

Diff_Inflam_Table_Sig_LPS_LPS$"mutations" <- NA 
Diff_Inflam_Table_Sig_LPS_LPS$"unique_line_mutations" <- NA
Diff_Inflam_Table_Sig_LPS_LPS$"unique_group_mutations" <- NA

unique(Input_Data_Native$Line)

for(cell_line_counter in seq_along(unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)))
{
  if(unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)[cell_line_counter] %in% Input_Data_Native$Line)
  { 
    line_mutations <- Input_Data_Native$Mutation_name[Input_Data_Native$Line == unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)[cell_line_counter] & Input_Data_Native$Heteroplasia_value > 5.0]
    
    another_all_mutations <- unique(Input_Data_Native$Mutation_name[Input_Data_Native$Line != unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)[cell_line_counter] & Input_Data_Native$Heteroplasia_value > 5.0])
    another_group_mutations <- unique(Input_Data_Native$Mutation_name[!(Input_Data_Native$Line %in% unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)) & Input_Data_Native$Heteroplasia_value > 5.0])
    unique_line_Mutations <- setdiff(line_mutations, another_all_mutations)
    unique_group_Mutations <- setdiff(line_mutations, another_group_mutations)
    Diff_Inflam_Table_Sig_LPS_LPS$mutations[which(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line == unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)[cell_line_counter])] <- paste(line_mutations, collapse = " + ")
    Diff_Inflam_Table_Sig_LPS_LPS$unique_line_mutations[which(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line == unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)[cell_line_counter])] <- ifelse(length(unique_line_Mutations) >= 1, unique_line_Mutations, "уникальных мутаций в линии нет")
    Diff_Inflam_Table_Sig_LPS_LPS$unique_group_mutations[which(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line == unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)[cell_line_counter])] <- ifelse(length(unique_group_Mutations) >= 1, unique_group_Mutations, "уникальных мутаций в группе линий нет")
  } else
  {
    Diff_Inflam_Table_Sig_LPS_LPS$mutations[which(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line == unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)[cell_line_counter])] <- "нет данных"
    Diff_Inflam_Table_Sig_LPS_LPS$unique_line_mutations[which(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line == unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)[cell_line_counter])] <- "нет данных"
    Diff_Inflam_Table_Sig_LPS_LPS$unique_group_mutations[which(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line == unique(Diff_Inflam_Table_Sig_LPS_LPS$Cell_line)[cell_line_counter])] <- "нет данных" 
  }
}

Diff_Inflam_Table_Sig_LPS_LPS

####### Сохранение длинного формата таблицы результатов в Еxcel файл

List_Of_Result_Table <- list(Diff_Inflam_Table_LPS_LPS, Diff_Inflam_Table_Sig_LPS_LPS)
names(List_Of_Result_Table) <- c("Все линии и медиаторы", "Отличающиеся линии и медиаторы")
Table_ResultFileName <- file.path(PathToResults, "Сравнительный_анализ_K-LPS_LPS-LPS_14-10-2021.xlsx")
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)


####### тест Фридмана  для толерантности

Input_Inflammatory_Data$type <- str_replace_all(string = Input_Inflammatory_Data$type, pattern = "-", replacement = "_")

Input_Inflammatory_Data$type <- as.factor(Input_Inflammatory_Data$type)
Input_Inflammatory_Data$mediator <- as.factor(Input_Inflammatory_Data$mediator)

Input_Inflammatory_Data_Means <- aggregate(Input_Inflammatory_Data$concentration, by = list(type = Input_Inflammatory_Data$type,  mediator = Input_Inflammatory_Data$mediator), FUN = mean_without_NA)

Input_Inflammatory_Data_Means

names(Input_Inflammatory_Data_Means)[3] <- "mean_concentration"

Friedman_results <- friedman.test(mean_concentration ~ type | mediator, data = Input_Inflammatory_Data_Means)

class(Friedman_results)
str(Friedman_results)

Friedman_results

Friedman_results_Table <- data.frame(
  "Название теста" = Friedman_results$method,
  "Название набора данных" = Friedman_results$data.name,
  "Значение статистики" = Friedman_results$statistic,
  "р-значение" = Friedman_results$p.value,
  "Результат" = ifelse(Friedman_results$p.value > 0.05, "Различий нет", "Различия есть")
)
Friedman_results_Table

Friedman_post_hoc_Table <- friedmanmc(Input_Inflammatory_Data_Means$mean_concentration, groups = Input_Inflammatory_Data_Means$type, blocks = Input_Inflammatory_Data_Means$mediator)

names(Friedman_post_hoc_Table$dif.com) <- c("Уровень различий", "Критический уровень различий", "Достоверность различий")

Friedman_post_hoc_Table <- data.frame("Пары сравнения" = row.names(Friedman_post_hoc_Table$dif.com),
                                      Friedman_post_hoc_Table$dif.com[,1:3],
                                      "p-значение" = Friedman_post_hoc_Table$p.value
)

Friedman_post_hoc_Table


######## Тест Крускала-Уолиса



Kruskal_results_Table <- data.frame()
Kruskal_posthoc_Table <- data.frame()

for(mediator_counter in seq_along(unique(Input_Inflammatory_Data$mediator)))
{
  Input_Data_Mediator <- subset.data.frame(Input_Inflammatory_Data, Input_Inflammatory_Data$mediator == unique(Input_Inflammatory_Data$mediator)[mediator_counter])
  
  kruskal_results <- kruskal.test(concentration ~ type, data = Input_Data_Mediator)
  
  temp_Kruskal_results_Table <- data.frame(
    "Медиатор" = unique(Input_Inflammatory_Data$mediator)[mediator_counter],
    "Название теста" = kruskal_results$method,
    "Название набора данных" = kruskal_results$data.name,
    "Значение статистики" = kruskal_results$statistic,
    "р-значение" = kruskal_results$p.value,
    "Результат" = ifelse(kruskal_results$p.value > 0.05, "Различий нет", "Различия есть")
  )
  Kruskal_results_Table <- rbind.data.frame(temp_Kruskal_results_Table, Kruskal_results_Table)
  
  if(temp_Kruskal_results_Table$р.значение < 0.05)
  {
    kruskal_posthoc_results <- kruskalmc(concentration ~ type, data = Input_Data_Mediator, probs = 0.05)
    
    #kruskal_posthoc_results <- posthoc.kruskal.conover.test(concentration ~ type, data = Input_Data_Mediator, probs = 0.05)
    
    names(kruskal_posthoc_results$dif.com) <- c("Уровень различий", "Критический уровень различий", "Достоверность различий")
    
    temp_kruskal_posthoc_results_Table <- data.frame(
      "Медиатор" = unique(Input_Inflammatory_Data$mediator)[mediator_counter],
      "Пары сравнения" = row.names(kruskal_posthoc_results$dif.com),
      kruskal_posthoc_results$dif.com[,1:3],
      "p-значение" = kruskal_posthoc_results$signif.level
    )
    Kruskal_posthoc_Table <- rbind.data.frame(temp_kruskal_posthoc_results_Table, Kruskal_posthoc_Table)
  }
}

row.names(Kruskal_results_Table) <- NULL
row.names(Kruskal_posthoc_Table) <- NULL

Kruskal_results_Table
Kruskal_posthoc_Table

Kruskal_posthoc_Table[which(Kruskal_posthoc_Table$Достоверность.различий),]

##### kruskal.test; to reorder factor levels see relevel; for other functions about median multiple comparison see posthoc.kruskal.conover.test, kruskal

Kruskal_posthoc_Table$Пары.сравнения <- as.character(Kruskal_posthoc_Table$Пары.сравнения)

for (row_counter in 1:nrow(Kruskal_posthoc_Table))
  {
   Kruskal_posthoc_Table[row_counter, c("First_type", "Second_Type")] <- unlist(strsplit(Kruskal_posthoc_Table$Пары.сравнения[row_counter], split = "-", fixed = TRUE))
   
   Kruskal_posthoc_Table[row_counter, c("First_mean", "Second_mean")] <- c(mean(Input_Inflammatory_Data$concentration[Input_Inflammatory_Data$mediator == Kruskal_posthoc_Table$Медиатор[row_counter] & Input_Inflammatory_Data$type == Kruskal_posthoc_Table$First_type[row_counter]], na.rm = TRUE), 
                                                                           mean(Input_Inflammatory_Data$concentration[Input_Inflammatory_Data$mediator == Kruskal_posthoc_Table$Медиатор[row_counter] & Input_Inflammatory_Data$type == Kruskal_posthoc_Table$Second_Type[row_counter]], na.rm = TRUE)
                                                                           ) 
  }

Kruskal_posthoc_Table$"First_stimulation" <- NA
Kruskal_posthoc_Table$"Second_stimulation" <- NA

index <- match(Kruskal_posthoc_Table$First_type, Input_Inflammatory_Data$type)
Kruskal_posthoc_Table[!is.na(index), "First_stimulation"] <- as.character(Input_Inflammatory_Data[index[!is.na(index)], "stimulation"])

index <- match(Kruskal_posthoc_Table$Second_Type, Input_Inflammatory_Data$type)
Kruskal_posthoc_Table[!is.na(index), "Second_stimulation"] <- as.character(Input_Inflammatory_Data[index[!is.na(index)], "stimulation"])

Kruskal_posthoc_Table$"time_effect" <- NA

Kruskal_posthoc_Table$"time_effect"[which(Kruskal_posthoc_Table$Достоверность.различий & Kruskal_posthoc_Table$First_stimulation == Kruskal_posthoc_Table$Second_stimulation)] <- 1


Kruskal_posthoc_Table[which(Kruskal_posthoc_Table$time_effect == 1),]

####### Сохранение длинного формата таблицы результатов в Еxcel файл

List_Of_Result_Table <- list(Kruskal_posthoc_Table)
names(List_Of_Result_Table) <- c("Полное сравнение")
Table_ResultFileName <- file.path(PathToResults, "Полный_сравнительный_анализ_накопления_медиаторов_14-10-2021.xlsx")
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

############################################################### РАБОТА ПОСЛЕ 08-11-2021 ##########################################################################################

######### ФОРМУЛИРОВКА ЗАДАЧ:

######### 1. определить критические размеры эффектов
######### 2. на основании размеров эффектов определить критическое p-value и выполнить расчеты
######### 3. найти мутации характерные для трех групп (гипер, нон и толеранс)
######### корреляция вклада мутаций с гетероплазмией

######## ОКОНЧАНИЕ ФОРМУЛИРОВКИ

names(Input_Inflammatory_Data)
names(Ranking_Data)
names(Input_Data)

sort(unique(Input_Inflammatory_Data$cell_line))
sort(unique(Ranking_Data$Line))
sort(unique(Input_Data_Native$Line))
sort(unique(Input_Data_Mitophagia$cell_line))

####### Унификация наименования мутаций

Input_Inflammatory_Data$cell_line <- str_remove_all(Input_Inflammatory_Data$cell_line, "-")
Ranking_Data$Line <- str_remove_all(Ranking_Data$Line, "-")
Input_Data_Native$Line <- str_remove_all(Input_Data_Native$Line, "-")
Input_Data_Mitophagia$cell_line <- str_remove_all(Input_Data_Mitophagia$cell_line, "-")

Input_Inflammatory_Data$cell_line <- str_remove_all(Input_Inflammatory_Data$cell_line, "Native ")
Ranking_Data$Line <- str_remove_all(Ranking_Data$Line, "Native ")
Input_Data_Native$Line <- str_remove_all(Input_Data_Native$Line, "Native ")
Input_Data_Mitophagia$cell_line <- str_remove_all(Input_Data_Mitophagia$cell_line, "Native ")

setdiff(Input_Data_Native$Line, Input_Inflammatory_Data$cell_line)
setdiff(Input_Inflammatory_Data$cell_line, Input_Data_Native$Line)
setdiff(Input_Inflammatory_Data$cell_line, Ranking_Data$cell_line)
setdiff(Input_Data_Native$Line, Ranking_Data$Line)
setdiff(Ranking_Data$Line, Input_Data_Native$Line)
setdiff(Input_Data_Mitophagia$cell_line, Input_Data_Native$Line)
setdiff(Input_Data_Native$Line, Input_Data_Mitophagia$cell_line)

###### Добавление рейтинговых данных в таблицу с гетероплазмией

Input_Data_Native_Rank <- merge.data.frame(Input_Data_Native, Ranking_Data, by = "Line")

Input_Data_Native_Rank

unique(Input_Data_Native_Rank$Line)

Input_Data_Native[Input_Data_Native$Line %in% setdiff(Input_Data_Native$Line, Ranking_Data$Line), ]

###### Добавление рейтинговых данных в таблицу с митофагией

names(Ranking_Data)[which(names(Ranking_Data) == "cell_line")] <- "Line"
names(Input_Data_Mitophagia)[which(names(Input_Data_Mitophagia) == "cell_line")] <- "Line"

Input_Data_Native_Mitophagia <- merge.data.frame(Input_Data_Native, Input_Data_Mitophagia, by = "Line")

Input_Data_Native_Mitophagia

unique(Input_Data_Native_Mitophagia$Line)

Input_Data_Native[Input_Data_Native$Line %in% setdiff(Input_Data_Native$Line, Input_Data_Mitophagia$cell_line), ]

###### Добавление данных о концентрациях провосполительных медиаторах в таблицу с гетероплазмией

names(Input_Inflammatory_Data)[which(names(Input_Inflammatory_Data) == "cell_line")] <- "Line"

Input_Data_Native_Inflammatory <- merge.data.frame(Input_Data_Native, Input_Inflammatory_Data, by = "Line")

Input_Data_Native_Inflammatory

unique(Input_Data_Native_Inflammatory$Line)

Input_Data_Native_Inflammatory[Input_Data_Native_Inflammatory$Line %in% setdiff(Input_Data_Native$Line, Input_Data_Native_Inflammatory$cell_line), ]

##### Дальше работаем и строим систему линейных уравнений только с Input_Data_Native

Hetero_Table <- aggregate(Input_Data_Native_Rank$Heteroplasia_value, by = list(Input_Data_Native_Rank$Line, Input_Data_Native_Rank$Mutation_name), FUN = mean)

str(Hetero_Table)

names(Hetero_Table) <- c("Line", "Mutation_name", "Heteroplasmia_value")

Hetero_Matrix <- as.matrix(xtabs(Heteroplasmia_value ~ Line + Mutation_name, data = Hetero_Table))

Hetero_Matrix

######## работа с провоспалением для каждого медиатора в отдельности

names(Input_Data_Native_Inflammatory)

Input_Data_Native_Inflammatory_first <- subset.data.frame(Input_Data_Native_Inflammatory, Input_Data_Native_Inflammatory$type == "LPS")

Inf_Table <- aggregate(Input_Data_Native_Inflammatory_first[,c("concentration")], by = list(Input_Data_Native_Inflammatory_first$Line, Input_Data_Native_Inflammatory_first$mediator), FUN = mean)

Inf_Table_melt <- melt(Input_Data_Native_Inflammatory_first, measure.vars = c("concentration"), id.vars = c( "Line", "mediator"))

Inf_Table_cast <- cast(Input_Data_Native_Inflammatory_first, Line ~ mediator, mean_without_NA)


######## работа с рейтингами

Ranks_Table <- aggregate(Input_Data_Native_Rank[,c("Inflammation_Rank_First", "Inflammation_Rank_Second")], by = list(Input_Data_Native_Rank$Line), FUN = mean)
Ranks_Table
names(Ranks_Table)[1] <- "Line"

####### работа с митофагией

setdiff(Input_Data_Mitophagia$Line, Input_Data_Native$Line)
setdiff(Input_Data_Native$Line, Input_Data_Mitophagia$Line)

unique(Input_Data_Mitophagia$Line)

names(Input_Data_Mitophagia)

Mitos_Table <- aggregate(Input_Data_Native_Mitophagia[,c("Mitochondrial.mass.(FACS)__Median", "Mitophagy.rate_FCCP_Mean")], by = list(Input_Data_Native_Mitophagia$Line), FUN = mean)
Mitos_Table
names(Mitos_Table)[1] <- "Line"

####### Вычисление определителя матрицы гетероплазмии  и доказательство линейной независимости СЛАУ
dim(Hetero_Matrix)
det(Hetero_Matrix)

Hetero_SuperMatrix <- matrix(data = Hetero_Matrix, nrow = 12, ncol = 10, dimnames = list(row.names(Hetero_Matrix)))
Hetero_SuperMatrix
class(Hetero_SuperMatrix)

det(Hetero_SuperMatrix)

Null_Matrix <- matrix(data = rep(0, times = length(nrow(Hetero_Matrix))), nrow = nrow(Hetero_Matrix), ncol = 1)
Null_Matrix
Solve_Undependence <- solve(Hetero_Matrix, Null_Matrix)

Solve_Undependence


##### Работа с неквадратными матрицами

dim(Hetero_Matrix)
det(Hetero_Matrix)

Hetero_Matrix

##### Расчитать число комбинаций по 2 или по 3 строки из 12 , генерировать эти пары и последовательно их удалять, проверяя устойчивость решения

list_rows_pairs <- combn(sort(unique(as.character(row.names(Hetero_Matrix)))), 2, simplify = FALSE)

####### работа с рангами

Rank_Matrix <- as.matrix(xtabs(Inflammation_Rank_Second ~ Line, data = Ranks_Table))
colnames(Rank_Matrix) <- "Inflammation_Rank_Sum"
Rank_Matrix
Hetero_Matrix <- Hetero_Matrix[intersect(row.names(Hetero_Matrix), row.names(Rank_Matrix)), ]

Solve_Inflammation_Rank <- data.frame()

for (list_counter in 1:length(list_rows_pairs))
  {
   temp_Hetero_Matrix <- Hetero_Matrix[which(!(row.names(Hetero_Matrix) %in% unlist(list_rows_pairs[list_counter]))),]
   temp_Rank_Matrix <- Rank_Matrix[which(!(row.names(Rank_Matrix) %in% unlist(list_rows_pairs[list_counter]))),]
   temp_Solve_Inflammation_Rank_Sum <- solve(temp_Hetero_Matrix, temp_Rank_Matrix)
   temp_Solve_Inflammation_Rank_Sum <- as.data.frame(temp_Solve_Inflammation_Rank_Sum)
   temp_Solve_Inflammation_Rank_Sum$"Mutation" <- row.names(temp_Solve_Inflammation_Rank_Sum)
   temp_Solve_Inflammation_Rank_Sum$"Removed_Lines" <- paste(list_rows_pairs[[list_counter]], collapse = " and ") 
   row.names(temp_Solve_Inflammation_Rank_Sum) <- NULL
   Solve_Inflammation_Rank <- rbind.data.frame(Solve_Inflammation_Rank, temp_Solve_Inflammation_Rank_Sum)
  }

Solve_Inflammation_Rank

####### работа с митофагией

Mito_Matrix <- as.matrix(xtabs(Mitos_Table$`Mitophagy.rate_FCCP_Mean` ~ Line, data = Mitos_Table))
colnames(Mito_Matrix) <- "Mitophagy.rate_FCCP_Mean"
Mito_Matrix
Hetero_Matrix <- Hetero_Matrix[intersect(row.names(Hetero_Matrix), row.names(Mito_Matrix)), ]

Hetero_Matrix

Solve_Inflammation_Mito <- data.frame()

for (list_counter in 1:length(list_rows_pairs))
{
  temp_Hetero_Matrix <- Hetero_Matrix[which(!(row.names(Hetero_Matrix) %in% unlist(list_rows_pairs[list_counter]))),]
  temp_Mito_Matrix <- Mito_Matrix[which(!(row.names(Mito_Matrix) %in% unlist(list_rows_pairs[list_counter]))),]
  temp_Solve_Inflammation_Mito_Sum <- solve(temp_Hetero_Matrix, temp_Mito_Matrix)
  temp_Solve_Inflammation_Mito_Sum <- as.data.frame(temp_Solve_Inflammation_Mito_Sum)
  temp_Solve_Inflammation_Mito_Sum$"Mutation" <- row.names(temp_Solve_Inflammation_Mito_Sum)
  temp_Solve_Inflammation_Mito_Sum$"Removed_Lines" <- paste(list_rows_pairs[[list_counter]], collapse = " and ") 
  row.names(temp_Solve_Inflammation_Mito_Sum) <- NULL
  Solve_Inflammation_Mito <- rbind.data.frame(Solve_Inflammation_Mito, temp_Solve_Inflammation_Mito_Sum)
}

Solve_Inflammation_Mito

######## расcчитать p-value

P_value_DataSet <- data.frame()

for(Mutations_counter in 1: length(unique(Solve_Inflammation_Rank$Mutation)))
  {
   temp_DataSet <- subset.data.frame(Solve_Inflammation_Rank, Solve_Inflammation_Rank$Mutation == unique(Solve_Inflammation_Rank$Mutation)[Mutations_counter])
   
   temp_CI <- CI(temp_DataSet$temp_Solve_Inflammation_Rank_Sum)
   
   temp_result_DataSet <- data.frame("Mutation" = unique(Solve_Inflammation_Rank$Mutation)[Mutations_counter],
                                     "p-value" = p.value(temp_DataSet$temp_Solve_Inflammation_Rank_Sum),
                                     "CI_lower" = temp_CI["lower"],
                                     "CI_upper" = temp_CI["upper"],
                                     "Result" = ifelse(p.value(temp_DataSet$temp_Solve_Inflammation_Rank_Sum) < 0.05, "Вклад мутации статистически значимо отличается от нуля", "Вклад мутации статистически НЕ отличается от нуля")
                                    )
   P_value_DataSet <- rbind.data.frame(P_value_DataSet, temp_result_DataSet)
  }

row.names(P_value_DataSet) <- NULL

names(P_value_DataSet) <- c("Мутация", "p-значение", "нижняя_граница_ДИ_для_p-значения", "верхняя_граница_ДИ_для_p-значения", "Результат")

Solve_Mito_Sum_Desc_Stat <- aggregate(temp_Solve_Inflammation_Rank_Sum ~ Mutation, data = Solve_Inflammation_Rank, FUN = get_descr_stat_add)

Solve_Mito_Sum_Desc_Stat <- data.frame("Мутация" = Solve_Mito_Sum_Desc_Stat$Mutation, Solve_Mito_Sum_Desc_Stat$temp_Solve_Inflammation_Rank_Sum)

names(Solve_Mito_Sum_Desc_Stat) <- c(
  "Мутация",
  "Число случаев",
  "Минимальное значение",
  "Макскимальное значение",
  "Размах",
  "Cумма",
  "Среднее",
  "Медиана",
  "Квантиль 25%",
  "Квантиль 75%",
  "Mежквартильный интервал",
  "Нижняя граница выбросов",
  "Дисперсия",
  "Стандартное отклонение",
  "Стандартная ошибка среднего",
  "Нижняя граница 95% доверительного интервала_среднего",
  "Верхняя граница 95% доверительного интервала_среднего",
  "Асимметрия",
  "Эксцесс"
)

Solve_Mito_Sum_Desc_Stat

Solve_Mito_Sum_Desc_Stat <- merge(Solve_Mito_Sum_Desc_Stat, P_value_DataSet, by = "Мутация")

Solve_Mito_Sum_Desc_Stat

###### Собственно сохранение результатов работы с не квадратными матрицами

List_Of_Result_Table <- list(Solve_Mito_Sum_Desc_Stat)
names(List_Of_Result_Table) <- c("Описательная статистика")
Table_ResultFileName <- file.path(PathToResults, paste("Таблица вклада мутаций на основании второй стимуляции", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

################### Работа с суммарными рангами

Rank_Matrix <- as.matrix(xtabs(Inflammation_Rank_Sum ~ Line, data = Ranks_Table))
colnames(Rank_Matrix) <- "Inflammation_Rank_Sum"
Rank_Matrix
Hetero_Matrix <- Hetero_Matrix[intersect(row.names(Hetero_Matrix), row.names(Rank_Matrix)), ]
Solve_Inflammation_Rank_Sum <- solve(Hetero_Matrix, Rank_Matrix)
Solve_Inflammation_Rank_Sum <- as.data.frame(Solve_Inflammation_Rank_Sum)
Solve_Inflammation_Rank_Sum

################### Работа с рангами первой стимуляции

Rank_Matrix <- as.matrix(xtabs(Inflammation_Rank_First ~ Line, data = Ranks_Table))
colnames(Rank_Matrix) <- "Inflammation_Rank_First"
Rank_Matrix
Hetero_Matrix <- Hetero_Matrix[intersect(row.names(Hetero_Matrix), row.names(Rank_Matrix)), ]
Solve_Inflammation_Rank_First <- solve(Hetero_Matrix, Rank_Matrix)
Solve_Inflammation_Rank_First <- as.data.frame(Solve_Inflammation_Rank_First)
Solve_Inflammation_Rank_First

################### Работа с рангами второй стимуляции

Rank_Matrix <- as.matrix(xtabs(Inflammation_Rank_Second ~ Line, data = Ranks_Table))
colnames(Rank_Matrix) <- "Inflammation_Rank_Second"
Rank_Matrix
Hetero_Matrix <- Hetero_Matrix[intersect(row.names(Hetero_Matrix), row.names(Rank_Matrix)), ]
Solve_Inflammation_Rank_Second <- solve(Hetero_Matrix, Rank_Matrix)
Solve_Inflammation_Rank_Second

class(Solve_Inflammation_Rank_Second)

Solve_Inflammation_Rank_Second <- as.data.frame(Solve_Inflammation_Rank_Second)
Solve_Inflammation_Rank_Second

######## объединение полученных результатов

Solve_Inflammation_Rank <- merge(Solve_Inflammation_Rank_Sum, Solve_Inflammation_Rank_First, by = "row.names", all = TRUE)
row.names(Solve_Inflammation_Rank) <- Solve_Inflammation_Rank[,"Row.names"]
Solve_Inflammation_Rank$Row.names <- NULL
Solve_Inflammation_Rank <- merge(Solve_Inflammation_Rank, Solve_Inflammation_Rank_Second, by = "row.names", all = TRUE)
Solve_Inflammation_Rank
names(Solve_Inflammation_Rank)[which(names(Solve_Inflammation_Rank) == "Row.names")] <- "Mutation_name"

###### Сохранение полученных результатов в Еxcel файл

###### Преобразование длинной формы таблицы гетерозиготности в широкую

reshape(Hetero_Table, idvar = "Line", timevar = "Mutation_name", direction = "wide")

xtabs(Heteroplasmia_value ~ Mutation_name + Line, data = Hetero_Table)

Hetero_data_frame <- spread(Hetero_Table, key = Mutation_name, value = Heteroplasmia_value)

Hetero_data_frame[is.na(Hetero_data_frame)] <- 0

class(Hetero_data_frame)

###### Собственно сохранение

List_Of_Result_Table <- list(Hetero_data_frame, Ranks_Table[,-c(which(names(Ranks_Table) == "OK_value"))], Solve_Inflammation_Rank)
names(List_Of_Result_Table) <- c("Матрица_гетероплазмии", "Векторы_рангов", "Векторы_решений_СЛАУ")
Table_ResultFileName <- file.path(PathToResults, "Системы_линейных уравнений_гетероплазмия-провоспаление_с_решениями_12-11-2021.xlsx")
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

############################# ОТВЕТВЛЕНИЕ К СЛУ ГЕТЕРОПЛАЗМИЯ-ДЫХАНИЕ ###############################################

Ranks_Table$"OK_Rank" <- rank(Ranks_Table$OK_value)

##### расчеты по скорости поглощения кислорода

Rank_Matrix <- as.matrix(xtabs(OK_value ~ Line, data = Ranks_Table))
colnames(Rank_Matrix) <- "OK_value"
Rank_Matrix
Hetero_Matrix <- Hetero_Matrix[intersect(row.names(Hetero_Matrix), row.names(Rank_Matrix)), ]
Solve_OK_Rank <- solve(Hetero_Matrix, Rank_Matrix)
Solve_OK_Rank <- as.data.frame(Solve_OK_Rank)
Solve_OK_Rank

##### расчеты по рангам скорости поглощения кислорода

Rank_Matrix <- as.matrix(xtabs(OK_Rank ~ Line, data = Ranks_Table))
colnames(Rank_Matrix) <- "OK_Rank"
Rank_Matrix
Hetero_Matrix <- Hetero_Matrix[intersect(row.names(Hetero_Matrix), row.names(Rank_Matrix)), ]
Solve_Rank_OK_Rank <- solve(Hetero_Matrix, Rank_Matrix)
Solve_Rank_OK_Rank <- as.data.frame(Solve_Rank_OK_Rank)
Solve_Rank_OK_Rank

Solve_OK_Rank <- merge(Solve_OK_Rank, Solve_Rank_OK_Rank, by = "row.names", all = TRUE)
Solve_OK_Rank
names(Solve_OK_Rank)[which(names(Solve_OK_Rank) == "Row.names")] <- "Mutation_name"

###### Сохранение полученных результатов в Еxcel файл

List_Of_Result_Table <- list(Solve_OK_Rank)
names(List_Of_Result_Table) <- c("Векторы_решений_СЛАУ")
PathToResults_OK <- "C:\\Users\\Андрей\\OneDrive\\Documents\\Work\\НИИМЧ\\mtDNA мутации и клеточное дыхание\\Results"
Table_ResultFileName <- file.path(PathToResults_OK, "Системы_линейных уравнений_гетероплазмия-дыхание_с_решениями_12-11-2021.xlsx")
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)


####### работа с провоспалительными медиаторами в отдельности

##### Расчитать число комбинаций по 2 или по 3 строки из 12 , генерировать эти пары и последовательно их удалять, проверяя устойчивость решения

list_rows_pairs <- combn(sort(unique(as.character(row.names(Hetero_Matrix)))), 2, simplify = FALSE)

List_Solve_Inflammation <- list()
List_Solve_Sum_Desc_Stat <- list()

for (mediator_counter in 2:length(names(Inf_Table_cast)))
{

temp_Inf_Matrix <- matrix(data = Inf_Table_cast[, mediator_counter], nrow = nrow(Inf_Table_cast), ncol = 1, dimnames = list(Inf_Table_cast$Line))

Hetero_Matrix <- Hetero_Matrix[intersect(row.names(Hetero_Matrix), row.names(temp_Inf_Matrix)), ]

Hetero_Matrix

Solve_Inflammation <- data.frame()

for (list_counter in 1:length(list_rows_pairs))
{
  temp_Hetero_Matrix <- Hetero_Matrix[which(!(row.names(Hetero_Matrix) %in% unlist(list_rows_pairs[list_counter]))),]
  temp_Matrix <- temp_Inf_Matrix[which(!(row.names(temp_Inf_Matrix) %in% unlist(list_rows_pairs[list_counter]))),]
  temp_Solve_Inflammation_Sum <- solve(temp_Hetero_Matrix, temp_Matrix)
  temp_Solve_Inflammation_Sum <- as.data.frame(temp_Solve_Inflammation_Sum)
  temp_Solve_Inflammation_Sum$"Mutation" <- row.names(temp_Solve_Inflammation_Sum)
  temp_Solve_Inflammation_Sum$"Removed_Lines" <- paste(list_rows_pairs[[list_counter]], collapse = " and ") 
  row.names(temp_Solve_Inflammation_Sum) <- NULL
  Solve_Inflammation <- rbind.data.frame(Solve_Inflammation, temp_Solve_Inflammation_Sum)
}

Solve_Inflammation_Table <- xtabs(temp_Solve_Inflammation_Sum ~ Mutation + Removed_Lines, data = Solve_Inflammation)

List_Solve_Inflammation[[names(Inf_Table_cast)[mediator_counter]]] <- Solve_Inflammation_Table


######## расcчитать p-value

P_value_DataSet <- data.frame()

for(Mutations_counter in 1: length(unique(Solve_Inflammation$Mutation)))
{
  temp_DataSet <- subset.data.frame(Solve_Inflammation, Solve_Inflammation$Mutation == unique(Solve_Inflammation$Mutation)[Mutations_counter])
  
  temp_CI <- CI(temp_DataSet$temp_Solve_Inflammation_Sum)
  
  temp_result_DataSet <- data.frame("Mutation" = unique(Solve_Inflammation$Mutation)[Mutations_counter],
                                    "p-value" = p.value(temp_DataSet$temp_Solve_Inflammation_Sum),
                                    "CI_lower" = temp_CI["lower"],
                                    "CI_upper" = temp_CI["upper"],
                                    "Result" = ifelse(p.value(temp_DataSet$temp_Solve_Inflammation_Sum) < 0.05, "Вклад мутации статистически значимо отличается от нуля", "Вклад мутации статистически НЕ отличается от нуля")
  )
  P_value_DataSet <- rbind.data.frame(P_value_DataSet, temp_result_DataSet)
}

row.names(P_value_DataSet) <- NULL

names(P_value_DataSet) <- c("Мутация", "p-значение", "нижняя_граница_ДИ", "верхняя_граница_ДИ", "Результат")

Solve_Sum_Desc_Stat <- aggregate(temp_Solve_Inflammation_Sum ~ Mutation, data = Solve_Inflammation, FUN = get_descr_stat_add)

Solve_Sum_Desc_Stat <- data.frame("Мутация" = Solve_Sum_Desc_Stat$Mutation, Solve_Sum_Desc_Stat$temp_Solve_Inflammation_Sum)

names(Solve_Sum_Desc_Stat) <- c(
  "Мутация",
  "Число случаев",
  "Минимальное значение",
  "Макскимальное значение",
  "Размах",
  "Cумма",
  "Среднее",
  "Медиана",
  "Квантиль 25%",
  "Квантиль 75%",
  "Mежквартильный интервал",
  "Нижняя граница выбросов",
  "Дисперсия",
  "Стандартное отклонение",
  "Стандартная ошибка среднего",
  "Асимметрия",
  "Эксцесс"
)

Solve_Sum_Desc_Stat

Solve_Sum_Desc_Stat <- merge(Solve_Sum_Desc_Stat, P_value_DataSet, by = "Мутация")

List_Solve_Sum_Desc_Stat[[names(Inf_Table_cast)[mediator_counter]]] <- Solve_Sum_Desc_Stat

message("Медиатор ", names(Inf_Table_cast)[mediator_counter], " обработан.\n")

}

###### Собственно сохранение результатов работы с не квадратными матрицами

Table_ResultFileName <- file.path(PathToResults, paste("Таблица вклада мутаций для каждого медиатора на основании первой стимуляции", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Solve_Sum_Desc_Stat)

########### АНАЛИЗ ОТНОШЕНИЯ РИСКОВ И ШАНСОВ ДЛЯ ПРОВОСПАЛИТЕЛЬНОГО ОТВЕТА проведен в файле s_mtDNA_mutations_28-12-2021.R ################

######################## РАБОТА С ДАННЫМИ С МИТОФАГИЕЙ И С ПРОВОСПАЛЕНИЕМ ######################################

Input_Data
str(Input_Data)

Line_Mitophagy_Inflamation_DataSet <- Input_Data

names(Line_Mitophagy_Inflamation_DataSet) <- c("Line", "Mitophagy", "Inflamantion_type_number", "Inflamantion_rank", "Inflamantion_type")

Line_Mitophagy_Inflamation_DataSet

####### расчет средних

Line_Mitophagy_Inflamation_DataSet$Inflamantion_type <- as.factor(Line_Mitophagy_Inflamation_DataSet$Inflamantion_type)

Input_Data_Means <- aggregate(Line_Mitophagy_Inflamation_DataSet$Mitophagy, by = list(Inflamantion_type = Line_Mitophagy_Inflamation_DataSet$Inflamantion_type), FUN = mean)

Input_Data_Means

names(Input_Data_Means)[2] <- "mean_Mitophagy"

###### Подготовка средних к сохранению и печати

Input_Data_Means$mean_Mitophagy <- round(Input_Data_Means$mean_Mitophagy, digits = 3)

Input_Data_Means

###### конец подготовки средних к сохранению и печати

Kruskal_results_Table <- data.frame()
Kruskal_posthoc_Table <- data.frame()

  kruskal_results <- kruskal.test(Mitophagy ~ Inflamantion_type, data = Line_Mitophagy_Inflamation_DataSet)
  
  temp_Kruskal_results_Table <- data.frame(
    "Название теста" = kruskal_results$method,
    "Название набора данных" = kruskal_results$data.name,
    "Значение статистики" = kruskal_results$statistic,
    "р-значение" = kruskal_results$p.value,
    "Результат" = ifelse(kruskal_results$p.value > 0.05, "Различий нет", "Различия есть")
  )
  Kruskal_results_Table <- rbind.data.frame(temp_Kruskal_results_Table, Kruskal_results_Table)
  
  if(temp_Kruskal_results_Table$р.значение < 0.05)
  {
    kruskal_posthoc_results <- kruskalmc(Mitophagy ~ Inflamantion_type, data = Line_Mitophagy_Inflamation_DataSet, probs = 0.05)
    
    names(kruskal_posthoc_results$dif.com) <- c("Уровень различий", "Критический уровень различий", "Достоверность различий")
    
    temp_kruskal_posthoc_results_Table <- data.frame(
      "Пары сравнения" = row.names(kruskal_posthoc_results$dif.com),
      kruskal_posthoc_results$dif.com[,1:3],
      "p-значение" = kruskal_posthoc_results$signif.level
    )
    Kruskal_posthoc_Table <- rbind.data.frame(temp_kruskal_posthoc_results_Table, Kruskal_posthoc_Table)
  }

row.names(Kruskal_results_Table) <- NULL
row.names(Kruskal_posthoc_Table) <- NULL

Kruskal_results_Table
Kruskal_posthoc_Table


###### Сохранение таблиц в файл Excel

List_Of_Result_Table <- list(Input_Data_Means, Kruskal_results_Table, Kruskal_posthoc_Table)
names(List_Of_Result_Table) <- c("Средний уровень холестерина", "Тест Краскела", "Post-hoc tests")
Table_ResultFileName <- file.path(PathToResults, paste("Сравнительный анализ Средние и тест Краскела-Уолиса для митофагии и типа провоспаления", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

###### Отрисовка графика митофагии относительно типа провоспаления

my_comparisons <- combn(sort(unique(as.character(Line_Mitophagy_Inflamation_DataSet$Inflamantion_type))),2, simplify = FALSE)



my_comparisons <- rev(list(
  c("non-responder", "non-tolerance"),
  c("non-responder", "tolerance")
)
)

Plot_File_Name <- file.path(PathToResults, paste("Box plots of mitophagy by inflamation type of cell line", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

ggplot(Line_Mitophagy_Inflamation_DataSet, aes(x = Inflamantion_type, y = Mitophagy)) +
  geom_boxplot(aes(fill = Inflamantion_type), position = position_dodge(width = 1), outlier.size = 3) +
  stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 1), show.legend = FALSE,) +
  #stat_pvalue_manual(data = df_p_val, label = "p.adj.signif") +
  #stat_compare_means(method =  "wilcox.test", comparisons = my_comparisons, label = "p.format") +
  stat_compare_means(method =  "wilcox.test", comparisons = my_comparisons, label = "p.signif",  bracket.size = 0.5, size = 6) + 
  scale_fill_brewer(palette = "Set3", name = "Inflamantion type") +
  #scale_x_discrete(limits = c("control", "ldl_only", "ANXA1_ldl_knockdown", "F2RL1_ldl_knockdown", "IL15_ldl_knockdown", "PERK_ldl_knockdown", "TSPYL2_ldl_knockdown"), labels = c( "All genes", "All genes", "ANXA1", "F2RL1", "IL15", "PERK", "TSPYL2")) +
  scale_y_continuous(breaks=seq(0, max(Line_Mitophagy_Inflamation_DataSet$Mitophagy), 1)) +
  xlab(label = "Inflamantion type in cell line") +
  ylab(label = "Mitophagy change") +
  geom_point(aes(shape = "mean"),  alpha = 0) +
  geom_point(aes(shape = "median"),  alpha = 0) +
  geom_point(aes(shape = "outlier"),  alpha = 0) +
  guides(fill = FALSE, shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("dodgerblue","black", "black"), shape = c(18, 95, 20), size = 5))) +
  #ggtitle("mitophagy by inflamation type of cell line") +
  theme_AO +
  theme(legend.key.size = unit(1, "cm"))

dev.off()

############## Объединение наборов данных с митофагией и мутациями

#### Унификация названий клеточных линий

unique(Input_Data_Lines_Mutations_Heteroplasmy$Line)
unique(Line_Mitophagy_Inflamation_DataSet$Line)

Input_Data_Lines_Mutations_Heteroplasmy$Line <- str_remove_all(Input_Data_Lines_Mutations_Heteroplasmy$Line, "Native ")
Input_Data_Lines_Mutations_Heteroplasmy$Line <- str_remove_all(Input_Data_Lines_Mutations_Heteroplasmy$Line, "-")
Input_Data_Lines_Mutations_Heteroplasmy$Line[which(Input_Data_Lines_Mutations_Heteroplasmy$Line == "HSMAM1")] <- "MAM1"
Input_Data_Lines_Mutations_Heteroplasmy$Line[which(Input_Data_Lines_Mutations_Heteroplasmy$Line == "HSMAM2")] <- "MAM2"
Input_Data_Lines_Mutations_Heteroplasmy$Line[which(Input_Data_Lines_Mutations_Heteroplasmy$Line == "HSMAM3")] <- "MAM3"
Input_Data_Lines_Mutations_Heteroplasmy$Line[which(Input_Data_Lines_Mutations_Heteroplasmy$Line == "THP1")] <- "THP1"

Line_Mitophagy_Inflamation_DataSet$Line <- str_remove_all(Line_Mitophagy_Inflamation_DataSet$Line, "Native ")
Line_Mitophagy_Inflamation_DataSet$Line <- str_remove_all(Line_Mitophagy_Inflamation_DataSet$Line, "-")
Line_Mitophagy_Inflamation_DataSet$Line[which(Line_Mitophagy_Inflamation_DataSet$Line == "ТНР1")] <- "THP1"

setdiff(Input_Data_Lines_Mutations_Heteroplasmy$Line, Line_Mitophagy_Inflamation_DataSet$Line)
setdiff(Line_Mitophagy_Inflamation_DataSet$Line, Input_Data_Lines_Mutations_Heteroplasmy$Line)

#### Подготовка широкой таблицы с мутациями и гетероплазмией по клеточным линиям

######## подсчет частот мутаций в линиях

Count_Table <- xtabs(Heteroplasia_value ~ Line + Mutation_name, data = Input_Data_Lines_Mutations_Heteroplasmy_base_cons)
Count_Table
class(Count_Table)
str(Count_Table)

Heteroplasmy_dataframe <- as.data.frame.matrix(Count_Table)

Heteroplasmy_dataframe$"Line" <- row.names(Heteroplasmy_dataframe)

str(Heteroplasmy_dataframe)

row.names(Line_Mitophagy_Inflamation_DataSet) <- Line_Mitophagy_Inflamation_DataSet$Line

Line_Mitophagy_Inflamation_whole <- merge(Line_Mitophagy_Inflamation_DataSet, Heteroplasmy_dataframe, by = "Line", all = TRUE) 

Line_Mitophagy_Inflamation_long <- gather(Line_Mitophagy_Inflamation_whole, Mutation_name, Heteroplasmy, a1555g:t3336c, factor_key = TRUE)

####### Сохранение длинной таблицы данных о митофагии, мутациях и гетероплазмии в файл Rdata

Save_FileName <- file.path(PathToSources, paste("Митофагия_мутации_гетероплазмия_13-04-2022.Rdata", sep = ''))
saveRDS(Line_Mitophagy_Inflamation_long, file = Save_FileName)

list.files(PathToSources)

###### Отрисовка графика митофагии относительно названий мутаций

Line_Mitophagy_Inflamation_long <- subset.data.frame(Line_Mitophagy_Inflamation_long, Line_Mitophagy_Inflamation_long$Heteroplasmy > 5)

my_comparisons <- combn(sort(unique(as.character(Line_Mitophagy_Inflamation_long$Mutation_name))),2, simplify = FALSE)

my_comparisons <- rev(list(
  c("a1555g", "g15059a"),
  c("c3256t", "g15059a"),
  c("c5178a", "g15059a"),
  c("g12315a", "g15059a"),
  c("g13513a", "g15059a"),
  c("g14846a", "g15059a"),
  c("t3336c", "g15059a")
)
)

Plot_File_Name <- file.path(PathToResults, paste("Box plots of mitophagy by mt mutation", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

ggplot(Line_Mitophagy_Inflamation_long, aes(x = Mutation_name, y = Mitophagy)) +
  geom_boxplot(aes(fill = Mutation_name), position = position_dodge(width = 1), outlier.size = 3) +
  stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 1), show.legend = FALSE,) +
  #stat_pvalue_manual(data = df_p_val, label = "p.adj.signif") +
  #stat_compare_means(method =  "wilcox.test", comparisons = my_comparisons, label = "p.format") +
  stat_compare_means(method =  "wilcox.test", comparisons = my_comparisons, label = "p.signif",  bracket.size = 0.5, size = 6) + 
  scale_fill_brewer(palette = "Set3", name = "Mutations") +
  #scale_x_discrete(limits = c("control", "ldl_only", "ANXA1_ldl_knockdown", "F2RL1_ldl_knockdown", "IL15_ldl_knockdown", "PERK_ldl_knockdown", "TSPYL2_ldl_knockdown"), labels = c( "All genes", "All genes", "ANXA1", "F2RL1", "IL15", "PERK", "TSPYL2")) +
  scale_y_continuous(breaks=seq(0, max(Line_Mitophagy_Inflamation_DataSet$Mitophagy), 0.1)) +
  xlab(label = "Mutations") +
  ylab(label = "Mitophagy change") +
  geom_point(aes(shape = "mean"),  alpha = 0) +
  geom_point(aes(shape = "median"),  alpha = 0) +
  geom_point(aes(shape = "outlier"),  alpha = 0) +
  guides(fill = FALSE, shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("dodgerblue","black", "black"), shape = c(18, 95, 20), size = 5))) +
  #ggtitle("mitophagy by inflamation type of cell line") +
  theme_AO +
  theme(legend.key.size = unit(1, "cm"))

dev.off()

####### расчет средних

Line_Mitophagy_Inflamation_long$Mutation_name <- as.factor(Line_Mitophagy_Inflamation_long$Mutation_name)

Input_Data_Means <- aggregate(Line_Mitophagy_Inflamation_long$Mitophagy, by = list(Inflamantion_type = Line_Mitophagy_Inflamation_long$Mutation_name), FUN = mean_without_NA)

Input_Data_Means

names(Input_Data_Means)[2] <- "mean_Mitophagy"

###### Подготовка средних к сохранению и печати

Input_Data_Means$mean_Mitophagy <- round(Input_Data_Means$mean_Mitophagy, digits = 3)

Input_Data_Means

###### конец подготовки средних к сохранению и печати

Kruskal_results_Table <- data.frame()
Kruskal_posthoc_Table <- data.frame()

kruskal_results <- kruskal.test(Mitophagy ~ Mutation_name, data = Line_Mitophagy_Inflamation_long)

temp_Kruskal_results_Table <- data.frame(
  "Название теста" = kruskal_results$method,
  "Название набора данных" = kruskal_results$data.name,
  "Значение статистики" = kruskal_results$statistic,
  "р-значение" = kruskal_results$p.value,
  "Результат" = ifelse(kruskal_results$p.value > 0.05, "Различий нет", "Различия есть")
)
Kruskal_results_Table <- rbind.data.frame(temp_Kruskal_results_Table, Kruskal_results_Table)

if(temp_Kruskal_results_Table$р.значение < 0.05)
{
  kruskal_posthoc_results <- kruskalmc(Mitophagy ~ Mutation_name, data = Line_Mitophagy_Inflamation_long, probs = 0.05)
  
  names(kruskal_posthoc_results$dif.com) <- c("Уровень различий", "Критический уровень различий", "Достоверность различий")
  
  temp_kruskal_posthoc_results_Table <- data.frame(
    "Пары сравнения" = row.names(kruskal_posthoc_results$dif.com),
    kruskal_posthoc_results$dif.com[,1:3],
    "p-значение" = kruskal_posthoc_results$signif.level
  )
  Kruskal_posthoc_Table <- rbind.data.frame(temp_kruskal_posthoc_results_Table, Kruskal_posthoc_Table)
}

row.names(Kruskal_results_Table) <- NULL
row.names(Kruskal_posthoc_Table) <- NULL

Kruskal_results_Table
Kruskal_posthoc_Table

Kruskal_posthoc_Table[Kruskal_posthoc_Table$Достоверность.различий, ]


###### Сохранение таблиц в файл Excel

List_Of_Result_Table <- list(Input_Data_Means, Kruskal_results_Table, Kruskal_posthoc_Table)
names(List_Of_Result_Table) <- c("Средний уровень холестерина", "Тест Краскела", "Post-hoc tests")
Table_ResultFileName <- file.path(PathToResults, paste("Сравнительный анализ Средние и тест Краскела-Уолиса для митофагии и митмутаций", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

########### Средние гетероплазмии по мутациям

aggregate(Heteroplasmy ~ Mutation_name, data = Line_Mitophagy_Inflamation_long[Line_Mitophagy_Inflamation_long$Heteroplasmy > 5,], FUN = mean_without_NA)

aggregate(Heteroplasmy ~ Mutation_name, data = Line_Mitophagy_Inflamation_long[Line_Mitophagy_Inflamation_long$Heteroplasmy > 5,], FUN = get_descr_stat_add)

aggregate(Heteroplasmy ~ Mutation_name + Inflamantion_type, data = Line_Mitophagy_Inflamation_long, FUN = mean_without_NA)

###### Открытие файла с таблицей данных о митофагии, мутациях и гетероплазмии

list.files(PathToSources)

InputData_FileName <- file.path(PathToSources, "Митофагия_мутации_гетероплазмия_13-04-2022.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
Input_Data_Mitophagia_Line <- ListOfInputData_FileSheets

names(Input_Data_Mitophagia_Line)
str(Input_Data_Mitophagia_Line)

###### Отрисовка графика митофагии по клеточным линиям для каждой интересующей нас мутации

Input_Data_Mitophagia_Line <- subset.data.frame(Input_Data_Mitophagia_Line, Input_Data_Mitophagia_Line$Heteroplasmy > 5)

##### Добавление поля с наличием или отсутствием мутации

names(Input_Data_Mitophagia_Line)

Input_Data_Mitophagia_Line$"There_is_mutation" <- NA 

Input_Data_Mitophagia_Line$There_is_mutation <- ifelse(Input_Data_Mitophagia_Line$Mutation_name == "g15059a", "g15059a", "without g15059a")


my_comparisons <- combn(sort(unique(as.character(Input_Data_Mitophagia_Line$Line))),2, simplify = FALSE)

my_comparisons <- rev(list(
  c("a1555g", "g15059a"),
  c("c3256t", "g15059a"),
  c("c5178a", "g15059a"),
  c("g12315a", "g15059a"),
  c("g13513a", "g15059a"),
  c("g14846a", "g15059a"),
  c("t3336c", "g15059a")
)
)

Plot_File_Name <- file.path(PathToResults, paste("Box plots of mitophagy by cell lines for g15059a mutation", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

ggplot(Input_Data_Mitophagia_Line, aes(x = Line, y = Mitophagy)) +
  geom_boxplot(aes(fill = There_is_mutation), position = position_dodge(width = 1), outlier.size = 3) +
  stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 1), show.legend = FALSE,) +
  #stat_pvalue_manual(data = df_p_val, label = "p.adj.signif") +
  #stat_compare_means(method =  "wilcox.test", comparisons = my_comparisons, label = "p.format") +
  stat_compare_means(method =  "wilcox.test", comparisons = my_comparisons, label = "p.signif",  bracket.size = 0.5, size = 6) + 
  scale_fill_brewer(palette = "Set3", name = "Mutation") +
  #scale_x_discrete(limits = c("control", "ldl_only", "ANXA1_ldl_knockdown", "F2RL1_ldl_knockdown", "IL15_ldl_knockdown", "PERK_ldl_knockdown", "TSPYL2_ldl_knockdown"), labels = c( "All genes", "All genes", "ANXA1", "F2RL1", "IL15", "PERK", "TSPYL2")) +
  scale_y_continuous(breaks=seq(0, max(Line_Mitophagy_Inflamation_DataSet$Mitophagy), 0.1)) +
  xlab(label = "Cell line") +
  ylab(label = "Mitophagy change") +
  geom_point(aes(shape = "mean"),  alpha = 0) +
  geom_point(aes(shape = "median"),  alpha = 0) +
  geom_point(aes(shape = "outlier"),  alpha = 0) +
  guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("dodgerblue","black", "black"), shape = c(18, 95, 20), size = 5))) +
  #ggtitle("mitophagy by inflamation type of cell line") +
  theme_AO +
  theme(legend.key.size = unit(1, "cm"))

dev.off()

############################## Расчеты для объединенных данных Никиты и Собенина ###########################################

####### Открытие преобразованного входного файла данных

PathToSources_Sobenin_Data <- "C:\\Users\\Андрей\\OneDrive\\Documents\\Work\\НИИМЧ\\База данных по проекту толерантность\\Sources"

InputData_FileName <- file.path(PathToSources_Sobenin_Data, "all data 12_04_22.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
DataSet <- ListOfInputData_FileSheets

names(DataSet)
str(DataSet)

### Исправление некоторый полей таблицы исходных данных

DataSet$ATS <- factor(DataSet$ATS, labels = c("Норма", "Предрасположенность", "Атеросклероз"), levels = c(1, 2, 3), ordered = TRUE)
unique(DataSet$ATS)

write_json(x = DataSet, path = file.path(PathToSources_Sobenin_Data, "Clinical study tolerance_ALL_DATA_02_06_2022.json"))

write_json(x = Input_Data_Mitophagia_Line, path = file.path(PathToSources_Sobenin_Data, "Митофагия_мутации_гетероплазмия_04_06_2022.json"))

#DataSet_Json <- toJSON(x = DataSet, factor = "string", pretty = TRUE)
#write_json(x = DataSet_Json, path = file.path(PathToSources_Sobenin_Data, "Clinical study tolerance_ALL_DATA_02_06_2022_2.json"))

####### Добавление типа иммуногенности на основании классификации от Никиты

#grep(pattern = "_control_2_LPS|_1_LPS_2_LPS", x = names(DataSet), ignore.case = FALSE)
#temp_row <- DataSet[rows_counter, c(1, grep(pattern = "_control_2_LPS|_1_LPS_2_LPS", x = names(DataSet), ignore.case = FALSE))]

df_Immunogenicity <- data.frame()

for(rows_counter in 1:nrow(DataSet))
  {
   temp_row <- DataSet[rows_counter, c(1, grep(pattern = "_control_2_LPS|_1_LPS_2_LPS", x = names(DataSet), ignore.case = FALSE))]
   temp_row_long <- gather(temp_row, parametr, value, TNF_1_control_2_LPS:IL10_1_LPS_2_LPS, factor_key = TRUE)
   temp_row_long$parametr <- as.character(temp_row_long$parametr)
   temp_row_long$"mediator" <- str_extract(string = temp_row_long$parametr, pattern = "^[[A-Z,a-z]|[0-9]]+")
   
   for(temp_rows_counter in  1:nrow(temp_row_long))
     {
     if(grepl(pattern = "_control_2_LPS", x = temp_row_long$parametr[temp_rows_counter], ignore.case = FALSE)) 
      {
       temp_row_long$parametr[temp_rows_counter] <- "K_LPS"
      } else if(grepl(pattern = "_1_LPS_2_LPS", x = temp_row_long$parametr[temp_rows_counter], ignore.case = FALSE)) 
              {
               temp_row_long$parametr[temp_rows_counter] <- "LPS_LPS"
              }
     }
   
   for(temp_rows_counter in seq_along(unique(temp_row_long$mediator)))
   {
    temp_Data_mediator <- subset.data.frame(x = temp_row_long, temp_row_long$mediator == unique(temp_row_long$mediator)[temp_rows_counter])
    temp_mediator <- unique(temp_row_long$mediator)[temp_rows_counter]
    temp_K_LPS <- temp_Data_mediator$value[temp_Data_mediator$parametr == "K_LPS"]
    temp_LPS_LPS <- temp_Data_mediator$value[temp_Data_mediator$parametr == "LPS_LPS"]
    
    temp_Immunogenicity_type <- get_Immunogenicity_type(mediator = temp_mediator, KLPS = temp_K_LPS, LPSLPS = temp_LPS_LPS)
    
    temp_df_Immunogenicity <- data.frame("ID" = temp_Data_mediator$ID[1], "Mediator" = unique(temp_row_long$mediator)[temp_rows_counter], "Immunogenicity_type" = temp_Immunogenicity_type)
    
    df_Immunogenicity <- rbind.data.frame(df_Immunogenicity, temp_df_Immunogenicity)
   }
   
  }

df_Immunogenicity

df_Immunogenicity$Mediator <- paste(df_Immunogenicity$Mediator, "Immunogenicity_type", sep = "_")

unique(df_Immunogenicity$Immunogenicity_type)

####### Преобразование таблицы иммунного ответа в широкую форму

df_Immunogenicity_wide <- spread(data = df_Immunogenicity, Mediator, Immunogenicity_type)

####### Объединение таблиц

DataSet_new <- merge.data.frame(DataSet, df_Immunogenicity_wide, by = "ID")

DataSet <- DataSet_new

DataSet$CCL2_Immunogenicity_type <- as.factor(DataSet$CCL2_Immunogenicity_type)
DataSet$IL10_Immunogenicity_type <- as.factor(DataSet$IL10_Immunogenicity_type)
DataSet$IL1b_Immunogenicity_type <- as.factor(DataSet$IL1b_Immunogenicity_type)
DataSet$IL6_Immunogenicity_type <- as.factor(DataSet$IL6_Immunogenicity_type)
DataSet$IL8_Immunogenicity_type <- as.factor(DataSet$IL8_Immunogenicity_type)
DataSet$TNF_Immunogenicity_type <- as.factor(DataSet$TNF_Immunogenicity_type)

#################### Работа с новыми данными из Орла от 14-09-2022 ##################

###### Открытие файла с таблицей данных о митофагии, мутациях и гетероплазмии

list.files(PathToSources)

InputData_FileName <- file.path(PathToSources, "Митофагия_мутации_гетероплазмия_13-04-2022.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
Input_Data_Mitophagia_Line <- ListOfInputData_FileSheets

names(Input_Data_Mitophagia_Line)
str(Input_Data_Mitophagia_Line)

##### открытие файла Excel с данными о митофагии из Орла

FileName <- file.path(PathToSources, "Митофагия_3 повторности.xlsx")
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$R_данные

names(InputData_Excel_Sheet)
str(InputData_Excel_Sheet)

Data_Orel_long <- gather(InputData_Excel_Sheet, experiment, Value, "кол-во.клеток_базовая_повторность.1":"медиана_стимулированная.пируватом_повторность.3")

names(Data_Orel_long) <- c("Line", "experiment", "Mitophagy")

names(Data_Orel_long)

unique(Data_Orel_long$experiment)

Data_Orel_long$Value_type <- str_remove(str_extract(Data_Orel_long$experiment, "^[[:alpha:]|[\\.]|[-]]+_"), "_")
Data_Orel_long$experiment_type <- str_remove_all(str_extract(Data_Orel_long$experiment, "_[[:alpha:]|[\\.]|[-]]+_"), "_")
Data_Orel_long$experiment_number <- as.integer(str_extract(Data_Orel_long$experiment, "[:digit:]"))

str(Data_Orel_long)

###### Унификация названий клеточных линий

unique(Data_Orel_long$Line)[order(unique(Data_Orel_long$Line))]
unique(Input_Data_Mitophagia_Line$Line)[order(unique(Input_Data_Mitophagia_Line$Line))]

length(unique(Data_Orel_long$Line))
length(unique(Input_Data_Mitophagia_Line$Line))

Data_Orel_long$Line <- str_remove(Data_Orel_long$Line, "-")

Data_Orel_long$Line[which(Data_Orel_long$Line == "TCHSM1")] <- "HSM1"
Data_Orel_long$Line[which(Data_Orel_long$Line == "TCHSM2")] <- "HSM2"
Data_Orel_long$Line[which(Data_Orel_long$Line == "TCHSMAM1")] <- "MAM1"
Data_Orel_long$Line[which(Data_Orel_long$Line == "TCHSMAM2")] <- "MAM2"
Data_Orel_long$Line[which(Data_Orel_long$Line == "TCHSMAM2")] <- "MAM2"
Data_Orel_long$Line[which(Data_Orel_long$Line == "TCHSMAM3")] <- "MAM3"
Data_Orel_long$Line[which(Data_Orel_long$Line == "TCLSM1")] <- "LSM1"
Data_Orel_long$Line[which(Data_Orel_long$Line == "TCLSM2")] <- "LSM2"
Data_Orel_long$Line[which(Data_Orel_long$Line == "TCLSM3")] <- "LSM3"

setdiff(unique(Input_Data_Mitophagia_Line$Line), unique(Data_Orel_long$Line))
setdiff(unique(Data_Orel_long$Line), unique(Input_Data_Mitophagia_Line$Line))

####### Сохранение длинной таблицы данных о митофагии из Орла в файл .Rdata

Save_FileName <- file.path(PathToSources, paste("Митофагия_Орел_17-09-2022.Rdata", sep = ''))
saveRDS(Data_Orel_long, file = Save_FileName)

list.files(PathToSources)

###### Открытие файла с таблицей данных о митофагии, мутациях и гетероплазмии из Орла 

list.files(PathToSources)

InputData_FileName <- file.path(PathToSources, "Митофагия_Орел_17-09-2022.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
Data_Orel_long <- ListOfInputData_FileSheets

names(Data_Orel_long)
str(Data_Orel_long)

Data_Orel_long

##### открытие файла Excel с данными о гетероплазмии от Жигмитовой

FileName <- file.path(PathToSources, "Сводная таблица 19.09.22 К+ Лены Ж.xlsx")
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$R_данные

names(InputData_Excel_Sheet)
str(InputData_Excel_Sheet)

Data_Heteroplasmy_Zhigmitova_long <- gather(InputData_Excel_Sheet, Mutation_name, Heteroplasmy, t3336c:del652)

setdiff(unique(Data_Heteroplasmy_Zhigmitova_long$Line), unique(Data_Orel_long$Line))
setdiff(unique(Data_Orel_long$Line), unique(Data_Heteroplasmy_Zhigmitova_long$Line))

unique(Data_Orel_long$Line)[order(unique(Data_Orel_long$Line))]
unique(Data_Heteroplasmy_Zhigmitova_long$Line)[order(unique(Data_Heteroplasmy_Zhigmitova_long$Line))]

Data_Orel_Mitophagy_Heteroplasmy_Zhigmitova_long <- merge.data.frame(Data_Orel_long, Data_Heteroplasmy_Zhigmitova_long, by = "Line")

Data_Orel_Mitophagy_Heteroplasmy_Zhigmitova_long

####### Сохранение длинной таблицы данных о митофагии из Орла и гетероплазмии от Жигмитовой в файл .Rdata

Save_FileName <- file.path(PathToSources, paste("Митофагия_Орел_гетероплазмия_Жигмитова_23-09-2022.Rdata", sep = ''))
saveRDS(Data_Orel_Mitophagy_Heteroplasmy_Zhigmitova_long, file = Save_FileName)

list.files(PathToSources)

################## Открытие файла с таблицей данных о митофагии, мутациях и гетероплазмии из Орла и Жигмитовой
##################
##################
##################
##################

list.files(PathToSources)

InputData_FileName <- file.path(PathToSources, "Митофагия_Орел_гетероплазмия_Жигмитова_05-10-2022.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
Input_Data_Mitophagia_Line_Heteroplasmy <- ListOfInputData_FileSheets

names(Input_Data_Mitophagia_Line_Heteroplasmy)
str(Input_Data_Mitophagia_Line_Heteroplasmy)

Input_Data_Mitophagia_Line_Heteroplasmy <- Input_Data_Mitophagia_Line_Heteroplasmy[Input_Data_Mitophagia_Line_Heteroplasmy$Value_type == "среднее.значение",]

#### Вычисление среднего уровня базовой митофагии

Mitophagy_Table <- aggregate(Mitophagy ~ Line, data = Input_Data_Mitophagia_Line_Heteroplasmy[Input_Data_Mitophagia_Line_Heteroplasmy$experiment_type == "базовая",], FUN = mean)
Input_Data_Mitophagia_Line_Heteroplasmy$basal_level_mitophagy <- mean(Mitophagy_Table$Mitophagy) 

####### Вычисленение нормализованной митофагии (относительной митофагии)

Input_Data_Mitophagia_Line_Heteroplasmy$Mitophagy_relative <- Input_Data_Mitophagia_Line_Heteroplasmy$Mitophagy / Input_Data_Mitophagia_Line_Heteroplasmy$basal_level_mitophagy

####### Расчет индивидуальной базовой митофагии

Input_Data_Mitophagia_Line_Heteroplasmy$basal_level_mitophagy_individual <- NA
Input_Data_Mitophagia_Line_Heteroplasmy$Mitophagy_relative_individual <- NA

for (rows_counter in 1:nrow(Input_Data_Mitophagia_Line_Heteroplasmy)) 
  {
   temp_Line <- Input_Data_Mitophagia_Line_Heteroplasmy$Line[rows_counter] 
   temp_experiment_number <- Input_Data_Mitophagia_Line_Heteroplasmy$experiment_number[rows_counter]
   temp_basal_level_mitophagy <- unique(Input_Data_Mitophagia_Line_Heteroplasmy$Mitophagy[(Input_Data_Mitophagia_Line_Heteroplasmy$Line == temp_Line)&(Input_Data_Mitophagia_Line_Heteroplasmy$experiment_number == temp_experiment_number)&(Input_Data_Mitophagia_Line_Heteroplasmy$experiment_type == "базовая")&(Input_Data_Mitophagia_Line_Heteroplasmy$Value_type == "среднее.значение")])
   Input_Data_Mitophagia_Line_Heteroplasmy$basal_level_mitophagy_individual[rows_counter] <- temp_basal_level_mitophagy
   Input_Data_Mitophagia_Line_Heteroplasmy$Mitophagy_relative_individual[rows_counter] <- Input_Data_Mitophagia_Line_Heteroplasmy$Mitophagy[rows_counter] / temp_basal_level_mitophagy 
   message("Строка № ", rows_counter, " из ", nrow(Input_Data_Mitophagia_Line_Heteroplasmy))
  }


DataSet_Zhigmitova <- Input_Data_Mitophagia_Line_Heteroplasmy

DataSet_Zhigmitova

####################

DataSet <- Input_Data_Mitophagia_Line_Heteroplasmy

unique(DataSet$Mutation_name)

names(DataSet)
str(DataSet)

######## Переименование названий мутаций в правильный формат

unique(DataSet$Mutation_name)

DataSet$Mutation_name[which(DataSet$Mutation_name == "a1555g")] <- "m.1555A>G"
DataSet$Mutation_name[which(DataSet$Mutation_name == "del652")] <- "del652G"
DataSet$Mutation_name[which(DataSet$Mutation_name == "c5178a")] <- "m.5178C>A"
DataSet$Mutation_name[which(DataSet$Mutation_name == "g15059a")] <- "m.15059G>A"
DataSet$Mutation_name[which(DataSet$Mutation_name == "g14459a")] <- "m.14459G>A"
DataSet$Mutation_name[which(DataSet$Mutation_name == "g14846a")] <- "m.14846G>A"
DataSet$Mutation_name[which(DataSet$Mutation_name == "g13513a")] <- "m.13513G>A"
DataSet$Mutation_name[which(DataSet$Mutation_name == "c3256t")] <- "m.3256C>T"
DataSet$Mutation_name[which(DataSet$Mutation_name == "t3336c")] <- "m.3336T>C"

####### Сохранение длинной таблицы данных о митофагии из Орла и гетероплазмии от Жигмитовой с исправленными мазваниями мутаций в файл .Rdata

Save_FileName <- file.path(PathToSources, paste("Митофагия_Орел_гетероплазмия_Жигмитова_05-10-2022.Rdata", sep = ''))
saveRDS(DataSet, file = Save_FileName)

list.files(PathToSources)

Friedman_results_Table <- data.frame()
Friedman_post_hoc_Table <- data.frame()

for(value_type_counter in seq_along(unique(DataSet$Value_type)))
  {
   for (experiment_type_counter in seq_along(unique(DataSet$experiment_type)))
     {
       message("Расчет ", unique(DataSet$Value_type)[value_type_counter], " для ", unique(DataSet$experiment_type)[experiment_type_counter])
     
       temp_DataSet <- subset.data.frame(DataSet, (DataSet$Value_type == unique(DataSet$Value_type)[value_type_counter]) & (DataSet$experiment_type == unique(DataSet$experiment_type)[experiment_type_counter]), select = c("Line", "Mitophagy", "experiment_number"))
       temp_DataSet$experiment_number <- as.factor(temp_DataSet$experiment_number)
       temp_DataSet$Line <- as.factor(temp_DataSet$Line)
       
       temp_DataSet <- reshape2::dcast(temp_DataSet, Line ~ experiment_number, value.var="Mitophagy", fun.aggregate = mean_without_NA)
       row.names(temp_DataSet) <- temp_DataSet$Line
       temp_DataSet_matrix <- as.matrix(temp_DataSet[complete.cases(temp_DataSet), c("1", "2", "3")])
       
       Friedman_results <- friedman.test(temp_DataSet_matrix)
       
       temp_Friedman_results_Table <- data.frame(
       "Тип_эксперимента" = unique(DataSet$experiment_type)[experiment_type_counter],
       "Тип_значения" = unique(DataSet$Value_type)[value_type_counter],
       "Название теста" = Friedman_results$method,
       "Название набора данных" = Friedman_results$data.name,
       "Значение статистики" = Friedman_results$statistic,
       "р-значение" = Friedman_results$p.value,
       "Результат" = ifelse(Friedman_results$p.value > 0.05, "Различий нет", "Различия есть")
       )

       Friedman_results_Table <- rbind.data.frame(Friedman_results_Table, temp_Friedman_results_Table)


       temp_Friedman_post_hoc_Table <- friedmanmc(temp_DataSet_matrix)


       names(temp_Friedman_post_hoc_Table$dif.com) <- c("Уровень различий", "Критический уровень различий", "Достоверность различий")


       temp_temp_Friedman_post_hoc_Table <- data.frame("Пары сравнения" = row.names(temp_Friedman_post_hoc_Table$dif.com),
                                      temp_Friedman_post_hoc_Table$dif.com[,1:3],
                                      "p-значение" = temp_Friedman_post_hoc_Table$p.value,
                                      "Тип_эксперимента" = unique(DataSet$experiment_type)[experiment_type_counter],
                                      "Тип_значения" = unique(DataSet$Value_type)[value_type_counter]
                                      )

Friedman_post_hoc_Table <- rbind.data.frame(Friedman_post_hoc_Table, temp_temp_Friedman_post_hoc_Table)

   }
  }

Friedman_results_Table
Friedman_post_hoc_Table

###### Сохранение таблиц в файл Excel

List_Of_Result_Table <- list(Friedman_results_Table, Friedman_post_hoc_Table)
names(List_Of_Result_Table) <- c("Тест_Фридмана", "post_hoc_тесты_Фридмана")
Table_ResultFileName <- file.path(PathToResults, paste("Исследование различий повторностей в митофагии из Орла", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

####### расчет коэффициентов корреляции для митофагии и гетероплазмии для каждого типа эксперимента в отдельности

unique(DataSet$Value_type)
names(DataSet)

Input_Data <- subset.data.frame(x = DataSet, DataSet$Value_type == "среднее.значение", select = c("Mitophagy", "experiment_type", "Heteroplasmy"))

for(experiment_type_counter in seq_along(unique(Input_Data$experiment_type)))
  {
   
  temp_Input_Data <- subset.data.frame(Input_Data, Input_Data$experiment_type == unique(Input_Data$experiment_type)[experiment_type_counter])

InputData_Cor <- melt(data = temp_Input_Data)

InputData_Cor <- InputData_Cor[, c("variable", "value")]

#str(InputData_Cor)

InputData_Cor$variable <- as.character(InputData_Cor$variable)

#######

Cor_Method <- c("pearson", "spearman", "kendall")

###### для отрисовки матриц корреляции

col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow",
                           "#7FFF7F", "cyan", "#007FFF", "blue","#00007F"))

# col4 <- colorRampPalette(c("#00007F","blue","#007FFF","cyan",
#                            "#7FFF7F", "yellow", "#FF7F00", "red","#7F0000"))
# 
# col4 <- colorRampPalette(c("Red","Yellow", "DarkCyan"))



for (Cor_Method_Counter in seq_along(Cor_Method))
{
  CorSet_List <- list()
  
  Corr_Results <- corr.test(x = temp_Input_Data[, c("Mitophagy", "Heteroplasmy")], use = "complete", method = Cor_Method[Cor_Method_Counter])
  print(Corr_Results, short = FALSE)
  
  r_table <- as.data.frame(Corr_Results$r)
  r_table$"Показатель" <- row.names(Corr_Results$r)
  r_table <- r_table[c(length(names(r_table)), 1:(length(names(r_table)))-1)]
  row.names(r_table) <- NULL
  CorSet_List[["Correlation matrix"]] <- r_table
  
  p_table <- as.data.frame(Corr_Results$p)
  p_table$"Показатель" <- row.names(Corr_Results$p)
  p_table <- p_table[c(length(names(p_table)), 1:(length(names(p_table)))-1)]
  row.names(p_table) <- NULL
  CorSet_List[["p-values"]] <- p_table 
  
  ci_table <- as.data.frame(Corr_Results$ci)
  ci_table$"Показатель" <- row.names(Corr_Results$ci)
  ci_table <- ci_table[c(length(names(ci_table)), 1:(length(names(ci_table)))-1)]
  row.names(ci_table) <- NULL
  CorSet_List[["Confidence intervals"]] <- ci_table
  
  #### Сохранение результатов  в Excel файл 
  
  ResultFileName <- paste("Корреляция по ", Cor_Method[Cor_Method_Counter], "_", unique(Input_Data$experiment_type)[experiment_type_counter], "_", Sys.Date(), ".xlsx", sep = "")
  Table_ResultFileName <- file.path(PathToResults, ResultFileName)
  save_to_excel_file_by_XLConnect(Table_ResultFileName, CorSet_List)
  
  temp_CorSet <- cor(temp_Input_Data[, c("Mitophagy", "Heteroplasmy")], method = Cor_Method[Cor_Method_Counter], use = "complete")
  temp_CorSet[is.na(temp_CorSet)] <- 0
  
  #CorSet_List[[Cor_Method_Counter]] <- as.data.frame(temp_CorSet)
  
  ###### отрисовка матриц корреляции
  
  File_Name <- paste(Cor_Method[Cor_Method_Counter], "_", unique(Input_Data$experiment_type)[experiment_type_counter], "_Матрица_корреляции", "_", Sys.Date(), ".jpg", sep = "")
  Saved_Fugure_File_Name <- file.path(PathToResults, File_Name)
  jpeg(filename = Saved_Fugure_File_Name, width = 800, height = 600, units = "px", pointsize =14, quality = 200, bg = "white", res = NA, restoreConsole = TRUE)
  
  Figure <- corrplot(temp_CorSet, method = "color", col = col4(20), cl.length = 21,
                     order = "AOE", addCoef.col="black",
                     type = "upper",
                     mar = c(1, 1, 1, 1),
                     diag = FALSE
  )
  print(Figure)
  dev.off()
  
}
}

####### Подготовка набора данных

names(DataSet)

new_DataSet <- subset.data.frame(DataSet, DataSet$Value_type == "медиана")




####### Диаграммы рассеяния с кондиционными переменными ############

names(DataSet)

temp_DataSet <- subset.data.frame(DataSet, DataSet$experiment_type == "стимулированная.пируватом" & DataSet$Value_type == "медиана") 

Mutation_name_Table <- aggregate(Heteroplasmy ~ Mutation_name, data = temp_DataSet, FUN = median)

Mutation_name_Table <- Mutation_name_Table[order(Mutation_name_Table$Heteroplasmy, decreasing = FALSE),]

Plot_File_Name <- file.path(PathToResults, paste("Scatter plots of mitophagy by heteroplasmy for mutation", unique(temp_DataSet$experiment_type), Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")

Figure_facet <- ggplot(temp_DataSet, aes(x = Heteroplasmy, y = Mitophagy)) +
geom_point(aes(colour = as.factor(temp_DataSet$experiment_type)), size = 2) +
facet_wrap( ~ factor(Mutation_name, levels = Mutation_name_Table$Mutation_name), nrow = 2) +
labs(col = "Experiment\nnumber") +
scale_y_continuous(breaks=seq(0, max(temp_DataSet$Mitophagy, na.rm = TRUE), 5)) +
scale_x_continuous(breaks=seq(0, 100, 20)) +
theme_AO

print(Figure_facet)

dev.off()

####### Диаграммы рассеяния с кондиционными переменными для всех типов экспериментов с объединением по линиям ############

names(DataSet)

DataSet$experiment_number_factor <- factor(DataSet$experiment_number, levels = c(1, 2, 3), labels = c("1st", "2nd", "3rd"))
DataSet$experiment_type <- factor(DataSet$experiment_type, levels = c("базовая", "стимулированная.пируватом", "стимулированная.FCCP"), labels = c("Базовая", "Пируват", "FCCP"))

for(Line_counter in seq_along(unique(DataSet$Line)))
  {
   temp_DataSet <- subset.data.frame(DataSet, DataSet$Line == unique(DataSet$Line)[Line_counter] & DataSet$Value_type == "медиана") 

   unique(temp_DataSet$Mutation_name)

#temp_DataSet[temp_DataSet$Mutation_name == "a1555g", c("Mitophagy", "experiment_type", "Heteroplasmy", "experiment_number_factor")]

   Mutation_name_Table <- aggregate(Heteroplasmy ~ Mutation_name, data = temp_DataSet, FUN = median)

   Mutation_name_Table <- Mutation_name_Table[order(Mutation_name_Table$Heteroplasmy, decreasing = FALSE),]

   Plot_File_Name <- file.path(PathToResults, paste("Scatter plots of mitophagy by heteroplasmy for mutation", unique(temp_DataSet$Value_type), unique(DataSet$Line)[Line_counter], Sys.Date(), ".jpg", sep = "_"))
   jpeg(filename = Plot_File_Name, width = 1500, height = 900, units = "px", pointsize = 18, quality = 200, bg = "white")

   Figure_facet <- ggplot(temp_DataSet, aes(x = Heteroplasmy, y = Mitophagy)) +
   geom_point(aes(colour = experiment_type), size = 5) +
   geom_line(aes(group = experiment_number_factor, linetype = experiment_number_factor), alpha = 0.5, size = 1) +
   geom_line(aes(x = 5), linetype = "dashed", col = "red", size = 0.5) +
   facet_wrap( ~ factor(Mutation_name, levels = Mutation_name_Table$Mutation_name), nrow = 3) +
   labs(colour = "Митофагия", linetype = "Номер\nповтора") +
  #guides(col = "none") + 
   scale_y_continuous(breaks=seq(0, max(temp_DataSet$Mitophagy, na.rm = TRUE), 5)) +
   scale_x_continuous(breaks=seq(0, 100, 20)) +
   theme_AO

   print(Figure_facet)

   dev.off()

   message("Работа с линией ", unique(DataSet$Line)[Line_counter], " завершена.")
  }

################ Определение типа митофагии - дефектная или индуцируемая ###########


Kruskal_results_Table <- data.frame()
Kruskal_posthoc_Table <- data.frame()


for(Line_counter in seq_along(unique(DataSet$Line)))
{
  temp_DataSet <- subset.data.frame(DataSet, DataSet$Line == unique(DataSet$Line)[Line_counter] & DataSet$Value_type == "медиана")
  kruskal_results <- kruskal.test(Mitophagy ~ experiment_type, data = temp_DataSet)
  
  Base_mean_median <- mean_without_NA(temp_DataSet$Mitophagy[temp_DataSet$experiment_type == "Базовая"])
  Pyruvate_mean_median <- mean_without_NA(temp_DataSet$Mitophagy[temp_DataSet$experiment_type == "Пируват"])
  FCCP_mean_median <- mean_without_NA(temp_DataSet$Mitophagy[temp_DataSet$experiment_type == "FCCP"])
  
  
  temp_Kruskal_results_Table <- data.frame(
      "Line" = unique(DataSet$Line)[Line_counter],
      "Название теста" = kruskal_results$method,
      "Название набора данных" = kruskal_results$data.name,
      "Значение статистики" = kruskal_results$statistic,
      "Базовая митофагия" = Base_mean_median,
      "Пируват" = Pyruvate_mean_median,
      "FCCP" = FCCP_mean_median,
      "р-значение" = kruskal_results$p.value,
      "Результат" = ifelse(kruskal_results$p.value >= 0.05, "дефектная", "индуцируемая"),
      "Результат_уточненный" = NA
    )
    
    if(temp_Kruskal_results_Table$р.значение < 0.05)
    {
      kruskal_posthoc_results <- kruskalmc(Mitophagy ~ experiment_type, data = temp_DataSet, probs = 0.05)
      
      names(kruskal_posthoc_results$dif.com) <- c("Уровень различий", "Критический уровень различий", "Достоверность различий")
      
      # Base_mean_median <- mean_without_NA(temp_DataSet$Mitophagy[temp_DataSet$experiment_type == "Базовая"])
      # Pyruvate_mean_median <- mean_without_NA(temp_DataSet$Mitophagy[temp_DataSet$experiment_type == "Пируват"])
      # FCCP_mean_median <- mean_without_NA(temp_DataSet$Mitophagy[temp_DataSet$experiment_type == "FCCP"])
      # 
      temp_kruskal_posthoc_results_Table <- data.frame(
        "Line" = unique(DataSet$Line)[Line_counter],
        "Пары сравнения" = row.names(kruskal_posthoc_results$dif.com),
        kruskal_posthoc_results$dif.com[,1:3],
        # "Базовая митофагия" = Base_mean_median,
        # "Пируват" = Pyruvate_mean_median,
        # "FCCP" = FCCP_mean_median,
        "p-значение" = kruskal_posthoc_results$signif.level
      )
      Kruskal_posthoc_Table <- rbind.data.frame(temp_kruskal_posthoc_results_Table, Kruskal_posthoc_Table)
      
      temp_Kruskal_results_Table$Результат_уточненный <- ifelse((Base_mean_median >= Pyruvate_mean_median | Base_mean_median >= FCCP_mean_median), "дефектная", "индуцируемая") 
    } else
             {
              temp_Kruskal_results_Table$Результат_уточненный  <- "дефектная"  
             }
    
    Kruskal_results_Table <- rbind.data.frame(temp_Kruskal_results_Table, Kruskal_results_Table)
    
    }

Kruskal_results_Table
Kruskal_posthoc_Table

####### Сохранение таблицы результатов в Еxcel файл

List_Of_Result_Table <- list(Kruskal_results_Table, Kruskal_posthoc_Table)
names(List_Of_Result_Table) <- c("Сводная_таблица", "Парные_сравнения")
Table_ResultFileName <- file.path(PathToResults, paste("Определение типа митофагии из Орла", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

Mitophagy_type_Table <- Kruskal_results_Table[, c("Line", "Результат_уточненный")]

row.names(Mitophagy_type_Table) <- NULL

names(Mitophagy_type_Table) <- c("Line", "Тип митофагии")

DataSet_Mutations_Mitophagy_type <- merge.data.frame(DataSet, Mitophagy_type_Table, by = "Line")

########## Изучение связи мутаций с типом митофагии

unique(DataSet_Mutations_Mitophagy_type$Mutation_name)
unique(DataSet_Mutations_Mitophagy_type$`Тип митофагии`)

DataSet_Mutations_Mitophagy_type <- DataSet_Mutations_Mitophagy_type[DataSet_Mutations_Mitophagy_type$Heteroplasmy > 5,]

Test_result <- data.frame()

for (Mutations_counter in seq_along(unique(DataSet_Mutations_Mitophagy_type$Mutation_name)))
{
  tempData <- subset.data.frame(DataSet_Mutations_Mitophagy_type, DataSet_Mutations_Mitophagy_type$Mutation_name == unique(DataSet_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter] & DataSet_Mutations_Mitophagy_type$Value_type == "медиана")
  
  
  tempData$"experiment_type_mitophagy_type" <- paste(tempData$experiment_type, tempData$'Тип митофагии', sep = "_")
  
  tempData_PowerPoint <-  data.frame("Тип митофагии" = tempData$`Тип митофагии`,
                                    "Ряд1" = as.numeric(tempData$Heteroplasmy)
                                    )
  tempData_PowerPoint$Тип.митофагии <- as.character(tempData_PowerPoint$Тип.митофагии)
  
  #str(tempData_PowerPoint)
  
  tempData_PowerPoint$Тип.митофагии[which(tempData_PowerPoint$Тип.митофагии == "индуцируемая")] <- "induced" 
  tempData_PowerPoint$Тип.митофагии[which(tempData_PowerPoint$Тип.митофагии == "дефектная")] <- "defective"
    
  ##### сохранение данных в файл Excel для графика PowerPoint
  
  List_Of_Result_Table <- list(tempData_PowerPoint)
  names(List_Of_Result_Table) <- c(unique(DataSet_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter])
  Table_ResultFileName <- file.path(PathToResults, paste("Гетероплазмия_мутации", str_extract(unique(DataSet_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter], "[0-9]+"), "для_PowerPoint", Sys.Date(), ".xlsx", sep = "_"))
  save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)
  
  #####
  
  my_comparisons1 <- rev(list(
    c("дефектная", "индуцируемая")))
  
  Plot_File_Name <- file.path(PathToResults, paste("Boxplots of heteroplasmy by mitophagy type for mutation", str_extract(unique(DataSet_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter], "[0-9]+"), Sys.Date(), ".jpg", sep = "_"))
  jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")
  
  
  Figure_boxplot <- ggplot(tempData, aes(x = tempData$`Тип митофагии`, y = Heteroplasmy)) +
    geom_boxplot(outlier.size = 3, position = position_dodge(width = 0.9)) +
    stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE,) +
    stat_compare_means(method =  "wilcox.test", comparisons = my_comparisons1, label = "p.signif", bracket.size = 0.5, size = 6) +
    xlab(label = "Type of mitophagy") +
    ylab(label = "Heteroplasmy, %") +
    scale_x_discrete(limits = c("дефектная", "индуцируемая"), labels = c( "defective", "induced")) +
    # geom_point(aes(shape = "mean"),  alpha = 0) +
    # geom_point(aes(shape = "median"),  alpha = 0) +
    # geom_point(aes(shape = "outlier"),  alpha = 0) +
    # guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
    # #ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
    theme_AO
  
  print(Figure_boxplot)
  dev.off()
  
  #combn(sort(unique(as.character(tempData$experiment_type_mitophagy_type))), 2, simplify = FALSE)
  
  
  my_comparisons2 <- rev(list(
    c("FCCP_дефектная", "Базовая_дефектная"),
    c("FCCP_дефектная", "Пируват_дефектная"),
    c("Базовая_дефектная", "Пируват_дефектная"),
    c("Базовая_индуцируемая", "Пируват_индуцируемая"),
    c("FCCP_индуцируемая",  "Пируват_индуцируемая"),
    c("FCCP_индуцируемая",  "Базовая_индуцируемая")
  )
  )
  
  tempData_PowerPoint <- dcast(tempData, `Тип митофагии` + Line + experiment_number ~ experiment_type, value.var = "Mitophagy", fun.aggregate = mean_without_NA)
  tempData_PowerPoint <- tempData_PowerPoint[c("Тип митофагии", "Базовая", "Пируват",  "FCCP")]
  tempData_PowerPoint$`Тип митофагии` <- as.character(tempData_PowerPoint$`Тип митофагии`)
  
  #str(tempData_PowerPoint)
  
  tempData_PowerPoint$`Тип митофагии`[which(tempData_PowerPoint$`Тип митофагии` == "индуцируемая")] <- "induced" 
  tempData_PowerPoint$`Тип митофагии`[which(tempData_PowerPoint$`Тип митофагии` == "дефектная")] <- "defective"
  
  ##### сохранение данных в файл Excel для графика PowerPoint
  
  List_Of_Result_Table <- list(tempData_PowerPoint)
  names(List_Of_Result_Table) <- c(unique(DataSet_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter])
  Table_ResultFileName <- file.path(PathToResults, paste("Митофагия_мутации", str_extract(unique(DataSet_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter], "[0-9]+"), "для_PowerPoint", Sys.Date(), ".xlsx", sep = "_"))
  save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)
  
  #####
  
  
  Plot_File_Name <- file.path(PathToResults, paste("Boxplots of mitophagy by mitophagy type for mutation", str_extract(unique(DataSet_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter], "[0-9]+"), Sys.Date(), ".jpg", sep = "_"))
  jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")
  
  
  Figure_boxplot <- ggplot(tempData, aes(x = tempData$experiment_type_mitophagy_type, y = Mitophagy, fill = experiment_type)) +
    geom_boxplot(outlier.size = 3, position = position_dodge(width = 0.9)) +
    stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE) +
    stat_compare_means(method =  "wilcox.test", comparisons = my_comparisons2, label = "p.signif", bracket.size = 0.5, size = 6) +
    xlab(label = "Type of mitophagy") +
    ylab(label = "Mitophagy, %") +
    labs(fill = "Type of mitophagy") +
    scale_x_discrete(limits = c("Базовая_дефектная", "Пируват_дефектная", "FCCP_дефектная",
                                "Базовая_индуцируемая", "Пируват_индуцируемая", "FCCP_индуцируемая"), 
                     labels = c("", "defective", "", "", "induced", "")) +
    scale_fill_discrete(limits = c("Базовая", "Пируват", "FCCP"), labels = c("Base", "Pyruvate", "FCCP")) +
    # geom_point(aes(shape = "mean"),  alpha = 0) +
    # geom_point(aes(shape = "median"),  alpha = 0) +
    # geom_point(aes(shape = "outlier"),  alpha = 0) +
    # guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
    guides(fill = guide_legend(title = NULL))  +
    # #ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
    theme_AO +
    theme(axis.ticks.x = element_blank())
  
  print(Figure_boxplot)
  dev.off()
  
  pairwise_wilcox_test_result <- pairwise.wilcox.test(tempData$Mitophagy, tempData$experiment_type_mitophagy_type , p.adjust.method = "holm")
  
  temp_pairwise_wilcox_test_result <- data.frame("Мутация" = unique(DataSet_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter], 
                                                 "Название_группы" = row.names(pairwise_wilcox_test_result$p.value),
                                                 pairwise_wilcox_test_result$p.value,
                                                 "Название теста" = pairwise_wilcox_test_result$method,
                                                 "Метод_коррекции" = pairwise_wilcox_test_result$p.adjust.method
                                                 )
  row.names(temp_pairwise_wilcox_test_result) <- NULL
                                                 
  ##### сохранение данных в файл Excel для графика множественных сравнений
  
  List_Of_Result_Table <- list(temp_pairwise_wilcox_test_result)
  names(List_Of_Result_Table) <- c(unique(DataSet_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter])
  Table_ResultFileName <- file.path(PathToResults, paste("Митофагия_мутации_множественные_сравнения", str_extract(unique(DataSet_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter], "[0-9]+"), "для_PowerPoint", Sys.Date(), ".xlsx", sep = "_"))
  save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)
  
  #####
  
  test_result <- wilcox.test(tempData$Heteroplasmy ~ tempData$`Тип митофагии`, paired = FALSE)
  
  Defect_mean <- mean_without_NA(tempData$Heteroplasmy[tempData$`Тип митофагии` == "дефектная"])
  Norm_mean <- mean_without_NA(tempData$Heteroplasmy[tempData$`Тип митофагии` == "индуцируемая"])
  
  temp_Test_result <- data.frame("Мутация" = unique(DataSet_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter], 
                                 "Название теста" = test_result$method,
                                 "Значение теста" = test_result$statistic,
                                 "Гетерозизотность_дефектной" = Defect_mean,
                                 "Гетерозизотность_индуцированной" = Norm_mean,
                                 "р-значение" = test_result$p.value,
                                 "Результат теста" = ifelse(test_result$p.value < 0.05, "Мутация связана", "Мутация не связана"))
  
  Test_result <- rbind.data.frame(Test_result, temp_Test_result)
  
}

Test_result

######### сохранение результатов проверки на тип распределения

List_Of_Result_Table <- list(Test_result)
names(List_Of_Result_Table) <-  c("Сводная_таблица")
Table_ResultFileName <- file.path(PathToResults, paste("Взаимосвязь гетероплазмии мутаций и типа митофагии из Орла", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

unique(DataSet$Line)

######### Работа с данными по митофагии и гетероплазмии от Журавлёва и Никифорова от мая 2022 года

##### открытие файла Excel с данными о гетероплазмии от Журавлёва и Никифорова от мая 2022 года

FileName <- file.path(PathToSources, "Цибриды_Митофагия_Гетероплазмия_Никифоров_в_единицах_Винокурова.xlsx")
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$`Heteroplasmy_05-2022_R`

names(InputData_Excel_Sheet)
str(InputData_Excel_Sheet)

Data_Heteroplasmy_Zhuravlev_long <- gather(InputData_Excel_Sheet, Mutation_name, Heteroplasmy, del652G:"m.14846G>A")
Data_Heteroplasmy_Zhuravlev_long

InputData_Excel_Sheet <- ListOfFileSheets$`Mitophagy_05-2022_R`

###### Открытие файла с новыми данными по митофагии от Никифорова от 16-10-2022

FileName <- file.path(PathToSources, "Митофагия (4-10%). Ноябрь 2021. FCCP. Первый прогон всех линий.xlsx")
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$Mitophagy_Nikiforov_R

######

###### Открытие файла с новыми данными по митофагии от Никифорова от 17-10-2022

FileName <- file.path(PathToSources, "Цибриды_Митофагия_Гетероплазмия_Никифоров_в_единицах_Винокурова.xlsx")
ListOfFileSheets <- get_excel_file_content_by_openxlsx(FileName)
InputData_Excel_Sheet <- ListOfFileSheets$`Mitophagy_05-2022_%_R`

######

Data_Mitophagy_Zhuravlev_long <- gather(InputData_Excel_Sheet, Experiment, Mitophagy, "базовая_среднее.значение":"стимулированная.FCCP_среднее.значение")
Data_Mitophagy_Zhuravlev_long

setequal(unique(Data_Mitophagy_Zhuravlev_long$Line), unique(Data_Heteroplasmy_Zhuravlev_long$Line))

Data_Mitophagy_Heteroplasmy_Zhuravlev_long <- merge.data.frame(Data_Mitophagy_Zhuravlev_long, Data_Heteroplasmy_Zhuravlev_long, by = "Line")

Data_Mitophagy_Heteroplasmy_Zhuravlev_long

Data_Mitophagy_Heteroplasmy_Zhuravlev_long$experiment_type <- str_remove(str_extract(Data_Mitophagy_Heteroplasmy_Zhuravlev_long$Experiment, "^[[[:alpha:]|[\\.]]]+_"), "_")
Data_Mitophagy_Heteroplasmy_Zhuravlev_long$Value_type <- str_remove_all(str_extract(Data_Mitophagy_Heteroplasmy_Zhuravlev_long$Experiment, "_[[[:alpha:]|[\\.]]]+"), "_")

#### Вычисление среднего уровня базовой митофагии

Mitophagy_Table <- aggregate(Mitophagy ~ Line, data = Data_Mitophagy_Heteroplasmy_Zhuravlev_long[Data_Mitophagy_Heteroplasmy_Zhuravlev_long$experiment_type == "базовая",], FUN = mean)
Data_Mitophagy_Heteroplasmy_Zhuravlev_long$basal_level_mitophagy <- mean(Mitophagy_Table$Mitophagy) 

####### Вычисленение нормализованной митофагии (относительной митофагии)

Data_Mitophagy_Heteroplasmy_Zhuravlev_long$Mitophagy_relative <- Data_Mitophagy_Heteroplasmy_Zhuravlev_long$Mitophagy / Data_Mitophagy_Heteroplasmy_Zhuravlev_long$basal_level_mitophagy

DataSet_Zhigmitova <- Input_Data_Mitophagia_Line_Heteroplasmy

DataSet_Zhigmitova

####### Сохранение длинной таблицы с данными по митофагии и гетероплазмии от Журавлёва и Никифорова от мая 2022 года в файл .Rdata

Save_FileName <- file.path(PathToSources, paste("Митофагия_гетероплазмия_Журавлев_Никифоров_17-10-2022.Rdata", sep = ''))
saveRDS(Data_Mitophagy_Heteroplasmy_Zhuravlev_long, file = Save_FileName)

list.files(PathToSources)

###### Открытие файла с таблицей с данными по митофагии и гетероплазмии от Журавлёва и Никифорова от мая 2022 года

list.files(PathToSources)

InputData_FileName <- file.path(PathToSources, "Митофагия_гетероплазмия_Журавлев_Никифоров_17-10-2022.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
Input_Data_Mitophagia_Line_Heteroplasmy <- ListOfInputData_FileSheets

names(Input_Data_Mitophagia_Line_Heteroplasmy)
str(Input_Data_Mitophagia_Line_Heteroplasmy)

####### Расчет индивидуальной базовой митофагии

Input_Data_Mitophagia_Line_Heteroplasmy$basal_level_mitophagy_individual <- NA
Input_Data_Mitophagia_Line_Heteroplasmy$Mitophagy_relative_individual <- NA

for (rows_counter in 1:nrow(Input_Data_Mitophagia_Line_Heteroplasmy)) 
{
  temp_Line <- Input_Data_Mitophagia_Line_Heteroplasmy$Line[rows_counter] 
  temp_experiment_number <- Input_Data_Mitophagia_Line_Heteroplasmy$experiment_number[rows_counter]
  temp_basal_level_mitophagy <- unique(Input_Data_Mitophagia_Line_Heteroplasmy$Mitophagy[(Input_Data_Mitophagia_Line_Heteroplasmy$Line == temp_Line)&(Input_Data_Mitophagia_Line_Heteroplasmy$experiment_number == temp_experiment_number)&(Input_Data_Mitophagia_Line_Heteroplasmy$experiment_type == "базовая")&(Input_Data_Mitophagia_Line_Heteroplasmy$Value_type == "среднее.значение")])
  Input_Data_Mitophagia_Line_Heteroplasmy$basal_level_mitophagy_individual[rows_counter] <- temp_basal_level_mitophagy
  Input_Data_Mitophagia_Line_Heteroplasmy$Mitophagy_relative_individual[rows_counter] <- Input_Data_Mitophagia_Line_Heteroplasmy$Mitophagy[rows_counter] / temp_basal_level_mitophagy 
  message("Строка № ", rows_counter, " из ", nrow(Input_Data_Mitophagia_Line_Heteroplasmy))
}

DataSet_Zuravlev <- Input_Data_Mitophagia_Line_Heteroplasmy

unique(DataSet_Zuravlev$Mutation_name)

DataSet_Zuravlev

############# Сравнение гетероплазмии и митофагии Журавлева и Жегмитовой

names(DataSet_Zhigmitova)
names(DataSet_Zuravlev)

setequal(names(DataSet_Zhigmitova), names(DataSet_Zuravlev))

names(DataSet_Zuravlev)[which(names(DataSet_Zuravlev) == "Experiment")] <- "experiment"

DataSet_Zhigmitova$Author <- "Zhigmitova"
DataSet_Zuravlev$Author <- "Zuravlev"
#DataSet_Zuravlev$experiment_number <- 1

DataSet_Zhigmitova_Zuravlev <- merge(DataSet_Zhigmitova, DataSet_Zuravlev, all = TRUE)

names(DataSet_Zhigmitova_Zuravlev)

unique(DataSet_Zhigmitova_Zuravlev$Line)
unique(DataSet_Zhigmitova_Zuravlev$experiment)
unique(DataSet_Zhigmitova_Zuravlev$Value_type)
unique(DataSet_Zhigmitova_Zuravlev$experiment_type)
unique(DataSet_Zhigmitova_Zuravlev$experiment_number)
unique(DataSet_Zhigmitova_Zuravlev$Mutation_name)
unique(DataSet_Zhigmitova_Zuravlev$Author)

####### Сохранение длинной объединенной таблицы с данными по митофагии и гетероплазмии от Журавлёва и Никифорова и Винокурова и Жегмитовой в файл .Rdata

Save_FileName <- file.path(PathToSources, paste("Митофагия_гетероплазмия_объединенный_файл_19-10-2022.Rdata", sep = ''))
saveRDS(DataSet_Zhigmitova_Zuravlev, file = Save_FileName)

list.files(PathToSources)

###### Открытие файла с объединенной таблицы с данными по митофагии и гетероплазмии от Журавлёва и Никифорова и Винокурова и Жегмитовой в файл .Rdata

list.files(PathToSources)

InputData_FileName <- file.path(PathToSources, "Митофагия_гетероплазмия_объединенный_файл_19-10-2022.Rdata")
ListOfInputData_FileSheets <- readRDS(file = InputData_FileName)
Input_Data_Mitophagia_Line_Heteroplasmy <- ListOfInputData_FileSheets

names(Input_Data_Mitophagia_Line_Heteroplasmy)
str(Input_Data_Mitophagia_Line_Heteroplasmy)

DataSet_Zhigmitova_Zuravlev <- Input_Data_Mitophagia_Line_Heteroplasmy
DataSet_Zhigmitova_Zuravlev
unique(DataSet_Zhigmitova_Zuravlev$Value_type)

#DataSet_Zhigmitova_Zuravlev$Mitophagy001 <- ifelse(DataSet_Zhigmitova_Zuravlev$Author == "Zhigmitova", DataSet_Zhigmitova_Zuravlev$Mitophagy/100, DataSet_Zhigmitova_Zuravlev$Mitophagy)

wilcox.test(DataSet_Zhigmitova_Zuravlev$Heteroplasmy ~ DataSet_Zhigmitova_Zuravlev$Author, paired = FALSE)
wilcox.test(DataSet_Zhigmitova_Zuravlev$Mitophagy ~ DataSet_Zhigmitova_Zuravlev$Author, paired = FALSE)

Plot_File_Name <- file.path(PathToResults, paste("Boxplots of comparation analysis heteroplasmy Zhigmitova vs Zuravlev_new ", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")


Figure_boxplot <- ggplot(DataSet_Zhigmitova_Zuravlev, aes(x = Mutation_name, y = Heteroplasmy, fill = Author)) +
  geom_boxplot(outlier.size = 3, position = position_dodge(width = 0.9)) +
  stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE,) +
  xlab(label = "Мутация") +
  ylab(label = "Гетероплазмия") +
  labs(fill = "Исследователь") +
  scale_fill_discrete(limits = c("Zhigmitova", "Zuravlev"), labels = c( "Жигмитова", "Никифоров")) +
  scale_y_continuous(breaks = seq(from = 0, to = max(DataSet_Zhigmitova_Zuravlev$Heteroplasmy) + 5, by = 10)) +
  # geom_point(aes(shape = "mean"),  alpha = 0) +
  # geom_point(aes(shape = "median"),  alpha = 0) +
  # geom_point(aes(shape = "outlier"),  alpha = 0) +
  # guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
  # #ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
  theme_AO +
  theme(axis.text.x = element_text(size = 12))

print(Figure_boxplot)
dev.off()

Plot_File_Name <- file.path(PathToResults, paste("Boxplots of comparation analysis basal mitophagy Vinokurov vs Nikiforov_new ", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")


Figure_boxplot <- ggplot(data = DataSet_Zhigmitova_Zuravlev[DataSet_Zhigmitova_Zuravlev$experiment_type == "базовая",], aes(x = Line, y = Mitophagy, fill = Author)) +
  geom_boxplot(outlier.size = 3, position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = basal_level_mitophagy, col = Author), size = 1) +
  stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE) +
  xlab(label = "Клеточная линия") +
  ylab(label = "Митофагия") +
  labs(fill = "Базальная\nмитофагия", col = "Средняя\nбазальная\nмитофагия" ) +
  scale_fill_discrete(limits = c("Zhigmitova", "Zuravlev"), labels = c( "Винокуров", "Никифоров")) +
  scale_color_discrete(limits = c("Zhigmitova", "Zuravlev"), labels = c( "Винокуров", "Никифоров")) +
  scale_y_continuous(breaks = seq(from = 0, to = max(DataSet_Zhigmitova_Zuravlev$Mitophagy, na.rm = TRUE) + 5, by = 5)) +
  # geom_point(aes(shape = "mean"),  alpha = 0) +
  # geom_point(aes(shape = "median"),  alpha = 0) +
  # geom_point(aes(shape = "outlier"),  alpha = 0) +
  # guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
  # #ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
  theme_AO +
  theme(axis.text.x = element_text(size = 8))

print(Figure_boxplot)
dev.off()

Plot_File_Name <- file.path(PathToResults, paste("Boxplots of comparation analysis FCCP mitophagy Vinokurov vs Nikiforov_new ", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")


Figure_boxplot <- ggplot(data = DataSet_Zhigmitova_Zuravlev[DataSet_Zhigmitova_Zuravlev$experiment_type == "стимулированная.FCCP",], aes(x = Line, y = Mitophagy, fill = Author)) +
  geom_boxplot(outlier.size = 3, position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = basal_level_mitophagy, col = Author), size = 1) +
  stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE) +
  xlab(label = "Клеточная линия") +
  ylab(label = "Митофагия") +
  labs(fill = "FCCP", col = "Средняя\nбазальная\nмитофагия" ) +
  scale_fill_discrete(limits = c("Zhigmitova", "Zuravlev"), labels = c( "Винокуров", "Никифоров")) +
  scale_color_discrete(limits = c("Zhigmitova", "Zuravlev"), labels = c( "Винокуров", "Никифоров")) +
  # geom_point(aes(shape = "mean"),  alpha = 0) +
  # geom_point(aes(shape = "median"),  alpha = 0) +
  # geom_point(aes(shape = "outlier"),  alpha = 0) +
  # guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
  # #ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
  theme_AO +
  theme(axis.text.x = element_text(size = 8))

print(Figure_boxplot)
dev.off()

##### Отрисовка сравнительных графиков с относительной митофагией

DataSet_basal <- aggregate(Mitophagy_relative_individual ~ Line + Author, data = DataSet_Zhigmitova_Zuravlev[DataSet_Zhigmitova_Zuravlev$experiment_type == "базовая",], FUN = mean)

Plot_File_Name <- file.path(PathToResults, paste("Barplot of comparation analysis relative basal mitophagy Vinokurov vs Nikiforov_new ", Sys.Date()+1, ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")


Figure_barplot <- ggplot(data = DataSet_basal, aes(x = Line, y = Mitophagy_relative_individual, fill = Author)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.8) +
  geom_hline(aes(yintercept = 1, col = Author), size = 1) +
  #stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE) +
  xlab(label = "Клеточная линия") +
  ylab(label = "Митофагия") +
  scale_y_continuous(breaks = seq(from = 0, to = max(DataSet_basal$Mitophagy_relative_) + 0.5, by = 0.2)) +
  labs(fill = "Относительная\nбазальная\nмитофагия", col = "Средняя\nбазальная\nмитофагия" ) +
  scale_fill_discrete(limits = c("Zhigmitova", "Zuravlev"), labels = c( "Винокуров", "Никифоров")) +
  scale_color_manual(limits = c("Zhigmitova"), labels = c(" "), values = c("darkgreen")) +
  # geom_point(aes(shape = "mean"),  alpha = 0) +
  # geom_point(aes(shape = "median"),  alpha = 0) +
  # geom_point(aes(shape = "outlier"),  alpha = 0) +
  # guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
  # #ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
  theme_AO +
  theme(axis.text.x = element_text(size = 8))

print(Figure_barplot)
dev.off()

DataSet_FCCP <- aggregate(Mitophagy_relative_individual ~ Line + Author, data = DataSet_Zhigmitova_Zuravlev[DataSet_Zhigmitova_Zuravlev$experiment_type == "стимулированная.FCCP",], FUN = mean)
DataSet_FCCP

Plot_File_Name <- file.path(PathToResults, paste("Barplot of comparation analysis relative FCCP mitophagy Vinokurov vs Nikiforov ", Sys.Date(), ".jpg", sep = "_"))
jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")


Figure_barplot <- ggplot(data = DataSet_FCCP, aes(x = Line, y = Mitophagy_relative_individual, fill = Author)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.8) +
  geom_hline(aes(yintercept = 1, col = Author), size = 1) +
  #stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE) +
  xlab(label = "Клеточная линия") +
  ylab(label = "Митофагия") +
  scale_y_continuous(breaks = seq(from = 0, to = max(DataSet_FCCP$Mitophagy_relative_individual) + 0.5, by = 0.2)) +
  labs(fill = "Относительная\nмитофагия\nпосле FCCP", col = "Средняя\nбазальная\nмитофагия" ) +
  scale_fill_discrete(limits = c("Zhigmitova", "Zuravlev"), labels = c( "Винокуров", "Никифоров")) +
  scale_color_manual(limits = c("Zhigmitova"), labels = c(" "), values = c("darkgreen")) +
  # geom_point(aes(shape = "mean"),  alpha = 0) +
  # geom_point(aes(shape = "median"),  alpha = 0) +
  # geom_point(aes(shape = "outlier"),  alpha = 0) +
  # guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
  # #ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
  theme_AO +
  theme(axis.text.x = element_text(size = 10))

print(Figure_barplot)
dev.off()

##################### РАБОТА С ДАННЫМИ ПО МИТОФАГИИ И ГЕТЕРОПЛАЗМИИ ОТ ЖУРАВЛЕВА И НИКИФОРОВА ОТ МАЯ 2022 ГОДА #########

################ Определение типа митофагии - дефектная или индуцируемая ###########



Kruskal_results_Table <- data.frame()
Kruskal_posthoc_Table <- data.frame()


for(Line_counter in seq_along(unique(DataSet_Zuravlev$Line)))
{
  temp_DataSet_Zuravlev <- subset.data.frame(DataSet_Zuravlev, DataSet_Zuravlev$Line == unique(DataSet_Zuravlev$Line)[Line_counter])
  kruskal_results <- kruskal.test(Mitophagy ~ experiment_type, data = temp_DataSet_Zuravlev)
  
  Base_mean_median <- mean_without_NA(temp_DataSet_Zuravlev$Mitophagy[temp_DataSet_Zuravlev$experiment_type == "базовая"])
 
  FCCP_mean_median <- mean_without_NA(temp_DataSet_Zuravlev$Mitophagy[temp_DataSet_Zuravlev$experiment_type == "стимулированная.FCCP"])
  
  
  temp_Kruskal_results_Table <- data.frame(
    "Line" = unique(DataSet_Zuravlev$Line)[Line_counter],
    "Название теста" = kruskal_results$method,
    "Название набора данных" = kruskal_results$data.name,
    "Значение статистики" = kruskal_results$statistic,
    "Базовая митофагия" = Base_mean_median,
    #"Пируват" = Pyruvate_mean_median,
    "FCCP" = FCCP_mean_median,
    "р-значение" = kruskal_results$p.value,
    "Результат" = ifelse(kruskal_results$p.value >= 0.05, "дефектная", "индуцируемая"),
    "Результат_уточненный" = NA
  )
  
  if(temp_Kruskal_results_Table$р.значение < 0.05)
  {
    kruskal_posthoc_results <- kruskalmc(Mitophagy ~ experiment_type, data = temp_DataSet_Zuravlev, probs = 0.05)
    
    names(kruskal_posthoc_results$dif.com) <- c("Уровень различий", "Критический уровень различий", "Достоверность различий")
    
    # Base_mean_median <- mean_without_NA(temp_DataSet_Zuravlev$Mitophagy[temp_DataSet_Zuravlev$experiment_type == "Базовая"])
    # Pyruvate_mean_median <- mean_without_NA(temp_DataSet_Zuravlev$Mitophagy[temp_DataSet_Zuravlev$experiment_type == "Пируват"])
    # FCCP_mean_median <- mean_without_NA(temp_DataSet_Zuravlev$Mitophagy[temp_DataSet_Zuravlev$experiment_type == "FCCP"])
    # 
    temp_kruskal_posthoc_results_Table <- data.frame(
      "Line" = unique(DataSet_Zuravlev$Line)[Line_counter],
      "Пары сравнения" = row.names(kruskal_posthoc_results$dif.com),
      kruskal_posthoc_results$dif.com[,1:3],
      # "Базовая митофагия" = Base_mean_median,
      # "Пируват" = Pyruvate_mean_median,
      # "FCCP" = FCCP_mean_median,
      "p-значение" = kruskal_posthoc_results$signif.level
    )
    Kruskal_posthoc_Table <- rbind.data.frame(temp_kruskal_posthoc_results_Table, Kruskal_posthoc_Table)
    
    temp_Kruskal_results_Table$Результат_уточненный <- ifelse((Base_mean_median >= FCCP_mean_median), "дефектная", "индуцируемая") 
  } else
  {
    temp_Kruskal_results_Table$Результат_уточненный  <- "дефектная"  
  }
  
  Kruskal_results_Table <- rbind.data.frame(temp_Kruskal_results_Table, Kruskal_results_Table)
  
}

Kruskal_results_Table
Kruskal_posthoc_Table

####### Сохранение таблицы результатов в Еxcel файл

List_Of_Result_Table <- list(Kruskal_results_Table, Kruskal_posthoc_Table)
names(List_Of_Result_Table) <- c("Сводная_таблица", "Парные_сравнения")
Table_ResultFileName <- file.path(PathToResults, paste("Определение типа митофагии по данным Журавлева и Никифорова", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

Mitophagy_type_Table <- Kruskal_results_Table[, c("Line", "Результат_уточненный")]

row.names(Mitophagy_type_Table) <- NULL

names(Mitophagy_type_Table) <- c("Line", "Тип митофагии")

DataSet_Zuravlev_Mutations_Mitophagy_type <- merge.data.frame(DataSet_Zuravlev, Mitophagy_type_Table, by = "Line")

########## Изучение связи мутаций с типом митофагии

unique(DataSet_Zuravlev_Mutations_Mitophagy_type$Mutation_name)

DataSet_Zuravlev_Mutations_Mitophagy_type <- DataSet_Zuravlev_Mutations_Mitophagy_type[DataSet_Zuravlev_Mutations_Mitophagy_type$Heteroplasmy > 5,]

Test_result <- data.frame()

for (Mutations_counter in seq_along(unique(DataSet_Zuravlev_Mutations_Mitophagy_type$Mutation_name)))
{
  tempData <- subset.data.frame(DataSet_Zuravlev_Mutations_Mitophagy_type, DataSet_Zuravlev_Mutations_Mitophagy_type$Mutation_name == unique(DataSet_Zuravlev_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter] & DataSet_Zuravlev_Mutations_Mitophagy_type$experiment_number == 1)
  
  Plot_File_Name <- file.path(PathToResults, paste("Boxplots of heteroplasmy by mitophagy type for mutation", str_extract(unique(DataSet_Zuravlev_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter], "[0-9]+"), Sys.Date(), ".jpg", sep = "_"))
  jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")
  
  
  Figure_boxplot <- ggplot(tempData, aes(x = tempData$`Тип митофагии`, y = Heteroplasmy)) +
    geom_boxplot(outlier.size = 3, position = position_dodge(width = 0.9)) +
    stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE,) +
    xlab(label = "Тип митофагии") +
    ylab(label = "Гетероплазмия") +
    # geom_point(aes(shape = "mean"),  alpha = 0) +
    # geom_point(aes(shape = "median"),  alpha = 0) +
    # geom_point(aes(shape = "outlier"),  alpha = 0) +
    # guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
    # #ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
    theme_AO
  
  print(Figure_boxplot)
  dev.off()
  
  Plot_File_Name <- file.path(PathToResults, paste("Boxplots of mitophagy by mitophagy type for mutation", str_extract(unique(DataSet_Zuravlev_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter], "[0-9]+"), Sys.Date(), ".jpg", sep = "_"))
  jpeg(filename = Plot_File_Name, width = 1000, height = 800, units = "px", pointsize =18, quality = 200, bg = "white")
  
  
  Figure_boxplot <- ggplot(tempData, aes(x = tempData$`Тип митофагии`, y = Mitophagy, fill = experiment_type)) +
    geom_boxplot(outlier.size = 3, position = position_dodge(width = 0.9)) +
    stat_summary(fun = mean, colour = "dodgerblue", geom = "point", shape = 18, size = 5, position = position_dodge(width = 0.9), show.legend = FALSE,) +
    xlab(label = "Тип митофагии") +
    ylab(label = "Митофагия") +
    labs(fill = "Тип митофагии") +
    # scale_fill_discrete(limits = c("Zhigmitova", "Zuravlev"), labels = c( "Винокуров", "Никифоров")) +
    # geom_point(aes(shape = "mean"),  alpha = 0) +
    # geom_point(aes(shape = "median"),  alpha = 0) +
    # geom_point(aes(shape = "outlier"),  alpha = 0) +
    # guides(shape = guide_legend(title = NULL, override.aes = list(alpha = 1, colour = c("black", "black", "dodgerblue"), shape = c(20, 95, 18), size = 5))) +
    # #ggtitle("Влияние LPS на накопление провоспалительных медиаторов") +
    theme_AO
  
  print(Figure_boxplot)
  dev.off()
  
  
  test_result <- wilcox.test(tempData$Heteroplasmy ~ tempData$`Тип митофагии`, paired = FALSE)
  
  Defect_mean <- mean_without_NA(tempData$Heteroplasmy[tempData$`Тип митофагии` == "дефектная"])
  Norm_mean <- mean_without_NA(tempData$Heteroplasmy[tempData$`Тип митофагии` == "индуцируемая"])
  
  temp_Test_result <- data.frame("Мутация" = unique(DataSet_Zuravlev_Mutations_Mitophagy_type$Mutation_name)[Mutations_counter], 
                                 "Название теста" = test_result$method,
                                 "Значение теста" = test_result$statistic,
                                 "Гетерозизотность_дефектной" = Defect_mean,
                                 "Гетерозизотность_индуцированной" = Norm_mean,
                                 "р-значение" = test_result$p.value,
                                 "Результат теста" = ifelse(test_result$p.value < 0.05, "Мутация связана", "Мутация не связана"))
  
  Test_result <- rbind.data.frame(Test_result, temp_Test_result)
  
}

Test_result

######### сохранение результатов проверки на тип распределения

List_Of_Result_Table <- list(Test_result)
names(List_Of_Result_Table) <-  c("Сводная_таблица")
Table_ResultFileName <- file.path(PathToResults, paste("Взаимосвязь гетероплазмии мутаций и типа митофагии по данным Журавлева и Никифорова", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, List_Of_Result_Table)

######### Исследование митофагия ~ гетерозиготность с помощью СЛАУ ################

unique(DataSet_Zhigmitova$Value_type)
names(DataSet_Zhigmitova)
DataSet_SOLE_Zhigmitova <- DataSet_Zhigmitova[, c("Line", "Mitophagy", "experiment_type", "Mutation_name", "Heteroplasmy")]

DataSet_SOLE_Zhigmitova <- DataSet_SOLE_Zhigmitova[DataSet_SOLE_Zhigmitova$experiment_type == "базовая",]

temp_Heteroplasmy_Table <- aggregate(Heteroplasmy ~ Line + Mutation_name, data = DataSet_SOLE_Zhigmitova, FUN = mean)

Heteroplasmy_Table <- spread(temp_Heteroplasmy_Table, Mutation_name, Heteroplasmy)

row.names(Heteroplasmy_Table) <- Heteroplasmy_Table$Line
Heteroplasmy_Table$Line <- NULL

Mitophagy_Table <- aggregate(Mitophagy ~ Line, data = DataSet_SOLE_Zhigmitova, FUN = mean)

row.names(Mitophagy_Table) <- Mitophagy_Table$Line
Mitophagy_Table$Line <- NULL

Heteroplasmy_Table
Mitophagy_Table

Modeling_Results_List <- get_SOLE_model_results(A_matrix_dataframe = Heteroplasmy_Table, B_vector_dataframe = Mitophagy_Table)

if(!is.na(Modeling_Results_List))
  {
   message("Моделирование выполнено.")
  } else {
          message("Моделирование не прошло!")
  }

Modeling_Results_List

###### Собственно сохранение результатов работы моделирования

Table_ResultFileName <- file.path(PathToResults, paste("Решение СЛАУ для митофагии по гетероплазмии по данным Винокурова и Жегмитовой", Sys.Date(), ".xlsx", sep = "_"))
save_to_excel_file_by_XLConnect(Table_ResultFileName, Modeling_Results_List)
