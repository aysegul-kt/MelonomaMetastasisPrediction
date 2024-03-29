library("BBmisc") #, lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library(readr)
library("caret")
library("stringr", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library(stringi)
library("e1071")
library("rlist")
library(data.table)
library(uuid)
library(Rmisc)

source('~/Google Drive/Studies/Thesis/AysegulPhdThesis/MelonomaMetastasisPrediction/src/helpers/libs.R')
source('~/Google Drive/Studies/Thesis/AysegulPhdThesis/MelonomaMetastasisPrediction/src/preprocess/dataMergeService.R')
source('~/Google Drive/Studies/Thesis/AysegulPhdThesis/MelonomaMetastasisPrediction/src/helpers/logger.R')
source('~/Google Drive/Studies/Thesis/AysegulPhdThesis/MelonomaMetastasisPrediction/src/preprocess/normalizationService.R')
source('~/Google Drive/Studies/Thesis/AysegulPhdThesis/MelonomaMetastasisPrediction/src/helpers/SignificancyDataCollection.R')
source('~/Google Drive/Studies/Thesis/AysegulPhdThesis/MelonomaMetastasisPrediction/src/training/trainingController.R')