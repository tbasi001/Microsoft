# clearing environment
rm(list=ls())

# set the working directory
setwd("/Users/csuftitan/Desktop/Prolific Study/Data Final")
getwd()

# loading packages
library(psych)
library(GPArotation)
library(mice)
library(rmarkdown)
library(careless)
library(nnet)
library(dominanceanalysis)
library(profileR)
library(ggpubr)
library(ggplot2)
library(EFAtools)

# importing data
data <- read.csv("Final.csv", header=TRUE, na.strings="", 
                 stringsAsFactors = TRUE)

write.csv(data, file = "Final.csv", row.names = F)

# missing data imputation 
impute <- mice(data, m = 1)
data <- complete(impute, 1)
sum(is.na(data))

# attention checks
data$ac1 <- abs((data$fsc.8.ac1a - (8-data$fsc.20.ac1b)))/7
data$ac2 <- abs((data$fsc.28.ac2a - (8-data$fsc.2.ac2b)))/7
data$ac3 <- abs((data$fsc.1.ac3a - (8-data$fsc.18.ac3b)))/7
data$ac4 <- abs((data$fsc.23.ac4a - (8-data$fsc.14.ac4b)))/7
data$ac5 <- abs((data$ccso.3.ac1a - (6-data$ccs.ac1b)))/5
data$ac6 <- abs((data$ccsi.52.ac2a - (6-data$ccs.ac2b)))/5
data$ac7 <- abs((data$dg.pe2.ac1a - (6-data$dg.ac1b)))/5
data$ac8 <- abs((data$hd.5.ac1a - (6-data$hd.ac1b)))/5
data$ac9 <- abs((data$bfico.33.ac1a - (6-data$bfi.ac1b)))/5
data$ac10 <- abs((data$bfiee.11r.ac2a - (6-data$bfi.ac2b)))/5
data$ac11 <- abs((data$bfiat.12r.ac3a - (6-data$bfi.ac3b)))/5
data$ac12 <- abs((data$bfind.24r.ac4a - (6-data$bfi.ac4b)))/5
data$ac13 <- abs((data$bfioa.5r.ac5a - (6-data$bfi.ac5b)))/5
data$ac14 <- abs((data$hexo.1.ac1a - (6-data$hex.ac1b)))/5
data$ac15 <- abs((data$hexg.3.ac1a - (6-data$hex.ac1b.1)))/5
data$ac16 <- abs((data$ipip.2.ac1a - (6-data$ipip.ac1b)))/5

# data$ac1pomp <- (data$ac1 / 7)
# data$ac2pomp <- (data$ac2 / 7)
# data$ac3pomp <- (data$ac3 / 7)
# data$ac4pomp <- (data$ac4 / 5)
# data$ac5pomp <- (data$ac5 / 5)
# data$ac6pomp <- (data$ac6 / 5)
# data$ac7pomp <- (data$ac7 / 5)
# data$ac8pomp <- (data$ac8 / 5)
# data$ac9pomp <- (data$ac9 / 5)
# data$ac10pomp <- (data$ac10 / 5)
# data$ac11pomp <- (data$ac11 / 5)
# data$ac12pomp <- (data$ac12 / 5)
# data$ac13pomp <- (data$ac13 / 5)
# data$ac14pomp <- (data$ac14 / 5)
# data$ac15pomp <- (data$ac15 / 5)


data$acpompsum <- data$ac1 + data$ac2 + data$ac3 + data$ac4 + data$ac5 + data$ac6 + 
  data$ac7 + data$ac8 + data$ac9 + data$ac10 + data$ac11 + data$ac12 + data$ac13 + 
  data$ac14 + data$ac15 + data$ac16

data$acpompsumz <- scale(data$acpompsum)

data$ac1f <- ifelse((data$fsc.8.ac1a > 4 & data$fsc.20.ac1b > 4) | (data$fsc.8.ac1a < 4 & data$fsc.20.ac1b < 4), 1, 0)
data$ac2f <- ifelse((data$fsc.28.ac2a > 4 & data$fsc.2.ac2b > 4) | (data$fsc.28.ac2a < 4 & data$fsc.2.ac2b < 4), 1, 0)
data$ac3f <- ifelse((data$fsc.1.ac3a > 4 & data$fsc.18.ac3b > 4) | (data$fsc.1.ac3a < 4 & data$fsc.18.ac3b < 4), 1, 0)
data$ac4f <- ifelse((data$fsc.23.ac4a > 4 & data$fsc.14.ac4b > 4) | (data$fsc.23.ac4a < 4 & data$fsc.14.ac4b < 4), 1, 0)
data$ac5f <- ifelse((data$ccso.3.ac1a > 3 & data$ccs.ac1b > 3) | (data$ccso.3.ac1a < 3 & data$ccs.ac1b < 3), 1, 0)
data$ac6f <- ifelse((data$ccsi.52.ac2a > 3 & data$ccs.ac2b > 3) | (data$ccsi.52.ac2a < 3 & data$ccs.ac2b < 3), 1, 0)
data$ac7f <- ifelse((data$dg.pe2.ac1a > 3 & data$dg.ac1b > 3) | (data$dg.pe2.ac1a < 3 & data$dg.ac1b < 3), 1, 0)
data$ac8f <- ifelse((data$hd.5.ac1a > 3 & data$hd.ac1b > 3) | (data$hd.5.ac1a < 3 & data$hd.ac1b < 3), 1, 0)
data$ac9f <- ifelse((data$bfico.33.ac1a > 3 & data$bfi.ac1b > 3) | (data$bfico.33.ac1a < 3 & data$bfi.ac1b < 3), 1, 0)
data$ac10f <- ifelse((data$bfiee.11r.ac2a > 3 & data$bfi.ac2b > 3) | (data$bfiee.11r.ac2a < 3 & data$bfi.ac2b < 3), 1, 0)
data$ac11f <- ifelse((data$bfiat.12r.ac3a > 3 & data$bfi.ac3b > 3) | (data$bfiat.12r.ac3a < 3 & data$bfi.ac3b < 3), 1, 0)
data$ac12f <- ifelse((data$bfind.24r.ac4a > 3 & data$bfi.ac4b > 3) | (data$bfind.24r.ac4a < 3 & data$bfi.ac4b < 3), 1, 0)
data$ac13f <- ifelse((data$bfioa.5r.ac5a > 3 & data$bfi.ac5b > 3) | (data$bfioa.5r.ac5a < 3 & data$bfi.ac5b < 3), 1, 0)
data$ac14f <- ifelse((data$hexo.1.ac1a > 3 & data$hex.ac1b > 3) | (data$hexo.1.ac1a < 3 & data$hex.ac1b < 3), 1, 0)
data$ac15f <- ifelse((data$hexg.3.ac1a > 3 & data$hex.ac1b.1 > 3) | (data$hexg.3.ac1a < 3 & data$hex.ac1b.1 < 3), 1, 0)
data$ac16f <- ifelse((data$ipip.2.ac1a > 3 & data$ipip.ac1b > 3) | (data$ipip.2.ac1a < 3 & data$ipip.ac1b < 3), 1, 0)

data$acsumf <- data$ac1f + data$ac2f + data$ac3f + data$ac4f + data$ac5f + data$ac6f + data$ac7f + data$ac8f + 
  data$ac9f + data$ac10f + data$ac11f + data$ac12f + data$ac13f + data$ac14f + data$ac15f + data$ac16f

data$acsumfz <- scale(data$acsumf)

# identifying individuals who answer with long strings of consecutive numbers
data$ls <- longstring(data)
data$lsz <- scale(data$ls)

# reverse coding
data[, c("ccso.1r", "ccso.5r", "ccso.6r", "ccso.7r", "ccso.9r", "ccsv.11r", "ccsv.12r", "ccsv.13r", "ccsv.14r", 
         "ccst.22r", "ccst.24r", "ccst.27r", "ccst.29r", "ccss.31r", "ccss.33r", "ccss.35r", "ccss.36r",
         "ccss.37r", "ccsr.44r", "ccsr.46r", "ccsr.47r", "ccsr.50r", "ccsi.53r", "ccsi.54r", "ccsi.58r", 
         "ccsi.59r", "ccsi.60r")] <-
  6-data[, c("ccso.1r", "ccso.5r", "ccso.6r", "ccso.7r", "ccso.9r", "ccsv.11r", "ccsv.12r", "ccsv.13r", "ccsv.14r", 
             "ccst.22r", "ccst.24r", "ccst.27r", "ccst.29r", "ccss.31r", "ccss.33r", "ccss.35r", "ccss.36r",
             "ccss.37r", "ccsr.44r", "ccsr.46r", "ccsr.47r", "ccsr.50r", "ccsi.53r", "ccsi.54r", "ccsi.58r", 
             "ccsi.59r", "ccsi.60r")]

data[, c("dg.ci1r", "dg.ci2r", "dg.ci3r", "dg.ci4r")] <-
  6-data[, c("dg.ci1r", "dg.ci2r", "dg.ci3r", "dg.ci4r")]

data[, c("bfico.3r", "bfina.4r", "bfioa.5r.ac5a", "bficp.8r", "bfind.9r", "bfiee.11r.ac2a", "bfiat.12r.ac3a", 
         "bfies.16r", "bfiac.17r", "bfiar.22r", "bficp.23r", "bfind.24r.ac4a", "bfioi.25r", "bfiee.26r", "bficr.28r", 
         "bfine.29r", "bfioc.30r", "bfies.31r", "bfiea.36r", "bfiar.37r", "bfiat.42r", "bfine.44r", "bfioc.45r", 
         "bfiac.47r", "bfico.48r", "bfina.49r", "bfioa.50r", "bfiea.51r", "bfioi.55r", "bficr.58r")] <-
  6-data[, c("bfico.3r", "bfina.4r", "bfioa.5r.ac5a", "bficp.8r", "bfind.9r", "bfiee.11r.ac2a", "bfiat.12r.ac3a", 
             "bfies.16r", "bfiac.17r", "bfiar.22r", "bficp.23r", "bfind.24r.ac4a", "bfioi.25r", "bfiee.26r", "bficr.28r", 
             "bfine.29r", "bfioc.30r", "bfies.31r", "bfiea.36r", "bfiar.37r", "bfiat.42r", "bfine.44r", "bfioc.45r", 
             "bfiac.47r", "bfico.48r", "bfina.49r", "bfioa.50r", "bfiea.51r", "bfioi.55r", "bficr.58r")]

data[, c("hexr.4r", "hexe.7r", "hexr.8r", "hexo.9r", "hexd.10r", "hexo.13r", "hexd.14r", "hexr.16r")] <-
  6-data[, c("hexr.4r", "hexe.7r", "hexr.8r", "hexo.9r", "hexd.10r", "hexo.13r", "hexd.14r", "hexr.16r")]

data[, c("hexs.1r", "hexf.2r", "hexf.6r", "hexg.7r", "hexs.9r", "hexg.11r", "hexm.12r", "hexf.14r", 
         "hexg.15r", "hexm.16r")] <-
  6-data[, c("hexs.1r", "hexf.2r", "hexf.6r", "hexg.7r", "hexs.9r", "hexg.11r", "hexm.12r", "hexf.14r", 
             "hexg.15r", "hexm.16r")]

data[, c("ipip.1r", "ipip.3r", "ipip.5r", "ipip.8r", "ipip.9r")] <-
  6-data[, c("ipip.1r", "ipip.3r", "ipip.5r", "ipip.8r", "ipip.9r")]

# composite variables
data$ccso <- (data$ccso.1r + data$ccso.2 + data$ccso.3.ac1a + data$ccso.4 + data$ccso.5r + data$ccso.6r + 
                data$ccso.7r + data$ccso.8 + data$ccso.9r + data$ccso.10)/10

ccso <- data[,c("ccso.1r", "ccso.2", "ccso.3.ac1a", "ccso.4", "ccso.5r", "ccso.6r", "ccso.7r", "ccso.8", 
                "ccso.9r", "ccso.10")]
psych::alpha(ccso)

data$ccsv <- (data$ccsv.11r + data$ccsv.12r + data$ccsv.13r + data$ccsv.14r + data$ccsv.15 + data$ccsv.16 + 
                data$ccsv.17 + data$ccsv.18 + data$ccsv.19 + data$ccsv.20)/10

ccsv <- data[,c("ccsv.11r", "ccsv.12r", "ccsv.13r", "ccsv.14r", "ccsv.15", "ccsv.16", "ccsv.17", "ccsv.18", 
                "ccsv.19", "ccsv.20")]
psych::alpha(ccsv)

data$ccst <- (data$ccst.21 + data$ccst.22r + data$ccst.23 + data$ccst.24r + data$ccst.25 + data$ccst.26 + 
                data$ccst.27r + data$ccst.28 + data$ccst.29r + data$ccst.30)/10

ccst <- data[,c("ccst.21", "ccst.22r", "ccst.23", "ccst.24r", "ccst.25", "ccst.26", "ccst.27r", "ccst.28", 
                "ccst.29r", "ccst.30")]
psych::alpha(ccst)

data$ccss <- (data$ccss.31r + data$ccss.32 + data$ccss.33r + data$ccss.34 + data$ccss.35r + data$ccss.36r + 
                data$ccss.37r + data$ccss.38 + data$ccss.39 + data$ccss.40)/10

ccss <- data[,c("ccss.31r", "ccss.32", "ccss.33r", "ccss.34", "ccss.35r", "ccss.36r", "ccss.37r", "ccss.38", 
                "ccss.39", "ccss.40")]
psych::alpha(ccss)

data$ccsr <- (data$ccsr.41 + data$ccsr.42 + data$ccsr.43 + data$ccsr.44r + data$ccsr.45 + data$ccsr.46r + 
                data$ccsr.47r + data$ccsr.48 + data$ccsr.49 + data$ccsr.50r)/10

ccsr <- data[,c("ccsr.41", "ccsr.42", "ccsr.43", "ccsr.44r", "ccsr.45", "ccsr.46r", "ccsr.47r", "ccsr.48", 
                "ccsr.49", "ccsr.50r")]
psych::alpha(ccsr)

data$ccsi <- (data$ccsi.51 + data$ccsi.52.ac2a + data$ccsi.53r + data$ccsi.54r + data$ccsi.55 + data$ccsi.56 + 
                data$ccsi.57 + data$ccsi.58r + data$ccsi.59r + data$ccsi.60r)/10

ccsi <- data[,c("ccsi.51", "ccsi.52.ac2a", "ccsi.53r", "ccsi.54r", "ccsi.55", "ccsi.56", "ccsi.57", "ccsi.58r", 
                "ccsi.59r", "ccsi.60r")]
psych::alpha(ccsi)

data$dg <- (data$dg.pe1 + data$dg.pe2.ac1a + data$dg.pe3 + data$dg.pe4 + data$dg.ci1r + data$dg.ci2r + 
              data$dg.ci3r + data$dg.ci4r)/8
data$dgpe <- (data$dg.pe1 + data$dg.pe2.ac1a + data$dg.pe3 + data$dg.pe4)/4

dgpe <- data[,c("dg.pe1", "dg.pe2.ac1a", "dg.pe3", "dg.pe4")]
psych::alpha(dgpe)

data$dgci <- (data$dg.ci1r + data$dg.ci2r + data$dg.ci3r + data$dg.ci4r)/4

dgci <- data[,c("dg.ci1r", "dg.ci2r", "dg.ci3r", "dg.ci4r")]
psych::alpha(dgci)

data$hd <- (data$hd.1 + data$hd.2 + data$hd.3 + data$hd.4 + data$hd.5.ac1a + data$hd.6 + data$hd.7 + data$hd.8 + 
              data$hd.9 + data$hd.10)/10

hd <- data[,c("hd.1", "hd.2", "hd.3", "hd.4", "hd.5.ac1a", "hd.6", "hd.7", "hd.8", "hd.9", "hd.10")]
psych::alpha(hd)

data$bfie <- (data$bfies.1 + data$bfiea.6 + data$bfiee.11r.ac2a + data$bfies.16r + data$bfiea.21 + 
                data$bfiee.26r + data$bfies.31r + data$bfiea.36r + data$bfiee.41 + data$bfies.46 + 
                data$bfiea.51r + data$bfiee.56)/12
data$bfies <- (data$bfies.1 + data$bfies.16r + data$bfies.31r + data$bfies.46)/4

bfies <- data[,c("bfies.1", "bfies.16r", "bfies.31r", "bfies.46")]
psych::alpha(bfies)

data$bfiea <- (data$bfiea.6 + data$bfiea.21 + data$bfiea.36r + data$bfiea.51r)/4

bfiea <- data[,c("bfiea.6", "bfiea.21", "bfiea.36r", "bfiea.51r")]
psych::alpha(bfiea)

data$bfiee <- (data$bfiee.11r.ac2a + data$bfiee.26r + data$bfiee.41 + data$bfiee.56)/4

bfiee <- data[,c("bfiee.11r.ac2a", "bfiee.26r", "bfiee.41", "bfiee.56")]
psych::alpha(bfiee)

data$bfia <- (data$bfiac.2 + data$bfiar.7 + data$bfiat.12r.ac3a + data$bfiac.17r + data$bfiar.22r + 
                data$bfiat.27 + data$bfiac.32 + data$bfiar.37r + data$bfiat.42r + data$bfiac.47r + 
                data$bfiar.52 + data$bfiat.57)/12
data$bfiac <- (data$bfiac.2 + data$bfiac.17r + data$bfiac.32 + data$bfiac.47r)/4

bfiac <- data[,c("bfiac.2", "bfiac.17r", "bfiac.32", "bfiac.47r")]
psych::alpha(bfiac)

data$bfiar <- (data$bfiar.7 + data$bfiar.22r + data$bfiar.37r + data$bfiar.52)/4

bfiar <- data[,c("bfiar.7", "bfiar.22r", "bfiar.37r", "bfiar.52")]
psych::alpha(bfiar)

data$bfiat <- (data$bfiat.12r.ac3a + data$bfiat.27 + data$bfiat.42r + data$bfiat.57)/4

bfiat <- data[,c("bfiat.12r.ac3a", "bfiat.27", "bfiat.42r", "bfiat.57")]
psych::alpha(bfiat)

data$bfic <- (data$bfico.3r + data$bficp.8r + data$bficr.13 + data$bfico.18 + data$bficp.23r + 
                data$bficr.28r + data$bfico.33.ac1a + data$bficp.38 + data$bficr.43 + data$bfico.48r + 
                data$bficp.53 + data$bficr.58r)/12
data$bfico <- (data$bfico.3r + data$bfico.18 + data$bfico.33.ac1a + data$bfico.48r)/4

bfico <- data[,c("bfico.3r", "bfico.18", "bfico.33.ac1a", "bfico.48r")]
psych::alpha(bfico)

data$bficp <- (data$bficp.8r + data$bficp.23r + data$bficp.38 + data$bficp.53)/4

bficp <- data[,c("bficp.8r", "bficp.23r", "bficp.38", "bficp.53")]
psych::alpha(bficp)

data$bficr <- (data$bficr.13 + data$bficr.28r + data$bficr.43 + data$bficr.58r)/4

bficr <- data[,c("bficr.13", "bficr.28r", "bficr.43", "bficr.58r")]
psych::alpha(bficr)

data$bfin <- (data$bfina.4r + data$bfind.9r + data$bfine.14 + data$bfina.19 + data$bfind.24r.ac4a + 
                data$bfine.29r + data$bfina.34 + data$bfind.39 + data$bfine.44r + data$bfina.49r + 
                data$bfind.54 + data$bfine.59)/12
data$bfina <- (data$bfina.4r + data$bfina.19 + data$bfina.34 + data$bfina.49r)/4

bfina <- data[,c("bfina.4r", "bfina.19", "bfina.34", "bfina.49r")]
psych::alpha(bfina)

data$bfind <- (data$bfind.9r + data$bfind.24r.ac4a + data$bfind.39 + data$bfind.54)/4

bfind <- data[,c("bfind.9r", "bfind.24r.ac4a", "bfind.39", "bfind.54")]
psych::alpha(bfind)

data$bfine <- (data$bfine.14 + data$bfine.29r + data$bfine.44r + data$bfine.59)/4

bfine <- data[,c("bfine.14", "bfine.29r", "bfine.44r", "bfine.59")]
psych::alpha(bfine)

data$bfio <- (data$bfioa.5r.ac5a + data$bfioi.10 + data$bfioc.15 + data$bfioa.20 + data$bfioi.25r + 
                data$bfioc.30r + data$bfioa.35 + data$bfioi.40 + data$bfioc.45r + data$bfioa.50r + 
                data$bfioi.55r + data$bfioc.60)/12
data$bfioa <- (data$bfioa.5r.ac5a + data$bfioa.20 + data$bfioa.35 + data$bfioa.50r)/4

bfioa <- data[,c("bfioa.5r.ac5a", "bfioa.20", "bfioa.35", "bfioa.50r")]
psych::alpha(bfioa)

data$bfioi <- (data$bfioi.10 + data$bfioi.25r + data$bfioi.40 + data$bfioi.55r)/4

bfioi <- data[,c("bfioi.10", "bfioi.25r", "bfioi.40", "bfioi.55r")]
psych::alpha(bfioi)

data$bfioc <- (data$bfioc.15 + data$bfioc.30r + data$bfioc.45r + data$bfioc.60)/4

bfioc <- data[,c("bfioc.15", "bfioc.30r", "bfioc.45r", "bfioc.60")]
psych::alpha(bfioc)

data$hexc <- (data$hexo.1.ac1a + data$hexd.2 + data$hexe.3 + data$hexr.4r + data$hexo.5 + 
                data$hexd.6 + data$hexe.7r + data$hexr.8r + data$hexo.9r + data$hexd.10r + 
                data$hexe.11 + data$hexr.12 + data$hexo.13r + data$hexd.14r + data$hexe.15 + 
                data$hexr.16r)/16
data$hexo <- (data$hexo.1.ac1a + data$hexo.5 + data$hexo.9r + data$hexo.13r)/4

hexo <- data[,c("hexo.1.ac1a", "hexo.5", "hexo.9r", "hexo.13r")]
psych::alpha(hexo)

data$hexd <- (data$hexd.2 + data$hexd.6 + data$hexd.10r + data$hexd.14r)/4

hexd <- data[,c("hexd.2", "hexd.6", "hexd.10r", "hexd.14r")]
psych::alpha(hexd)

data$hexe <- (data$hexe.3 + data$hexe.7r + data$hexe.11 + data$hexe.15)/4

hexe <- data[,c("hexe.3", "hexe.7r", "hexe.11", "hexe.15")]
psych::alpha(hexe)

data$hexr <- (data$hexr.4r + data$hexr.8r + data$hexr.12 + data$hexr.16r)/4

hexr <- data[,c("hexr.4r", "hexr.8r", "hexr.12", "hexr.16r")]
psych::alpha(hexr)

data$hexh <- (data$hexs.1r + data$hexf.2r + data$hexg.3.ac1a + data$hexm.4 + data$hexs.5 + 
                data$hexf.6r + data$hexg.7r + data$hexm.8 + data$hexs.9r + data$hexf.10 + 
                data$hexg.11r + data$hexm.12r + data$hexs.13 + data$hexf.14r + data$hexg.15r + 
                data$hexm.16r)/16
data$hexs <- (data$hexs.1r + data$hexs.5 + data$hexs.9r + data$hexs.13)/4

hexs <- data[,c("hexs.1r", "hexs.5", "hexs.9r", "hexs.13")]
psych::alpha(hexs)

data$hexf <- (data$hexf.2r + data$hexf.6r + data$hexf.10 + data$hexf.14r)/4

hexf <- data[,c("hexf.2r", "hexf.6r", "hexf.10", "hexf.14r")]
psych::alpha(hexf)

data$hexg <- (data$hexg.3.ac1a + data$hexg.7r + data$hexg.11r + data$hexg.15r)/4

hexg <- data[,c("hexg.3.ac1a", "hexg.7r", "hexg.11r", "hexg.15r")]
psych::alpha(hexg)

data$hexm <- (data$hexm.4 + data$hexm.8 + data$hexm.12r + data$hexm.16r)/4

hexm <- data[,c("hexm.4", "hexm.8", "hexm.12r", "hexm.16r")]
psych::alpha(hexm)

data$ipip <- (data$ipip.1r + data$ipip.2.ac1a + data$ipip.3r + data$ipip.4 + data$ipip.5r + 
                data$ipip.6 + data$ipip.7 + data$ipip.8r + data$ipip.9r + data$ipip.10)/10

ipip <- data[,c("ipip.1r", "ipip.2.ac1a", "ipip.3r", "ipip.4", "ipip.5r", "ipip.6", "ipip.7", 
                "ipip.8r", "ipip.9r", "ipip.10")]
psych::alpha(ipip)

# removing poor attention partcipants
data <- subset(data, data$trust != 1)
data <- subset(data, data$ls < 50)
data <- subset(data, data$acsumf < 8)

# EFA of FSC

fscfa <- data[, c("fsc.1.ac3a", "fsc.3", "fsc.4", "fsc.5", "fsc.6", "fsc.7", "fsc.8.ac1a","fsc.9", "fsc.10", 
                  "fsc.11", "fsc.12", "fsc.13", "fsc.15", "fsc.16", "fsc.17", "fsc.19", "fsc.21", "fsc.22", 
                  "fsc.23.ac4a", "fsc.24", "fsc.25", "fsc.26", "fsc.27", "fsc.28.ac2a", "fsc.29", "fsc.30")]

# "fsc.31", "fsc.32", "fsc.33", "fsc.34", "fsc.35", "fsc.36", Stop Signal items
# "fsc.37", "fsc.38", "fsc.39", "fsc.40", "fsc.41", "fsc.42")] Affective forecasting items

eigen.fscfa <- eigen(cor(fscfa))$values
eigen(cor(fscfa))$values

write.csv(eigen.fscfa, file = "eigenvalues fsc fa study 4.csv", row.names = T)

plot(eigen(cor(fscfa))$values)
scree(fscfa, factors = FALSE)

fscfa3 <- fa(fscfa, 3, rotate = "oblimin", fm = "ml")
print(fscfa3, sort = TRUE, cutoff = 0)

fscfasort <- fa.sort(fscfa3)
write.csv(fscfasort$loadings, file = "sorted fsc fa.csv", row.names = T)

fscfa4 <- fa(fscfa, 4, rotate = "oblimin", fm = "ml")
print(fscfa4, sort = TRUE, cutoff = 0)

fscfa5 <- fa(fscfa, 5, rotate = "oblimin", fm = "ml")
print(fscfa5, sort = TRUE, cutoff = 0)

## No Affective Forecasting 
fscfanoaf <- data[, c("fsc.1.ac3a", "fsc.3", "fsc.4", "fsc.5", "fsc.6", "fsc.7", "fsc.8.ac1a","fsc.9", "fsc.10", 
                      "fsc.11", "fsc.12", "fsc.13", "fsc.15", "fsc.16", "fsc.17", "fsc.19", "fsc.21", "fsc.22", 
                      "fsc.23.ac4a", "fsc.24", "fsc.25", "fsc.26", "fsc.27", "fsc.28.ac2a", "fsc.29", "fsc.30", 
                      "fsc.31", "fsc.32", "fsc.33", "fsc.34", "fsc.35", "fsc.36")]

eigen(cor(fscfanoaf))$values
plot(eigen(cor(fscfanoaf))$values)
scree(fscfanoaf, factors = FALSE)

fscfanoaf4 <- fa(fscfanoaf, 4, rotate = "oblimin", fm = "ml")
print(fscfanoaf4, sort = TRUE, cutoff = 0)

fscfanoaf3 <- fa(fscfanoaf, 3, rotate = "oblimin", fm = "ml")
print(fscfanoaf3, sort = TRUE, cutoff = 0)

## No Affective Forecasting or Stop Signal Sensitivity 
fscfanoafss <- data[, c("fsc.1.ac3a", "fsc.3", "fsc.4", "fsc.5", "fsc.6", "fsc.7", "fsc.8.ac1a","fsc.9", "fsc.10", 
                        "fsc.11", "fsc.12", "fsc.13", "fsc.15", "fsc.16", "fsc.17", "fsc.19", "fsc.21", "fsc.22", 
                        "fsc.23.ac4a", "fsc.24", "fsc.25", "fsc.26", "fsc.27", "fsc.28.ac2a", "fsc.29", "fsc.30")]

eigen(cor(fscfanoafss))$values
plot(eigen(cor(fscfanoafss))$values)
scree(fscfanoafss, factors = FALSE)

fscfanoafss4 <- fa(fscfanoafss, 4, rotate = "oblimin", fm = "ml")
print(fscfanoafss4, sort = TRUE, cutoff = 0)

fscfanoafss3 <- fa(fscfanoafss, 3, rotate = "oblimin", fm = "ml")
print(fscfanoafss3, sort = TRUE, cutoff = 0)
fscefa3 <- print(fscfanoafss3, sort = TRUE, cutoff = 0)

fsctest <- factanal(fscfanoafss, factors = 3, rotation = "oblimin")
fsctestprint <- print(fsctest, sort = TRUE, cutoff = .4)

fasorted <- fa.sort(fsctest, polar = FALSE)

write.csv(fscfanoafss3$loadings, file = "EFA of 30 Q FSC.csv", row.names = T)
write.csv(fsctestprint$loadings, file = "EFA of 30 Q FSC Print.csv", row.names = T)

N_FACTORS(as.data.frame(fscfa3))

# FSC Factors
data$fsc1 <- (data$fsc.11 + data$fsc.7 + data$fsc.24 + data$fsc.26 + data$fsc.29 + data$fsc.22 + 
                data$fsc.27 + data$fsc.21 + data$fsc.4 + data$fsc.8.ac1a + data$fsc.30 + 
                data$fsc.1.ac3a + data$fsc.9)/13

data$fsc2 <- (data$fsc.19 + data$fsc.17 + data$fsc.23.ac4a + data$fsc.12 + data$fsc.3)/5

data$fsc3 <- (data$fsc.10 + data$fsc.13 + data$fsc.28.ac2a)/3

# Stop Signal data$fsc4 <- (data$fsc.31 + data$fsc.32 + data$fsc.33 + data$fsc.34 + data$fsc.36)/5

# FSC Factors with No Affective Forecasting or Stop Signal Sensitivity 
data$fsc1b <- (data$fsc.11 + data$fsc.7 + data$fsc.24 + data$fsc.26 + data$fsc.29 + data$fsc.22 +
                data$fsc.21 + data$fsc.27 + data$fsc.4 + data$fsc.8.ac1a + data$fsc.1.ac3a + 
                data$fsc.9 + data$fsc.30 + data$fsc.16 + data$fsc.6)/15
data$fsc2b <- (data$fsc.23.ac4a + data$fsc.17 + data$fsc.19 + data$fsc.12 + data$fsc.3)/5
data$fsc3b <- (data$fsc.28.ac2a + data$fsc.13 + data$fsc.10)/3

correlation <- cor(data.frame(data$fsc1b, data$fsc2b, data$fsc3b, data$ccss), use = "pairwise")
print(correlation)

# correlation 
correlation <- cor(data.frame(data$ccso, data$ccsv, data$ccst, data$ccss, data$ccsr, data$ccsi, data$dg, data$dgpe, 
                              data$dgci, data$hd, data$bfie, data$bfies, data$bfiea, data$bfiee, data$bfia, 
                              data$bfiac, data$bfiar, data$bfiat, data$bfic, data$bfico, data$bficp, 
                              data$bficr, data$bfin, data$bfina, data$bfind, data$bfine, data$bfio, 
                              data$bfioa, data$bfioi, data$bfioc, data$hexc, data$hexo, data$hexd, 
                              data$hexe, data$hexr, data$hexh, data$hexs, data$hexf, data$hexg, data$hexm), use="pairwise")

write.csv(round(correlation, 2), file = "Correlation CCSS and FSC.csv", row.names = T)

# Conscientiousness EFA

confa <- data[, c("fsc1", "fsc2", "fsc3", "ccso", "ccsv", "ccst", "ccss", "ccsr", "ccsi", 
                  "dgpe", "dgci", "hd", "bfico", "bficp", "bficr", "hexo", "hexd", "hexe", "hexr", "ipip")]

eigen(cor(confa))$values
plot(eigen(cor(confa))$values)

confa3 <- fa(confa, 3, rotate = "oblimin", fm = "ml")
print(confa3, sort = TRUE, cutoff = 0)

confa4 <- fa(confa, 4, rotate = "oblimin", fm = "ml")
print(confa4, sort = TRUE, cutoff = 0)

confa5 <- fa(confa, 5, rotate = "oblimin", fm = "ml")
print(confa5, sort = TRUE, cutoff = 0)

confasort <- fa.sort(confa5)
write.csv(confasort$loadings, file = "sorted con fa 5.csv", row.names = T)

confa6 <- fa(confa, 6, rotate = "oblimin", fm = "ml")
print(confa6, sort = TRUE, cutoff = 0)

bfi <- data[, c("fsc1", "fsc2", "fsc3", "bfiea", "bfiee", "bfies", "bfiac", "bfiar", "bfiat", 
                "bfico", "bficp", "bficr", "bfina", "bfine", "bfind", "bfioa", "bfioi", "bfioc")]
bfifa <- fa(bfi, 5, rotate = "oblimin", fm = "ml")
print(bfifa, sort = TRUE, cutoff = 0)

fscfbifasort <- fa.sort(bfifa)
write.csv(fscfbifasort$loadings, file = "sorted fsc bfi fa.csv", row.names = T)

write.csv(bfifa$loadings, file = "EFA FSC BFI.csv", row.names = T)

# Stop Signal Sensitivity Item Response by Conscientiousness Scores and Alpha Values
# data$fsc4 <- (data$fsc.31 + data$fsc.32 + data$fsc.33 + data$fsc.34 + data$fsc.36)/5
# function doBy and specific function is summaryBy

library(ltm)
cronbach.alpha(data.frame(data$fsc.31, data$fsc.32, data$fsc.33, data$fsc.34, data$fsc.36))

# ss1 fsc 31

table(data$fsc.31)

ss1.1 <- subset(data, data$fsc.31==1)
ss1.2 <- subset(data, data$fsc.31==2)
ss1.3 <- subset(data, data$fsc.31==3)
ss1.4 <- subset(data, data$fsc.31==4)
ss1.5 <- subset(data, data$fsc.31==5)
ss1.6 <- subset(data, data$fsc.31==6)
ss1.7 <- subset(data, data$fsc.31==7)

print(mean(ss1.1$bfic))
print(mean(ss1.2$bfic))
print(mean(ss1.3$bfic))
print(mean(ss1.4$bfic))
print(mean(ss1.5$bfic))
print(mean(ss1.6$bfic))
print(mean(ss1.7$bfic))

ss1 <- c((mean(ss1.1$bfic)), (mean(ss1.2$bfic)), (mean(ss1.3$bfic)), (mean(ss1.4$bfic)), (mean(ss1.5$bfic)), 
         (mean(ss1.6$bfic)), (mean(ss1.7$bfic)))
level <- 1:7
ss1.df <- cbind.data.frame(level, ss1)

# ggplot(data = ss1.df, aes(x = level, y = ss1))
# ggplot(data = ss1.df, mapping = aes(level, ss1))
# qplot(level, ss1, data = ss1.df)       

ggplot(data = ss1.df, aes(x = level)) + 
  geom_line(aes(y = ss1)) +
  ylim(1, 5)

# Amanda method with reduced calculations

library(Rmisc)
ss1.summary <- summarySE(data = data, "bfic", "fsc.31", conf.interval = .95)

ggss1 <- ggplot(data = ss1.summary, aes(x = fsc.31, y = bfic)) +  geom_line() + ylim(1, 5)

ggsave("ggss1.jpeg", units="in", width=8, height=8, dpi=300)
dev.off() # always need to run this after ggsave - Bible according to Amanda

# ss2 fsc 32

table(data$fsc.32)

ss2.1 <- subset(data, data$fsc.32==1)
ss2.2 <- subset(data, data$fsc.32==2)
ss2.3 <- subset(data, data$fsc.32==3)
ss2.4 <- subset(data, data$fsc.32==4)
ss2.5 <- subset(data, data$fsc.32==5)
ss2.6 <- subset(data, data$fsc.32==6)
ss2.7 <- subset(data, data$fsc.32==7)

ss2 <- c((mean(ss2.1$bfic)), (mean(ss2.2$bfic)), (mean(ss2.3$bfic)), (mean(ss2.4$bfic)), (mean(ss2.5$bfic)), 
         (mean(ss2.6$bfic)), (mean(ss2.7$bfic)))
level <- 1:7
ss2.df <- cbind.data.frame(level, ss2)

ggplot(data = ss2.df, aes(x = level)) + 
  geom_line(aes(y = ss2)) + 
  ylim(1, 5)

# ss3 fsc 33

table(data$fsc.33)

ss3.1 <- subset(data, data$fsc.33==1)
ss3.2 <- subset(data, data$fsc.33==2)
ss3.3 <- subset(data, data$fsc.33==3)
ss3.4 <- subset(data, data$fsc.33==4)
ss3.5 <- subset(data, data$fsc.33==5)
ss3.6 <- subset(data, data$fsc.33==6)
ss3.7 <- subset(data, data$fsc.33==7)

ss3 <- c((mean(ss3.1$bfic)), (mean(ss3.2$bfic)), (mean(ss3.3$bfic)), (mean(ss3.4$bfic)), (mean(ss3.5$bfic)), 
         (mean(ss3.6$bfic)), (mean(ss3.7$bfic)))
level <- 1:7
ss3.df <- cbind.data.frame(level, ss3)

ggplot(data = ss3.df, aes(x = level)) + 
  geom_line(aes(y = ss3)) + 
  ylim(1, 5)

# ss4 fsc 34

table(data$fsc.34)

ss4.1 <- subset(data, data$fsc.34==1)
ss4.2 <- subset(data, data$fsc.34==2)
ss4.3 <- subset(data, data$fsc.34==3)
ss4.4 <- subset(data, data$fsc.34==4)
ss4.5 <- subset(data, data$fsc.34==5)
ss4.6 <- subset(data, data$fsc.34==6)
ss4.7 <- subset(data, data$fsc.34==7)

ss4 <- c((mean(ss4.1$bfic)), (mean(ss4.2$bfic)), (mean(ss4.3$bfic)), (mean(ss4.4$bfic)), (mean(ss4.5$bfic)), 
         (mean(ss4.6$bfic)), (mean(ss4.7$bfic)))
level <- 1:7
ss4.df <- cbind.data.frame(level, ss4)

ggplot(data = ss4.df, aes(x = level)) + 
  geom_line(aes(y = ss4)) +
  ylim(1, 5)

# ss5 fsc 36

table(data$fsc.36)

ss5.1 <- subset(data, data$fsc.36==1)
ss5.2 <- subset(data, data$fsc.36==2)
ss5.3 <- subset(data, data$fsc.36==3)
ss5.4 <- subset(data, data$fsc.36==4)
ss5.5 <- subset(data, data$fsc.36==5)
ss5.6 <- subset(data, data$fsc.36==6)
ss5.7 <- subset(data, data$fsc.36==7)

ss5 <- c((mean(ss5.1$bfic)), (mean(ss5.2$bfic)), (mean(ss5.3$bfic)), (mean(ss5.4$bfic)), (mean(ss5.5$bfic)), 
         (mean(ss5.6$bfic)), (mean(ss5.7$bfic)))
level <- 1:7
ss5.df <- cbind.data.frame(level, ss5)

ggplot(data = ss5.df, aes(x = level)) + 
  geom_line(aes(y = ss5)) +
  ylim(1, 5)

library(Rmisc)
ss5.summary <- summarySE(data = data, "bfic", "fsc.36", conf.interval = .95)

ggss1 <- ggplot(data = ss5.summary, aes(x = fsc.36, y = bfic)) +  geom_line() + ylim(1, 5)

ggsave("ggss1.jpeg", units="in", width=8, height=8, dpi=300)
dev.off() # always need to run this after ggsave - Bible according to Amanda

# Stop Signal BFI Correlation
ss.bfi.cor <- cor(data.frame(data$fsc4, data$bfie, data$bfies, data$bfiea, data$bfiee, data$bfia, 
                             data$bfiac, data$bfiar, data$bfiat, data$bfic, data$bfico, data$bficp, 
                             data$bficr, data$bfin, data$bfina, data$bfind, data$bfine, data$bfio, 
                             data$bfioa, data$bfioi, data$bfioc), use="pairwise")
print(ss.bfi.cor)
write.csv(round(ss.bfi.cor, 2), file = "SS BFI COR.csv", row.names = T)

# FSC BFI Correlation
fsc.bfi.cor <- cor(data.frame(data$fsc2, data$fsc1, data$fsc3, data$bfico, data$bficp, data$bficr), use="pairwise")

print(fsc.bfi.cor)
write.csv(round(fsc.bfi.cor, 2), file = "FSC BFI COR.csv", row.names = T)
