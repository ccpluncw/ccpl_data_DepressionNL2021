
library(tidyr)
library(dplyr)
library(nlme)
library(lme4)
#require(Hmisc)
library(multcomp)
library(psych)


# Read in desired dataset
sdThresh <- 5
sdThresh2 <- 5

tmp <- read.table("dat.ALLNumlineRT.txt", header=T, sep="\t")
rawSNs <- length(unique(tmp$sn))

tmpa <- tmp[!is.na(tmp$bdi),]
SNnoBDI <- length(unique(tmpa$sn))

snWithoutBetas <- unique(tmpa[is.na(tmpa$finalB),"sn"])
numline.bdi.tmp <- tmpa[!(tmpa$sn %in% snWithoutBetas),]
SNnoBeta <- length(unique(numline.bdi.tmp$sn))
# create column that calculates finalB SD values according to condition and event valence

numline.bdi.out <- numline.bdi.tmp %>% group_by(rtPressure, questionKey) %>% transform(finalBs = scale(finalB))

# Eliminate outliers beyond specified sdThres (Threshold: 5 SDs above or below the mean)

numline.bdi <- numline.bdi.out[abs(numline.bdi.out$finalBs) < sdThresh ,]

# Report eliminated observations (outliers outside sdThresh)

sdThresh.outliers <- numline.bdi.out[abs(numline.bdi.out$finalBs) > sdThresh ,]

SNwithinSDthresh <- length(unique(numline.bdi$sn))

sink("sdThresh.outliers.txt")
  cat("***** REMOVED BECAUSE BETA WAS MORE THAN 5 SDs FROM THE MEAN ******\n\n")
  print(sdThresh.outliers)
sink()


# Use filtered data (filters: sdThresh filter) and choose only desired columns

numline.bdi.tmp1 <- numline.bdi[,c(1,2,11,12,13,14,15,16,18,20:33)]

# Reshape data to get Participant Demographics

preFilteredParticipant.data <- tmpa[,c("sn", "rtPressure", "gender", "bdi",  "PANAS.PA", "PANAS.NA", "year", "age", "Hispanic_Latino", "American_Indian.Alaska_Native", "Asian" , "Black.African_American", "Native_Hawaiian.Pacific_Islander", "White", "Unknown")]
preFilteredParticipant.data <- preFilteredParticipant.data[!duplicated(preFilteredParticipant.data),]

# Reshape the data so that we can calculate finalB difference scores

numline.bdi.tmp2 <- reshape(numline.bdi.tmp1, idvar = c("sn", "model", "session", "bound", "estTask", "rtPressure", "estStimTime", "gender", "bdi",  "PANAS.PA", "PANAS.NA", "year", "age", "Hispanic_Latino", "American_Indian.Alaska_Native", "Asian" , "Black.African_American", "Native_Hawaiian.Pacific_Islander", "White", "Unknown"), timevar="questionKey", direction="wide")

#############################################################################
#create ethnicity variable

numline.bdi.tmp2$Hispanic_Latino <- ifelse(is.na(numline.bdi.tmp2$Hispanic_Latino), "unknown", as.character(numline.bdi.tmp2$Hispanic_Latino))
numline.bdi.tmp2$gender <- ifelse(is.na(numline.bdi.tmp2$gender), "unknown", as.character(numline.bdi.tmp2$gender))

numline.bdi.tmp2$ethnicity <- ifelse(numline.bdi.tmp2$American_Indian.Alaska_Native == "yes", "AIAN", NA)
numline.bdi.tmp2$ethnicity <- ifelse(numline.bdi.tmp2$Asian == "yes" & is.na(numline.bdi.tmp2$ethnicity), "Asian", ifelse(numline.bdi.tmp2$Asian == "yes" & !is.na(numline.bdi.tmp2$ethnicity), "Mult", numline.bdi.tmp2$ethnicity))
numline.bdi.tmp2$ethnicity <- ifelse(numline.bdi.tmp2$Black.African_American == "yes" & is.na(numline.bdi.tmp2$ethnicity), "Black", ifelse(numline.bdi.tmp2$Black.African_American == "yes" & !is.na(numline.bdi.tmp2$ethnicity), "Mult", numline.bdi.tmp2$ethnicity))
numline.bdi.tmp2$ethnicity <- ifelse(numline.bdi.tmp2$Native_Hawaiian.Pacific_Islander == "yes" & is.na(numline.bdi.tmp2$ethnicity), "NAPI", ifelse(numline.bdi.tmp2$Native_Hawaiian.Pacific_Islander == "yes" & !is.na(numline.bdi.tmp2$ethnicity), "Mult", numline.bdi.tmp2$ethnicity))
numline.bdi.tmp2$ethnicity <- ifelse(numline.bdi.tmp2$White == "yes" & is.na(numline.bdi.tmp2$ethnicity), "White", ifelse(numline.bdi.tmp2$White == "yes" & !is.na(numline.bdi.tmp2$ethnicity), "Mult", numline.bdi.tmp2$ethnicity))
numline.bdi.tmp2$ethnicity <- ifelse(numline.bdi.tmp2$Unknown == "yes" & is.na(numline.bdi.tmp2$ethnicity), "Unknown", numline.bdi.tmp2$ethnicity)


preFilteredParticipant.data$Hispanic_Latino <- ifelse(is.na(preFilteredParticipant.data$Hispanic_Latino), "unknown", as.character(preFilteredParticipant.data$Hispanic_Latino))
preFilteredParticipant.data$gender <- ifelse(is.na(preFilteredParticipant.data$gender), "unknown", as.character(preFilteredParticipant.data$gender))

preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$American_Indian.Alaska_Native == "yes", "AIAN", NA)
preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$Asian == "yes" & is.na(preFilteredParticipant.data$ethnicity), "Asian", ifelse(preFilteredParticipant.data$Asian == "yes" & !is.na(preFilteredParticipant.data$ethnicity), "Mult", preFilteredParticipant.data$ethnicity))
preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$Black.African_American == "yes" & is.na(preFilteredParticipant.data$ethnicity), "Black", ifelse(preFilteredParticipant.data$Black.African_American == "yes" & !is.na(preFilteredParticipant.data$ethnicity), "Mult", preFilteredParticipant.data$ethnicity))
preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$Native_Hawaiian.Pacific_Islander == "yes" & is.na(preFilteredParticipant.data$ethnicity), "NAPI", ifelse(preFilteredParticipant.data$Native_Hawaiian.Pacific_Islander == "yes" & !is.na(preFilteredParticipant.data$ethnicity), "Mult", preFilteredParticipant.data$ethnicity))
preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$White == "yes" & is.na(preFilteredParticipant.data$ethnicity), "White", ifelse(preFilteredParticipant.data$White == "yes" & !is.na(preFilteredParticipant.data$ethnicity), "Mult", preFilteredParticipant.data$ethnicity))
preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$Unknown == "yes" & is.na(preFilteredParticipant.data$ethnicity), "Unknown", preFilteredParticipant.data$ethnicity)

#############################################################################


# Calculate a new variable called "nnFinalB": neutralFinalB - negativeFinalB

numline.bdi.tmp2$nnFinalB <- numline.bdi.tmp2$finalB.negative - numline.bdi.tmp2$finalB.neutral

# Calculate a new variable called "npFinalB": neutralFinalB - positiveFinalB

numline.bdi.tmp2$npFinalB <- numline.bdi.tmp2$finalB.positive - numline.bdi.tmp2$finalB.neutral


############################################################################

tmp2 <- numline.bdi.tmp2

numline.bdi.tmp <- tmp2[!is.na(tmp2$nnFinalB) & !is.na(tmp2$npFinalB) ,]

# create columns that calculate nnfinalB & npfinalB SD values according to condition

numline.bdi.tmp2 <- numline.bdi.tmp %>% group_by(rtPressure) %>% transform ( nnfinalBs = scale(nnFinalB), npfinalBs = scale(npFinalB))

# Eliminate outliers beyond specified sdThres (Threshold: 3 SDs above or below the mean)

numline.bdi.tmp.keep <- numline.bdi.tmp2[abs(numline.bdi.tmp2$nnfinalBs) < sdThresh2 & abs(numline.bdi.tmp2$npfinalBs) < sdThresh2 , ]

# Report eliminated observations (outliers outside sdThresh)

sdThresh.outliers <- numline.bdi.tmp2[abs(numline.bdi.tmp2$nnfinalBs) > sdThresh2 | abs(numline.bdi.tmp2$npfinalBs) > sdThresh2, ]

snFinal <- length(unique(numline.bdi.tmp.keep$sn))
sink("sdThresh.outliers.txt", append=TRUE)
  cat("***** REMOVED BECAUSE BETA DIFF WAS MORE THAN 5 SDs FROM THE MEAN ******\n\n")
  print(sdThresh.outliers)
sink()


sink("Demographics.txt")
cat("\nRaw subjects: ", rawSNs)
cat("\nN removed because of missing BDI scores: ", rawSNs - SNnoBDI, " Resulting in a new total of: ", SNnoBDI)
cat("\nN removed because of missing Betas: ", SNnoBDI - SNnoBeta, " Resulting in a new total of: ", SNnoBeta)
cat("\nN removed because their were Betas more than ", sdThresh, " SDs from the mean: ", SNnoBeta - SNwithinSDthresh, " Resulting in a new total of: ", SNwithinSDthresh)
cat("\nN removed because their were Beta Diffs more than ", sdThresh2, " SDs from the mean: ", SNwithinSDthresh - snFinal, " Resulting in a new total of: ", snFinal)
cat("\nThus, the final participant total N: ", snFinal)

#Participant Demographics
#Gender
cat("\n\n**** Participant Demographics before filtering *****\n\n")
print(with(preFilteredParticipant.data, table(gender)))
#Mean Age
print(with(preFilteredParticipant.data, mean(age, na.rm=T)))
#SD Age
print(with(preFilteredParticipant.data, sd(age, na.rm=T)))
#Ethnicity
print(with(preFilteredParticipant.data, table(Hispanic_Latino)))
print(with(preFilteredParticipant.data, table(ethnicity)))

  #Participant Demographics
  cat("\n\n**** Participant Demographics after filtering *****\n\n")
  #Gender
  print(with(numline.bdi.tmp.keep, table(gender)))
  #Mean Age
  print(with(numline.bdi.tmp.keep, mean(age)))
  #SD Age
  print(with(numline.bdi.tmp.keep, sd(age)))
  #Ethnicity
  print(with(numline.bdi.tmp.keep, table(Hispanic_Latino)))
  print(with(numline.bdi.tmp.keep, table(ethnicity)))
sink(NULL)
#############################################################################

numline.bdi.tmp.keep$rtPressureF <- as.factor(numline.bdi.tmp.keep$rtPressure)

#############################################################################
# Reshape the data so that finalB difference scores are organized in a single column

numline.bdi.tmp3 <- reshape(numline.bdi.tmp.keep, varying = c("nnFinalB", "npFinalB"), idvar='sn', direction="long", v.names=c("change.score"), times=c("nn", "np"), timevar='diffB')

#############################################################################

# Recode diffB variable to create new variable sctype (score type)

numline.bdi.tmp3$sctype <- ifelse(numline.bdi.tmp3$diffB == "nn", -1, 1)

# Data for NRT (Control) Condition only

numline.NRT <- numline.bdi.tmp3[numline.bdi.tmp3$rtPressure==FALSE,]

# Data for RT Condition only

numline.RT <- numline.bdi.tmp3[numline.bdi.tmp3$rtPressure==TRUE,]
#############################################################################

# Test Beck Model for BDI score
# NRT Condition

Beck.BDI.NRT <- numline.NRT %>% group_by(sctype,bdi) %>% summarise (change.score.mean=mean(change.score, na.rm=TRUE))
Beck.BDI.NRT$sctype.bdi <- Beck.BDI.NRT$sctype*Beck.BDI.NRT$bdi
Beck.BDI.NRT$sctype.neg <- Beck.BDI.NRT$sctype*-1

Beck.BDI.NRT.lm <- with(Beck.BDI.NRT,  lm(change.score.mean ~ 0 + sctype.bdi))

# RT Condition

Beck.BDI.RT <- numline.RT %>% group_by(sctype,bdi) %>% summarise (change.score.mean=mean(change.score, na.rm=TRUE))
Beck.BDI.RT$sctype.bdi <- Beck.BDI.RT$sctype*Beck.BDI.RT$bdi
Beck.BDI.RT$sctype.neg <- Beck.BDI.RT$sctype*-1

Beck.BDI.RT.lm <- with(Beck.BDI.RT,  lm(change.score.mean ~ 0 + sctype.bdi))

# Test DR Model for BDI score
# NRT Condition

DR.BDI.NRT.lm <- with(Beck.BDI.NRT,  lm(change.score.mean ~ 0 + sctype.neg + sctype.bdi))

# RT Condition

DR.BDI.RT.lm <- with(Beck.BDI.RT,  lm(change.score.mean ~ 0 + sctype.neg + sctype.bdi))

#############################################################################

sink("Beck and DR Models.out.txt")

print("BECK")
print("NO RT Pressure Condition")
print(summary(Beck.BDI.NRT.lm))
BIC(Beck.BDI.NRT.lm)
print("RT Pressure Condition")
print(summary(Beck.BDI.RT.lm))
BIC(Beck.BDI.RT.lm)

print("DR")
print("NO RT Pressure Condition")
print(summary(DR.BDI.NRT.lm))
BIC(DR.BDI.NRT.lm)
print("RT Pressure Condition")
print(summary(DR.BDI.RT.lm))
BIC(DR.BDI.RT.lm)

sink()

#############################################################################
# t-test comparing RTs for the reaction time pressure condition vs. control condition

numline.RTPvNP <- numline.bdi.tmp1 %>% group_by(sn,rtPressure) %>% summarise (rt.mean=mean(mRT, na.rm=TRUE), beta.mean=mean(finalB, na.rm=TRUE))

RTt.test <- t.test(numline.RTPvNP$rt.mean ~ numline.RTPvNP$rtPressure)
RT.cohensD <- cohen.d(numline.RTPvNP$rt.mean, numline.RTPvNP$rtPressure,alpha=.05,std=TRUE)

Betat.test <- t.test(numline.RTPvNP$beta.mean ~ numline.RTPvNP$rtPressure)
Beta.cohensD <- cohen.d(numline.RTPvNP$beta.mean, numline.RTPvNP$rtPressure,alpha=.05,std=TRUE)

numline.sum <- data.frame(numline.bdi.tmp1 %>% group_by(rtPressure) %>% summarise(rt.mean=mean(mRT, na.rm=TRUE), beta.mean=mean(finalB, na.rm=TRUE), rt.sd=sd(mRT, na.rm=TRUE), beta.sd=sd(finalB, na.rm=TRUE)))

sink("t-test.out.txt")
  cat("\n**** ttest comparing the RT across rt pressure conditions.\n\n")
  print(RTt.test)
  cat("\nCohen's D: ", RT.cohensD$cohen.d[2])
  cat("\n\n**** ttest comparing the Betas across rt pressure conditions.\n\n")
  print(Betat.test)
  cat("\nCohen's D: ", Beta.cohensD$cohen.d[2])
  cat("\n A table of the means and betas for the rt conditions \n")
  print(numline.sum)
sink()

#############################################################################
# BDI score Histogram numline.bdi.tmp3
#hist.data <- tapply(numline.bdi.tmp3$bdi, numline.bdi.tmp3$sn, mean)

#par(mar=c(5.1, 6.1, 4.1, 2.1))
hist(numline.bdi.tmp2$bdi, las=1, main=NULL, xlab="BDI Score", ylab=NULL, xlim=c(0,63), ylim=c(0,40), col=NULL)
mtext("N", side=2, las=1, adj=1, line=3)
dev.copy(pdf, 'BDI score Final Subjects.pdf')
dev.off()

hist(preFilteredParticipant.data$bdi, las=1, main=NULL, xlab="BDI Score", ylab=NULL, xlim=c(0,63), ylim=c(0,40), col=NULL)
mtext("N", side=2, las=1, adj=1, line=3)
dev.copy(pdf, 'BDI score ALL Subjects.pdf')
dev.off()

pdf("Change Score by BDI Score for the 500 ms Condition.pdf")
plot(Beck.BDI.NRT$sctype.bdi, Beck.BDI.NRT$change.score.mean, xlab="BDI Score", ylab="change score", main="Change Score by BDI Score * sctype for the NRT Condition")
abline(Beck.BDI.NRT.lm)
plot(Beck.BDI.RT$sctype.bdi, Beck.BDI.RT$change.score.mean, xlab="BDI Score", ylab="change score", main="Change Score by BDI Score * sctype for the RT Condition")
abline(Beck.BDI.RT.lm)
dev.off()

pdf("Predicte Functions.pdf")
plot(NULL, NULL, xlab="BDI Score", ylab="change score", main="Beck Model", ylim = c(-0.1, 0.1), xlim = c(-40,40), bty="n")
abline(0,0, col="grey", lty=3)
abline(0,0.001, col="black")

xBDIn <- c(-40, 0)
xBDIp <- c(0, 40)
yBDIn <- 0.04 + .001*xBDIn
yBDIp <- -0.04 + .001*xBDIp
plot(NULL, NULL, xlab="BDI Score", ylab="change score", main="Depressive Realism Model", ylim = c(-0.1, 0.1), xlim = c(-40,40), bty="n")
lines(yBDIn~xBDIn)
lines(yBDIp~xBDIp)
abline(0,0, col="grey", lty=3)
dev.off()
