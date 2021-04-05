
library(tidyr)
library(dplyr)
library(nlme)
library(lme4)
#require(Hmisc)
library(multcomp)
library(sjstats)

# Read in desired dataset
sdThresh <- 5
sdThresh2 <- 5

tmp <- read.table("dat.ALLNumlineLPT.txt", header=T, sep="\t")
rawSNs <- length(unique(tmp$sn))

tmpa <- tmp[!is.na(tmp$bdi),]
SNnoBDI <- length(unique(tmpa$sn))

snWithoutBetas <- unique(tmpa[is.na(tmpa$finalB),"sn"])
numline.bdi.tmp <- tmpa[!(tmpa$sn %in% snWithoutBetas),]
SNnoBeta <- length(unique(numline.bdi.tmp$sn))

# create column that calculates finalB SD values according to condition and event valence

numline.bdi.out <- numline.bdi.tmp %>% group_by(estStimTime, questionKey) %>% transform(finalBs = scale(finalB))

# Eliminate outliers beyond specified sdThres (Threshold: 5 SDs above or below the mean)

numline.bdi <- numline.bdi.out[abs(numline.bdi.out$finalBs) <= sdThresh ,]

# Report eliminated observations (outliers outside sdThresh)

sdThresh.outliers <- numline.bdi.out[abs(numline.bdi.out$finalBs) > sdThresh ,]

SNwithinSDthresh <- length(unique(numline.bdi$sn))

sink("sdThresh.outliers.txt")
  cat("***** REMOVED BECAUSE BETA WAS MORE THAN 5 SDs FROM THE MEAN ******\n\n")
  print(sdThresh.outliers)
sink(NULL)

#############################################################################

# Use filtered data (filters: sdThresh filter) and choose only desired columns

numline.bdi.tmp1 <- numline.bdi[,c(1,2,11,12,13,14,15,16, 18, 20:33)]

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
numline.bdi.tmp2$ethnicity <- ifelse(is.na(numline.bdi.tmp2$ethnicity), "Unknown", numline.bdi.tmp2$ethnicity)

preFilteredParticipant.data$Hispanic_Latino <- ifelse(is.na(preFilteredParticipant.data$Hispanic_Latino), "unknown", as.character(preFilteredParticipant.data$Hispanic_Latino))
preFilteredParticipant.data$gender <- ifelse(is.na(preFilteredParticipant.data$gender), "unknown", as.character(preFilteredParticipant.data$gender))

preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$American_Indian.Alaska_Native == "yes", "AIAN", NA)
preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$Asian == "yes" & is.na(preFilteredParticipant.data$ethnicity), "Asian", ifelse(preFilteredParticipant.data$Asian == "yes" & !is.na(preFilteredParticipant.data$ethnicity), "Mult", preFilteredParticipant.data$ethnicity))
preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$Black.African_American == "yes" & is.na(preFilteredParticipant.data$ethnicity), "Black", ifelse(preFilteredParticipant.data$Black.African_American == "yes" & !is.na(preFilteredParticipant.data$ethnicity), "Mult", preFilteredParticipant.data$ethnicity))
preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$Native_Hawaiian.Pacific_Islander == "yes" & is.na(preFilteredParticipant.data$ethnicity), "NAPI", ifelse(preFilteredParticipant.data$Native_Hawaiian.Pacific_Islander == "yes" & !is.na(preFilteredParticipant.data$ethnicity), "Mult", preFilteredParticipant.data$ethnicity))
preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$White == "yes" & is.na(preFilteredParticipant.data$ethnicity), "White", ifelse(preFilteredParticipant.data$White == "yes" & !is.na(preFilteredParticipant.data$ethnicity), "Mult", preFilteredParticipant.data$ethnicity))
preFilteredParticipant.data$ethnicity <- ifelse(preFilteredParticipant.data$Unknown == "yes" & is.na(preFilteredParticipant.data$ethnicity), "Unknown", preFilteredParticipant.data$ethnicity)
preFilteredParticipant.data$ethnicity <- ifelse(is.na(preFilteredParticipant.data$ethnicity), "Unknown", preFilteredParticipant.data$ethnicity)

#############################################################################


# Calculate a new variable called "nnFinalB": neutralFinalB - negativeFinalB

numline.bdi.tmp2$nnFinalB <- numline.bdi.tmp2$finalB.negative - numline.bdi.tmp2$finalB.neutral

# Calculate a new variable called "npFinalB": neutralFinalB - positiveFinalB

numline.bdi.tmp2$npFinalB <- numline.bdi.tmp2$finalB.positive - numline.bdi.tmp2$finalB.neutral

#############################################################################
# Data for 0 (Control) Condition only

numline.0.t <- numline.bdi.tmp1[numline.bdi.tmp1$estStimTime==0,]

# t-test comparing RTs for the reaction time pressure condition vs. control condition

LPTt.test <- t.test(numline.0.t$mRT, mu=1000)
LPTt.cohensD <- (mean(numline.0.t$mRT) - 1000)/sd(numline.0.t$mRT)


sink("LPT effects.txt")
  cat("\n**** ttest comparing the RT of the control condition to the 1000 ms time pressure condition.\n\n")
  print(LPTt.test)
  print(LPTt.cohensD)
sink(NULL)
#############################################################################

tmp2 <- numline.bdi.tmp2

numline.bdi.tmp <- tmp2[!is.na(tmp2$nnFinalB) & !is.na(tmp2$npFinalB) ,]

# create columns that calculate nnfinalB & npfinalB SD values according to condition

numline.bdi.tmp2 <- numline.bdi.tmp %>% group_by(estStimTime) %>% transform ( nnfinalBs = scale(nnFinalB), npfinalBs = scale(npFinalB))

# Eliminate outliers beyond specified sdThres (Threshold: 5 SDs above or below the mean)

numline.bdi.tmp.keep <- numline.bdi.tmp2[abs(numline.bdi.tmp2$nnfinalBs) < sdThresh2 & abs(numline.bdi.tmp2$npfinalBs) < sdThresh2 , ]

# Report eliminated observations (outliers outside sdThresh)

sdThresh.outliers <- numline.bdi.tmp2[abs(numline.bdi.tmp2$nnfinalBs) > sdThresh2 | abs(numline.bdi.tmp2$npfinalBs) > sdThresh2, ]
snFinal <- length(unique(numline.bdi.tmp.keep$sn))

sink("sdThresh.outliers.txt", append=TRUE)
  cat("***** REMOVED BECAUSE BETA DIFF WAS MORE THAN 5 SDs FROM THE MEAN ******\n\n")
  print(sdThresh.outliers)
sink(NULL)

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
  #Gender
  cat("\n\n**** Participant Demographics after filtering *****\n\n")
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

numline.bdi.tmp.keep$estStimTimeF <- as.factor(numline.bdi.tmp.keep$estStimTime)

#############################################################################
# Reshape the data so that finalB difference scores are organized in a single column

numline.bdi.tmp3 <- reshape(numline.bdi.tmp.keep, varying = c("nnFinalB", "npFinalB"), idvar='sn', direction="long", v.names=c("change.score"), times=c("nn", "np"), timevar='diffB')

#############################################################################

# Recode diffB variable to create new variable sctype (score type)

numline.bdi.tmp3$sctype <- ifelse(numline.bdi.tmp3$diffB == "nn", -1, 1)

# Isolate the 0 (control) LPT condition

numline.0 <- numline.bdi.tmp3[numline.bdi.tmp3$estStimTime==0,]

# Isolate the 500 LPT condition

numline.500 <- numline.bdi.tmp3[numline.bdi.tmp3$estStimTime==500,]

# Isolate the 1000 LPT condition

numline.1000 <- numline.bdi.tmp3[numline.bdi.tmp3$estStimTime==1000,]

# Isolate the 500 & 1000 LPT conditions

numline.limit <- numline.bdi.tmp3[numline.bdi.tmp.keep$estStimTime!=0,]

#############################################################################

# Test Beck Model for BDI score
# 0 Condition

Beck.BDI.0 <- numline.0 %>% group_by(sctype,bdi) %>% summarise (change.score.mean=mean(change.score, na.rm=TRUE))

Beck.BDI.0$sctype.bdi <- Beck.BDI.0$sctype*Beck.BDI.0$bdi
Beck.BDI.0$sctype.neg <- Beck.BDI.0$sctype*-1

Beck.BDI.0.lm <- with(Beck.BDI.0,  lm(change.score.mean ~ 0 + sctype.bdi))

# 500 Condition

Beck.BDI.500 <- numline.500 %>% group_by(sctype,bdi) %>% summarise (change.score.mean=mean(change.score, na.rm=TRUE))
Beck.BDI.500$sctype.bdi <- Beck.BDI.500$sctype*Beck.BDI.500$bdi
Beck.BDI.500$sctype.neg <- Beck.BDI.500$sctype*-1

Beck.BDI.500.lm <- with(Beck.BDI.500,  lm(change.score.mean ~ 0 + sctype.bdi))

# 1000 Condition

Beck.BDI.1000 <- numline.1000 %>% group_by(sctype,bdi) %>% summarise (change.score.mean=mean(change.score, na.rm=TRUE))
Beck.BDI.1000$sctype.bdi <- Beck.BDI.1000$sctype*Beck.BDI.1000$bdi
Beck.BDI.1000$sctype.neg <- Beck.BDI.1000$sctype*-1

Beck.BDI.1000.lm <- with(Beck.BDI.1000,  lm(change.score.mean ~ 0 + sctype.bdi))

# limit Conditions

Beck.BDI.limit <- numline.limit %>% group_by(sctype,bdi) %>% summarise (change.score.mean=mean(change.score, na.rm=TRUE))
Beck.BDI.limit$sctype.bdi <- Beck.BDI.limit$sctype*Beck.BDI.limit$bdi
Beck.BDI.limit$sctype.neg <- Beck.BDI.limit$sctype*-1

Beck.BDI.limit.lm <- with(Beck.BDI.limit,  lm(change.score.mean ~ 0 + sctype.bdi))

# Test DR Model for BDI score
# 0 Condition

DR.BDI.0.lm <- with(Beck.BDI.0,  lm(change.score.mean ~ 0 + sctype.neg + sctype.bdi))

# 500 Condition

DR.BDI.500.lm <- with(Beck.BDI.500,  lm(change.score.mean ~ 0 + sctype.neg + sctype.bdi))

# 1000 Condition

DR.BDI.1000.lm <- with(Beck.BDI.1000,  lm(change.score.mean ~ 0 + sctype.neg + sctype.bdi))

# limit Conditions

DR.BDI.limit.lm <- with(Beck.BDI.limit,  lm(change.score.mean ~ 0 + sctype.neg + sctype.bdi))

#############################################################################

sink("Beck and DR Models.out.txt")

  print("Beck control")
  print(summary(Beck.BDI.0.lm))
  BIC(Beck.BDI.0.lm)
  print("Beck 500ms")
  print(summary(Beck.BDI.500.lm))
  BIC(Beck.BDI.500.lm)
  print("Beck 1000ms")
  print(summary(Beck.BDI.1000.lm))
  BIC(Beck.BDI.1000.lm)

  print("DR control")
  print(summary(DR.BDI.0.lm))
  BIC(DR.BDI.0.lm)
  print("DR 500")
  print(summary(DR.BDI.500.lm))
  BIC(DR.BDI.500.lm)
  print("DR 1000")
  print(summary(DR.BDI.1000.lm))
  BIC(DR.BDI.1000.lm)

sink(NULL)

#############################################################################
#ANOVA (Mixed-effects Model)
numline.bdi.tmp1$estStimTimeF <- as.factor(numline.bdi.tmp1$estStimTime)

beta.means.LPT <- numline.bdi.tmp1 %>% group_by(sn, estStimTimeF) %>% summarise (beta.mean=mean(finalB, na.rm=TRUE))

beta.means.LPT.aov <- aov(beta.mean ~ estStimTimeF,data = beta.means.LPT)

mRT.means.LPT <- numline.bdi.tmp1 %>% group_by(estStimTime) %>% summarise (mRT.mean=mean(mRT, na.rm=TRUE), mRT.sd=sd(mRT, na.rm=TRUE), beta.mean=mean(finalB, na.rm=TRUE), beta.sd=sd(finalB, na.rm=TRUE))


sink("LPT effects.txt", append = T)
  cat("\n**** ANOVA comparing the Betas across rt pressure conditions (including control condition).\n\n")
  print(anova(beta.means.LPT.aov))
  cat("\n ETA SQ for the anova \n")
  print(eta_sq(beta.means.LPT.aov))
  cat("\n A table of the means and betas for the rt conditions \n")
  print(mRT.means.LPT)

sink(NULL)

#############################################################################
pdf("Change Score by BDI Score for the 500 ms Condition.pdf")
plot(Beck.BDI.0$sctype.bdi, Beck.BDI.0$change.score.mean, xlab="sctype*BDI Score", ylab="Beta Diff", main="No Time Pressure", ylim=c(-0.15, 0.15), xlim=c(-40, 40), bty="n")
abline(0,0, col="grey", lty=3)
abline(Beck.BDI.0.lm)
plot(Beck.BDI.500$sctype.bdi, Beck.BDI.500$change.score.mean, xlab="sctype*BDI Score", ylab="Beta Diff", main="500 ms", ylim=c(-0.15, 0.15), xlim=c(-40, 40), bty="n")
abline(0,0, col="grey", lty=3)
abline(Beck.BDI.500.lm)
plot(Beck.BDI.1000$sctype.bdi, Beck.BDI.1000$change.score.mean, xlab="sctype*BDI Score", ylab="Beta Diff", main="1000 ms", ylim=c(-0.15, 0.15), xlim=c(-40, 40), bty="n")
abline(0,0, col="grey", lty=3)
abline(Beck.BDI.1000.lm)

dev.off()

#par(mar=c(5.1, 6.1, 4.1, 2.1))
hist(numline.bdi.tmp2$bdi, las=1, main=NULL, xlab="BDI Score", ylab=NULL, xlim=c(0,63), ylim=c(0,80), col=NULL)
mtext("N", side=2, las=1, adj=1, line=3)
dev.copy(pdf, 'BDI score Final Subjects.pdf')
dev.off()

hist(preFilteredParticipant.data$bdi, las=1, main=NULL, xlab="BDI Score", ylab=NULL, xlim=c(0,63), ylim=c(0,80), col=NULL)
mtext("N", side=2, las=1, adj=1, line=3)
dev.copy(pdf, 'BDI score ALL Subjects.pdf')
dev.off()
