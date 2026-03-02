# replicationCCES.r
# This program reproduces the CCES results from:
#     Rice, Douglas, Brian Schaffner, and David Barney (2020) "Political Ideology and 
#     Issue Importance." Political Research Quarterly
#
# dr 9.29.2020

library(foreign)
library(stringr)
library(ggplot2)
library(rstan)
library(readr)
library(haven)
library(plyr)
library(gdata)
library(parallel)
library(edstan)
library(ltm)

# 
set.seed(176)

# some stan preliminaries
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ----------------------------------------------
# ----------------------------------------------
# CCES Data
# ----------------------------------------------
# ----------------------------------------------

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Estimate basic IRT w/ just questions from CCES
# =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# set working directory
setwd("c:/users/drrice/dropbox/Schaffner Rice Barney/Data - 2016")

# read in data
myData <- read_dta("CCES16_Recoded_Relabeled.dta")

# select questions
myQuestions <- c("a_guncontrol_1","a_guncontrol_2","a_guncontrol_3","a_guncontrol_4",
				"a_abortion_1", "a_abortion_2", "a_abortion_3", "a_abortion_4", "a_abortion_5", "a_abortion_6",
				"a_taxes_1", "a_taxes_2",
				"a_immigration_1", "a_immigration_2", "a_immigration_3", "a_immigration_4", "a_immigration_5", "a_immigration_6", "a_immigration_7", "a_immigration_8", "a_immigration_9",
				"a_defensespending_1", "a_defensespending_2",
				"a_environment_1", "a_environment_2", "a_environment_3", "a_environment_4",
				"a_jobs_1", "a_jobs_2",
				"a_crime_1", "a_crime_2", "a_crime_3", "a_crime_4",
				"a_natsec_1", "a_natsec_2", "a_natsec_3", "a_natsec_4", "a_natsec_5", "a_natsec_6", "a_natsec_7", "a_natsec_8", 
				"a_healthcare_1","a_healthcare_2",
				"a_gaymarriage_1")
				
# create data for IRT 
myData$myIds <- myData$V101
J <- length(myData$myIds)
I <- length(myQuestions)
cols <- which(colnames(myData) %in% c(myQuestions))
y <- as.matrix(myData[,cols])

# delete missing data
N <- length(y)
jj <- rep(1:J, times=I)
ii <- rep(1:I, each=J)
miss <- which(is.na(y))
N <- N - length(miss)
jj <- jj[-miss]
ii <- ii[-miss]
y <- y[-miss]

# create indicator for 5 & 2
x <- rep(0,J)
x[5] <- 1
x[2] <- -1

# base IRT stan code
stanCode <- "
data {
  int<lower=1> J; // respondents
  int<lower=1> I; // questions
  int<lower=1> N; //no. of observations
  int<lower=1, upper=J> jj[N]; // respondent for observation n
  int<lower=1, upper=I> ii[N]; // question for observation n
  int<lower=0, upper=1> y[N]; // response of observation n
  real x[J]; // covariate for respondents
}
parameters {
  vector<lower=0>[I] alpha;
  vector[I] beta;
  vector[J] theta;	
  real gamma;
}
model {
  vector[N] eta;
  alpha ~ lognormal(0.5,1);
  beta ~ normal(0,3);
  for (j in 1:J)
	  theta[j] ~ normal(gamma * x[j], 1);
  for (n in 1:N)
	  eta[n] = beta[ii[n]] + (theta[jj[n]] * alpha[ii[n]]);
  y ~ bernoulli_logit(eta);
}"

# estimate
stanData <- list(J=J, I=I, N=N, jj=jj, ii=ii, y=y)

stanFit <- stan(model_code=stanCode, data=stanData, iter=5000, warmup=2500, chains=4, verbose=TRUE, cores=4, control = list(max_treedepth = 10))

# =-=-=-=-=-=-=-=-=-=-=-=-=-=
#  Extract Thetas
# =-=-=-=-=-=-=-=-=-=-=-=-=-=

sumFit <- summary(stanFit, pars=c("theta"))$summary
ips <- sumFit[,1]
myData$ips <- ips

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Estimate Weighted IRT w/ questions from CCES
# =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# drop everything except stanFit and ips
keep(stanFit, ips, sure=T)

# read in data
myData <- read_dta("CCES16_Recoded_Relabeled.dta")

# select questions
myQuestions <- c("a_guncontrol_1","a_guncontrol_2","a_guncontrol_3","a_guncontrol_4",
				"a_abortion_1", "a_abortion_2", "a_abortion_3", "a_abortion_4", "a_abortion_5", "a_abortion_6",
				"a_taxes_1", "a_taxes_2",
				"a_immigration_1", "a_immigration_2", "a_immigration_3", "a_immigration_4", "a_immigration_5", "a_immigration_6", "a_immigration_7", "a_immigration_8", "a_immigration_9",
				"a_defensespending_1", "a_defensespending_2",
				"a_environment_1", "a_environment_2", "a_environment_3", "a_environment_4",
				"a_jobs_1", "a_jobs_2",
				"a_crime_1", "a_crime_2", "a_crime_3", "a_crime_4",
				"a_natsec_1", "a_natsec_2", "a_natsec_3", "a_natsec_4", "a_natsec_5", "a_natsec_6", "a_natsec_7", "a_natsec_8", 
				"a_healthcare_1","a_healthcare_2",
				"a_gaymarriage_1")

# select scales
myScales <- c("i_guncontrol",
			  "i_abortion",
			  "i_taxes",
			  "i_immigration",
			  "i_defensespending",
			  "i_environment", 
			  "i_jobs",
			  "i_crime",
			  "i_natsec",
			  "i_healthcare",
			  "i_gaymarriage")

# create data for IRT 
myData$myIds <- myData$V101
J <- length(myData$myIds)
I <- length(myQuestions)
cols <- which(colnames(myData) %in% c(myQuestions))
y <- as.matrix(myData[,cols])

# create matrix of weights
cols <- which(colnames(myData) %in% c(myScales))
spots <- as.matrix(myData[,cols])
topics <- sub("i_", "", myScales)
spotsMat <- matrix(NA, J, I)
for (i in 1:I){
	for (k in 1:length(myScales)){
		if(length(grep(topics[k], myQuestions[i])) > 0){
			spotsMat[,i] <- spots[,k]
		}
	}
}

#spotsMat <- prop.table(spotsMat, 1)
spots <- spotsMat

# delete missing data
N <- length(y)
jj <- rep(1:J, times=I)
ii <- rep(1:I, each=J)
miss <- which(is.na(y))
N <- N - length(miss)
jj <- jj[-miss]
ii <- ii[-miss]
y <- y[-miss]
spots <- spots[-miss]
miss <- which(is.na(spots))
N <- N - length(miss)
jj <- jj[-miss]
ii <- ii[-miss]
y <- y[-miss]
spots <- spots[-miss]

# create indicator for 5 & 2
x <- rep(0,J)
x[5] <- 1
x[2] <- -1

# weighted IRT stan code
stanCode <- "
data {
  int<lower=1> J; // respondents
  int<lower=1> I; // questions, or number of columns
  int<lower=1> N; // number of observations
  int<lower=1, upper=J> jj[N]; // respondent for observation n
  int<lower=1, upper=I> ii[N]; // question for observation n
  int<lower=0, upper=1> y[N]; // response of observation n
  real<lower=1, upper=5> spots[N]; // scale of observation n
  real x[J]; // covariate for respondents
}
parameters {
  vector<lower=0>[I] alpha;
  vector[J] theta;	
  vector[I] beta;
  real gamma;
}
model {
  vector[N] eta;
  alpha ~ lognormal(0.5,1);
  beta ~ normal(0,3);
  for (j in 1:J)
	  theta[j] ~ normal(gamma * x[j],1);
  for (n in 1:N)
	  eta[n] = beta[ii[n]] + theta[jj[n]] * (alpha[ii[n]] + spots[n]);
  y ~ bernoulli_logit(eta);
}"

# create data
stanData <- list(J=J, I=I, N=N, jj=jj, ii=ii, y=y, spots=spots)

# estimate
weightedFit <- stan(model_code=stanCode, data=stanData, iter=5000, warmup=2500, chains=4, verbose=TRUE, cores=4)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=
#  Extract Estimates
# =-=-=-=-=-=-=-=-=-=-=-=-=-=

sumW <- summary(weightedFit, pars=c("theta"))$summary
wip <- sumW[,1]
myData$wip <- wip

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Estimate rating scale model of CCES
# =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# drop everything except stanFit and ips
keep(stanFit, ips, weightedFit, wip, sure=T)

# read in data
myData <- read_dta("CCES16_Recoded_Relabeled.dta")

# select questions
myQuestions <- c("a_guncontrol_1","a_guncontrol_2","a_guncontrol_3","a_guncontrol_4",
				"a_abortion_1", "a_abortion_2", "a_abortion_3", "a_abortion_4", "a_abortion_5", "a_abortion_6",
				"a_taxes_1", "a_taxes_2",
				"a_immigration_1", "a_immigration_2", "a_immigration_3", "a_immigration_4", "a_immigration_5", "a_immigration_6", "a_immigration_7", "a_immigration_8", "a_immigration_9",
				"a_defensespending_1", "a_defensespending_2",
				"a_environment_1", "a_environment_2", "a_environment_3", "a_environment_4",
				"a_jobs_1", "a_jobs_2",
				"a_crime_1", "a_crime_2", "a_crime_3", "a_crime_4",
				"a_natsec_1", "a_natsec_2", "a_natsec_3", "a_natsec_4", "a_natsec_5", "a_natsec_6", "a_natsec_7", "a_natsec_8", 
				"a_healthcare_1","a_healthcare_2",
				"a_gaymarriage_1")

# select scales
myScales <- c("i_guncontrol",
			  "i_abortion",
			  "i_taxes",
			  "i_immigration",
			  "i_defensespending",
			  "i_environment", 
			  "i_jobs",
			  "i_crime",
			  "i_natsec",
			  "i_healthcare",
			  "i_gaymarriage")

# create data for IRT 
myData$myIds <- myData$V101
J <- length(myData$myIds)
I <- length(myQuestions)
cols <- which(colnames(myData) %in% c(myQuestions))
y <- as.matrix(myData[,cols])

# create matrix of weights
cols <- which(colnames(myData) %in% c(myScales))
spots <- as.matrix(myData[,cols])
topics <- sub("i_", "", myScales)
spotsMat <- matrix(NA, J, I)
for (i in 1:I){
	for (k in 1:length(myScales)){
		if(length(grep(topics[k], myQuestions[i])) > 0){
			spotsMat[,i] <- spots[,k]
		}
	}
}
spots <- spotsMat

y[y==0] <- -1
y <- y * spots
y <- y + 5

# format for Furr GRM
myList <- irt_data(y)

# estimate
grmFit <- irt_stan(myList, "grsm_latent_reg.stan", chains=4, iter=1000, control = list(max_treedepth = 20))

# =-=-=-=-=-=-=-=-=-=-=-=-=-=
#  Extract Thetas
# =-=-=-=-=-=-=-=-=-=-=-=-=-=

sumFit <- summary(grmFit, pars=c("theta"))$summary
grips <- sumFit[,1]
myData$grips <- grips
 
 
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Effort Moderated IRT w/ questions from CCES
# =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# drop everything 
keep(stanFit, weightedFit, grmFit, ips, wips, grips, sure=T)

# read in data
myData <- read_dta("CCES16_Recoded_Relabeled.dta")

# select questions
myQuestions <- c("a_guncontrol_1","a_guncontrol_2","a_guncontrol_3","a_guncontrol_4",
				"a_abortion_1", "a_abortion_2", "a_abortion_3", "a_abortion_4", "a_abortion_5", "a_abortion_6",
				"a_taxes_1", "a_taxes_2",
				"a_immigration_1", "a_immigration_2", "a_immigration_3", "a_immigration_4", "a_immigration_5", "a_immigration_6", "a_immigration_7", "a_immigration_8", "a_immigration_9",
				"a_defensespending_1", "a_defensespending_2",
				"a_environment_1", "a_environment_2", "a_environment_3", "a_environment_4",
				"a_jobs_1", "a_jobs_2",
				"a_crime_1", "a_crime_2", "a_crime_3", "a_crime_4",
				"a_natsec_1", "a_natsec_2", "a_natsec_3", "a_natsec_4", "a_natsec_5", "a_natsec_6", "a_natsec_7", "a_natsec_8", 
				"a_healthcare_1","a_healthcare_2",
				"a_gaymarriage_1")

# select scales
myScales <- c("i_guncontrol",
			  "i_abortion",
			  "i_taxes",
			  "i_immigration",
			  "i_defensespending",
			  "i_environment", 
			  "i_jobs",
			  "i_crime",
			  "i_natsec",
			  "i_healthcare",
			  "i_gaymarriage")

# create data for IRT 
myData$myIds <- myData$V101
J <- length(myData$myIds)
I <- length(myQuestions)
cols <- which(colnames(myData) %in% c(myQuestions))
y <- as.matrix(myData[,cols])

# create matrix of weights
cols <- which(colnames(myData) %in% c(myScales))
spots <- as.matrix(myData[,cols])
topics <- sub("i_", "", myScales)
spotsMat <- matrix(NA, J, I)
for (i in 1:I){
	for (k in 1:length(myScales)){
		if(length(grep(topics[k], myQuestions[i])) > 0){
			spotsMat[,i] <- spots[,k]
		}
	}
}

# convert to dummy by item
#spots <- spotsMat
#for (s in 1:ncol(spots)){
#	spots[,s] <- ifelse(spots[,s] >= median(spots[,s], na.rm=T), 1, 0)
#}

# convert to dummy by respondent
spots <- spotsMat
for (s in 1:nrow(spots)){
	spots[s,] <- ifelse(spots[s,] >= median(spots[s,], na.rm=T), 1, 0)
}


# delete missing data
N <- length(y)
jj <- rep(1:J, times=I)
ii <- rep(1:I, each=J)
miss <- which(is.na(y))
N <- N - length(miss)
jj <- jj[-miss]
ii <- ii[-miss]
y <- y[-miss]
spots <- spots[-miss]
miss <- which(is.na(spots))
N <- N - length(miss)
jj <- jj[-miss]
ii <- ii[-miss]
y <- y[-miss]
spots <- spots[-miss]

# create indicator for 26 & 8
x <- rep(0,J)
x[26] <- 1
x[8] <- -1

# effort moderated IRT stan code
stanCode <- "
data {
  int<lower=1> J; // respondents
  int<lower=1> I; // questions, or number of columns
  int<lower=1> N; // number of observations
  int<lower=1, upper=J> jj[N]; // respondent for observation n
  int<lower=1, upper=I> ii[N]; // question for observation n
  int<lower=0, upper=1> y[N]; // response of observation n
  real<lower=0, upper=1> spots[N]; // dummy of observation n
  real x[J]; // covariate for respondents
}
parameters {
  vector<lower=0>[I] alpha;
  vector[J] theta;	
  vector[I] beta;
  real gamma;
}
model {
  vector[N] eta;
  alpha ~ lognormal(0.5,1);
  beta ~ normal(0,3);
  for (j in 1:J)
	  theta[j] ~ normal(gamma * x[j],1);
  for (n in 1:N)
	  eta[n] = (spots[n] * (beta[ii[n]] + (alpha[ii[n]] * theta[jj[n]]))) + ((1 - spots[n]) * 0.5);
  y ~ bernoulli_logit(eta);
}"

# create data
stanData <- list(J=J, I=I, N=N, jj=jj, ii=ii, y=y, spots=spots)

# estimate
emFit <- stan(model_code=stanCode, data=stanData, iter=5000, warmup=2500, chains=4, verbose=TRUE, cores=2)

sumFit <- summary(emFit, pars=c("theta"))$summary
remips <- sumFit[,1]
myData$remips <- remips

# =-=-=-=-=-=-=-=-=-=-=--=-=-
#     Create dataset
# =-=-=-=-=-=-=-=--=-=-=-=--=

# recode direction
ips <- ips * -1
wip <- wip * -1
grips <- grips * -1
remips <- remips * -1

# add to dataset
myData$ips <- ips
myData$grips <- grips
myData$remips <- remips
myData$wip <- wip

# =-=-=-=-=-=-=-=-=-=-=-=-=-=
# Scatterplot Comparisons
# =-=-=-=-=-=-=-=-=-=-=-=-=-=
library(GGally)

# comparisons within groups
myData$myPid <- NA
myData$myPid[which(myData$pid7 < 4)] <- "Democrat"
myData$myPid[which(myData$pid7 > 4 & myData$pid7 < 8)] <- "Republican"
myData$myPid[which(myData$pid7 == 4 | myData$pid7 > 7)] <- "Other"


my_dens <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_density(..., alpha = 0.6, bw=.25) 
}

my_points <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_point(..., alpha = 0.1) 
}

# scatterplot
p <- ggpairs(myData, mapping = aes(color = myPid), lower = list(continuous = my_points), diag = list(continuous = my_dens), 
			columns = which(colnames(myData) %in% c("ips", "wip", "grips", "remips")), columnLabels = c("Standard", "Scale", "Moderated", "Weighted")) +
			theme_bw()

# update colors
for (i in 1:p$nrow){
	for(j in 1:p$ncol){
	p[i,j] <- p[i,j] +  
		scale_colour_manual(breaks = c("Democrat", "Other", "Republican"), 
                      values=c("blue", "grey", "red"))  
#		scale_fill_manual(breaks = c("Democrat", "Other", "Republican"), 
#                      values=c("blue", "grey", "red"))
	}
}
 

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#  Compare Thetas to Obama approval
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

myData$bo <- NA
myData$bo[myData$CC16_320a == 4] <- 1
myData$bo[myData$CC16_320a == 3] <- 2
myData$bo[myData$CC16_320a == 5] <- 3
myData$bo[myData$CC16_320a == 2] <- 4
myData$bo[myData$CC16_320a == 1] <- 5

# estimate models; each of these questions is on a 5 point scale
i1 <- polr(as.factor(bo) ~ ips, data=myData)
w1 <- polr(as.factor(bo) ~ wip, data=myData)
g1 <- polr(as.factor(bo) ~ grips, data=myData)
r1 <- polr(as.factor(bo) ~ remips, data=myData)

# compute accuracy for Obama approval
ipreds1 <- table(predict(i1, type="class"), myData$bo[-which(is.na(myData$bo))])
iacc1 <- (ipreds1[1,1] + ipreds1[2,2] + ipreds1[3,3] + ipreds1[4,4] + ipreds1[5,5])/(sum(ipreds1))
wpreds1 <- table(predict(w1, type="class"), myData$bo[-which(is.na(myData$bo))])
wacc1 <- (wpreds1[1,1] + wpreds1[2,2] + wpreds1[3,3] + wpreds1[4,4] + wpreds1[5,5])/(sum(wpreds1))
gpreds1 <- table(predict(g1, type="class"), myData$bo[-which(is.na(myData$bo))])
gacc1 <- (gpreds1[1,1] + gpreds1[2,2] + gpreds1[3,3] + gpreds1[4,4] + gpreds1[5,5])/(sum(gpreds1))
rpreds1 <- table(predict(r1, type="class"), myData$bo[-which(is.na(myData$bo))])
racc1 <- (rpreds1[1,1] + rpreds1[2,2] + rpreds1[3,3] + rpreds1[4,4] + rpreds1[5,5])/(sum(rpreds1))

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#  Compare Thetas to Spending Scales
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

myData$deficit <- myData$CC16_415r
myData$deficit[myData$deficit > 100] <- NA

# four plots: deficit against each, labels = R^2 and correlation
pdf("cces-deficit-ips.pdf", height=5, width=5)
ggplot(myData, aes(x=ips, y=deficit, group=myPid)) +
	geom_point(aes(color=myPid), alpha=.2) + 
	scale_colour_manual(breaks = c("Democrat", "Other", "Republican"), 
                      values=c("blue", "grey", "red")) + 
	scale_fill_manual(breaks = c("Democrat", "Other", "Republican"), 
                   values=c("blue", "grey", "red")) +
	xlab(expression(paste("Correlation =  0.66, ", R^{2}," = 0.43"))) + 
	ylab("Deficit Preference") + 
	ggtitle("Standard") + 
	theme_bw() + 
	theme(legend.position="none")
dev.off()


# four plots: deficit against each, labels = R^2 and correlation
pdf("cces-deficit-wip.pdf", height=5, width=5)
ggplot(myData, aes(x=wip, y=deficit, group=myPid)) +
	geom_point(aes(color=myPid), alpha=.2) + 
	scale_colour_manual(breaks = c("Democrat", "Other", "Republican"), 
                      values=c("blue", "grey", "red")) + 
	scale_fill_manual(breaks = c("Democrat", "Other", "Republican"), 
                   values=c("blue", "grey", "red")) +
	xlab(expression(paste("Correlation =  0.66, ", R^{2}," = 0.43"))) + 
	ylab("Deficit Preference") + 
	ggtitle("Weighted") + 
	theme_bw() + 
	theme(legend.position="none")
dev.off()

	
# four plots: deficit against each, labels = R^2 and correlation
pdf("cces-deficit-grips.pdf", height=5, width=5)
ggplot(myData, aes(x=grips, y=deficit, group=myPid)) +
	geom_point(aes(color=myPid), alpha=.2) + 
	scale_colour_manual(breaks = c("Democrat", "Other", "Republican"), 
                      values=c("blue", "grey", "red")) + 
	scale_fill_manual(breaks = c("Democrat", "Other", "Republican"), 
                   values=c("blue", "grey", "red")) +
	xlab(expression(paste("Correlation =  0.63, ", R^{2}," = 0.40"))) + 
	ylab("Deficit Preference") + 
	ggtitle("Scale") + 
	theme_bw() +
	theme(legend.position="none")
dev.off()

	
# four plots: deficit against each, labels = R^2 and correlation
pdf("cces-deficit-remips.pdf", height=5, width=5)
ggplot(myData, aes(x=remips, y=deficit, group=myPid)) +
	geom_point(aes(color=myPid), alpha=.2) + 
	scale_colour_manual(breaks = c("Democrat", "Other", "Republican"), 
                      values=c("blue", "grey", "red")) + 
	scale_fill_manual(breaks = c("Democrat", "Other", "Republican"), 
                   values=c("blue", "grey", "red")) +
	ylab("Deficit Preference") + 
	xlab(expression(paste("Correlation =  0.65, ", R^{2}," = 0.42"))) + 
	ggtitle("Moderated") + 
	theme_bw() + 
	theme(legend.position="none")
dev.off()

# =-=-=-=-=-=-=-=-=-=-=-=-=-=
# Predict 2016 Primary Votes
# =-=-=-=-=-=-=-=-=-=-=-=-=-=

library(nnet)

# recode; excluded category are those who voted for another Democrat
demData <- subset(myData, CC16_328 < 4)
demData$vote <- 0
demData$vote[demData$CC16_328==1] <- 1
demData$vote[demData$CC16_328==2] <- 2

demIps <- multinom(as.factor(vote) ~ ips, data=demData)
demWip <- multinom(as.factor(vote) ~ wip, data=demData)
demGrips <- multinom(as.factor(vote) ~ grips, data=demData)
demRemips <- multinom(as.factor(vote) ~ remips, data=demData)
table(predict(demIps) == demData$vote) #2641, all Clinton
table(predict(demWip) == demData$vote) #2744, 234 correct Sanders, 131 incorrectly predicted Sanders
table(predict(demGrips) == demData$vote) #2641, all Clinton but two Sanders votes incorrectly predicted as 0s
table(predict(demRemips) == demData$vote) # 2641, all Clinton but one Sanders vote incorrectly predicted as 0

gopData <- subset(myData, CC16_328 > 3 & CC16_328 < 9)
gopData$vote <- 0
gopData$vote[gopData$CC16_328==4] <- 1
gopData$vote[gopData$CC16_328==5] <- 2
gopData$vote[gopData$CC16_328==6] <- 3
gopData$vote[gopData$CC16_328==7] <- 4

gopIps <- multinom(as.factor(vote) ~ ips, data=gopData)
gopWip <- multinom(as.factor(vote) ~ wip, data=gopData)
gopGrips <- multinom(as.factor(vote) ~ grips, data=gopData)
gopRemips <- multinom(as.factor(vote) ~ remips, data=gopData)
table(predict(gopIps) == gopData$vote) # 2134
table(predict(gopWip) == gopData$vote) # 2130, predicted more Trump, fewer Kasich
table(predict(gopGrips) == gopData$vote) # 2116, predicted more Trump, fewer Cruz
table(predict(gopRemips) == gopData$vote) # 2118, predicted fewer Trump, more Cruz


# create predicted probability data for Democrats
demIpsData <- data.frame(ips = seq(from = round(mean(demData$ips)-2*sd(demData$ips),1),
							 to = round(mean(demData$ips)+2*sd(demData$ips),1),
							 by = .01))
demWipData <- data.frame(wip = seq(from = round(mean(demData$wip)-2*sd(demData$wip),1),
							 to = round(mean(demData$wip)+2*sd(demData$wip),1),
							 by = .01))	
demGripsData <- data.frame(grips = seq(from = round(mean(demData$grips)-2*sd(demData$grips),1),
							 to = round(mean(demData$grips)+2*sd(demData$grips),1),
							 by = .01))
demRemipsData <- data.frame(remips = seq(from = round(mean(demData$remips)-2*sd(demData$remips),1),
							 to = round(mean(demData$remips)+2*sd(demData$remips),1),
							 by = .01))			 

# create predicted probability data for Republicans							 
gopIpsData <- data.frame(ips = seq(from = round(mean(gopData$ips)-2*sd(gopData$ips),1),
							 to = round(mean(gopData$ips)+2*sd(gopData$ips),1),
							 by = .01))
gopWipData <- data.frame(wip = seq(from = round(mean(gopData$wip)-2*sd(gopData$wip),1),
							 to = round(mean(gopData$wip)+2*sd(gopData$wip),1),
							 by = .01))
gopGripsData <- data.frame(grips = seq(from = round(mean(gopData$grips)-2*sd(gopData$grips),1),
							 to = round(mean(gopData$grips)+2*sd(gopData$grips),1),
							 by = .01))
gopRemipsData <- data.frame(remips = seq(from = round(mean(gopData$remips)-2*sd(gopData$remips),1),
							 to = round(mean(gopData$remips)+2*sd(gopData$remips),1),
							 by = .01))

library(reshape2)
library(directlabels)
							 
# predicted probabilities and plots for democrats
demIpsPP <- cbind(demIpsData, predict(demIps, newdata=demIpsData, type="probs"))
colnames(demIpsPP) <- c("ips", "Other", "Clinton", "Sanders") 
demIpsPP <- melt(demIpsPP, id.vars = c("ips"), value.name = "probability")
demWipPP <- cbind(demWipData, predict(demWip, newdata=demWipData, type="probs"))
colnames(demWipPP) <- c("wip", "Other", "Clinton", "Sanders") 
demWipPP <- melt(demWipPP, id.vars = c("wip"), value.name = "probability")
demGripsPP <- cbind(demGripsData, predict(demGrips, newdata=demGripsData, type="probs"))
colnames(demGripsPP) <- c("grips", "Other", "Clinton", "Sanders") 
demGripsPP <- melt(demGripsPP, id.vars = c("grips"), value.name = "probability")
demRemipsPP <- cbind(demRemipsData, predict(demRemips, newdata=demRemipsData, type="probs"))
colnames(demRemipsPP) <- c("remips", "Other", "Clinton", "Sanders") 
demRemipsPP <- melt(demRemipsPP, id.vars = c("remips"), value.name = "probability")


demIpsPlot <- ggplot(demIpsPP, aes(x=ips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.65) +
	xlim(-2.75,0.75) +
	xlab("Standard (Accuracy = 55.9)") + 
	ylab("Probability of Vote") +
	scale_color_grey(start=.8, end=.2) + 
	theme_bw() +
	theme(legend.position="none")
demWipPlot <- ggplot(demWipPP, aes(x=wip, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.65) +
	xlim(-.8,.25) +
	xlab("Weighted (Accuracy = 58.1)") + 
	ylab("Probability of Vote") +
	scale_color_grey(start=.8, end=.2) + 
	theme_bw() +
	theme(legend.position="none")
demGripsPlot <- ggplot(demGripsPP, aes(x=grips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.65) +
	xlim(-.5,0.25) +
	xlab("Scale (Accuracy = 55.9)") + 
	ylab("Probability of Vote") +
	scale_color_grey(start=.8, end=.2) + 
	theme_bw() +
	theme(legend.position="none")
demRemipsPlot <- ggplot(demRemipsPP, aes(x=remips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.65) +
	xlim(-2.75,0.75) +
	xlab("Moderated (Accuracy = 55.9)") + 
	ylab("Probability of Vote") +
	scale_color_grey(start=.8, end=.2) + 
	theme_bw() +
	theme(legend.position="none")
	
pdf("cces-demPrimary-ips.pdf", width=4, height=4)
direct.label(demIpsPlot, list("first.points", fontface="bold", cex=.8))
dev.off()
pdf("cces-demPrimary-wip.pdf", width=4, height=4)
direct.label(demWipPlot, list("first.points", fontface="bold", cex=.8))
dev.off()
pdf("cces-demPrimary-grips.pdf", width=4, height=4)
direct.label(demGripsPlot, list("first.points", fontface="bold", cex=.8))
dev.off()
pdf("cces-demPrimary-remips.pdf", width=4, height=4)
direct.label(demRemipsPlot, list("first.points", fontface="bold", cex=.8))
dev.off()

# predicted probabilities and plots for republicans
gopIpsPP <- cbind(gopIpsData, predict(gopIps, newdata=gopIpsData, type="probs"))
colnames(gopIpsPP) <- c("ips", "Other", "Trump", "Cruz", "Kasich", "Rubio") 
gopIpsPP <- melt(gopIpsPP, id.vars = c("ips"), value.name = "probability")
gopWipPP <- cbind(gopWipData, predict(gopWip, newdata=gopWipData, type="probs"))
colnames(gopWipPP) <- c("wip", "Other", "Trump", "Cruz", "Kasich", "Rubio")
gopWipPP <- melt(gopWipPP, id.vars = c("wip"), value.name = "probability")
gopGripsPP <- cbind(gopGripsData, predict(gopGrips, newdata=gopGripsData, type="probs"))
colnames(gopGripsPP) <- c("grips", "Other", "Trump", "Cruz", "Kasich", "Rubio")
gopGripsPP <- melt(gopGripsPP, id.vars = c("grips"), value.name = "probability")
gopRemipsPP <- cbind(gopRemipsData, predict(gopRemips, newdata=gopRemipsData, type="probs"))
colnames(gopRemipsPP) <- c("remips", "Other", "Trump", "Cruz", "Kasich", "Rubio")
gopRemipsPP <- melt(gopRemipsPP, id.vars = c("remips"), value.name = "probability")


gopIpsPlot <- ggplot(gopIpsPP, aes(x=ips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.65) +
	xlim(-1.25,2.5) +
	xlab("Standard (Accuracy = 49.8)") + 
	ylab("Probability of Vote") +
	scale_color_grey(start=.8, end=.2) + 
	theme_bw() +
	theme(legend.position="none")
gopWipPlot <- ggplot(gopWipPP, aes(x=wip, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.65) +
	xlim(-.4,.5) +
	xlab("Weighted (Accuracy = 49.7)") + 
	ylab("Probability of Vote") +
	scale_color_grey(start=.8, end=.2) + 
	theme_bw() +
	theme(legend.position="none")
gopGripsPlot <- ggplot(gopGripsPP, aes(x=grips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.65) +
	xlim(-.25,.35) +
	xlab("Scale (Accuracy = 49.3)") + 
	ylab("Probability of Vote") +
	scale_color_grey(start=.8, end=.2) + 
	theme_bw() +
	theme(legend.position="none")
gopRemipsPlot <- ggplot(gopRemipsPP, aes(x=remips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.6) +
	xlim(-1,2.25) +
	xlab("Moderated (Accuracy = 49.4)") + 
	ylab("Probability of Vote") +
	scale_color_grey(start=.8, end=.2) + 
	theme_bw() +
	theme(legend.position="none")
	
pdf("cces-gopPrimary-ips.pdf", width=4, height=4)
direct.label(gopIpsPlot, list("first.points", fontface="bold", cex=.8))
dev.off()
pdf("cces-gopPrimary-wip.pdf", width=4, height=4)
direct.label(gopWipPlot, list("first.points", fontface="bold", cex=.8))
dev.off()
pdf("cces-gopPrimary-grips.pdf", width=4, height=4)
direct.label(gopGripsPlot, list("first.points", fontface="bold", cex=.8))
dev.off()
pdf("cces-gopPrimary-remips.pdf", width=4, height=4)
direct.label(gopRemipsPlot, list("first.points", fontface="bold", cex=.6))
dev.off()			

# =-=-=-=-=-=-=-=-=-
# Bimodality
# =-=-=-=-=-=-=-=-=-

library(mousetrap)
bimodality_coefficient(myData$ips)
bimodality_coefficient(myData$wip)
bimodality_coefficient(myData$remips)
bimodality_coefficient(myData$grips)

