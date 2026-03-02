# replicationLucid.r
# This program reproduces the Lucid results from:
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
# Lucid Data
# ----------------------------------------------
# ----------------------------------------------

# read in data
myData <- read_dta("Schaffner_Barney_Lucid_Module.dta")

# recode a few of these based on factor loadings.
myData$immig3 <- -1 * (myData$immig3 - 1)
myData$abort3 <- -1 * (myData$abort3 - 1)
myData$guns1 <- -1 * (myData$guns1 - 1)
myData$guns2 <- -1 * (myData$guns2 - 1)

# flip everything so conservative is positive
myData$immig1 <- -1 * (myData$immig1 - 1)
myData$immig2 <- -1 * (myData$immig2 - 1)
myData$immig3 <- -1 * (myData$immig3 - 1)
myData$guns1 <- -1 * (myData$guns1 - 1)
myData$guns2 <- -1 * (myData$guns2 - 1)
myData$guns3 <- -1 * (myData$guns3 - 1)
myData$abort1 <- -1 * (myData$abort1 - 1)
myData$abort2 <- -1 * (myData$abort2 - 1)
myData$abort3 <- -1 * (myData$abort3 - 1)

# select questions
myQuestions <- c("immig1", "immig2", "immig3", "abort1", "abort2", "abort3", "guns1", "guns2", "guns3")
												
# create data for IRT 
myData$myIds <- c(1:nrow(myData))
J <- length(myData$myIds)
I <- length(myQuestions)
cols <- which(colnames(myData) %in% c(myQuestions))
y <- as.vector(as.matrix(myData[,cols]))
N <- length(y)
jj <- rep(1:J, times=I)
ii <- rep(1:I, each=J)

# create indicator for 26 & 8
x <- rep(0,J)
x[26] <- 1
x[8] <- -1

# base IRT stan code
stanCode <- "
data {
  int<lower=1> I; // # questions
  int<lower=1> J; // # respondents
  int<lower=1> N; // # observations
  int<lower=1, upper=I> ii[N]; // question for observation n
  int<lower=1, upper=J> jj[N]; // respondent for observation n
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

stanFit <- stan(model_code=stanCode, data=stanData, iter=7000, warmup=3500, chains=4, verbose=TRUE, cores=4, control = list(max_treedepth = 20))

# =-=-=-=-=-=-=-=-=-=-=-=-=-=
#  Extract Thetas
# =-=-=-=-=-=-=-=-=-=-=-=-=-=

sumFit <- summary(stanFit, pars=c("theta"))$summary
ips <- sumFit[,1]

# =-=-=-=-=-=-=-=-=-=-=-=-=-=
#   Estimate rating scale model
# =-=-=-=-=-=-=-=-=-=-=-=-=-=

# drop everything except stanFit and ips
keep(stanFit, ips, sure=T)

# read in data
myData <- read_dta("Schaffner_Barney_Lucid_Module.dta")

# recode a few of these based on factor loadings.
myData$immig3 <- -1 * (myData$immig3 - 1)
myData$abort3 <- -1 * (myData$abort3 - 1)
myData$guns1 <- -1 * (myData$guns1 - 1)
myData$guns2 <- -1 * (myData$guns2 - 1)

# flip everything so conservative is positive
myData$immig1 <- -1 * (myData$immig1 - 1)
myData$immig2 <- -1 * (myData$immig2 - 1)
myData$immig3 <- -1 * (myData$immig3 - 1)
myData$guns1 <- -1 * (myData$guns1 - 1)
myData$guns2 <- -1 * (myData$guns2 - 1)
myData$guns3 <- -1 * (myData$guns3 - 1)
myData$abort1 <- -1 * (myData$abort1 - 1)
myData$abort2 <- -1 * (myData$abort2 - 1)
myData$abort3 <- -1 * (myData$abort3 - 1)

# select questions
myQuestions <- c("immig1", "immig2", "immig3", "abort1", "abort2", "abort3", "guns1", "guns2", "guns3")

# select scales
myScales <- c("imp_immig1",
			  "imp_immig2",
			  "imp_immig3",
			  "imp_abort1",
			  "imp_abort2",
			  "imp_abort3", 
			  "imp_guns1",
			  "imp_guns2",
			  "imp_guns3"
			  )
												
# create data for IRT 
myData$myIds <- c(1:nrow(myData))
J <- length(myData$myIds)
I <- length(myQuestions)
cols <- which(colnames(myData) %in% c(myQuestions))
y <- as.vector(as.matrix(myData[,cols]))
N <- length(y)
jj <- rep(1:J, times=I)
ii <- rep(1:I, each=J)
			  
# create matrix of weights
cols <- which(colnames(myData) %in% c(myScales))
spots <- as.matrix(myData[,cols])
spots <- 1 + round(spots/25)

# rescale data
y[y==0] <- -1
y <- y * spots
y <- y + 5

# format for Furr GRM
myList <- irt_data(y)

# estimate GRM
grmFit <- irt_stan(myList, "grsm_latent_reg.stan", chains=4, iter=300)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=
#  Extract Estimates
# =-=-=-=-=-=-=-=-=-=-=-=-=-=

sumFit <- summary(grmFit, pars=c("theta"))$summary
grips <- sumFit[,1]

keep(ips, grips, stanFit, grmFit, sure=T)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Effort Moderated IRT w/ questions from Lucid
# =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# read in data
myData <- read_dta("Schaffner_Barney_Lucid_Module.dta")

# recode a few of these based on factor loadings.
myData$immig3 <- -1 * (myData$immig3 - 1)
myData$abort3 <- -1 * (myData$abort3 - 1)
myData$guns1 <- -1 * (myData$guns1 - 1)
myData$guns2 <- -1 * (myData$guns2 - 1)

# flip everything so conservative is positive
myData$immig1 <- -1 * (myData$immig1 - 1)
myData$immig2 <- -1 * (myData$immig2 - 1)
myData$immig3 <- -1 * (myData$immig3 - 1)
myData$guns1 <- -1 * (myData$guns1 - 1)
myData$guns2 <- -1 * (myData$guns2 - 1)
myData$guns3 <- -1 * (myData$guns3 - 1)
myData$abort1 <- -1 * (myData$abort1 - 1)
myData$abort2 <- -1 * (myData$abort2 - 1)
myData$abort3 <- -1 * (myData$abort3 - 1)

# select questions
myQuestions <- c("immig1", "immig2", "immig3", "abort1", "abort2", "abort3", "guns1", "guns2", "guns3")

# select scales
myScales <- c("imp_immig1",
			  "imp_immig2",
			  "imp_immig3",
			  "imp_abort1",
			  "imp_abort2",
			  "imp_abort3", 
			  "imp_guns1",
			  "imp_guns2",
			  "imp_guns3"
			  )
												
# create data for IRT 
myData$myIds <- c(1:nrow(myData))
J <- length(myData$myIds)
I <- length(myQuestions)
cols <- which(colnames(myData) %in% c(myQuestions))
y <- as.vector(as.matrix(myData[,cols]))
N <- length(y)
jj <- rep(1:J, times=I)
ii <- rep(1:I, each=J)
			  
# create matrix of weights
cols <- which(colnames(myData) %in% c(myScales))
spots <- as.matrix(myData[,cols]/100)

# convert to dummy by respondent
for (s in 1:nrow(spots)){
	spots[s,] <- ifelse(spots[s,] >= median(spots[s,], na.rm=T), 1, 0)
}
spots <- as.vector(spots)


# create indicator for 26 & 8
x <- rep(0,J)
x[26] <- 1
x[8] <- -1

# effort moderated IRT stan code
stanCode <- "
data {
  int<lower=1> I; // # questions
  int<lower=1> J; // # respondents
  int<lower=1> N; // # observations
  int<lower=1, upper=I> ii[N]; // question for observation n
  int<lower=1, upper=J> jj[N]; // respondent for observation n
  int<lower=0, upper=1> y[N]; // response of observation n
  real<lower=0, upper=1> spots[N]; // differential item
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

# estimate
stanData <- list(J=J, I=I, N=N, jj=jj, ii=ii, y=y)

emFit <- stan(model_code=stanCode, data=stanData, iter=7000, warmup=3500, chains=4, verbose=TRUE, cores=4, control = list(max_treedepth = 20))

remFit <- emFit
sumFit <- summary(emFit, pars=c("theta"))$summary
remips <- sumFit[,1]

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Weighted IRT w/ questions from Lucid
# =-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# drop everything except stanFit and ips
keep(stanFit, grmFit, ips, grips, remFit, remips, sure=T)

# read in data
myData <- read_dta("Schaffner_Barney_Lucid_Module.dta")

# recode a few of these based on factor loadings.
myData$immig3 <- -1 * (myData$immig3 - 1)
myData$abort3 <- -1 * (myData$abort3 - 1)
myData$guns1 <- -1 * (myData$guns1 - 1)
myData$guns2 <- -1 * (myData$guns2 - 1)

# flip everything so conservative is positive
myData$immig1 <- -1 * (myData$immig1 - 1)
myData$immig2 <- -1 * (myData$immig2 - 1)
myData$immig3 <- -1 * (myData$immig3 - 1)
myData$guns1 <- -1 * (myData$guns1 - 1)
myData$guns2 <- -1 * (myData$guns2 - 1)
myData$guns3 <- -1 * (myData$guns3 - 1)
myData$abort1 <- -1 * (myData$abort1 - 1)
myData$abort2 <- -1 * (myData$abort2 - 1)
myData$abort3 <- -1 * (myData$abort3 - 1)

# select questions
myQuestions <- c("immig1", "immig2", "immig3", "abort1", "abort2", "abort3", "guns1", "guns2", "guns3")

# select scales
myScales <- c("imp_immig1",
			  "imp_immig2",
			  "imp_immig3",
			  "imp_abort1",
			  "imp_abort2",
			  "imp_abort3", 
			  "imp_guns1",
			  "imp_guns2",
			  "imp_guns3"
			  )
												
# create data for IRT 
myData$myIds <- c(1:nrow(myData))
J <- length(myData$myIds)
I <- length(myQuestions)
cols <- which(colnames(myData) %in% c(myQuestions))
y <- as.vector(as.matrix(myData[,cols]))
N <- length(y)
jj <- rep(1:J, times=I)
ii <- rep(1:I, each=J)
			  
# create matrix of weights
cols <- which(colnames(myData) %in% c(myScales))
spots <- as.matrix(myData[,cols]/100)
spots <- as.vector(spots)

# create indicator for 26 & 8
x <- rep(0,J)
x[26] <- 1
x[8] <- -1

# weighted IRT stan code
stanCode <- "
data {
  int<lower=1> I; // # questions
  int<lower=1> J; // # respondents
  int<lower=1> N; // # observations
  int<lower=1, upper=I> ii[N]; // question for observation n
  int<lower=1, upper=J> jj[N]; // respondent for observation n
  int<lower=0, upper=1> y[N]; // response of observation n
  real<lower=0, upper=1> spots[N]; // differential item
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

# estimate
stanData <- list(J=J, I=I, N=N, jj=jj, ii=ii, y=y)

wFit <- stan(model_code=stanCode, data=stanData, iter=7000, warmup=3500, chains=4, verbose=TRUE, cores=4, control = list(max_treedepth = 20))

sumFit <- summary(wFit, pars=c("theta"))$summary
wips <- sumFit[,1]


# drop everything
keep(stanFit, grmFit, wFit, ips, wips, grips, remFit, remips, sure=T)



# =-=-=-=-=-=-=-=-=-=-=--=-=-
#     Create dataset
# =-=-=-=-=-=-=-=--=-=-=-=--=

# recode direction
ips <- ips * -1
wip <- wips * -1
grips <- grips * -1
remips <- remips * -1


# add to dataset
myData$ips <- ips
myData$wip <- wips
myData$grips <- grips
myData$remips <- remips


# =-=-=-=-=-=-=-=-=-=-=-=-=-=
# Scatterplot Comparisons
# =-=-=-=-=-=-=-=-=-=-=-=-=-=

library(GGally)

# comparisons within groups
myData$myPid <- NA
myData$myPid[which(myData$Q49 < 4)] <- "Democrat"
myData$myPid[which(myData$Q49 > 4 & myData$Q49 < 8)] <- "Republican"
myData$myPid[which(myData$Q49 == 4 | myData$Q49 > 7)] <- "Other"

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
			columns = which(colnames(myData) %in% c("ips", "wip", "grips", "emips")), columnLabels = c("Standard", "Scale", "Moderated", "Weighted")) +
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
 
# adjust x axes
for (i in 1:p$nrow){
	p[i,1] <- p[i,1] + 
		xlim(-2,2.5)
	}
for (i in 2:p$nrow){
	p[i,2] <- p[i,2] + 
		xlim(-2,2.5)
	}
# adjust y axes
for (i in 1:(p$ncol-1)){
	p[4,i] <- p[4,i] + 
		ylim(-2.1,2.1)
	}
p[2,1] <- p[2,1] +
	ylim(-2.1, 2.1)
	

# =-=-=-=-=-=-=-=-=-=-=-=-=-=
# Predict policy choices
# =-=-=-=-=-=-=-=-=-=-=-=-=-=

library(MASS)

# estimate models; each of these questions is on a 5 point scale
i1 <- polr(as.factor(Q44_1) ~ ips, data=myData)
w1 <- polr(as.factor(Q44_1) ~ wip, data=myData)
g1 <- polr(as.factor(Q44_1) ~ grips, data=myData)
r1 <- polr(as.factor(Q44_1) ~ remips, data=myData)

i2 <- polr(as.factor(Q44_2) ~ ips, data=myData)
w2 <- polr(as.factor(Q44_2) ~ wip, data=myData)
g2 <- polr(as.factor(Q44_2) ~ grips, data=myData)
r2 <- polr(as.factor(Q44_2) ~ remips, data=myData)

i3 <- polr(as.factor(Q44_3) ~ ips, data=myData)
w3 <- polr(as.factor(Q44_3) ~ wip, data=myData)
g3 <- polr(as.factor(Q44_3) ~ grips, data=myData)
r3 <- polr(as.factor(Q44_3) ~ remips, data=myData)

# compute accuracy for 1st question
ipreds1 <- table(predict(i1, type="class"), myData$Q44_1)
iacc1 <- (ipreds1[1,1] + ipreds1[2,2] + ipreds1[3,3] + ipreds1[4,4] + ipreds1[5,5])/(sum(ipreds1))
# 42.4
wpreds1 <- table(predict(w1, type="class"), myData$Q44_1)
wacc1 <- (wpreds1[1,1] + wpreds1[2,2] + wpreds1[3,3] + wpreds1[4,4] + wpreds1[5,5])/(sum(wpreds1))
# 42.2 
gpreds1 <- table(predict(g1, type="class"), myData$Q44_1)
gacc1 <- (gpreds1[1,1] + gpreds1[2,2] + gpreds1[3,3] + gpreds1[4,4] + gpreds1[5,5])/(sum(gpreds1))
# 42.3
rpreds1 <- table(predict(r1, type="class"), myData$Q44_1)
racc1 <- (rpreds1[1,1] + rpreds1[2,2] + rpreds1[3,3] + rpreds1[4,4] + rpreds1[5,5])/(sum(rpreds1))
# 42.2

# compute accuracy for second question
ipreds2 <- table(predict(i2, type="class"), myData$Q44_2)
iacc2 <- (ipreds2[1,1] + ipreds2[2,2] + ipreds2[3,3] + ipreds2[4,4] + ipreds2[5,5])/(sum(ipreds2))
# 29.9
wpreds2 <- table(predict(w2, type="class"), myData$Q44_2)
wacc2 <- (wpreds2[1,1] + wpreds2[2,2] + wpreds2[3,3] + wpreds2[4,4] + wpreds2[5,5])/(sum(wpreds2))
# 30.0 (2 more)
gpreds2 <- table(predict(g2, type="class"), myData$Q44_2)
gacc2 <- (gpreds2[1,1] + gpreds2[2,2] + gpreds2[3,3] + gpreds2[4,4] + gpreds2[5,5])/(sum(gpreds2))
# 30.2 (5 more)
rpreds2 <- table(predict(r2, type="class"), myData$Q44_2)
racc2 <- (rpreds2[1,1] + rpreds2[2,2] + rpreds2[3,3] + rpreds2[4,4] + rpreds2[5,5])/(sum(rpreds2))
# 30.1

# compute accuracy for third question
ipreds3 <- table(predict(i3, type="class"), myData$Q44_3)
iacc3 <- (ipreds3[1,1] + ipreds3[2,2] + ipreds3[3,3] + ipreds3[4,4] + ipreds3[5,5])/(sum(ipreds3))
# 29.4
wpreds3 <- table(predict(w3, type="class"), myData$Q44_3)
wacc3 <- (wpreds3[1,1] + wpreds3[2,2] + wpreds3[3,3] + wpreds3[4,4] + wpreds3[5,5])/(sum(wpreds3))
# 30.4 (14 more)
gpreds3 <- table(predict(g3, type="class"), myData$Q44_3)
gacc3 <- (gpreds3[1,1] + gpreds3[2,2] + gpreds3[3,3] + gpreds3[4,4] + gpreds3[5,5])/(sum(gpreds3))
# 29.4
rpreds3 <- table(predict(r3, type="class"), myData$Q44_2)
racc3 <- (rpreds3[1,1] + rpreds3[2,2] + rpreds3[3,3] + rpreds3[4,4] + rpreds3[5,5])/(sum(rpreds3))
# 24.1

# create predicted probability plots
ipsData <- data.frame(ips = seq(from = round(mean(myData$ips)-2*sd(myData$ips),1),
							 to = round(mean(myData$ips)+2*sd(myData$ips),1),
							 by = .01))
wipData <- data.frame(wip = seq(from = round(mean(myData$wip)-2*sd(myData$wip),1),
							 to = round(mean(myData$wip)+2*sd(myData$wip),1),
							 by = .01))							 
gripsData <- data.frame(grips = seq(from = round(mean(myData$grips)-2*sd(myData$grips),1),
							 to = round(mean(myData$grips)+2*sd(myData$grips),1),
							 by = .01))
remipsData <- data.frame(remips = seq(from = round(mean(myData$remips)-2*sd(myData$remips),1),
							 to = round(mean(myData$remips)+2*sd(myData$remips),1),
							 by = .01))
							 
# predicted probabilities and plots for democrats
ips1PP <- cbind(ipsData, predict(i1, newdata=ipsData, type="probs"))
ips2PP <- cbind(ipsData, predict(i2, newdata=ipsData, type="probs"))
ips3PP <- cbind(ipsData, predict(i3, newdata=ipsData, type="probs"))
wip1PP <- cbind(wipData, predict(w1, newdata=wipData, type="probs"))
wip2PP <- cbind(wipData, predict(w2, newdata=wipData, type="probs"))
wip3PP <- cbind(wipData, predict(w3, newdata=wipData, type="probs"))
grips1PP <- cbind(gripsData, predict(g1, newdata=gripsData, type="probs"))
grips2PP <- cbind(gripsData, predict(g2, newdata=gripsData, type="probs"))
grips3PP <- cbind(gripsData, predict(g3, newdata=gripsData, type="probs"))
remips1PP <- cbind(remipsData, predict(r1, newdata=remipsData, type="probs"))
remips2PP <- cbind(remipsData, predict(r2, newdata=remipsData, type="probs"))
remips3PP <- cbind(remipsData, predict(r3, newdata=remipsData, type="probs"))

# convert to plotting forms
ips1PP <- melt(ips1PP, id.vars = c("ips"), value.name = "probability")
ips2PP <- melt(ips2PP, id.vars = c("ips"), value.name = "probability")
ips3PP <- melt(ips3PP, id.vars = c("ips"), value.name = "probability")
wip1PP <- melt(wip1PP, id.vars = c("wip"), value.name = "probability")
wip2PP <- melt(wip2PP, id.vars = c("wip"), value.name = "probability")
wip3PP <- melt(wip3PP, id.vars = c("wip"), value.name = "probability")
grips1PP <- melt(grips1PP, id.vars = c("grips"), value.name = "probability")
grips2PP <- melt(grips2PP, id.vars = c("grips"), value.name = "probability")
grips3PP <- melt(grips3PP, id.vars = c("grips"), value.name = "probability")
remips1PP <- melt(remips1PP, id.vars = c("remips"), value.name = "probability")
remips2PP <- melt(remips2PP, id.vars = c("remips"), value.name = "probability")
remips3PP <- melt(remips3PP, id.vars = c("remips"), value.name = "probability")


# plot for response 1
ips1Plot <- ggplot(ips1PP, aes(x=ips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Standard (Accuracy = 42.4)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")

wip1Plot <- ggplot(wip1PP, aes(x=wip, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Weighted (Accuracy = 42.2)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")
	
grips1Plot <- ggplot(grips1PP, aes(x=grips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Scale (Accuracy = 42.3)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")
	
remips1Plot <- ggplot(remips1PP, aes(x=remips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Moderated (Accuracy = 42.2)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")
	
grid.arrange(ips1Plot, wip1Plot, grips1Plot, remips1Plot, nrow=1)

# plot for response 2
ips2Plot <- ggplot(ips2PP, aes(x=ips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Standard (Accuracy = 29.9)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")
	
wip2Plot <- ggplot(wip2PP, aes(x=wip, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Weighted (Accuracy = 30.0)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")
	
grips2Plot <- ggplot(grips2PP, aes(x=grips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Scale (Accuracy = 30.2)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")
	
remips2Plot <- ggplot(remips2PP, aes(x=remips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Moderated (Accuracy = 30.1)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")
	
grid.arrange(ips2Plot, wip2Plot, grips2Plot, emips2Plot, nrow=1)
	
# plot for response 3
ips3Plot <- ggplot(ips3PP, aes(x=ips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Standard (Accuracy = 29.4)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")

wip3Plot <- ggplot(wip3PP, aes(x=wip, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Weighted (Accuracy = 30.4)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")

grips3Plot <- ggplot(grips3PP, aes(x=grips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Scale (Accuracy = 29.4)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")
	
remips3Plot <- ggplot(remips3PP, aes(x=remips, y=probability, group=variable)) + 
	geom_line(aes(color=variable), size=1.5, alpha=.8) +
	ylim(0,0.5) +
	xlab("Moderated (Accuracy = 24.1)") +
	ylab("Probability of Vote") +
	scale_color_grey() +
	theme_bw() +
	theme(legend.position="none")
	
grid.arrange(ips3Plot, wip3Plot, grips3Plot, emips3Plot, nrow=1)
		

# =-=-=-=-=-=-=-=-=-=-=-=-=
# Bimodality coefficient
# =-=-=-=-=-=-=-=-=-=-=-=-=

library(mousetrap)
bimodality_coefficient(myData$ips) # 40.2
bimodality_coefficient(myData$wip) # 40.8
bimodality_coefficient(myData$remips) # 41.3
bimodality_coefficient(myData$grips) # 51.1