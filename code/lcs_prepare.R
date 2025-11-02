# Application of the Outcome Weights Framework for Double Machine Learning to the Lalonde Study
# Laura Kreisel
# 2025

## Female Sample Constructed by Calonico and Smith (2017)
# @imbens&xu (2024)

library(haven)
library(labelled)

############################ Data preparation ############################

d <- haven::read_dta("data/cs/NSW_AFDC_CS.dta")
d <- as.data.frame(zap_label(d))
for (i in 1:ncol(d)) {
  if (is.numeric(d[, i])) {
    d[, i] <- as.vector(d[, i])
  }
}

# keep reconstructed Lalonde Sample & PSID-1
d <- subset(d, lalonde == 1 | psid1 == 1)
names(d)
d$u75 <- ifelse(d$re75 == 0, 1, 0)

# keep relevant variables
covar <- c("age", "educ", "nodegree", "married", "black", "hisp", "nchildren75", "re75", "u75", "afdc75")
Y <- "re79"

var <- c(Y, "treated", covar, "psid1")
d <- d[, var]

colnames(d)[2] <- "treat"
names(d)
d$treat[which(d$psid1 == 1)] <- 0 # nonexperimental controls

# experimental sample
lcs <- subset(d, psid1 == 0)
lcs$psid1 <- NULL
table(lcs$treat)
names(lcs)
sum(!complete.cases(lcs)) # check missingness

# nonexperimental sample
lcs_psid <- subset(d, treat == 1 | psid1 == 1)
lcs_psid$treat[which(lcs_psid$psid1 == 1)] <- 0
lcs_psid$psid1 <- NULL
table(lcs_psid$treat)
names(lcs_psid)
sum(!complete.cases(lcs_psid)) # check missingness

# nonexperimental sample (PSID2)
lcs_psid2 <- subset(lcs_psid, afdc75 == 1)
lcs_psid2$treat[which(lcs_psid2$treat == 0)] <- 0 
table(lcs_psid2$treat)
sum(!complete.cases(lcs_psid2)) # check missingness
dim(lcs_psid2)

# extended sample
lcs_psid.plus <- d
lcs_psid.plus$exp <- 1 - lcs_psid.plus$psid1

############################ Summary Statistics ############################

X <- c("age", "educ", "nodegree", "married", "black", "nchildren75", "hisp", "re75", "u75", "afdc75")

# mean
stat.lcs.tr <- apply(lcs[which(lcs$treat == 1), X], 2, mean)
stat.lcs.co <- apply(lcs[which(lcs$treat == 0), X], 2, mean)
stat.psid <- apply(lcs_psid[which(lcs_psid$treat == 0), X], 2, mean)
stat.psid2 <- apply(lcs_psid2[which(lcs_psid2$treat == 0), X], 2, mean)
out <- cbind.data.frame(stat.lcs.tr, stat.lcs.co, stat.psid, stat.psid2)
out[c("re75"), ] <- out[c("re75"), ]/1000

# sd
stat.lcs.tr <- apply(lcs[which(lcs$treat == 1), X], 2, sd)
stat.lcs.co <- apply(lcs[which(lcs$treat == 0), X], 2, sd)
stat.psid <- apply(lcs_psid[which(lcs_psid$treat == 0), X], 2, sd)
stat.psid2 <- apply(lcs_psid2[which(lcs_psid2$treat == 0), X], 2, sd)
out2 <- cbind.data.frame(stat.lcs.tr, stat.lcs.co, stat.psid, stat.psid2)
out2[c("re75"), ] <- out2[c("re75"), ]/1000

#count
count.lcs.tr <- nrow(lcs[which(lcs$treat == 1), ])
count.lcs.co <- nrow(lcs[which(lcs$treat == 0), ])
count.psid <- nrow(lcs_psid[which(lcs_psid$treat == 0), ])
count.psid2 <- nrow(lcs_psid2[which(lcs_psid2$treat == 0), ])
out <- rbind(out, Counts = c(count.lcs.tr, count.lcs.co, count.psid, count.psid2))
out2 <- rbind(out2, Counts = c(NA, NA, NA, NA))

# columns are samples
n <- nrow(out)   
sav <- matrix("", n*2, ncol(out))
for (j in 1:n) {
  a <- out[j, ]
  b <- out2[j, ]
  if (rownames(out)[j] == "Counts") {
    sav[j*2-1, ] <- sprintf("%.0f", a)
    sav[j*2, ] <- rep("", ncol(out))  
  } else {
    sav[j*2-1, ] <- sprintf("%.2f", a)
    sav[j*2, ] <- paste0("(", sprintf("%.2f",b), ")")
  }
}
rownames(sav) <- rep(rownames(out), each = 2)
print(sav)

table(lcs$treat)
table(lcs_psid$treat)

write.csv(sav, file = "tables/stats_lcs.csv", row.names = FALSE)

############################ Trimming ############################

Y <- "re79"
treat <- "treat"

# redefine covariates: removing "nchildren75" to be used as placebo outcome
covar <- c("age", "educ", "nodegree", "married", "black", "hisp", "re75", "u75")

## trimming
p.forest <- probability_forest(X = lcs_psid.plus[, covar], 
                               Y = as.factor(lcs_psid.plus$exp), seed = 1234, num.trees = 4000) # experimental indicator as treated
lcs_psid.plus$ps <- p.forest$predictions[,2]
range(lcs_psid.plus$ps)
#hist(lcs_psid.plus$ps)

# histogram of propensity scores
round(quantile(lcs_psid.plus$ps[which(lcs_psid.plus$exp == 0 & lcs_psid.plus$ps > 0.6)], probs = seq(0, 1, 0.1)),3)
threshold <- 0.9 # only 50 PSID control units have ps > 0.9
nrow(subset(lcs_psid.plus, exp == 0 & ps > threshold)) # 27 control units
nrow(subset(lcs_psid.plus, exp == 1 & ps > threshold))

# trim experimental data
lcs_trim_psid <- subset(lcs_psid.plus, exp == 1 & ps <= threshold) 
lcs_trim_psid$ps <- NULL
lcs_trim_psid$exp <- NULL
lcs_trim_psid$psid1 <- NULL
table(lcs_trim_psid$treat) # 773 units: 393 treated vs 380 controls; 1185-773 = 412 experimental units dropped
head(lcs_trim_psid)
summary(lm(re79 ~ treat, data = lcs_trim_psid))

# propensity score matching without replacement
set.seed(1234) # need to set seed b/c tie-breaking is random
data <- subset(lcs_psid.plus, (treat == 1 | exp == 0) & ps <= threshold) # only use psid1 controls
table(data$treat)
data$ps <- probability_forest(X = data[, covar], 
                              Y = as.factor(data$treat), seed = 1234, num.trees = 4000)$predictions[,2]
mout <- Match(Y = data$re78, Tr = data$treat, X = data$ps, estimand = "ATT", M = 1,
              BiasAdjust = FALSE, replace=FALSE, ties = FALSE)
lcs_psid_trim <- data[c(mout$index.treated, mout$index.control), ]
table(lcs_psid_trim$treat) # 1064 units: 532 treated vs 532 controls

# estimate propensity score again
lcs_psid_trim$ps_new <- probability_forest(X = lcs_psid_trim[, covar], 
                                           Y = as.factor(lcs_psid_trim$treat), 
                                           seed = 1234, num.trees = 4000)$predictions[,2]
table(lcs_psid_trim$treat)

save(lcs, lcs_psid, lcs_trim_psid, lcs_psid_trim, file = "data/lcs.RData")


