# Application of the Outcome Weights Framework for Double Machine Learning to the LaLonde Study
# Laura Kreisel
# 2025

## Original data from LaLonde
# @imbens&xu (2024)

library(haven)
library(labelled)

############################ Data preparation ############################

process.data  <- function(filename) {
  d <- haven::read_dta(filename)
  d <- as.data.frame(zap_label(d))
  for (i in 1:ncol(d)) {
    if (is.numeric(d[, i])) {
      d[, i] <- as.vector(d[, i])
    }
  }
  if ("re74" %in% colnames(d)) {d$u74 <- ifelse(d$re74 == 0, 1, 0)}
  d$u75 <- ifelse(d$re75 == 0, 1, 0)
  return(d)
}

# experimental sample
nsw <- process.data("../data/lalonde/nsw.dta")
table(nsw$treat)
nsw$sample <- NA
nsw$sample[which(nsw$treat == 1)] <- 0
nsw$sample[which(nsw$treat == 0)] <- 0.5
table(nsw$treat, nsw$sample)
nsw_tr <- subset(nsw, treat == 1)
nsw_co <- subset(nsw, treat == 0)

ldw <- process.data("../data/lalonde/nsw_dw.dta")
ldw$sample <- NA
ldw$sample[which(ldw$treat == 1)] <- 1
ldw$sample[which(ldw$treat == 0)] <- 2
table(ldw$treat, ldw$sample)
ldw_tr <- ldw[which(ldw$treat == 1), ] # experimental treated 
ldw_co <- ldw[which(ldw$treat == 0), ] # experimental control
cps1 <- process.data("../data/lalonde/cps_controls.dta")
psid1 <- process.data("../data/lalonde/psid_controls.dta")

# nonexperimental samples
cps1$sample <- 3
ldw_cps <- rbind.data.frame(ldw_tr, cps1)
table(ldw_cps$treat)

psid1$sample <- 4
ldw_psid <- rbind.data.frame(ldw_tr, psid1)
table(ldw_psid$treat)

cps1a <- subset(cps1, select =  -c(re74, u74))
names(nsw)
names(cps1a)
nsw_cps <- rbind.data.frame(nsw_tr, cps1a)
table(nsw_cps$treat)

psid1a <- subset(psid1, select =  -c(re74, u74))
nsw_psid <- rbind.data.frame(nsw_tr, psid1a)
table(nsw_psid$treat)

############################ Summary Statistics ############################

X0 <- c("age", "education", "nodegree", "married", "black", "hispanic",  
        "re75", "u75")

X <- c("age", "education", "nodegree", "married", "black", "hispanic", 
       "re75", "u75", "re74", "u74")

# mean
nsw.tr <- c(apply(nsw_tr[, X0], 2, mean), rep(NA, 2))
nsw.co <- c(apply(nsw_co[, X0], 2, mean), rep(NA, 2))
ldw.tr <- c(apply(ldw_tr[, X], 2, mean))
ldw.co <- apply(ldw_co[, X], 2, mean)
cps <- apply(cps1[, X], 2, mean)
psid <- apply(psid1[, X], 2, mean)
out <- cbind.data.frame(nsw.tr, nsw.co, cps, psid, ldw.tr, ldw.co)
out[c("re75", "re74"), ] <- out[c("re75", "re74"), ]/1000

# sd
nsw.tr <- c(apply(nsw_tr[, X0], 2, sd), rep(NA, 2))
nsw.co <- c(apply(nsw_co[, X0], 2, sd), rep(NA, 2))
ldw.tr <- c(apply(ldw_tr[, X], 2, sd))
ldw.co <- apply(ldw_co[, X], 2, sd)
cps <- apply(cps1[, X], 2, sd)
psid <- apply(psid1[, X], 2, sd)
out2 <- cbind.data.frame(nsw.tr, nsw.co, cps, psid, ldw.tr, ldw.co)
out2[c("re75", "re74"), ] <- out2[c("re75", "re74"), ]/1000

# count
count.nsw.tr <- nrow(nsw[which(nsw$treat == 1), ])
count.nsw.co <- nrow(nsw[which(nsw$treat == 0), ])
count.ldw.tr <- nrow(ldw[which(ldw$treat == 1), ])
count.ldw.co <- nrow(ldw[which(ldw$treat == 0), ])
count.cps     <- nrow(cps1)                       
count.psid    <- nrow(psid1)                      
out <- rbind(out, Counts = c(count.nsw.tr, count.nsw.co, count.cps, count.psid, count.ldw.tr, count.ldw.co))
out2 <- rbind(out2, Counts = c(NA, NA, NA, NA, NA, NA))

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
write.csv(sav, file = "../tables/stats_nsw_ldw.csv", row.names = FALSE)
