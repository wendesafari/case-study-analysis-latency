###############################################################################
#
#  Title:    Emulating the results shown in Section 4 of the paper titled:
#            "Latency function estimation under the mixture cure model when the cure
#            status is available" 
#
#  Script:    case_study_analysis.R
#
#  Authors:   W. C. Safari | I.  Lopez-de-Ullibarri | M.A. Jacome
#
#  Date:      24 October 2022
#
#  Purpose:   This script contains the code for estimating the survival function 
#             of the individuals experiencing the event S0(t) and S0(t|x) 
###############################################################################

rm(list=ls()); while(!dev.cur())dev.off(); cat('\014')

# Necessary packages, if not installed the next code will install them and then load
PkgList <- c("data.table", "DescTools", "MASS", 
             "microbenchmark", "foreach", 
             "doParallel","doRNG","doSNOW", "readr", "gtsummary", "dplyr") 
new.packages <- PkgList[!(PkgList %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(PkgList, function(x)do.call("require", list(x)))

# Source functions 
invisible(source("functions/functions_definitions.R"))

# Read in the simulated COVID-19 data
sim_data <- read_csv("data/simulated_data.csv")

# The simulated_data.csv dataframe has 500 rowa and 6 columns, mimicking  
# the COVID-19 dataset in Safari et al (2022b) latency paper. The dataframe 
# contains the following columns:

# 1.- Sex - Dummy variable: 1 = Female, 0 = Male.

# 2.- Age (in years) - A numeric vector contains age of the patients.

# 3.- time (in days) - Time  of patients from hospital ward (HW) to either ICU or to being discharged
#       or to death or time spent in the hospital ward.

# 4.- delta - Dummy variable: 1 = ICU admission, 0 = otherwise.

# 5.- cure.known - Dummy variable: (0 = still in HW or ICU, 1 = otherwise).

# 6.- status -  the status of the patient at the end of the study 
#           (HW, ICU, discharged, Died).


sim_data <- sim_data %>% mutate_if(is.character, as.factor) 
glimpse(sim_data)

# Summary statistics
tab <- sim_data %>%
  dplyr::select(status, time) %>%
  mutate_if(is.character, as.factor) %>%
  tbl_summary(by = status) %>%
  bold_labels()
tab

# Extract the dataframe for the estimation
dfr <- sim_data[, c("age", "time", "delta", "cure.known")]

# Renaming the column 
# x -> age
# t -> time
# d -> delta
# xinu -> cure.known

colnames(dfr) <- c("x", "t", "d", "xinu")

# With the existence of tied cases at some of the survival times, we ordered in the
# following order: uncensored – censored unknown cure status – known cures.

ord.dfr <- as.data.frame(dfr[order(dfr$t, - dfr$d, dfr$xinu),])

# Estimate the latency for patients aged 54 year old

x0 <- 54

# Compute the bootstrap bandwidths
# Warning: it takes an average of ~2 seconds to run. 

latencyhboot <- latencybootfun(dfr = dfr, x0 = x0)

# Estimate the latency
# With the two bandwidths estimator

# We need to estimate the bootstrap bandwidths (h1, h2) from
# the estimators Shc(t|x) and 1-phc(x) to be used in the latency estimation.

# First estimate h1 from

Sch <- survfitcurePK(x, t, d, xinu, dfr, x0, latencyhboot[1])$S

# and h2 from
pch <-  survfitcurePK(x, t, d, xinu, dfr , x0, latencyhboot[2])$p

# Latency estimation S0(t|x) using the proposed estimator with two bandwidths.

if(pch == 0) {
  S0ch <- pch*Sch
  tmax <- quantile(ord.dfr[, 2], probs = 1 - 0.9)[[1]]
} else {
  S0ch <- ifelse(((Sch - (1 - pch))/pch) > 0 & ord.dfr$t > tail(ord.dfr$t[ord.dfr$d == 1], 1), 0,
                ifelse(((Sch - (1 - pch)) / pch) < 0, 0, ((Sch - (1 - pch)) / pch)))
  tmax <- data.table::first(ord.dfr$t[S0ch == data.table::first(DescTools::Closest(S0ch, 1 - 0.9, na.rm = T))])
  
}  

# Median time and its associated Intequartile range
t75 <- first(ord.dfr$t[S0ch == first(Closest(S0ch, 0.75, na.rm = T))])
tmed <- first(ord.dfr$t[S0ch == first(Closest(S0ch, 0.5, na.rm = T))])
t25 <- first(ord.dfr$t[S0ch == first(Closest(S0ch, 0.25, na.rm = T))])
c(t75, tmed, t25)

# Latency estimated plot for 54-year old patients
tgrid <- seq(0, tmax, length.out = nrow(dfr))
plot(c(0, tgrid), c(1, S0ch), type = "s", lwd=2,  ylab = "Latency", 
     xlab = "Time (days)", ylim = c(0, 1))

# The estimation of the latency conditioned on sex (for male and female separately)

# Extract the dataframe for the unconditional estimation
dfr <- sim_data[, c("sex", "time", "delta", "cure.known")]
colnames(dfr) <- c("sex", "t", "d", "xinu")

# Female
S_female <- survfitcurePK_un(dataset = dfr[dfr$sex == "Female", 2:4])

# p : probability of experiencing the final outcome
p_female <- S_female[[2]]

# Time
time <- S_female[[3]]

# S0(t): Survival function of the individuals experiencing the event
S0_female <- (S_female[[1]] - (1 - p_female))/p_female

tmed <- first(time[S0_female == first(Closest(S0_female, 0.5, na.rm = T))])
t25 <- first(time[S0_female == first(Closest(S0_female, 0.25, na.rm = T))])
t75 <- first(time[S0_female == first(Closest(S0_female, 0.75, na.rm = T))])
c(t75, tmed, t25)


# male
S_male <- survfitcurePK_un(dataset = dfr[dfr$sex == "Male", 2:4])

# p : probability of experiencing the final outcome
p_male <- S_male[[2]]

# Time
time <- S_male[[3]]

# S0(t): Survival function of the individuals experiencing the event
S0_male <- (S_male[[1]] - (1 - p_male))/p_male

tmed <- first(time[S0_male == first(Closest(S0_male, 0.5, na.rm = T))])
t25 <- first(time[S0_male == first(Closest(S0_male, 0.25, na.rm = T))])
t75 <- first(time[S0_male == first(Closest(S0_male, 0.75, na.rm = T))])
c(t75, tmed, t25)


# Latency by sex
plot(c(0,S_female[[3]]), c(1,S0_female), type = 's', lwd=2, 
     ylab =  "Latency", xlab = "Time (days)", col="black",cex.lab=1.2, main = "Gender")
lines(c(0,S_male[[3]]), c(1,S0_male), type = 's', lty=1, lwd=2, col="blue")
legend("topright", 
       legend = c("Female", "Male"), 
       col = c("black",  "blue"), 
       bty = "n", 
       cex = 1.2, 
       text.col = "black",
       lty = c(1,1), lwd = 2)
