# MLE method for Individual size distributions from Edwards et al. 

# this is an old copied script from a previous version of the NEON size spectra analysis


# load packages
# packages
library(sizeSpectra)
library(tidyverse)
library(brms)
library(janitor)
library(ggridges)

### Custom Functions ###
# the following function is modified from the eightMethodsMEE() function in the sizeSpectra package
# it is modified to only return the MLE estimate and 95 %CI for the bounded power law exponent, b, and is modified to work with list-columns in the tidyverse.  

# df is the data frame with every observation in it's own row
# rsp_var is the column to perform MLE on
# i.e. in this analysis the column is "dw" for estimated dry weight
MLE_tidy <- function(df, rsp_var){
  # define variables
  x <- df[[rsp_var]]
  xmin = min(x)
  xmax = max(x)
  log.x = log(x)
  sum.log.x = sum(log.x)
  
  # initial starting point for parameter estimate
  PL.bMLE = 1/(log(min(x)) - sum.log.x/length(x)) - 1
  
  # non-linear minimization  
  PLB.minLL = nlm(negLL.PLB, 
                  p = PL.bMLE, x = x, n = length(x), 
                  xmin = xmin, xmax = xmax,
                  sumlogx = sum.log.x)
  
  # estimate for b
  PLB.bMLE = PLB.minLL$estimate
  # minimum estimate of b
  PLB.minNegLL = PLB.minLL$minimum
  
  ## 95% CI calculation
  bvec = seq(PLB.bMLE - 1, PLB.bMLE + 1, 1e-05) # original =-0.5
  PLB.LLvals = vector(length = length(bvec))
  for (i in 1:length(bvec)) {
    PLB.LLvals[i] = negLL.PLB(bvec[i],
                              x = x,
                              n = length(x), 
                              xmin = xmin,
                              xmax = xmax,
                              sumlogx = sum.log.x)
  }
  critVal = PLB.minNegLL + qchisq(0.95, 1)/2
  bIn95 = bvec[PLB.LLvals < critVal]
  # confidence interval
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))
  if (PLB.MLE.bConf[1] == min(bvec) | 
      PLB.MLE.bConf[2] == max(bvec)) {
    dev.new()
    plot(bvec, PLB.LLvals)
    abline(h = critVal, col = "red")
    stop("Need to make bvec larger - see R window")
  }
  # return b estimate and min/max 95% CI
  return(data.frame(b = PLB.bMLE,
                    minCI = min(bIn95),
                    maxCI = max(bIn95)))
}

# modified function to plot MLE estimate of size spectra (ISD)
isd_plot <- function (x, b, confVals = NULL,
                      panel = "b", log.xy = "xy",
                      mgpVals = c(1.6, 0.5, 0),
                      inset = c(0, -0.04),
                      xlim_global = NA,
                      ylim_global = NA, ...) 
{
  if (is.na(xlim_global[1])) {
    xlim_global = c(min(x), max(x))
  }
  if (is.na(ylim_global[1])) {
    ylim_global = c(1, length(x))
  }
  plot(sort(x, decreasing = TRUE), 1:length(x), log = log.xy, 
       xlab = expression(paste("Values, ", italic(x))), 
       ylab = expression(
         paste("Number of ", values >= x), sep = ""),
       mgp = mgpVals, xlim = xlim_global, 
       ylim = ylim_global, axes = FALSE, ...)
  xLim = 10^par("usr")[1:2]
  yLim = 10^par("usr")[3:4]
  if (log.xy == "xy") {
    logTicks(xLim, yLim, xLabelSmall = c(5, 50, 500))
  }
  if (log.xy == "x") {
    mgpVal = c(2, 0.5, 0)
    logTicks(xLim, yLim = NULL, xLabelSmall = c(5, 50, 500), 
             mgpVal = mgpVal)
    yBig = c(0, 500, 1000)
    axis(2, at = yBig, labels = yBig, mgp = mgpVal)
    axis(2, seq(yBig[1], yBig[length(yBig)], by = 100),
         labels = rep("", 11), tcl = -0.2, mgp = mgpVal)
  }
  x.PLB = seq(min(x), max(x), length = 1000)
  y.PLB = (1 - pPLB(x = x.PLB,
                    b = b,
                    xmin = min(x.PLB),
                    xmax = max(x.PLB))) * 
    length(x)
  lines(x.PLB, y.PLB, col = "red")
  if (panel == "b") {
    for (i in c(1, length(confVals))) {
      lines(x.PLB,
            (1 - pPLB(x = x.PLB, b = confVals[i],
                      xmin = min(x.PLB),
                      xmax = max(x.PLB))) * length(x), 
            col = "red", lty = 2)
    }
    #legend("topright", "(b)", bty = "n", inset = inset)
  }
  if (panel == "h") {
    legJust(c("(h) MLE",
              paste("b=", signif(b, 3), sep = "")),
            inset = inset, logxy = TRUE)
  }
}

# helper function to plot across lists of results
plot_b_est <- function(dat, b, grp_var, ...){
  isd_plot(sample(dat$dw, 10000),
           b = b$b,
           confVals = c(b$minCI, b$maxCI))
  # add labels
  mtext(paste0('Site = ',
               dat$siteID[1],
               as.character(dat[[grp_var]][1])),
        side = 3, line = -6, adj = 0.01)
  mtext(paste0("b = ", round(b, digits = 2)),
        side = 3, line = -7, adj = 0.05)
  # some of the plots have a text error in b-estimate, not sure why
}

# load and wrangle data--------------------------------------------------


# read in estimated dry weight data
dw <- readRDS("data/ark_dw.RDS") 
# remove lake and large river sites


dw %>% group_by(site) %>%
  summarize(min_dw = min(dw),
             max_dw = max(dw))


# MLE ---------------------------------------------------------------------


# estimate MLE for each site and collection 
# note, this is running across 18+ million rows of data, it can take a while
mle.dw <- dw %>%
  #filter(!is.na(dw)) %>%
  #mutate(date = as.Date(collectDate)) %>%
  group_by(site) %>%
  # create list-column
  nest() %>% 
  # estimate b and 95% CI
  mutate(mle = map(data, MLE_tidy, rsp_var = "dw")) %>%
  unnest(cols = mle)

mle.dw.rep <- dw %>%
  group_by(site, rep) %>%
  # create list-column
  nest() %>% 
  # estimate b and 95% CI
  mutate(mle = map(data, MLE_tidy, rsp_var = "dw")) %>%
  unnest(cols = mle)

mle.dw.rep %>% group_by(site) %>%
  summarize(mean(b),
            sd(b),
            mean(minCI),
            mean(maxCI)) 