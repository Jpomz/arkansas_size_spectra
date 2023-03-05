# mle tidy

library(tidyverse)
library(sizeSpectra)
library(ggpubr)

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
  bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 1e-05)
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


dat <- readRDS("data/ark_dw.RDS")

dat <- dat %>%
  mutate(y_fact = 
           factor(year,
                  levels = c("1990", "2012", "2019")))

# how many individuals by site and year?
dat %>%
  filter(dw >0.0026) %>%
  group_by(site, y_fact) %>%
  count()

# how many individuals per sample?
dat %>%
  filter(dw >0.0026) %>%
  group_by(site, rep, y_fact) %>%
  count() 

dat %>%
  filter(dw >0.0026) %>%
  group_by(site, rep, y_fact) %>%
  count() %>%
  ggplot(aes(x = y_fact,
             y = n, 
             color = site,
             shape = rep)) +
  geom_point(
    position = position_jitter(
      width = 0.1,
      height = 0
    ))

mle_lambda <- dat %>%
  filter(dw >0.0026) %>%
  group_by(site, y_fact) %>%
  nest() %>%
  mutate(lambda = map(data,
                   MLE_tidy,
                   "dw")) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup()

mle_lambda %>%
  arrange(site, y_fact)

mle_lambda %>%
  ggplot(aes(y = b,
             ymin = minCI,
             ymax = maxCI,
             x = y_fact,
             color = site)) +
  geom_pointrange(
    size = 1,
    position = position_dodge(
      width = 0.25
    )) +
  theme_bw() +
  labs(y = expression(lambda),
       x = "year")

ggline(mle_lambda, 
       x = "y_fact",
       y = "b",
       color = "site", 
       size = 1)

# double the minimum filter
dat %>%
  filter(dw >0.0052) %>%
  group_by(site, y_fact) %>%
  count() 
# only a handful fewer individuals  

dat %>%
  filter(dw >0.0052) %>%
  group_by(site, y_fact) %>%
  nest() %>%
  mutate(lambda = map(data,
                      MLE_tidy,
                      "dw")) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup()%>%
  arrange(site, y_fact) %>%
  ggplot(aes(y = b,
             ymin = minCI,
             ymax = maxCI,
             x = y_fact,
             color = site)) +
  geom_pointrange(size = 1) +
  theme_bw() +
  labs(y = expression(lambda),
       x = "year")

dat %>%
  filter(dw >0.0026) %>%
  group_by(site, rep, y_fact) %>%
  count()

dat %>%
  filter(dw >0.0026) %>%
  group_by(site, y_fact, rep) %>%
  nest() %>%
  mutate(lambda = map(data,
                      MLE_tidy,
                      "dw")) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup() %>%
  group_by(site, y_fact) %>%
  summarise(mean_b = mean(b),
            sd(b))
