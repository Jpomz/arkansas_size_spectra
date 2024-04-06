# lambda MLE by rep

library(sizeSpectra)
library(tidyverse)


# custom function 
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

# read in dry weight data
dat <- readRDS("data/ark_dw.RDS")
dat %>%
  distinct(year) %>%
  arrange(year)

# set years as ordered factors
dat <- dat %>%
  mutate(y_fact = 
           factor(year,
                  levels = c("1990",
                             "1999",
                             "2012",
                             "2019",
                             "2021"))) %>%
  as_tibble()

dat %>%
  filter(dw <=0.0026)
# 0 rows, already filtered


# Per replicate -----------------------------------------------------------


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


mle_lambda_rep <- dat %>%
  group_by(year, site, rep) %>%
  nest() %>%
  mutate(lambda = map(data,
                      MLE_tidy,
                      "dw")) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup()


mle_lambda_rep %>%
  ggplot(aes(x = year, 
             y = b,
             ymin = minCI, 
             ymax = maxCI,
             color = site, 
             group = rep)) +
  geom_pointrange(
    position = position_dodge(width = 1)
  ) +
  scale_color_manual(values = c("#019AFF", "#FF1984")) +
  # stat_smooth(method = "lm", 
  #             inherit.aes = FALSE,
  #             aes(x = year, y = b, color = site)) +
  NULL


ols <- lm(b~site*year, dat = mle_lambda_rep)
summary(ols)
ggplot(ols,
       aes(x = .fitted,
           y = .resid)) +
  geom_point()

# weighted regression 
# gamma = SE
# gamma = (b_high - b_low) / (2 * 1.96)

mle_lambda_rep %>%
  ggplot(aes(sample = b)) +
  stat_qq()+ 
  stat_qq_line()

mle_lambda_rep <- mle_lambda_rep %>% 
  mutate(se = (maxCI - minCI) / 2 * 1.96,
         var = se**2)

summary(lm(b~site*year, dat = mle_lambda_rep, weights = var))
ggplot(mle_lambda_rep,
       aes(x = year - mean(year),
           y = b,
           ymin = minCI,
           ymax = maxCI,
           color = site)) +
  geom_pointrange() +
  geom_smooth(method = "lm", 
              aes(weight = var))


mle_lambda_rep %>%
  group_by(site, year) %>%
  summarize(mean_b = mean(b)) %>%
  ggplot(aes(x = year,
             y = mean_b,
             color = site)) +
  geom_point() +
  geom_smooth(method = "lm")

mle_lambda_rep %>%
  group_by(site, year) %>%
  summarize(mean_b = mean(b)) %>%
  mutate(year_c = year - mean(year)) %>%
  lm(mean_b ~ site*year_c, data = .) %>%
  summary()


# combine reps ------------------------------------------------------------



# estimate lambda exponent of a power law
# N ~ M^b
# N = number of individuals
# M = individual body size
# b = lambda, exponent describing power law
## This is unknown parameter that we are estimating 
mle_lambda <- dat %>%
  #filter(dw >0.0026) %>%
  group_by(site, year) %>%
  nest() %>%
  mutate(lambda = map(data,
                      MLE_tidy,
                      "dw")) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup()

mle_lambda <- mle_lambda %>% 
  mutate(se = (maxCI - minCI) / 2 * 1.96,
         var = se**2,
         year_0 = year - min(year))

weighted_ols <- lm(b~site*year_0, dat = mle_lambda, weights = var)
summary(weighted_ols)

ggplot(mle_lambda,
       aes(x = year,
           y = b,
           ymin = minCI,
           ymax = maxCI,
           color = site)) +
  geom_pointrange() +
  scale_color_manual(values = c("#019AFF", "#FF1984")) +
  geom_smooth(method = "lm", 
              aes(weight = var),
              alpha = 0.15) +
  labs(y = "Estimated \U03BB",
       x = "Year") +
  theme_bw()
ggsave("plots/mle_weighted_ols_Apr_2024.png",
       units = "in",
       width = 7,
       height = 6)


# Temperature -------------------------------------------------------------

# generated and downloaded from:
# https://wcc.sc.egov.usda.gov/reportGenerator/view/customWaterYearGroupByMonthReport/annual_calendar_year/start_of_period/369:CO:SNTL%7C485:CO:SNTL%7C547:CO:SNTL%7C938:CO:SNTL%7Cid=%22%22%7Cname/1990-01-01,2024-01-01/stationId,name,TAVG::value?fitToScreen=false

# 369	Brumley	CO	SNTL
# 938	Buckskin Joe	CO	SNTL
# 485	Fremont Pass	CO	SNTL
# 547	Ivanhoe	CO	SNTL

snotel <- read_csv("data/SNOTEL_temperature.csv", skip = 58)
names(snotel) <- c("year", "station_id", "station_name", "mean_air_temp")

# usgs gauge 07086000
usgs <- read_tsv("data/USGS_07086000_temperature.csv", skip = 34)
# remove first row and only keep last two columns
usgs <- usgs[-1,5:6]
names(usgs) <- c("year", "mean_temp")
usgs <- usgs %>%
  mutate(mean_temp = as.numeric(mean_temp))

# plots
snotel %>%
  ggplot(aes(x = year, y = mean_air_temp,
             color = station_name)) +
  geom_point() +
  stat_smooth(method = "lm") +
  labs(y = "Annual mean air temperature",
       x = "Year") +
  theme_bw()

snotel %>%
  ggplot(aes(x = year, y = mean_air_temp,
             color = station_name)) +
  geom_point() +
  stat_smooth(method = "lm",
              inherit.aes = FALSE,
              aes(x = year, y = mean_air_temp)) +
  labs(y = "Annual mean air temperature \U00B0 C",
       x = "Year") +
  theme_bw()
ggsave("plots/snotel_temp_Apr_2024.png",
       units = "in",
       width = 7,
       height = 6)

usgs %>%
  ggplot(aes(x = year, 
             y = mean_temp)) +
  geom_point() +
  stat_smooth(method = "lm")+
  labs(y = "Annual mean water temperature \U00B0 C",
       x = "Year") +
  theme_bw()
ggsave("plots/usgs_temp_Apr_2024.png",
       units = "in",
       width = 7,
       height = 6)

# lm stats
summary(lm(mean_air_temp ~ year, data = snotel))
# average air temp across snotel sites increasing by:
# 0.17 degrees C per year
summary(lm(mean_temp ~ year, data = usgs))
# average water temperature at usgs gauge downstream increasing by:
# 0.033 degrees per year. 
