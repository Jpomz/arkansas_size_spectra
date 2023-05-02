# mle tidy

# `sizeSpectra` is not hosted on CRAN
# to install the package directly from github, you need to have the `devtools` package installed. Run the following once if you don't already have it downloaded
# install.packages("devtools") 

#To install the latest version of sizeSpectra, run the following:
# devtools::install_github("andrew-edwards/sizeSpectra")

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
                             "2021")))

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

# plot of the count of individuals by sample
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

# estimate lambda exponent of a power law
# N ~ M^b
# N = number of individuals
# M = individual body size
# b = lambda, exponent describing power law
## This is unknown parameter that we are estimating 
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

# view estimated lambda values and 95% CI
mle_lambda %>%
  arrange(site, y_fact)

# plot point range of lambda values
# years on correct scale
mle_lambda %>%
  ggplot(aes(y = b,
             ymin = minCI,
             ymax = maxCI,
             x = year(as.Date(y_fact, 
                              format = "%Y")),
             color = site)) +
  geom_pointrange(
    size = 1,
    position = position_dodge(
      width = 0.25
    )) +
  theme_bw() +
  labs(y = expression(lambda),
       x = "year")
# plot point range of lambda values
# years as factor
(lambda_year_plot <- mle_lambda %>%
  ggplot(aes(y = b,
             ymin = minCI,
             ymax = maxCI,
             x = y_fact,
             color = site)) +
  geom_pointrange(
    size = 0.5,
    position = position_dodge(
      width = 0.25
    )) +
  theme_bw() +
    scale_color_viridis_d(option = "plasma", begin = 0.25, end = 0.75) +
  labs(y = expression(lambda),
       x = "year"))
ggsave(lambda_year_plot, 
       file = "plots/lambda_year.png",
       units = "in",
       width = 6,
       height = 6)


# adding line for visualization
# I think I like pointrange above better
# function from a different ggplot add-on??
ggline(mle_lambda, 
       x = "y_fact",
       y = "b",
       color = "site", 
       size = 1)

# plot isds
# join dw and MLE estimates
dat_plot <- dat %>%
  left_join(mle_lambda) %>%
  group_by(site, y_fact) %>% 
  mutate(sample_int = cur_group_id())

# select data and calculate rank order for dw observations
dat_plot = dat_plot %>%
  ungroup() %>%
  # global xmin
  mutate(xmin = min(dw, na.rm = TRUE)) %>%
  group_by(sample_int) %>%
  select(site, year, sample_int,
         dw, b, minCI, maxCI, xmin) %>%
  arrange(site, year, desc(dw)) %>% 
  # order of dw, and local xmax
  mutate(y_order = 1:n(),
         xmax = max(dw, na.rm = TRUE)) 

# split data by site and year
dat_split = dat_plot %>% 
  group_by(sample_int) %>% 
  group_split

nsamples = 1000

xy.PLB = NULL
for(i in 1:length(dat_split)) {
  sample_int = unique(dat_split[[i]]$sample_int)
  site = unique(dat_split[[i]]$site)
  year = unique(dat_split[[i]]$year)
  xmin = min(dat_split[[i]]$xmin)
  xmax = max(dat_split[[i]]$xmax)
  
  lambda = unique(dat_split[[i]]$b)
  .lower = unique(dat_split[[i]]$minCI)
  .upper = unique(dat_split[[i]]$maxCI)
  
  x.PLB = seq(min(dat_split[[i]]$dw),
              max(dat_split[[i]]$dw),
              length=nsamples) # x values to plot PLB
  
  y.PLB = (1 - (x.PLB^(lambda + 1) - (xmin^(lambda+1)))/(xmax^(lambda + 1) - (xmin^(lambda+1))))*nsamples
  ymin.PLB = (1 - (x.PLB^(.lower + 1) - (xmin^(.lower+1)))/(xmax^(.lower + 1) - (xmin^(.lower+1))))*nsamples
  ymax.PLB = (1 - (x.PLB^(.upper + 1) - (xmin^(.upper+1)))/(xmax^(.upper + 1) - (xmin^(.upper+1))))*nsamples
  
  xy.PLB[[i]] = tibble(dw = x.PLB, y_order = y.PLB,
                       ymin = ymin.PLB,
                       ymax = ymax.PLB,
                       xmax = xmax,
                       xmin = xmin) %>%
    mutate(sample_int = sample_int,
           site = site,
           year = year,
           lambda = lambda)
}

lines_toplot = bind_rows(xy.PLB)  %>% 
  # filter(sample_int <= 10) %>% 
   mutate(facet_name = paste(site, year)) # %>% 
  # left_join(dat %>% ungroup %>% distinct(site_id, mean, log_gpp_s, log_om_s))

# plots?
(isd_per_sample_plot = dat_plot %>% 
  mutate(facet_name = paste(site, sample_int)) %>% 
  ggplot(aes(x = dw, y = y_order, group = sample_int)) + 
  geom_point(shape = 21, size = 0.3, aes(color = site)) +
  geom_line(data = lines_toplot, aes(color = site)) +
  geom_ribbon(data = lines_toplot , aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~year) +
  theme_bw() +
    scale_color_viridis_d(option = "plasma", begin = 0.25, end = 0.75) +
  theme(
  #  strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text = element_blank()) +
  labs(y = "Number of values \u2265 x",
       x = "Individual dry mass (mg)",
       color = "Site") +
  #guides(color = "none") +
  coord_cartesian(ylim = c(limits = c(min(dat_plot$y_order), NA))))


# steepest and shallowest lines
dat_2 <- dat_plot %>%
  filter(site == "AR3",
         year == 1990 | year == 2021)
lines_2 <- lines_toplot%>%
  filter(site == "AR3",
         year == 1990 | year == 2021)
  
  
(two_plots <- dat_2 %>% 
  mutate(facet_name = paste(site, sample_int)) %>% 
  ggplot(aes(x = dw, y = y_order, group = sample_int)) + 
  geom_point(shape = 21,
             size = 1,
             alpha = 0.5,
             aes(color = as.factor(round(b, 2)))) +
  geom_line(data = lines_2, 
            aes(color = as.factor(round(lambda, 2)))) +
  geom_ribbon(data = lines_2 , 
              aes(ymin = ymin,
                  ymax = ymax), alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  #facet_wrap(~year) +
  theme_bw() +
  theme(
    #  strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.text = element_blank()) +
  labs(y = "Number of values \u2265 x",
       x = "Individual dry mass (mg)",
       color = "Lambda") +
  guides(fill = "none") +scale_color_viridis_d(option = "plasma", begin = 0.25, end = 0.75) +
  coord_cartesian(ylim = c(limits = c(min(dat_2$y_order), NA))))


ggsave(two_plots, 
       file = "plots/steep_shallow.png",
       units = "in",
       width = 6,
       height = 6)



# Below is just playing around and seeing how much results change


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

# 
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
            sd_b = sd(b)) %>%
  ggplot(aes(x = y_fact,
             y = mean_b,
             ymin = mean_b - sd_b,
             ymax = mean_b + sd_b)) +
  geom_pointrange()

dat %>%
  filter(dw >0.0026) %>%
  group_by(site, y_fact, rep) %>%
  nest() %>%
  mutate(lambda = map(data,
                      MLE_tidy,
                      "dw")) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ggplot(aes(x = y_fact,
             y = b,
             ymin = minCI,
             ymax = maxCI,
             color = site)) +
  geom_pointrange(
    position = position_dodge(
      width = 0.25
    )
  ) 
