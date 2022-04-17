# appends LW coeffs

# library ####
library(tidyverse)

# read in data made in first script
dat <- readRDS("data/ark_csv_stitched.R")

# read in csv with length weight equations
lw_coef <- read.csv("data/LW_coeffs.csv")

# which taxa don't have lw coefs?
no_coef <- setdiff(unique(dat$Label), unique(lw_coef$taxon))

# percent of data that has lw coeffs?
(nrow(dat[dat$Label %in% no_coef,])/ nrow(dat))*100
# ~ 9% no lw [4/17/22]

# make string vector of column names we care about
lw_cols <- c("taxon",
             "a",
             "b",
             "L_units",
             "dw_units",
             "formula_type",
             "formula")

# check formulas in LW for taxa in arkansas
equations <- lw_coef[lw_coef$taxon %in% unique(dat$Label),]


# I think this is old? 
# commenting out for now
# dat2 <- merge(dat2, lw_coef[,lw_cols], 
#       by.x = "Label",
#       by.y = "taxon", all.x = TRUE)

# join the ark data with lw coef values
fulldat <- merge(dat, lw_coef[,lw_cols], 
                 by.x = "Label",
                 by.y = "taxon", all.x = FALSE)

# dimensions of the ark and merged data
dim(dat)
dim(fulldat)

# different lw formulas
unique(fulldat$formula_type)

# estimate dw based on formula type
fulldat <- fulldat %>% 
  mutate(dw = case_when(
    formula_type == 1 ~ a * Value^b,
    formula_type == 2 ~ exp(a + b * log(Value))))

dim(fulldat)

# remove any dry weight values < or = to 0
dw <- fulldat %>%
  filter(dw >0) %>%
  select(file, dw)

dim(dw)
names(dw)


### need to also pull out year infor
# pull out info from "file" column to site, date, rep (surber sample), etc.
dw <- dw %>%
  separate(file, into = c("site", "rep", "date", "photo"),
           sep="_")
# separate date into day month and year
dw <- dw %>%
  separate(date, into = c("day", "month", "year"), sep="-")

# fix "1900" dates to 1990
dw <- dw %>%
  mutate(year = case_when(year == 1900 ~ 1990,
                          TRUE ~ as.numeric(year)))


# bin function ####

# this is a custom written function which estimates normalized size spectra using log2 width bins
bin_and_center <- function(data, var, breaks, ...){
  
  # data is a data frame
  # var is a string, and is the name of a column in data which you want to bin
  # breaks controls the number of bins as defined in hist(). Can be a vector giving breakpoints of bins [i.e. Log10 breaks = 10^seq(min, max)], a single number designating number of bins, etc. See ?hist for details. If not supplied by the user a default of log2 width bins are calculated 
  
  # are breaks supplied by the user?
  if(exists("breaks") == FALSE){
    
    # calculate log2 width bins which are inclusive of the range of data supplied
    breaks = 2^seq(floor(range(log2(data[[var]]))[1]),
                   ceiling(range(log2(data[[var]]))[2]))
    message("breaks not supplied, using log2 width bins")
  }
  
  # bin values using hist()
  binned_hist = hist(data[[var]], 
                     breaks = breaks,
                     include.lowest = TRUE, plot = FALSE)
  # calculate "left" and "right" edge of bins
  breaks_orig = binned_hist$breaks[1:(length(breaks)-1)]
  breaks_offset = binned_hist$breaks[2:length(breaks)]
  # total bin width = right edge - left edge
  bin_width = breaks_offset - breaks_orig
  count = binned_hist$counts
  log_mids = log10(binned_hist$mids)
  biomass = count * 10**log_mids
  nbiomass = log10(biomass / bin_width)
  dataout = data.frame(
    count = count,
    log_count = log10(count),
    # normalize counts =count / width (White et al 1997)
    log_count_corrected = log10(count / bin_width),
    # original midpoint of bin log10 transformed
    log_mids = log_mids,
    bin_width = bin_width,
    biomass = biomass,
    nbiomass = nbiomass)
  # remove bins with 0 counts
  # -Inf comes from log10(count / break_width) above
  dataout = dataout[dataout$log_count_corrected !=-Inf,]
  # recenter data at x=0
  mid_row = ceiling(nrow(dataout)/2)
  # subtract value of mid row from all mids
  dataout$log_mids_center =
    dataout[,"log_mids"] - dataout[mid_row,"log_mids"]
  dataout
}

# filter out data < or = 0.0026 mg
# See SI from Perkins et al. 2018 Ecology Letters
dw <- dw %>%
  filter(dw>0.0026)

min(dw$dw)

# save data with estimated dry weights
saveRDS(dw, "data/ark_dw.RDS")


# bin and center data #### 

# define break widths
# this code is for log2 width bins
breaks = 2^seq(floor(range(log2(dw$dw))[1]),
               ceiling(range(log2(dw$dw))[2]))

dw_bin <- dw %>%
  group_by(site, year) %>%
  select(dw) %>%
  nest(size_data = dw) %>%
  mutate(bin = map(size_data,
                   bin_and_center,
                   "dw",
                   breaks = breaks)) %>%
  unnest(cols = bin) %>%
  select(-size_data) %>%
  ungroup()

ggplot(dw_bin,
       aes(x = log_mids_center,
           y = log_count_corrected,
           color = as.factor(year), 
           shape = site)) +
  geom_point() +
  geom_smooth(
    #aes(group = interaction(site, year)),
    method = "lm",
    se = FALSE) +
  facet_wrap(.~site)

ggplot(dw_bin,
       aes(x = log_mids,
           y = log_count_corrected,
           color = as.factor(year), 
           shape = site)) +
  geom_point() +
  geom_smooth(
    #aes(group = interaction(site, year)),
    method = "lm",
    se = FALSE) +
  facet_wrap(.~site)

dw_bin %>%
  filter(log_mids_center >-1.50) %>%
ggplot(aes(x = log_mids_center,
           y = log_count_corrected,
           color = as.factor(year), 
           shape = site)) +
  geom_point() +
  geom_smooth(
    #aes(group = interaction(site, year)),
    method = "lm",
    se = FALSE) +
  facet_wrap(.~site)


# preliminary statistics
dw_bin <- dw_bin %>%
  mutate(y_fact = factor(year, 
                            levels = c("1990", "2012", "2019")))

# summary(lm(log_count_corrected ~ log_mids_center + site*y_fact, data = dw_bin))
# 
# summary(lm(log_count_corrected ~ log_mids_center*site, data = dw_bin))
# 
# summary(lm(log_count_corrected ~ log_mids_center+site, data = dw_bin))
# 
# summary(lm(log_count_corrected ~ log_mids_center*y_fact, data = dw_bin))
# 
# summary(lm(log_count_corrected ~ log_mids_center+y_fact, data = dw_bin))
# 
# summary(lm(log_count_corrected ~ log_mids_center, data = dw_bin))


dw_bin %>%
  group_by(site, y_fact) %>%
  mutate(id = cur_group_id()) %>%
  select(id, site, y_fact) %>%
  arrange(id) %>%
  unique()

dw_bin %>%
  group_by(site, year) %>%
  mutate(id = cur_group_id()) %>%
  #filter(log_mids_center >-1.50) %>%
  lm(log_count_corrected ~ log_mids_center:as.factor(id), data = .) %>% summary()

mod2 <- dw_bin %>%
  group_by(site, year) %>%
  mutate(id = cur_group_id()) %>%
  #filter(log_mids_center >-1.50) %>%
  lm(log_count_corrected ~ log_mids_center:as.factor(id), data = .) 

dw_bin %>%
  lm(log_count_corrected ~ log_mids_center+(log_mids_center:site:y_fact) , data = .) %>% summary()

mod1 <- dw_bin %>%
  filter(log_mids_center >-1.50) %>%
  lm(log_count_corrected ~ log_mids_center +
       #site + 
       #y_fact + 
       log_mids_center:site + 
       log_mids_center:y_fact + 
       log_mids_center:site:y_fact +
       NULL,
     data = .) #%>% summary()

# dw_bin %>%
#   filter(log_mids_center >-1.50) %>%
#   lm(log_count_corrected ~ log_mids_center*site, data = .)%>% summary()
# 
# dw_bin %>%
#   filter(log_mids_center >-1.50) %>%
#   lm(log_count_corrected ~ log_mids_center+site, data = .)%>% summary()
# 
# dw_bin %>%
#   filter(log_mids_center >-1.50) %>%
#   lm(log_count_corrected ~ log_mids_center*y_fact, data = .)%>% summary()
# 
# dw_bin %>%
#   filter(log_mids_center >-1.50) %>%
#   lm(log_count_corrected ~ log_mids_center+y_fact, data = .)%>% summary()
# 
# dw_bin %>%
#   filter(log_mids_center >-1.50) %>%
#   lm(log_count_corrected ~ log_mids_center, data = .)%>% summary()
  


# mod1 <- dw_bin %>%
#   filter(log_mids_center >-1.50) %>%
#   lm(log_count_corrected ~ log_mids_center + site*y_fact, data = .) 

# plot model fit
newdata <- dw_bin %>% 
  #filter(log_mids_center >-1.50) %>%
  select(site, y_fact, log_mids_center,
         log_count_corrected)

newdata <- cbind(newdata,
                 predict(mod1,
                         newdata,
                         interval = "predict"))

# figure with "predicted" data
dw_bin %>% 
  filter(log_mids_center >-1.50) %>%
  ggplot(
       aes(y = log_count_corrected,
           x = log_mids_center,
           color = y_fact)) + 
  labs(x = expression(Log[10]~M),
       y = expression(Log[10]~N)) +
  geom_point(size = 0.7) +
  geom_line(data = newdata,
            aes(y=fit,
                x = log_mids_center,
                group = interaction(site, y_fact),
                color = y_fact), 
            size = 0.5) +
  geom_ribbon(data = newdata,
              aes(ymin = lwr,
                  ymax = upr,
                  group = interaction(site, y_fact),
                  fill = y_fact),
              alpha = 0.1) +
  theme_classic() +
  theme(
        strip.text.x = element_text(size = 8, margin = margin(.09, 0, .09, 0, "cm"))) +
  facet_wrap(.~site)


newdata2 <- dw_bin %>% 
  group_by(site, y_fact) %>%
  mutate(id = cur_group_id(),
         id = factor(id, levels = c(1, 2, 3, 4, 5, 6))) %>%
  select(site, y_fact, log_mids_center,
         log_count_corrected, id) %>%
  ungroup()

newdata2 <- cbind(newdata,
                 predict(mod2,
                         newdata,
                         interval = "predict"))
dw_bin %>% 
  ggplot(
    aes(y = log_count_corrected,
        x = log_mids_center,
        color = y_fact)) + 
  labs(x = expression(Log[10]~M),
       y = expression(Log[10]~N)) +
  geom_point(size = 0.7) +
  geom_line(data = newdata2,
            aes(y=fit,
                x = log_mids_center,
                group = interaction(site, y_fact),
                color = y_fact), 
            size = 0.5) +
  geom_ribbon(data = newdata2,
              aes(ymin = lwr,
                  ymax = upr,
                  group = interaction(site, y_fact),
                  fill = y_fact),
              alpha = 0.1) +
  theme_classic() +
  theme(
    strip.text.x = element_text(size = 8, margin = margin(.09, 0, .09, 0, "cm"))) +
  facet_wrap(.~site)
