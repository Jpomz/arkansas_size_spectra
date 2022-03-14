# appends LW coeffs

library(tidyverse)

# read in data made in first script
dat <- readRDS("data/ark_csv_stitched.R")

# read in csv with length weight equations
lw_coef <- read.csv("data/LW_coeffs.csv")

# which taxa don't have lw coefs?
no_coef <- setdiff(unique(dat$Label), unique(lw_coef$taxon))

# percent of data that has lw coeffs?
nrow(dat[dat$Label %in% no_coef,])/ nrow(dat)

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
equations %>%
  select(lw_cols) %>% View


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

# pull out info from "file" column to site, date, rep (surber sample), etc.
dw <- dw %>%
  separate(file, into = "site", sep = 6, remove = TRUE)

dw <- dw %>%
  separate(site, into = c("site","rep"), sep = "_", remove = FALSE)

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

# define break widths
# this code is for log2 width bins
breaks = 2^seq(floor(range(log2(dw$dw))[1]),
               ceiling(range(log2(dw$dw))[2]))

min(dw$dw)

# save data with estimated dry weights
saveRDS(dw, "data/ark_dw.RDS")

# split data into 2 sites before estimating size spectra
# need to fix function to do this automatically
dw_ar1 <- dw %>%
  filter(site == "AR1")

dw_ar3 <- dw %>%
  filter(site == "AR3")

# estimate spectra for each site independently
bin_ar1 <- bin_and_center(dw_ar1, "dw", breaks = breaks)
bin_ar3 <- bin_and_center(dw_ar3, "dw", breaks = breaks)

bin_ar1$site <- "AR1"
bin_ar3$site <- "AR3"

# put two sites back together
bin <- bind_rows(bin_ar1, bin_ar3)

# plot size spectra results
ggplot(bin, aes(x = log_mids, y = log_count_corrected, color = site))+
  geom_point() +
  stat_smooth(method = "lm")

ggplot(bin, aes(x = log_mids_center, y = log_count_corrected, color = site))+
  geom_point() +
  stat_smooth(method = "lm")

# prelimanary statistics
spectra <- lm(log_count_corrected ~ log_mids_center * site, data = bin)
summary(spectra)

summary(lm(log_count_corrected ~ log_mids_center, data = bin_ar1))
summary(lm(log_count_corrected ~ log_mids_center, data = bin_ar3))


summary(lm(log_count_corrected ~ log_mids*site, data = bin))
summary(lm(log_count_corrected ~ log_mids+site, data = bin))


summary(lm(log_count_corrected ~ log_mids*site + I(log_mids^2), data = bin))

        
quad_mod <- lm(log_count_corrected ~ log_mids*site + I(log_mids^2), data = bin)


# code below this is preliminary and exploratory, I wouldn't worry about it too much for now. 

newdata <- bin %>% 
  select(site, log_mids_center,
         log_count_corrected, log_mids)

newdata <- cbind(newdata,
                 predict(quad_mod,
                         newdata,
                         interval = "predict"))

# figure with "predicted" data
ggplot(bin, aes(y = log_count_corrected, x = log_mids)) + 
  labs(x = expression(Log[10]~M), y = expression(Log[10]~N)) +
  geom_point(size = 0.7) +
  geom_line(data = newdata, aes(y=fit, x = log_mids, group = site), 
            size = 0.5) +
  geom_ribbon(data = newdata, aes(ymin = lwr, ymax = upr, group = site), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 8, margin = margin(.09, 0, .09, 0, "cm"))) 





rep.dw <- dw %>%
  group_by(site, rep) %>%
  # create list-column
  nest() %>% 
  # estimate b and 95% CI
  mutate(bin = map(data, bin_and_center, var = "dw", breaks = breaks)) %>%
  unnest(cols = bin)

ggplot(rep.dw, 
       aes(x = log_mids, y = log_count_corrected, color = site))+
  geom_point(aes(shape = rep)) +
  stat_smooth(method = "lm", aes())

summary(lm(log_count_corrected ~ log_mids_center *site, data = rep.dw))

summary(lm(log_count_corrected ~ log_mids_center *site +I(log_mids_center^2), data = rep.dw))

ggplot(rep.dw, 
       aes(x = log_mids, y = log_count_corrected, color = site))+
  geom_point(aes(shape = rep)) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2))


fulldat %>%
  separate(file, into = "site", sep = 3, remove = TRUE) %>%
  group_by(site, Label) %>%
  summarise(sp.dw = mean(dw),
            n.sp = n()) %>%
  ggplot(aes(x = sp.dw, y = n.sp, color = site)) +
  geom_point()+
  scale_x_log10() +
  scale_y_log10() +
  stat_smooth(method = "lm")



breaks_log_6 <- exp(seq(floor(log(min(dw$dw))),
                        ceiling(log(max(dw$dw))),
                        length.out = 7))

dw %>%
  group_by(site, rep) %>%
  # create list-column
  nest() %>% 
  # estimate b and 95% CI
  mutate(bin = map(data, bin_and_center, var = "dw", breaks = breaks_log_6)) %>%
  unnest(cols = bin) %>%
  ggplot(aes(x = log_mids_center, y = log_count_corrected, color = site))+
  geom_point(aes(shape = rep)) +
  stat_smooth(method = "lm")
