source("R/paretocounts.R")
library(brms)
library(tidyverse)
library(tidybayes)
library(modelr)

dat <- readRDS("ark_dw.RDS")

isd_data = dat %>%
  mutate(xmin = min(dw), # xmin
         xmax = max(dw),
         year = as.factor(year)) %>% #xmax
  group_by(site, year, dw) %>%
  add_tally(name = "counts")

# code for setting priors
get_prior(dw | vreal(counts, xmin, xmax) ~ site * year + (1|rep),
          data = isd_data,
          stanvars = stanvars,
          family = paretocounts())

fit = brm(dw | vreal(counts, xmin, xmax) ~ site * year + (1|rep) + (1|year),
          data = isd_data,
          stanvars = stanvars,
          family = paretocounts(),
          prior =
            c(prior(normal(0, 0.1),
                    class = "b"),
              prior(normal(-1.5, 0.25),
                    class = "Intercept"),
            prior(exponential(4),
                  class="sd")
              ),
          chains = 1,
          iter = 100)

fit

fit2 <- update(fit,
               iter = 4000)

#saveRDS(fit2, "R/model_fits/site-x-year-rand-rep_year-intercept.RDS")

fit2 <- readRDS("R/model_fits/site-x-year-rand-rep_year-intercept.RDS")

fit2
conditional_effects(fit2)

offsets = fit2$data %>%
  distinct(site, year, xmin, xmax, rep) %>%
  mutate(counts = 1,
         year = parse_number(as.character(year))) %>%  # this is just a placeholder, keep it as is
  tidybayes::add_epred_draws(fit2, re_formula = NA) %>% group_by(rep, year, site) %>% median_qi(.epred)

fit2$data %>%
  distinct(site, year, xmin, xmax) %>%
  mutate(counts = 1,
         year = parse_number(as.character(year))) %>%  # this is just a placeholder, keep it as is
  tidybayes::add_epred_draws(fit2, re_formula = NA) %>%
  ggplot(aes(x = year,
             y = .epred,
             color = site)) +
  stat_pointinterval() +
  # geom_pointrange(data = offsets,
  #            aes(y = .epred, ymin = .lower,
  #                ymax = .upper), color = "black",
  #            shape = 21, size = 1) +
  NULL


fit2$data %>%
  distinct(site, year, xmin, xmax) %>%
  mutate(counts = 1,
         year = parse_number(as.character(year))) %>%  # this is just a placeholder, keep it as is
  tidybayes::add_epred_draws(fit2, re_formula = NA) %>%
  ggplot(aes(x = year,
             y = .epred,
             color = site,
             fill = site)) +
  stat_halfeye(alpha = 0.5) +
  brms::theme_default()


plot_data <- fit2$data %>%
  distinct(site, year, xmin, xmax) %>%
  mutate(counts = 1,
         year = parse_number(as.character(year))) %>%  # this is just a placeholder, keep it as is
  tidybayes::add_epred_draws(fit2, re_formula = NA)


plot_data %>%
  ggplot(aes(x = year,
             y = .epred,
             color = site,
             fill = site)) +
  stat_halfeye(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("#019AFF", "#FF1984")) +
  scale_fill_manual(values = c("#019AFF", "#FF1984")) +
  labs(y = expression(lambda))


plot_data %>%
  filter(site == "AR3") %>%
  ggplot(aes(x = year,
             y = .epred,
             color = site,
             fill = site)) +
  stat_halfeye(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c( "deeppink4")) +
  scale_fill_manual(values = c( "#FF1984")) +
  labs(y = expression(lambda))

ggsave("R/figures/ar3.png")

plot_data %>%
  filter(site == "AR1") %>%
  ggplot(aes(x = year,
             y = .epred,
             color = site,
             fill = site)) +
  stat_halfeye(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("darkslategray", "deeppink4")) +
  scale_fill_manual(values = c("#019AFF", "#FF1984")) +
  labs(y = expression(lambda))

ggsave("R/figures/ar1.png")


plot_data %>%
  #filter(site == "AR3") %>%
  ggplot(aes(x = year,
             y = .epred,
             color = site,
             fill = site)) +
  stat_halfeye(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("darkslategray", "deeppink4")) +
  scale_fill_manual(values = c("#019AFF", "#FF1984")) +
  labs(y = expression(lambda))

ggsave("R/figures/ar3_ar1.png")


plot_data %>%
  #filter(site == "AR3") %>%
  ggplot(aes(x = year,
             y = .epred,
             color = site,
             fill = site)) +
  stat_halfeye(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("darkslategray", "deeppink4")) +
  scale_fill_manual(values = c("#019AFF", "#FF1984")) +
  labs(y = expression(lambda)) +
  stat_smooth(method = "lm")

ggsave("R/figures/ar3_ar1_lm-smooth.png")


plot_data %>%
  summarise(med = median(.epred),
            l95 = quantile(.epred, probs = (0.025)),
            u95 = quantile(.epred, probs = (0.975))) %>%
  ungroup() %>%
  select(site, year, med, l95, u95)
