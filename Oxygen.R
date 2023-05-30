###########################################
#### Project: Detrital thermal optimum ####
#### Purpose: Analysis of oxygen data  ####
#### Author: Luka Seamus Wright        ####
###########################################

#### Load data and packages ####
### Data ####
oxy <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Oxygen.csv")
str(oxy) # a plethora of  measurements taken by the oxygen metre, some of which
# are purely technical and not of interest in analysis
# The important variables are:
# ID = sample ID
# Type = "Blank" or "Sample"
# Group = 23 measurement groups
# PAR = 0 or 420 µmol photons m^-2 s^-1, which was derived empirically from
PAR <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/PAR.csv")
mean(PAR$PAR) # mean = 419.925 µmol photons m^-2 s^-1
sd(PAR$PAR) # sd = 18.85625 µmol photons m^-2 s^-1
# (measurement error in PAR will be incorporated later)
# delta_t = min
# Value = dissolved oxygen (µM)
# Temp = incubation temperature (°C)
mean(oxy$Temp) # mean = 17.59242°C
sd(oxy$Temp) # sd = 0.3823833°C
# Pressure = incubation pressure (hPa)
mean(oxy$Pressure) # mean = 1021.857 hPa
sd(oxy$Pressure) # sd = 4.274041 hPa
# Salinity = incubation salinity (‰)
mean(oxy$Salinity) # mean = 34.8625‰
sd(oxy$Salinity) # sd = 0.6846104‰
# these confounders will be incorporated later

### Packages ####
require(tidyverse) # this is my preferred data processing and visualisation framework 
require(magrittr) # this extends the pipe operator (%>%) of the tidyverse
require(rethinking) # this provides my Bayesian data analysis framework of choice

#### Test bulk modelling pipeline ####
# I first use simple frequentist linear models to see if the pipeline works
# to build 160 models (160 measurement series) simultaneously

# models need to be stratified by measurment, a combination of ID and PAR
oxy %<>% mutate(IDm = paste(ID, PAR, sep = ".")) # create variable measurement ID

# make trimmed version of dataframe with relevant variables only
oxy.trim <- oxy %$% data.frame(IDm = IDm, # measurement ID
                               IDs = ID, # sample ID
                               Group = Group, # measurement group
                               PAR = PAR, # µmol photons m^-2 s^-1
                               Temp = Temp, # °C
                               t = delta_t, # min
                               O2 = Value) # µM

# feed nested data into frequentist linear model (note that Temp cannot
# be included as a factor in the nesting because it varies within IDm)
fmod <- oxy.trim %>% nest_by(IDm, IDs, Group, PAR) %>%
        mutate(mod = list(lm(O2 ~ t, data = data)))

# add Temp to model meta-dataframe
fmod <- bind_cols(fmod, oxy %>% group_by(IDm) %>%
                  reframe(Tempsd = sd(Temp),
                          Temp = mean(Temp)) %>%
                  select(Temp, Tempsd))

# use broom package to extract coefficients
tidy.fmod <- fmod %>% reframe(broom::tidy(mod))
str(tidy.fmod) # Temp was knocked out again

# add Temp again
tidy.fmod <- bind_cols(tidy.fmod, oxy %>% group_by(IDm) %>%
                       reframe(Tempsd = sd(Temp),
                               Temp = mean(Temp)) %>%
                       select(Temp, Tempsd) %>%
                       slice(rep(1:n(), each = 2)))


#### Actual modelling ####
### Model structure ####

# O2 ~ dnorm(mean = µ, sd = σ)
# µ <- ɑ + β * t
# ɑ ~ dnorm(mean, sd)
# β ~ dnorm(mean, sd)
# σ ~ dexp(rate)

### Prior means ####
## Intercept ####
# before scrutinising the data, there is prior information to consider
# Rose et al. 2012 (doi: 10.5194/os-8-545-2012) provide annual baseline dissolved 
# oxygen levels of ~7 mg l^-1 for several coastal sites close to our study site 
# off Perth on the west coast of Australia
g.mol <- 15.9994 # g mol^-1 for oxygen
a.mu <- 7 / g.mol * 1e3
a.mu # 437.5164 µM

# this estimate seems extremely high and needs to be reconsidered in light of 
# other studies. Woo & Pattiaratchi 2008 (doi: 10.1016/j.dsr.2008.05.005)
# provide regional surface water dissolved oxygen estimates of 220-245 µM that
# do not need to be converted and are therefore more reliable
mean(c(220, 245)) # 232.5 µM

## Slopes ####
# here, three priors are required, one for blanks (which are expected to have a neutral
# slope), one for samples incubated in the dark (which are expected to have
# a negative slope) and one for samples incubated in the light (which are
# expected to have a positive slope)

# Staehr & Wernberg 2009 (doi: 10.1111/j.1529-8817.2008.00635.x) and 
# Wernberg et al. 2016 (doi: 10.1002/lno.10362) performed a similar experiment on the 
# same species (Ecklonia radiata), at the same site (Marmion Lagoon) and with similar 
# sized tissue (laterals of ~0.22 g dry mass)
# Priors for b can therefore be estimated by back-transforming their oxygen 
# evolution slopes (mg g^-1 h^-1) under light and dark conditions:

# get prior Ecklonia photosynthesis data
prior <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Prior.csv")
m <- 0.22 # g
V <- 0.12 # l
t <- 60 # min

# calculate light prior mean
b.mu.l <- prior %>% filter(Temperature >= 15, Temperature <= 20) %>% 
                           # these are closest to my experimental temperature
                    pull(NP) %>% mean() * m / V / t / g.mol * 1e3
b.mu.l # 3.056483 µM min^-1



# calculate dark prior mean
b.mu.d <- prior %>% filter(Temperature >= 15, Temperature <= 20) %>%
                    pull(R) %>% mean() * m / V / t / g.mol * 1e3 * -1
b.mu.d # -1.371941 µM min^-1

# blank prior mean is just assumed to be 0

### Prior standard deviations ####
## Intercept ####
# this process is iterative, so I choose a reasonable initial 
# standard deviation and test its implications
a.prior <- rnorm(n = 1e5 , mean = 232.5 , sd = 40)

ggplot() + 
  geom_density(aes(x = a.prior)) +
  theme_minimal()

## Slopes ####
# light prior
b.prior.l <- rnorm(n = 1e5, mean = 3.056483, sd = 4)

ggplot() + 
  geom_density(aes(x = b.prior.l)) +
  theme_minimal()

# dark prior
b.prior.d <- rnorm(n = 1e5, mean = -1.371941, sd = 2) 
# sd chosen to be lower for dark samples because there is higher confidence in
# maintenance of detrital respiration in dark than there is of detrital 
# photosynthesis in light, plus sd usually scales with mean which is smaller

ggplot() + 
  geom_density(aes(x = b.prior.d)) +
  theme_minimal()

# blank prior
b.prior.b <- rnorm(n = 1e5, mean = 0, sd = 1)
# there is reasonable confidence that the blank slope is in the vicinity of 0

ggplot() + 
  geom_density(aes(x = b.prior.b)) +
  theme_minimal()

## Cumulative effect of intercept and slopes ####
# light prior
predictor <- oxy.trim %$% seq(min(t), max(t), length.out = 100)

plot(NULL, xlim = c(0, 5), ylim = c(120, 410), # base plot is better with loops
     xlab = "min", ylab = "µM") # empty plot
oxy.trim %$% abline(h = min(O2), lty = 2) # data range
oxy.trim %$% abline(h = max(O2), lty = 2)
for(i in 1:1e3) curve(a.prior[i] + b.prior.l[i] * x, # loop
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = col.alpha("black", 0.1))
# intercept and slope vary too much


# dark prior
plot(NULL, xlim = c(0, 5), ylim = c(120, 410),
     xlab = "min", ylab = "µM")
oxy.trim %$% abline(h = min(O2), lty = 2)
oxy.trim %$% abline(h = max(O2), lty = 2)
for(i in 1:1e3) curve(a.prior[i] + b.prior.d[i] * x, 
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = col.alpha("black", 0.1))
# intercept and slope vary too much


# blank prior
plot(NULL, xlim = c(0, 5), ylim = c(120, 410),
     xlab = "min", ylab = "µM")
oxy.trim %$% abline(h = min(O2), lty = 2)
oxy.trim %$% abline(h = max(O2), lty = 2)
for(i in 1:1e3) curve(a.prior[i] + b.prior.b[i] * x, 
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = col.alpha("black", 0.1))
# intercept and slope vary too much


# narrower priors required
a.prior <- rnorm(n = 1e5 , mean = 232.5 , sd = 30)

ggplot() + 
  geom_density(aes(x = a.prior)) +
  geom_vline(aes(xintercept = c(232.5 - 2 * 30,  
                                232.5 + 2 * 30))) +
  theme_minimal()
# 95% of prior probability density still lies within 172.5 and 292.5 µM

b.prior.l <- rnorm(n = 1e5, mean = 3.056483, sd = 2)

ggplot() + 
  geom_density(aes(x = b.prior.l)) +
  geom_vline(aes(xintercept = c(3.056483 - 2 * 2,  
                                3.056483 + 2 * 2))) +
  theme_minimal()
# 95% of prior probability density still lies between -0.943517 and 7.056483 µM

b.prior.d <- rnorm(n = 1e5, mean = -1.371941, sd = 1)

ggplot() + 
  geom_density(aes(x = b.prior.d)) +
  geom_vline(aes(xintercept = c(-1.371941 - 2 * 1,  
                                -1.371941 + 2 * 1))) +
  theme_minimal()
# 95% of prior probability density still lies between -3.371941 and 0.628059 µM

b.prior.b <- rnorm(n = 1e5, mean = 0, sd = 0.5)

ggplot() + 
  geom_density(aes(x = b.prior.b)) +
  geom_vline(aes(xintercept = c(0 - 2 * 0.5,  
                                0 + 2 * 0.5))) +
  theme_minimal()
# 95% of prior probability density still lies between -1 and 1 µM

## Alternative test of cumulative effect ####
# trim dataframe down to only the variables used by the model
oxy.trim <- oxy %$% data.frame(t = delta_t, # min
                               O2 = Value) # µM

# the following models are meaningless because they are fit to the entire data,
# but I am only using it to get at the priors

# light prior
prior.test.l <- ulam( # ulam is a Hamiltonian Monte Carlo engine
                alist(
                  O2 ~ dnorm(mean = mu, sd = sigma),
                  mu <- a + b * t,
                  a ~ dnorm(mean = 232.5, sd = 30),
                  b ~ dnorm(mean = 3.056483, sd = 2),
                  sigma ~ dexp(rate = 1)
                ),
                data = oxy.trim, chains = 1, iter = 1e3
                )

ab.prior.l <- extract.prior(prior.test.l)

ab.mu.l <- link(prior.test.l, post = ab.prior.l,
                data = data.frame(t = predictor))

plot(NULL, xlim = c(0, 5), ylim = c(120, 410),
     xlab = "min", ylab = "µM")
oxy.trim %$% abline(h = min(O2), lty = 2)
oxy.trim %$% abline(h = max(O2), lty = 2)
for(i in 1:1e3) lines(predictor, ab.mu.l[i,],
                      col = col.alpha("black", 0.1))
# looks more reasonable

# dark prior
prior.test.d <- ulam(
                alist(
                  O2 ~ dnorm(mean = mu, sd = sigma),
                  mu <- a + b * t,
                  a ~ dnorm(mean = 232.5, sd = 30),
                  b ~ dnorm(mean = -1.371941, sd = 1),
                  sigma ~ dexp(rate = 1)
                ),
                data = oxy.trim, chains = 1, iter = 1e3
                )

ab.prior.d <- extract.prior(prior.test.d)

ab.mu.d <- link(prior.test.d, post = ab.prior.d,
                data = data.frame(t = predictor))

plot(NULL, xlim = c(0, 5), ylim = c(120, 410),
     xlab = "min", ylab = "µM")
oxy.trim %$% abline(h = min(O2), lty = 2)
oxy.trim %$% abline(h = max(O2), lty = 2)
for(i in 1:1e3) lines(predictor, ab.mu.d[i,],
                      col = col.alpha("black", 0.1))
# looks more reasonable


# blank prior
prior.test.b <- ulam(
                alist(
                  O2 ~ dnorm(mean = mu, sd = sigma),
                  mu <- a + b * t,
                  a ~ dnorm(mean = 232.5, sd = 30),
                  b ~ dnorm(mean = 0, sd = 0.5),
                  sigma ~ dexp(rate = 1)
                ),
                data = oxy.trim, chains = 1, iter = 1e3
                )

ab.prior.b <- extract.prior(prior.test.b)

ab.mu.b <- link(prior.test.b, post = ab.prior.b,
                data = data.frame(t = predictor))

plot(NULL, xlim = c(0, 5), ylim = c(120, 410),
     xlab = "min", ylab = "µM")
oxy.trim %$% abline(h = min(O2), lty = 2)
oxy.trim %$% abline(h = max(O2), lty = 2)
for(i in 1:1e3) lines(predictor, ab.mu.b[i,],
                      col = col.alpha("black", 0.1))
# looks more reasonable


### Diagnostic model ####
# the trimmed data need the variables Type, PAR and IDm to be startified by
oxy.trim <- oxy %$% data.frame(IDm = IDm, # measurement ID
                               IDs = ID, # sample ID
                               Group = Group, # measurment group
                               Type = Type, # "Blank" or "Sample"
                               PAR = PAR, # µmol photons m^-2 s^-1
                               t = delta_t, # min
                               O2 = Value) # µM

# run diagnostic model with one chain
test.mod <- oxy.trim %>% nest_by(IDm, IDs, Group, Type, PAR) %>%
            mutate(
            mod = ifelse(Type == "Blank",
                         list(
                           ulam(
                             alist(
                               O2 ~ dnorm(mean = mu, sd = sigma),
                               mu <- a + b * t,
                               a ~ dnorm(mean = 232.5, sd = 30),
                               b ~ dnorm(mean = 0, sd = 0.5),
                               sigma ~ dexp(rate = 1)
                             ),
                             data = data, chains = 1, iter = 1e3
                           )),
                         ifelse(Type == "Sample" & PAR == 0,
                                list(
                                  ulam(
                                    alist(
                                      O2 ~ dnorm(mean = mu, sd = sigma),
                                      mu <- a + b * t,
                                      a ~ dnorm(mean = 232.5, sd = 30),
                                      b ~ dnorm(mean = -1.371941, sd = 1),
                                      sigma ~ dexp(rate = 1)
                                    ),
                                    data = data, chains = 1, iter = 1e3
                                  )),
                                list(
                                  ulam(
                                    alist(
                                      O2 ~ dnorm(mean = mu, sd = sigma),
                                      mu <- a + b * t,
                                      a ~ dnorm(mean = 232.5, sd = 30),
                                      b ~ dnorm(mean = 3.056483, sd = 2),
                                      sigma ~ dexp(rate = 1)
                                    ),
                                    data = data, chains = 1, iter = 1e3
                                  ))
                                )
                         )
            )

lapply(test.mod$mod, precis, depth = 2) # quite low n_eff and some high Rhat4 but no warnings
lapply(test.mod$mod, pairs) # very strong negative correlation between a and b
# this is likely the reason for low n_eff and will be dealt with by centring
lapply(test.mod$mod, traceplot) # chains are healthy


# run diagnostic models with one chain on centred data
oxy.trim %<>% mutate(tc = t - mean(t)) # min (centred)

test.mod <- oxy.trim %>% nest_by(IDm, IDs, Group, Type, PAR) %>%
            mutate(
            mod = ifelse(Type == "Blank",
                         list(
                           ulam(
                             alist(
                               O2 ~ dnorm(mean = mu, sd = sigma),
                               mu <- a + b * tc,
                               a ~ dnorm(mean = 232.5, sd = 30),
                               b ~ dnorm(mean = 0, sd = 0.5),
                               sigma ~ dexp(rate = 1)
                             ),
                             data = data, chains = 1, iter = 1e3
                           )),
                         ifelse(Type == "Sample" & PAR == 0,
                                list(
                                  ulam(
                                    alist(
                                      O2 ~ dnorm(mean = mu, sd = sigma),
                                      mu <- a + b * tc,
                                      a ~ dnorm(mean = 232.5, sd = 30),
                                      b ~ dnorm(mean = -1.371941, sd = 1),
                                      sigma ~ dexp(rate = 1)
                                    ),
                                    data = data, chains = 1, iter = 1e3
                                  )),
                                list(
                                  ulam(
                                    alist(
                                      O2 ~ dnorm(mean = mu, sd = sigma),
                                      mu <- a + b * tc,
                                      a ~ dnorm(mean = 232.5, sd = 30),
                                      b ~ dnorm(mean = 3.056483, sd = 2),
                                      sigma ~ dexp(rate = 1)
                                    ),
                                    data = data, chains = 1, iter = 1e3
                                  ))
                                )
                         )
            )

lapply(test.mod$mod, precis, depth = 2) # better (~double) n_eff and Rhat4
lapply(test.mod$mod, pairs) # correlation between a and b is gone
lapply(test.mod$mod, traceplot) # chains are healthy


### Final model ####
# run final model with 8 parallel chains spread over the 8 cores of my computer
# and 1e4 iterations

oxy.mod <- oxy.trim %>% nest_by(IDm, IDs, Group, Type, PAR) %>%
            mutate(
            mod = ifelse(Type == "Blank",
                         list(
                           ulam(
                             alist(
                               O2 ~ dnorm(mean = mu, sd = sigma),
                               mu <- a + b * tc,
                               a ~ dnorm(mean = 232.5, sd = 30),
                               b ~ dnorm(mean = 0, sd = 0.5),
                               sigma ~ dexp(rate = 1)
                             ),
                             data = data, chains = 8, cores = 8, iter = 1e4
                           )),
                         ifelse(Type == "Sample" & PAR == 0,
                                list(
                                  ulam(
                                    alist(
                                      O2 ~ dnorm(mean = mu, sd = sigma),
                                      mu <- a + b * tc,
                                      a ~ dnorm(mean = 232.5, sd = 30),
                                      b ~ dnorm(mean = -1.371941, sd = 1),
                                      sigma ~ dexp(rate = 1)
                                    ),
                                    data = data, chains = 8, cores = 8, iter = 1e4
                                  )),
                                list(
                                  ulam(
                                    alist(
                                      O2 ~ dnorm(mean = mu, sd = sigma),
                                      mu <- a + b * tc,
                                      a ~ dnorm(mean = 232.5, sd = 30),
                                      b ~ dnorm(mean = 3.056483, sd = 2),
                                      sigma ~ dexp(rate = 1)
                                    ),
                                    data = data, chains = 8, cores = 8, iter = 1e4
                                  ))
                                )
                         )
            )


lapply(oxy.mod$mod, precis, depth = 2) # n_eff and Rhat4 are good
lapply(oxy.mod$mod, pairs) # no correlation between variables
lapply(oxy.mod$mod, traceplot) 
lapply(oxy.mod$mod, trankplot) # chains are healthy


# Prior-posterior comparison
oxy.mod.prior <- lapply(oxy.mod$mod, extract.prior)
oxy.mod.posterior <- lapply(oxy.mod$mod, extract.samples)
oxy.mod %<>% rownames_to_column("model") %>%
  mutate(model = as.integer(model)) # turn row names into variable model
oxy.mod %<>% mutate(Timepoint = as.integer(str_sub(IDs, -1))) %>%
  group_by(Group) %>%
  mutate(Position = match(IDs, unique(IDs)))

oxy.mod.prior.df <- tibble(sample = oxy.mod.prior) %>%
  mutate(model = row_number()) %>%
  unnest_longer(sample, indices_to = "coefficient") %>%
  unnest_longer(sample, indices_to = "n") %>%
  left_join(oxy.mod %>% select(-c(data, mod)), by = "model")

oxy.mod.posterior.df <- tibble(sample = oxy.mod.posterior) %>%
  mutate(model = row_number()) %>%
  unnest_longer(sample, indices_to = "coefficient") %>%
  unnest_longer(sample, indices_to = "n") %>%
  left_join(oxy.mod %>% select(-c(data, mod)), by = "model")

oxy.mod.prior.posterior.df <- data.frame(rbind(oxy.mod.prior.df, oxy.mod.posterior.df),
                                         distribution = c(rep("prior", 480000),
                                                          rep("posterior", 4800000)))



require(ggh4x) # package that enables nested facets among other things
ggplot(data = oxy.mod.prior.posterior.df %>% filter(coefficient == "b")) +
  geom_density(aes(x = sample, colour = distribution)) +
  # geom_text(aes(-4, 1, label = IDs)) +
  facet_nested(Timepoint + Group ~ Type + Position + PAR, scales = "free",
               independent = "all", nest_line = element_line()) +
  theme_minimal() + 
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")


# Posterior prediction intervals
predictorc <- oxy.trim %$% seq(min(tc), max(tc), length.out = 100) # centred predictor

oxy.mod.mu <- lapply(oxy.mod$mod, link, data = data.frame(tc = predictorc))
oxy.mod.mu.mean <- lapply(oxy.mod.mu, function(x) apply(x, 2, mean))
oxy.mod.mu.PI <- lapply(oxy.mod.mu, function(x) apply(x, 2, PI, prob = 0.98))

oxy.mod.sigma <- lapply(oxy.mod$mod, sim, data = data.frame(tc = predictorc))
oxy.mod.sigma.PI <- lapply(oxy.mod.sigma, function(x) apply(x, 2, PI, prob = 0.98))

oxy.mod.mu.mean.df <- tibble(mean = oxy.mod.mu.mean) %>%
  mutate(model = row_number()) %>%
  unnest_longer(mean, indices_to = "trow") %>%
  left_join(oxy.mod %>% select(-c(data, mod)), by = "model") %>%
  mutate(tc = rep(predictorc, 160), t = rep(predictor, 160))



oxy.mod.mu.PI.df <- tibble(sample = oxy.mod.mu.PI) %>%
  mutate(model = row_number())

oxy.mod.mu.PI.df$sample <- lapply(oxy.mod.mu.PI.df$sample, t)

oxy.mod.mu.PI.df <- oxy.mod.mu.PI.df %>% 
  unnest_longer(sample, indices_to = "trow")

oxy.mod.mu.PI.df <- do.call(data.frame, oxy.mod.mu.PI.df) %>%
  rename(PI1 = 1, PI99 = 2) %>%
  left_join(oxy.mod %>% select(-c(data, mod)), by = "model") %>%
  mutate(tc = rep(predictorc, 160), t = rep(predictor, 160))



oxy.mod.sigma.PI.df <- tibble(sample = oxy.mod.sigma.PI) %>%
  mutate(model = row_number())

oxy.mod.sigma.PI.df$sample <- lapply(oxy.mod.sigma.PI.df$sample, t)

oxy.mod.sigma.PI.df <- oxy.mod.sigma.PI.df %>% 
  unnest_longer(sample, indices_to = "trow")

oxy.mod.sigma.PI.df <- do.call(data.frame, oxy.mod.sigma.PI.df) %>%
  rename(PI1 = 1, PI99 = 2) %>%
  left_join(oxy.mod %>% select(-c(data, mod)), by = "model") %>%
  mutate(tc = rep(predictorc, 160), t = rep(predictor, 160))

oxy.mod.df <- cbind(oxy.mod.mu.mean.df,
                    oxy.mod.mu.PI.df %>% select(PI1, PI99),
                    oxy.mod.sigma.PI.df %>% select(PI1, PI99)) %>%
              rename(PI1.mu = 13, PI99.mu = 14)


oxy.trim %<>% 
  left_join(oxy.mod.df %>% select(IDm, Timepoint, Position),
            by = "IDm", multiple = "any")

ggplot() +
  geom_point(data = oxy.trim, aes(x = t, y = O2)) +
  geom_line(data = oxy.mod.mu.mean.df, aes(x = t, y = mean)) +
  geom_ribbon(data = oxy.mod.mu.PI.df, aes(x = t, ymin = PI1, ymax = PI99),
              alpha = 0.5) +
  geom_ribbon(data = oxy.mod.sigma.PI.df, aes(x = t, ymin = PI1, ymax = PI99),
              alpha = 0.1) +
  # geom_text(aes(-4, 1, label = IDs)) +
  facet_nested(Timepoint + Group ~ Type + Position + PAR, scales = "free",
               independent = "all", nest_line = element_line()) +
  # facet_wrap(~IDm, scales = "free") +
  theme_minimal()



# blank correction

oxy.mod.blank.slope.df <- oxy.mod.posterior.df %>%
  filter(coefficient == "b", Type == "Blank")
  
  
oxy.mod.slope.df <- oxy.mod.posterior.df %>%
  filter(coefficient == "b", Type == "Sample") %>%
  left_join(oxy.mod.blank.slope.df, by = c("Group", "Timepoint", "PAR", "n")) %>%
  mutate(diff = sample.x - sample.y) # this is where the action happens!
# the blank posterior probability distribution is subtracted from the sample
# posterior probability distributions from the same measurement group

# check correctness of merged dataframe by matching IDs that should be in
# the same measurement group
View(oxy.mod.slope.df %>%
       group_by(IDm.x) %>% 
       summarise(IDm.y = sample(IDm.y, 1), 
                 mx = mean(sample.x), 
                 my = mean(sample.y)))
  

ggplot(oxy.mod.slope.df) +
  geom_density(aes(sample.x), colour = "orange") +
  geom_density(aes(diff), colour = "blue") +
  geom_density(aes(sample.y)) +
  # annotate("text", aes(-0.1, 0, label = IDs.x)) +
  facet_nested(Timepoint + Group ~ Position.x + PAR, scales = "free",
               independent = "all", nest_line = element_line()) +
  # facet_wrap(~IDm, scales = "free") +
  theme_minimal() + 
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())



# put together posterior probability density for all repsonse variables

meta <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Meta.csv") %>%
  mutate(MR = Dry / Blotted, # calculate sample-specific dry to blotted mass ratio
         Massdry = Mass * MR,
         Treatment = ifelse(Temperature == 15 & PAR == 0, 1, # numbering follows increasing  
                            ifelse(Temperature == 15 & PAR == 8, 2, # PAR and temperature
                                   ifelse(Temperature == 20, 3, 4)))) %>% # calculate sample-specific dry mass (g)
  select(-c(Temperature, PAR, MR, Blotted, Dry, Initials))

# randomly assign the 5 initial samples a treatment index and a corresponding tank index
set.seed(2)
meta$Treatment[is.na(meta$Treatment)] <- sample(rep(1:4, 2), 5)
set.seed(10)
meta$Tank[1] <- sample(5:8, 1)
set.seed(20)
meta$Tank[2] <- sample(9:12, 1)
set.seed(10)
meta$Tank[3] <- sample(1:4, 1)
set.seed(20)
meta$Tank[4] <- sample(5:8, 1)
set.seed(10)
meta$Tank[5] <- sample(13:16, 1)

# this is where all photosynthesis variables are calculated!
photo.density.df <- oxy.mod.slope.df %>% 
  select(n, IDs.x, Group, Timepoint, PAR, Position.x, diff) %>% # select blank-corrected values
  filter(PAR == 420) %>% # filter for light measurements
  left_join(oxy.mod.slope.df %>% 
              select(n, IDs.x, Group, Timepoint, PAR, Position.x, diff) %>% # select blank-corrected values
              filter(PAR == 0), # filter for dark measurements
            by = c("n", "IDs.x", "Group", "Timepoint", "Position.x")) %>%
  select(n, IDs.x, Group, Timepoint, Position.x, diff.x, diff.y) %>%
  rename(ID = 2, Position = 5, NP = 6, R = 7) %>%
  left_join(meta %>%
              select(ID, Tank, Treatment, Day, Mass, Massdry),
            by = "ID") %>%
  mutate(NPblot = NP * 0.13 * 60 / Mass, # here all photosynthesis variables are calculated and/or corrected
         Rblot = R * 0.13 * 60 / Mass, # this converts µM min^-1 to µmol g^-1 blotted mass h^-1 
         GPblot = (NP - R) * 0.13 * 60 / Mass, 
         NPDblot = (NP + R) * 0.13 * 60 / Mass * 12, # this converts µM min^-1 to µmol g^-1 blotted mass d^-1 
         NPdry = NP * 0.13 * 60 / Massdry, # this converts µM min^-1 to µmol h^-1 g^-1 dry mass
         Rdry = R * 0.13 * 60 / Massdry,
         GPdry = (NP - R) * 0.13 * 60 / Massdry,
         NPDdry = (NP + R) * 0.13 * 60 / Massdry * 12) %>% # this converts µM min^-1 to µmol g^-1 dry mass d^-1 
  select(-c(NP, R)) %>%
  rename(NP = "NPblot", R = "Rblot", GP = "GPblot", NPD = "NPDblot")
# above the posterior probability distribution of GP is calculated from the distributions
# of NP and R from the same tissue sample; note the - instead of + because respiration 
# is expressed as negative values, so this calculation is effectively GP = NP + R


ggplot(photo.density.df) +
  geom_density(aes(NP), colour = "orange") +
  geom_density(aes(R), colour = "purple") +
  geom_density(aes(GP), colour = "green") +
  # annotate("text", aes(-0.1, 0, label = IDs.x)) +
  facet_nested(Timepoint + Group ~ Position, scales = "free",
               independent = "all", nest_line = element_line()) +
  # facet_wrap(~IDm, scales = "free") +
  theme_minimal() + 
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

# now the slopes are sorted, the intercepts are also of interest as a potential confounder
oxy.mod.intercept.df <- oxy.mod.posterior.df %>%
  filter(coefficient == "a", Type == "Sample", PAR == 420) %>%
  select(-c(coefficient, Type, PAR)) %>%
  left_join(oxy.mod.posterior.df %>%
              filter(coefficient == "a", Type == "Sample", PAR == 0) %>%
              select(-c(coefficient, Type, PAR)),
            by = c("n", "IDs", "Group", "Timepoint", "Position")) %>%
  rename(O2NP = "sample.x", O2R = "sample.y", ID = "IDs") %>%
  mutate(O2NPR = (O2NP + O2R)/2)


# merge intercepts and slopes into one dataframe
photo.density.df %<>% 
  left_join(oxy.mod.intercept.df %>%
              select(-c(IDm.x, IDm.y, model.x, model.y)),
                     by = c("n", "ID", "Group", "Timepoint", "Position"))

# there are also other potential confounders that are not contained in oxy.mod
tech.summary.df <- oxy %>% 
  filter(Type == "Sample") %>%
  select(ID, PAR, Temp, Pressure, Salinity) %>%
  group_by(ID, PAR) %>%
  summarise(Tmean = mean(Temp),
            Tsd = sd(Temp),
            Pmean = mean(Pressure),
            Psd = sd(Pressure),
            Smean = mean(Salinity),
            Ssd = sd(Salinity)) %>% 
  pivot_wider(names_from = PAR, 
              values_from = c(Tmean, Tsd, Pmean, Psd, Smean, Ssd)) %>%
  rename(TmeanR = "Tmean_0", TsdR = "Tsd_0", TmeanNP = "Tmean_420", TsdNP = "Tsd_420",
         PmeanR = "Pmean_0", PsdR = "Psd_0", PmeanNP = "Pmean_420", PsdNP = "Psd_420",
         SmeanR = "Smean_0", SsdR = "Ssd_0", SmeanNP = "Smean_420", SsdNP = "Ssd_420") 
# same salinity for R and NP; pressure and salinity contain zeros in standard deviation so
# variation cannot be modelled -> remove non-usable variables
tech.summary.df %<>% 
  select(-c(PsdR, PsdNP, SmeanR, SsdR, SsdNP)) %>%
  rename(Smean = "SmeanNP")

# incubation temperature is the only one of these confounders with standard deviation 
# that can be modelled, but the standard deviation of the mean of light and dark incubation
# temperature cannot be calculated without use of frequentist theory (variance sum and product laws)
# it will be necessary to run a simple Bayesian model to determine the probability distribution
# of the incubation temperature mean 

oxy.trim.T <- oxy %>% filter(Type == "Sample") %$% 
  data.frame(IDm = IDm, # measurement ID
             IDs = ID, # sample ID
             Group = Group, # measurment group
             PAR = PAR, # µmol photons m^-2 s^-1
             Temp = Temp) %>% # µM
  group_by(IDm) %>%
  mutate(Temp = standardize(Temp)) %>%
  ungroup()

T.mod <- oxy.trim.T %>% nest_by(IDm, IDs, Group, PAR) %>%
  mutate(mod = list(
    ulam(
      alist(
        Temp ~ dnorm(mean = mu, sd = sigma),
        mu <- T_,
        T_ ~ dnorm(mean = 0, sd = 1),
        sigma ~ dexp(rate = 1)
      ), data = data, chains = 8, cores = 8, iter = 1e4
    )
  )
  )


lapply(T.mod$mod, precis, depth = 2) # n_eff and Rhat4 are good
lapply(T.mod$mod, pairs) # no correlation between variables
lapply(T.mod$mod, traceplot) 
lapply(T.mod$mod, trankplot) # chains are healthy


# Prior-posterior comparison
T.mod.prior <- lapply(T.mod$mod, extract.prior)
T.mod.posterior <- lapply(T.mod$mod, extract.samples)

T.mod %<>% rownames_to_column("model") %>%
  mutate(model = as.integer(model)) # turn row names into variable model

T.mod %<>% mutate(Timepoint = as.integer(str_sub(IDs, -1))) %>%
  group_by(Group) %>%
  mutate(Position = match(IDs, unique(IDs)))

T.mod.prior.df <- tibble(sample = T.mod.prior) %>%
  mutate(model = row_number()) %>%
  unnest_longer(sample, indices_to = "coefficient") %>%
  unnest_longer(sample, indices_to = "n") %>%
  left_join(T.mod %>% select(-c(data, mod)), by = "model")

T.mod.posterior.df <- tibble(sample = T.mod.posterior) %>%
  mutate(model = row_number()) %>%
  unnest_longer(sample, indices_to = "coefficient") %>%
  unnest_longer(sample, indices_to = "n") %>%
  left_join(T.mod %>% select(-c(data, mod)), by = "model")

T.mod.prior.posterior.df <- data.frame(rbind(T.mod.prior.df, T.mod.posterior.df),
                                       distribution = c(rep("prior", 228000),
                                                        rep("posterior", 2280000)))

ggplot(data = T.mod.prior.posterior.df %>% filter(coefficient == "T_")) +
  geom_density(aes(x = sample, colour = distribution)) +
  # geom_text(aes(-4, 1, label = IDs)) +
  facet_nested(Timepoint + Group ~ Position + PAR, scales = "free",
               independent = "all", nest_line = element_line()) +
  theme_minimal() + 
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")


T.mod.density.df <- T.mod.posterior.df %>% 
  filter(coefficient == "T_") %>%
  select(n, IDs, Group, Timepoint, PAR, Position, sample) %>%
  filter(PAR == 420) %>% # filter for light measurements
  left_join(T.mod.posterior.df %>% 
              filter(coefficient == "T_") %>%
              select(n, IDs, Group, Timepoint, PAR, Position, sample) %>% # select blank-corrected values
              filter(PAR == 0), # filter for dark measurements
            by = c("n", "IDs", "Group", "Timepoint", "Position")) %>%
  rename(ID = "IDs", TNP = "sample.x", TR = "sample.y")

photo.density.df %<>% 
  left_join(T.mod.density.df %>%
              select(-c(PAR.x, PAR.y, Group, Timepoint, Position)),
            by = c("n", "ID"))


# save progress
write.csv(photo.density.df, "~/Desktop/photo.density.df.csv")

