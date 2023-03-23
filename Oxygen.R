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
# (measurement error in Temp will be incorporated later)

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

# Wernberg et al. 2016 (doi: 10.1002/lno.10362) performed a similar experiment 
# on the same species (Ecklonia radiata), at the same site (Marmion Lagoon) and
# with similar sized tissue (laterals of ~1.5 g blotted mass, ~0.22 g dry mass).
# Priors for b can therefore be estimated by back-transforming their oxygen 
# evolution slopes (mg g^-1 h^-1) under light and dark conditions:

# get prior Ecklonia photosynthesis data
prior <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Prior.csv")
m <- 0.22 # g
V <- 0.12 # l
t <- 60 # min

# calculate light prior mean
b.mu.l <- prior %>% filter(Reference == "Wernberg et al. 2016", 
                           Temperature >= 15, Temperature <= 20) %>% 
                           # these are closest to my experimental temperature
                    pull(NP) %>% mean() * m / V / t / g.mol * 1e3
b.mu.l # 2.376159 µM min^-1



# calculate dark prior mean
b.mu.d <- prior %>% filter(Reference == "Wernberg et al. 2016", 
                           Temperature >= 15, Temperature <= 20) %>%
                    pull(R) %>% mean() * m / V / t / g.mol * 1e3 * -1
b.mu.d # -1.03052 µM min^-1

# blank prior mean is just 0

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
b.prior.l <- rnorm(n = 1e5, mean = 2.376159, sd = 4)

ggplot() + 
  geom_density(aes(x = b.prior.l)) +
  theme_minimal()

# dark prior
b.prior.d <- rnorm(n = 1e5, mean = -1.03052, sd = 2) 
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

b.prior.l <- rnorm(n = 1e5, mean = 2.376159, sd = 2)

ggplot() + 
  geom_density(aes(x = b.prior.l)) +
  geom_vline(aes(xintercept = c(2.376159 - 2 * 2,  
                                2.376159 + 2 * 2))) +
  theme_minimal()
# 95% of prior probability density still lies between -1.623841 and 6.376159 µM

b.prior.d <- rnorm(n = 1e5, mean = -1.03052, sd = 1)

ggplot() + 
  geom_density(aes(x = b.prior.d)) +
  geom_vline(aes(xintercept = c(-1.03052 - 2 * 1,  
                                -1.03052 + 2 * 1))) +
  theme_minimal()
# 95% of prior probability density still lies between -3.03052 and 0.96948 µM

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

# the following model is meaningless because it is fit to the entire data,
# but I am only using it to get at the priors

# light prior
prior.test.l <- ulam( # ulam is a Hamiltonian Monte Carlo engine
                alist(
                  O2 ~ dnorm(mean = mu, sd = sigma),
                  mu <- a + b * t,
                  a ~ dnorm(mean = 232.5, sd = 30),
                  b ~ dnorm(mean = 2.376159, sd = 2),
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
                  b ~ dnorm(mean = -1.03052, sd = 1),
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
                                      b ~ dnorm(mean = -1.03052, sd = 1),
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
                                      b ~ dnorm(mean = 2.376159, sd = 2),
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
                                      b ~ dnorm(mean = -1.03052, sd = 1),
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
                                      b ~ dnorm(mean = 2.376159, sd = 2),
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
                                      b ~ dnorm(mean = -1.03052, sd = 1),
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
                                      b ~ dnorm(mean = 2.376159, sd = 2),
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

photo.density.df <- oxy.mod.slope.df %>% 
  select(n, IDs.x, Group, Timepoint, PAR, Position.x, diff) %>%
  filter(PAR == 420) %>%
  left_join(oxy.mod.slope.df %>% 
              select(n, IDs.x, Group, Timepoint, PAR, Position.x, diff) %>%
              filter(PAR == 0), 
            by = c("n", "IDs.x", "Group", "Timepoint", "Position.x")) %>%
  select(n, IDs.x, Group, Timepoint, Position.x, diff.x, diff.y) %>%
  rename(ID = 2, Position = 5, NP = 6, R = 7) %>%
  mutate(GP = NP - R) # this is where the action happens!
# the posterior probability distribution of GP is calculated from the distributions
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

# now the slopes are sorted, the intercepts are also of interest as a potential confound
oxy.mod.intercept.df <- oxy.mod.posterior.df %>%
  filter(coefficient == "a", Type == "Sample", PAR == 420) %>%
  select(-c(coefficient, Type, PAR)) %>%
  left_join(oxy.mod.posterior.df %>%
              filter(coefficient == "a", Type == "Sample", PAR == 0) %>%
              select(-c(coefficient, Type, PAR)),
            by = c("n", "IDs", "Group", "Timepoint", "Position")) %>%
  rename(O2NP = "sample.x", O2R = "sample.y", ID = "IDs")


# merge intercepts and slopes into one dataframe
photo.density.df %<>% 
  left_join(oxy.mod.intercept.df %>%
              select(-c(IDm.x, IDm.y, model.x, model.y)),
                     by = c("n", "ID", "Group", "Timepoint", "Position"))


photo.summary.df <- photo.density.df %>%
  group_by(ID) %>% summarise(Group = mean(Group),
                             Position = mean(Position),
                             Timepoint = mean(Timepoint),
                             NPmean = mean(NP),
                             NPmode = chainmode(NP),
                             NPmedian = median(NP),
                             NPsd = sd(NP),
                             NPpi1 = PI(NP, prob = 0.98)[[1]],
                             NPpi99 = PI(NP, prob = 0.98)[[2]],
                             NPhpdilo = HPDI(NP, prob = 0.98)[[1]],
                             NPhpdihi = HPDI(NP, prob = 0.98)[[2]],
                             Rmean = mean(R),
                             Rmode = chainmode(R),
                             Rmedian = median(R),
                             Rsd = sd(R),
                             Rpi1 = PI(R, prob = 0.98)[[1]],
                             Rpi99 = PI(R, prob = 0.98)[[2]],
                             Rhpdilo = HPDI(R, prob = 0.98)[[1]],
                             Rhpdihi = HPDI(R, prob = 0.98)[[2]],
                             GPmean = mean(GP),
                             GPmode = chainmode(GP),
                             GPmedian = median(GP),
                             GPsd = sd(GP),
                             GPpi1 = PI(GP, prob = 0.98)[[1]],
                             GPpi99 = PI(GP, prob = 0.98)[[2]],
                             GPhpdilo = HPDI(GP, prob = 0.98)[[1]],
                             GPhpdihi = HPDI(GP, prob = 0.98)[[2]],
                             O2NPmean = mean(O2NP),
                             O2NPmode = chainmode(O2NP),
                             O2NPmedian = median(O2NP),
                             O2NPsd = sd(O2NP),
                             O2NPpi1 = PI(O2NP, prob = 0.98)[[1]],
                             O2NPpi99 = PI(O2NP, prob = 0.98)[[2]],
                             O2NPhpdilo = HPDI(O2NP, prob = 0.98)[[1]],
                             O2NPhpdihi = HPDI(O2NP, prob = 0.98)[[2]],
                             O2Rmean = mean(O2R),
                             O2Rmode = chainmode(O2R),
                             O2Rmedian = median(O2R),
                             O2Rsd = sd(O2R),
                             O2Rpi1 = PI(O2R, prob = 0.98)[[1]],
                             O2Rpi99 = PI(O2R, prob = 0.98)[[2]],
                             O2Rhpdilo = HPDI(O2R, prob = 0.98)[[1]],
                             O2Rhpdihi = HPDI(O2R, prob = 0.98)[[2]],
                             Samples = length(n))

meta <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Meta.csv") %>%
  mutate(MR = Dry / Blotted) %>% # calculate sample-specific dry to blotted mass ratio
  mutate(MassD = Mass * MR) # calculate sample-specific dry mass (g)

meta %<>% left_join(photo.summary.df,
                    by = "ID")
  

tech.summary.df <- oxy %>% 
  select(ID, PAR, Temp, Pressure, Salinity) %>%
  group_by(ID, PAR) %>%
  summarise(Tmean = mean(Temp),
            Tsd = sd(Temp),
            Pmean = mean(Pressure),
            Psd = sd(Pressure),
            Salinity = mean(Salinity)) %>% # no standard deviation for salinity as 
# it is based on a single measurement per day so there is no variation at this level
  pivot_wider(names_from = PAR, 
              values_from = c(Tmean, Tsd, Pmean, Psd, Salinity)) %>%
  select(-11) %>%
  rename(TmeanR = "Tmean_0", TsdR = "Tsd_0", TmeanNP = "Tmean_420", TsdNP = "Tsd_420",
         PmeanR = "Pmean_0", PsdR = "Psd_0", PmeanNP = "Pmean_420", PsdNP = "Psd_420",
         Salinity = "Salinity_0") # same salinity for R and NP

meta %<>% left_join(tech.summary.df, by = "ID")

# mass correction and treatment redefinition
meta %<>% mutate(NPmean.c = NPmean * 0.13 * 60 / Mass, # convert µM min^-1 to µmol g^-1 h^-1
                 NPsd.c = NPsd * 0.13 * 60 / Mass, # (0.13 l is the incubation volume)
                 Rmean.c = Rmean * 0.13 * 60 / Mass,
                 Rsd.c = Rsd * 0.13 * 60 / Mass,
                 GPmean.c = GPmean * 0.13 * 60 / Mass,
                 GPsd.c = GPsd * 0.13 * 60 / Mass,
                 Treatment = ifelse(Temperature == 15 & PAR == 0, 1, # numbering follows increasing  
                                    ifelse(Temperature == 15 & PAR == 8, 2, # PAR and temperature
                                           ifelse(Temperature == 20, 3, 4))))


write.csv(meta, "~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Meta.update.csv",
          row.names = FALSE) # save to hard drive to continue work at a later date





# Net photosynthesis model based on blotted mass
meta.NP <- meta %>%
  select(Tank, Day, Treatment, # main predictor variables
         NPmean.c, NPsd.c, # response variable
         O2NPmean, O2NPsd, TmeanNP, TsdNP, PmeanNP, PsdNP, # potential confounds with measurement error
         Salinity, Mass) %>% # potential confounds without measurement error
  rename(NPmean = "NPmean.c", NPsd = "NPsd.c", O2mean = "O2NPmean", O2sd = "O2NPsd", 
         Tmean = "TmeanNP", Tsd = "TsdNP", Pmean = "PmeanNP", Psd = "PsdNP")



# Wright et al. 2022 (doi: 10.1111/gcb.16299) provide detrital photosynthesis over detrital age
# slopes for three Laminaria species: -0.001, -0.009 and -0.029 mg C g^-1 dry mass h^-1 d^-1
# Frontier et al. 2021 (doi: 10.1016/j.marenvres.2021.105277) is the only similarly useful prior
# study but unfortunately does not provide any slopes due to categorical analysis
r <- 0.127003814 # Wright et al. 2022 also provide this mean dry to wet mass ratio
g.mol <- 12.0107 # g mol^-1 for carbon

# note that somewhat counter-intuitively their rates per dry mass need to be multiplied by r to get 
# rates per wet mass because a given dry mass photosynthesises more than the same amount of wet mass
b.mu <- mean(c(-0.001, -0.009, -0.029)) * r / g.mol * 1e3
b.mu # -0.1374649 µmol O2 g^-1 wet mass h^-1 d^-1

# Staehr & Wernberg 2009 and Wernberg et al. 2016 provide baseline dry mass photosynthesis rates 
# but no mass conversion, but I can use our mean dry to wet mass ratio:
r <- meta %>% filter(!is.na(MR)) %>% pull(MR) %>% mean()
g.mol <- 15.9994 # g mol^-1 for oxygen

a.mu <- prior %>% filter(Temperature >= 15, Temperature <= 20) %>% 
                    # these are closest to my experimental temperature
                    pull(NP) %>% mean() * r / g.mol * 1e3
a.mu # 15.39553 µmol O2 g^-1 wet mass h^-1


b.prior <- rnorm(n = 1e5, mean = b.mu, sd = 1) # high uncertainty because there are no data for Ecklonia

ggplot() + 
  geom_density(aes(x = b.prior)) +
  theme_minimal()

a.prior <- rnorm(n = 1e5, mean = a.mu, sd = 4)

ggplot() + 
  geom_density(aes(x = a.prior)) +
  theme_minimal()


predictor <- meta.NP %$% seq(min(Day), max(Day), length.out = 100)
dev.off() 
plot(NULL, xlim = c(0, 120), ylim = c(-30, 30), # base plot is better with loops
     xlab = "d", ylab = "µmol g^-1 h^-1") # empty plot
meta.NP %>% filter(!is.na(NPmean)) %$% abline(h = min(NPmean), lty = 2) # data range
meta.NP %>% filter(!is.na(NPmean)) %$% abline(h = max(NPmean), lty = 2)
abline(h = 0)
for(i in 1:1e3) curve(a.prior[i] + b.prior[i] * x, # loop
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = col.alpha("black", 0.1))
# slope and intercept are way too variable

b.prior <- rnorm(n = 1e5, mean = b.mu, sd = 0.05) 
a.prior <- rnorm(n = 1e5, mean = a.mu, sd = 2)

plot(NULL, xlim = c(0, 120), ylim = c(-30, 30), # base plot is better with loops
     xlab = "d", ylab = "µmol g^-1 h^-1") # empty plot
meta.NP %>% filter(!is.na(NPmean)) %$% abline(h = min(NPmean), lty = 2) # data range
meta.NP %>% filter(!is.na(NPmean)) %$% abline(h = max(NPmean), lty = 2)
abline(h = 0)
for(i in 1:1e3) curve(a.prior[i] + b.prior[i] * x, # loop
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = col.alpha("black", 0.1))
# still looks variable but allows for net heterotrophy after 20 d (cf. Wright et al. 2022)
# while indicating the general negative trend evident from the literature

bc.prior <- rnorm(n = 1e5, mean = 0, sd = 0.1) # potential confound slope prior, 
# assumed to be neutral with high certainty because of controlled environment and procedure

ggplot() + 
  geom_density(aes(x = bc.prior)) +
  theme_minimal()

# potential confounds measured with error additionally need a prior for the intercept
# for initial incubation seawater O2, the same prior as in oxy.mod is used: dnorm(232.5, 30)
aT.prior <- rnorm(n = 1e5, mean = mean(c(15,20)), sd = 0.2) # incubation temperature 
# probably lies somewhere between controlled-temperature room (20°C) and chilled seawater (~15°C)

ggplot() + 
  geom_density(aes(x = aT.prior)) +
  theme_minimal()

aP.prior <- rnorm(n = 1e5, mean = 1000, sd = 20) # pressure should lie around atmospheric (~1 bar, ~1000 hPa)

ggplot() + 
  geom_density(aes(x = aP.prior)) +
  theme_minimal()


# the planned model requires complete cases and needs to be in list format

meta.NP.cc <- meta.NP %>%
  filter(!is.na(Tank) & !is.na(NPmean)) %>%
  select(-Psd) # it turns out Psd has zeros which cannot be used to model standard deviation

meta.NP.cc.l <- meta.NP.cc %>% as.list()


# meta.trim.cc %<>% mutate(Treatment = as.integer(Treatment))

meta.NP.cc.l$N <- nrow(meta.NP.cc)


NP.mod <- ulam(
  alist(
    # likelihood 
    NPmean ~ dnorm(mean = NP, sd = NPsd), # here the observed measurement error is introduced
    vector[N]:NP ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP
    
    # linear model
    mu <- a[Treatment] + b[Treatment] * Day + # treatment of interest (temp/PAR)
      at[Tank] + bt[Tank] * Day + # tank effect
      ac + bct * T_[i] + bco * O2[i] + bcp * Pmean + 
      bcs * Salinity + bcm * Mass, # potential confounds
    
    # predictor measurement error
    Tmean ~ dnorm(mean = T_, sd = Tsd),
    vector[N]:T_ ~ dnorm(mean = 17.5, sd = 0.2),
    O2mean ~ dnorm(mean = O2, sd = O2sd),
    vector[N]:O2 ~ dnorm(mean = 232.5, sd = 30),
    
    # higher level priors for partial pooling
    a[Treatment] ~ dnorm(mean = a_, sd = sa),
    b[Treatment] ~ dnorm(mean = b_, sd = sb),
    at[Tank] ~ dnorm(mean = at_, sd = sat),
    bt[Tank] ~ dnorm(mean = bt_, sd = sbt),
    
    # priors
    c(a_,at_,ac) ~ dnorm(mean = 15.39553, sd = 1), # intercepts
    c(b_,bt_) ~ dnorm(mean = -0.1374649, sd = 0.05), # slopes across detrital age
    c(bct,bco,bcp,bcs,bcm) ~ dnorm(mean = 0, sd = 0.1), # potential confound slopes
    c(sigma,sa,sb,sat,sbt) ~ dexp(rate = 1) # standard deviations
  ), 
  data = meta.NP.cc.l, chains = 8, cores = 8, iter = 1e3,
  control = list(adapt_delta = 0.99),
)

# note warnings!
traceplot(NP.mod)
trankplot(NP.mod) # very unhealthy chains
precis(NP.mod, depth = 2) # shocking Rhat4 and n_eff

# try to fix with non-centred priors, i.e. priors that are not hierarchical

NP.mod.nc <- ulam(
  alist(
    # likelihood
    NPmean ~ dnorm(mean = mu, sd = NPsd), # here the observed measurement error is introduced
    vector[N]:NP ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP
    
    # linear model
    mu <- a[Treatment] + b[Treatment] * Day + # treatment of interest (temp/PAR)
      at[Tank] + bt[Tank] * Day + # tank effect
      ac + bct * T_[i] + bco * O2[i] + bcp * Pmean + 
      bcs * Salinity + bcm * Mass, # potential confounds
    
    # predictor measurement error
    Tmean ~ dnorm(mean = T_, sd = Tsd),
    vector[N]:T_ ~ dnorm(mean = 17.5, sd = 0.2),
    O2mean ~ dnorm(mean = O2, sd = O2sd),
    vector[N]:O2 ~ dnorm(mean = 232.5, sd = 30),
    
    # define multilevel coefficients using z-scores that HMC can sample better
    save> vector[4]:a <<- a_ + za * sa,
    save> vector[4]:b <<- b_ + zb * sb,
    save> vector[16]:at <<- at_ + zat * sat,
    save> vector[16]:bt <<- bt_ + zbt * sbt,
    
    # priors
    vector[4]:za ~ dnorm(mean = 0, sd = 1), # z-score priors
    vector[4]:zb ~ dnorm(mean = 0, sd = 1),
    vector[16]:zat ~ dnorm(mean = 0, sd = 1),
    vector[16]:zbt ~ dnorm(mean = 0, sd = 1),
    
    c(a_,at_,ac) ~ dnorm(mean = 15.39553, sd = 1), # intercept priors
    c(b_,bt_) ~ dnorm(mean = -0.1374649, sd = 0.05), # slope priors
    c(bct,bco,bcp,bcs,bcm) ~ dnorm(mean = 0, sd = 0.1), # potential confound slope priors
    c(sigma,sa,sb,sat,sbt) ~ dexp(rate = 1) # standard deviation priors
  ), 
  data = meta.NP.cc.l, chains = 8, cores = 8, iter = 1e3,
  control = list(adapt_delta = 0.99),
)

traceplot(NP.mod.nc)
trankplot(NP.mod.nc) # chains are even worse! who would have thought that possible?!
precis(NP.mod.nc, depth = 2) # extraordinarily small n_eff and large Rhat4


# the only remaining option is to model the effects of interest and 
# the potential confounds separately


#########################################
NP.mod.main <- ulam(
  alist(
    # likelihood 
    NPmean ~ normal(mean = NP, sd = NPsd), # here the observed measurement error is introduced
    vector[N]:NP ~ normal(mean = mu, sd = sigma), # this describes the true, unobserved NP
    
    # linear model
    mu <- atreat[Treatment] + btreat[Treatment] * Day + # treatment of interest (temp/PAR)
      atank[Tank] + btank[Tank] * Day, # tank effect
    
    # adaptive priors
    c(atreat, btreat)[Treatment] ~ multi_normal(c(a_treat, b_treat), rhotreat, sigmatreat),
    c(atank, btank)[Tank] ~ multi_normal(Mu = c(a_tank, b_tank), rhotank, sigmatank),
    
    # hyper-priors
    c(a_treat,a_tank) ~ normal(mean = 15.39553, sd = 1), # intercepts
    c(b_treat,b_tank) ~ normal(mean = -0.1374649, sd = 0.05), # slopes across detrital age
    c(sigma,sigmatreat,sigmatank) ~ exponential(rate = 1), # standard deviations
    c(rhotreat, rhotank) ~ lkj_corr(eta = 2)
  ), 
  data = meta.NP.cc.l, chains = 8, cores = 8, iter = 1e3,
  control = list(adapt_delta = 0.99),
)

traceplot(NP.mod.main)
trankplot(NP.mod.main) # chains are quite healthy
precis(NP.mod.main, depth = 2) # ok n_eff and no extraordinarily large large Rhat4;
# but note that the model has a hard time sampling the nested, centred priors
# -> try to fix with non-centred priors

NP.mod.main.nc <- ulam(
  alist(
    # likelihood
    NPmean ~ dnorm(mean = mu, sd = NPsd), # here the observed measurement error is introduced
    vector[N]:NP ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP
    
    # linear model
    mu <- a[Treatment] + b[Treatment] * Day + # treatment of interest (temp/PAR)
      at[Tank] + bt[Tank] * Day, # tank effect
    
    # define multilevel coefficients using z-scores that HMC can sample better
    save> vector[4]:a <<- a_ + za * sa,
    save> vector[4]:b <<- b_ + zb * sb,
    save> vector[16]:at <<- at_ + zat * sat,
    save> vector[16]:bt <<- bt_ + zbt * sbt,
    
    # priors
    vector[4]:za ~ dnorm(mean = 0, sd = 1), # z-score priors
    vector[4]:zb ~ dnorm(mean = 0, sd = 1),
    vector[16]:zat ~ dnorm(mean = 0, sd = 1),
    vector[16]:zbt ~ dnorm(mean = 0, sd = 1),
    
    c(a_,at_) ~ dnorm(mean = 15.39553, sd = 1), # intercept priors
    c(b_,bt_) ~ dnorm(mean = -0.1374649, sd = 0.05), # slope priors
    c(sigma,sa,sb,sat,sbt) ~ dexp(rate = 1) # standard deviation priors
  ), 
  data = meta.NP.cc.l, chains = 8, cores = 8, iter = 1e3,
  control = list(adapt_delta = 0.99),
)

traceplot(NP.mod.main.nc)
trankplot(NP.mod.main.nc) # chains are ok
precis(NP.mod.main.nc, depth = 2) # less good n_eff and some large Rhat4;
# note that this model is better at some aspects but worse at others
# hard to tell so compare n_eff between models visually


# visual cross-validation
precis_c <- precis(NP.mod.main, depth = 2)
precis_nc <- precis(NP.mod.main.nc, depth = 2)

pars <- c(paste("NP[",1:52,"]",sep="") , paste("a[",1:4,"]",sep=""),
          paste("b[",1:4,"]",sep="") , paste("at[",1:16,"]",sep=""),
          paste("bt[",1:16,"]",sep=""), "at_" , "a_", "bt_", "b_",
          "sbt", "sat", "sb", "sa", "sigma")

neff_table <- cbind(precis_c[pars,"n_eff"] , precis_nc[pars,"n_eff"] )
dev.off() 
plot( neff_table , xlim=range(neff_table) , ylim=range(neff_table) ,
xlab="n_eff (centred)" , ylab="n_eff (non-centred)" , lwd=2 )
abline( a=0 , b=1 , lty=2 )

# the centred prior model clearly wins!
  

# Prior-posterior comparison
NP.mod.main.posterior <- extract.samples(NP.mod.main)
NP.mod.main.prior <- extract.prior(NP.mod.main)

NP.mod.main.posterior.df <- data.frame(NP.mod.main.posterior) %>%
  rownames_to_column(var = "n") %>%
  select(-c(2:53)) %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "treatment")) %>%
  mutate(treatment = as.integer(treatment)) %>%
  arrange(coefficient, treatment, n) %>%
  mutate(treatment = factor(treatment), coefficient = factor(coefficient))


NP.mod.main.prior.df <- data.frame(NP.mod.main.prior) %>%
  rownames_to_column(var = "n") %>%
  select(-c(2:53)) %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "treatment")) %>%
  mutate(treatment = as.integer(treatment)) %>%
  arrange(coefficient, treatment, n) %>%
  mutate(treatment = factor(treatment), coefficient = factor(coefficient))


NP.mod.main.prior.posterior.df <- data.frame(rbind(NP.mod.main.prior.df, NP.mod.main.posterior.df),
                                             distribution = c(rep("prior", 49000),
                                                              rep("posterior", 196000))) %>%
  mutate(distribution = factor(distribution))

ggplot(data = NP.mod.main.prior.posterior.df) +
  geom_density(aes(x = sample, colour = distribution)) +
  facet_nested_wrap(~ coefficient + treatment, scales = "free",
                    nest_line = element_line()) +
  theme_minimal() + 
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")


# Posterior prediction intervals
predictor.df <- meta.NP.cc %$% expand_grid(Day = seq(min(Day), max(Day), length.out = 50),
                                           Treatment = factor(Treatment),
                                           Tank = 16)

NP.mod.main.mu <- link(NP.mod.main, data = predictor.df)
predictor.df$fit <- apply(NP.mod.main.mu, 2, mean)







meta.density <- photo.density.df %>%
  group_by(ID) %>% slice_sample(n = 20000, replace = FALSE) %>%
  left_join(meta %>% select(ID, Tank, Day, Timepoint, Temperature, PAR, Mass, TmeanNP, PmeanNP, Salinity, ),
            by = c("ID", "Timepoint")) %>%
  mutate(NP = NP * 0.13 * 60 / Mass, # convert µM min^-1 to µmol g^-1 h^-1
         GP = GP * 0.13 * 60 / Mass, # (0.13 l is the incubation volume)
         R = R * 0.13 * 60 / Mass,
         Treatment = ifelse(Temperature == 15 & PAR == 0, 1,
                            ifelse(Temperature == 15 & PAR == 8, 2,
                                   ifelse(Temperature == 20, 3, 4))))

ggplot() +
  geom_violin(data = meta.density, aes(x = Day, y = NP, group = ID, colour = Treatment, fill = Treatment)) +
  facet_grid(~Treatment)

require(ggdensity)
ggplot() +
  geom_hline(yintercept = 0) +
  geom_hdr(data = meta.density, aes(x = O2NP, y = NP, group = ID)) +
  labs(x = expression("Initial dissolved O"[2]*" ("*mu*"M)"),
       y = expression("Net photosynthesis ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  theme_minimal() +
  theme(legend.position = "top")


ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = meta.density %>% 
                filter(!is.na(Treatment)) %>%
                mutate(Treatment = ifelse(Treatment == 1, "15°C / dark",
                                          ifelse(Treatment == 2, "15°C / light",
                                                 ifelse(Treatment == 3, "20°C / light", "25°C / light")))), 
              aes(x = Day, y = NP, group = ID)) +
  labs(x = "Detrital age (d)",
       y = expression("Net photosynthesis ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  facet_wrap(~Treatment) +
  theme_minimal() +
  theme(legend.position = "top")


ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = meta.density, aes(x = TmeanNP, y = NP, group = ID)) +
  labs(x = "Incubation temperature (°C)",
       y = expression("Net photosynthesis ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  theme_minimal() +
  theme(legend.position = "top")


ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = meta.density, aes(x = PmeanNP, y = NP, group = ID)) +
  labs(x = "Incubation pressure (hPa)",
       y = expression("Net photosynthesis ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  theme_minimal() +
  theme(legend.position = "top")

ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = meta.density, aes(x = Salinity, y = NP, group = ID)) +
  labs(x = "Incubation salinity (‰)",
       y = expression("Net photosynthesis ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  theme_minimal() +
  theme(legend.position = "top")


ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = meta.density, aes(x = Mass, y = NP, group = ID)) +
  labs(x = "Incubation mass (g)",
       y = expression("Net photosynthesis ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")")) +
  theme_minimal() +
  theme(legend.position = "top")



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

# post <- extract.samples(NP.mod)
# NPmean_impute_mu <- apply(post$NPmean_impute , 2 , mean )
# NPmean_impute_ci <- apply(post$NPmean_impute , 2 , PI )
# 
# 
# plot( meta.list.pcc$Day, meta.list.pcc$NPmean , pch=16 , col=rangi2 ,
# xlab="Day" , ylab="NPmean" )
# miss_idx <- which( is.na(meta.list.pcc$NPmean) )
# Di <- meta.list.pcc$Day[miss_idx]
# points(Di , NPmean_impute_mu)
# for ( i in 1:44 ) lines( rep(Di[i],2), NPmean_impute_ci[,i] )
# 
# 
# plot( meta.list.pcc$Treatment, meta.list.pcc$NPmean , pch=16 , col=rangi2 ,
#       xlab="Day" , ylab="NPmean" )
# miss_idx <- which( is.na(meta.list.pcc$NPmean) )
# Ti <- meta.list.pcc$Treatment[miss_idx]
# points(Ti , NPmean_impute_mu)
# for ( i in 1:44 ) lines( rep(Ti[i],2), NPmean_impute_ci[,i] )
















require(dagitty)
require(ggdag)
  
dag <- dagitty("dag {
               A -> P -> P_
               A -> D -> P_
               T -> P -> P_
               e -> P_
               }") %>%
  tidy_dagitty() %>%
  mutate(Observation = ifelse(name == "P" | name == "e", "Unobserved", "Observed"))
  

ggplot(dag, aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_dag_point(aes(colour = Observation), shape = 21) +
  geom_dag_edges() +
  geom_dag_text(colour = "black") +
  scale_colour_manual(values = c("white", "black")) +
  theme_dag() +
  theme(legend.position = "none")




  
  



# data(rugged)
# d <- rugged
# d$log_gdp <- log(d$rgdppc_2000)
# dd <- d[ complete.cases(d$rgdppc_2000) , ]
# dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
# dd$rugged_std <- dd$rugged / max(dd$rugged)
# dd$cid <- ifelse( dd$cont_africa==1 , 1 , 2 )
# 
# dat_slim <- list(
#   log_gdp_std = dd$log_gdp_std,
#   rugged_std = dd$rugged_std,
#   cid = as.integer( dd$cid )
# )
# str(dat_slim)
# 
# m9.1 <- ulam(
#   alist(
#     log_gdp_std ~ dnorm( mu , sigma ) ,
#     mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
#     a[cid] ~ dnorm( 1 , 0.1 ) ,
#     b[cid] ~ dnorm( 0 , 0.3 ) ,
#     sigma ~ dexp( 1 )
#   ) , data=dat_slim , chains=1 )
# 
# precis( m9.1 , depth=2 )
# 
# 
# m9.1 <- ulam(
#   alist(
#     log_gdp_std ~ dnorm( mean = mu , sd = sigma ) ,
#     mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
#     a[cid] ~ dnorm( mean = 1 , sd = 0.1 ) ,
#     b[cid] ~ dnorm( mean = 0 , sd = 0.3 ) ,
#     sigma ~ dexp( rate = 1 )
#   ) , data=dat_slim , chains=4 , cores=4 )
# 
# show(m9.1) # model details
# precis( m9.1 , depth = 2 ) # model results, similar to summary(m9.1)
# pairs(m9.1) # pairs plot of posterior samples and distribution
# traceplot(m9.1, chains = 1) # traceplot of estimates vs samples
# trankplot( m9.1 , n_cols=2 ) # trankplot of frequency vs rank
# stancode(m9.1) # show underlying Stan code
# 
# 
# y <- c(-1,1)
# set.seed(11)
# m9.2 <- ulam(
#     alist(
#         y ~ dnorm( mean = mu , sd = sigma ) ,
#         mu <- alpha ,
#         alpha ~ dnorm( mean = 0 , sd = 1000 ) ,
#         sigma ~ dexp( rate = 0.0001 )
#         ) , data=list(y=y) , chains=3 )
# 
# precis(m9.2)
# pairs( m9.2@stanfit ) 
# # Stan version of pairs plot with divergent transitions in red
# 
# traceplot(m9.2, chains = 2)
# trankplot(m9.2)
# 
# set.seed(11)
# m9.3 <- ulam(
#   alist(
#     y ~ dnorm( mean = mu , sd = sigma ) ,
#     mu <- alpha ,
#     alpha ~ dnorm( mean = 1 , sd = 10 ) ,
#     sigma ~ dexp( rate = 1 )
#   ) , data=list(y=y) , chains=3 )
# 
# precis(m9.3)
# pairs( m9.3@stanfit ) 
# traceplot(m9.3)
# trankplot(m9.3)
# 
# 
# plot(seq(-100, 100, 0.01), dnorm(x = seq(-100, 100, 0.01), mean = 0, sd = 1000),
#      type = "l")
# 
# 
# set.seed(41)
# y <- rnorm(100, mean = 0, sd = 1)
# set.seed(384)
# m9.4 <- ulam(
#   alist(
#     y ~ dnorm(mean = mu, sd = sigma),
#     mu <- a1 + a2,
#     a1 ~ dnorm(mean = 0, sd = 1000),
#     a2 ~ dnorm(mean = 0, sd = 1000),
#     sigma ~ dexp(rate = 1)
#   ),
#   data = list(y = y), chains = 3)
# precis(m9.4)
# traceplot(m9.4)
# trankplot(m9.4)
# 
# 
# m9.5 <- ulam(
#   alist(
#     y ~ dnorm(mean = mu, sd = sigma),
#     mu <- a1 + a2,
#     a1 ~ dnorm(mean = 0, sd = 10),
#     a2 ~ dnorm(mean = 0, sd = 10),
#     sigma ~ dexp(rate = 1)
#   ),
#   data = list(y = y), chains = 3)
# precis(m9.5)
# traceplot(m9.5)
# trankplot(m9.5)
# 


# data(WaffleDivorce)
# d <- WaffleDivorce
# dlist <- list(
# D_obs = standardize( d$Divorce ),
# D_sd = d$Divorce.SE / sd( d$Divorce ),
# M = standardize( d$Marriage ),
# A = standardize( d$MedianAgeMarriage ),
# N = nrow(d)
# )
# 
# m15.1 <- ulam(
# alist(
# D_obs ~ dnorm( D_true , D_sd ),
# vector[N]:D_true ~ dnorm( mu , sigma ),
# mu <- a + bA*A + bM*M,
# a ~ dnorm(0,0.2),
# bA ~ dnorm(0,0.5),
# bM ~ dnorm(0,0.5),
# sigma ~ dexp(1)
# ) , data=dlist , chains=4 , cores=4 )


# data("bangladesh")
# d <- bangladesh
# dat <- list(C = d$use.contraception,
#             D = as.integer(d$district),
#             U = ifelse(d$urban == 1, 1, 0))
# mCDUnc <- ulam(
#   alist(
#     C ~ dbern(prob = p),
#     logit(p) <- a[D] + b[D] * U,
#     save> vector[61]:a <<- a_ + za * as,
#     save> vector[61]:b <<- b_ + zb * bs,
#     vector[61]:za ~ dnorm(mean = 0, sd = 1),
#     vector[61]:zb ~ dnorm(mean = 0, sd = 1),
#     c(a_,b_) ~ dnorm(mean = 0, sd = 1),
#     c(as,bs) ~ dexp(rate = 1)
#   ), 
#   data = dat, chains = 4, cores = 4
# )



