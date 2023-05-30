### Packages ####
require(tidyverse) # this is my preferred data processing and visualisation framework 
require(magrittr) # this extends the pipe operator (%>%) of the tidyverse
require(rethinking) # this provides my Bayesian data analysis framework of choice

# load progress 
photo.density.df <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/photo.density.df.csv") %>%
                    select(-1)

# load original data
oxy <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Oxygen.csv")

oxy.summary <- oxy %>% filter(Type == "Sample", PAR == 420) %>%
  group_by(ID) %>% summarise(TNPmean = mean(Temp),
                             TNPsd = sd(Temp)) %>%
  left_join(oxy %>% filter(Type == "Sample", PAR == 0) %>%
              group_by(ID) %>% summarise(TRmean = mean(Temp),
                                         TRsd = sd(Temp)),
            by = "ID")

# unstandardise TNP and TR
photo.density.df %<>% left_join(oxy.summary, by = "ID") %>%
  group_by(ID) %>%
  mutate(TNP = TNP * TNPsd + TNPmean,
         TR = TR * TRsd + TRmean,
         TNPR = (TNP + TR)/2) %>%
  select(-24:-27)



# now summarise all posterior distributions
photo.summary.df <- photo.density.df %>%
  group_by(ID) %>% summarise(Tank = mean(Tank),
                             Treatment = mean(Treatment),
                             Day = mean(Day),
                             NPmean = mean(NP),
                             NPsd = sd(NP),
                             # NPpi1 = PI(NP, prob = 0.98)[[1]],
                             # NPpi99 = PI(NP, prob = 0.98)[[2]],
                             Rmean = mean(R),
                             Rsd = sd(R),
                             GPmean = mean(GP),
                             GPsd = sd(GP),
                             NPDmean = mean(NPD),
                             NPDsd = sd(NPD),
                             NPdrymean = mean(NPdry),
                             NPdrysd = sd(NPdry),
                             Rdrymean = mean(Rdry),
                             Rdrysd = sd(Rdry),
                             GPdrymean = mean(GPdry),
                             GPdrysd = sd(GPdry),
                             NPDdrymean = mean(NPDdry),
                             NPDdrysd = sd(NPDdry),
                             Mass = mean(Mass),
                             Massdry = mean(Massdry),
                             O2NPmean = mean(O2NP),
                             O2NPsd = sd(O2NP),
                             O2Rmean = mean(O2R),
                             O2Rsd = sd(O2R),
                             O2NPRmean = mean(O2NPR),
                             O2NPRsd = sd(O2NPR),
                             TNPmean = mean(TNP),
                             TNPsd = sd(TNP),
                             TRmean = mean(TR),
                             TRsd = sd(TR),
                             TNPRmean = mean(TNPR),
                             TNPRsd = sd(TNPR),
                             Samples = length(n))


# add confounders for which it doesn't make sense to calculate posterior distributions
# because of the lack of variability (see Oxygen.R)
tech.summary.df <- oxy %>% 
  filter(Type == "Sample") %>%
  select(ID, PAR, Pressure, Salinity) %>%
  group_by(ID, PAR) %>%
  summarise(Pmean = mean(Pressure),
            Psd = sd(Pressure),
            Smean = mean(Salinity),
            Ssd = sd(Salinity)) %>% 
  pivot_wider(names_from = PAR, 
              values_from = c(Pmean, Psd, Smean, Ssd)) %>%
  rename(PRmean = "Pmean_0", PRsd = "Psd_0", PNPmean = "Pmean_420", PNPsd = "Psd_420",
         SRmean = "Smean_0", SRsd = "Ssd_0", SNPmean = "Smean_420", SNPsd = "Ssd_420") 
# the zeros in the standard deviation demonstrate the lack of variability
# also, salinity is the same for light and dark measurements

tech.summary.df %<>% 
  select(ID, PNPmean, PRmean, SNPmean) %>%
  rename(Smean = "SNPmean") %>%
  mutate(PNPRmean = (PNPmean + PRmean)/2)

# join the two summary dataframes
photo.summary.df %<>% left_join(tech.summary.df, by = "ID")


# Net photosynthesis model based on blotted mass
NP.summary.df <- photo.summary.df %>%
  select(Tank, Treatment, Day, # main predictor variables
         NPmean, NPsd, # response variable
         O2NPmean, O2NPsd, TNPmean, TNPsd, # potential confounders with measurement error
         PNPmean, Smean, Mass) %>% # potential confounders without measurement error
  rename(O2mean = "O2NPmean", O2sd = "O2NPsd", Tmean = "TNPmean", 
         Tsd = "TNPsd", Pmean = "PNPmean")


# Prior selection
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
# but no mass conversion, but I can use our mean dry to blotted mass ratio:
r <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Meta.csv") %>%
  mutate(MR = Dry / Blotted) %>% filter(!is.na(MR)) %>% pull(MR) %>% mean()
g.mol <- 15.9994 # g mol^-1 for oxygen

a.mu <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Prior.csv") %>%
  filter(Temperature >= 15, Temperature <= 20) %>%
  # these are closest to my experimental temperature
  pull(NP) %>% mean() * r / g.mol * 1e3
a.mu # 15.39553 µmol O2 g^-1 wet mass h^-1


b.prior <- rnorm(n = 1e5, mean = b.mu, sd = 1)
# first guess: high uncertainty because there are no data for Ecklonia

ggplot() +
  geom_density(aes(x = b.prior)) +
  geom_vline(aes(xintercept = c(b.mu - 2 * 1,
                                b.mu + 2 * 1))) +
  theme_minimal()

a.prior <- rnorm(n = 1e5, mean = a.mu, sd = 5)

ggplot() +
  geom_density(aes(x = a.prior)) +
  geom_vline(aes(xintercept = c(a.mu - 2 * 5,
                                a.mu + 2 * 5))) +
  theme_minimal()


predictor <- NP.summary.df %$% seq(min(Day), max(Day), length.out = 100)
plot(NULL, xlim = c(0, 120), ylim = c(-30, 30), # base plot is better with loops
     xlab = "d", ylab = "µmol g^-1 h^-1") # empty plot
NP.summary.df %>% filter(!is.na(NPmean)) %$% abline(h = min(NPmean), lty = 2) # data range
NP.summary.df %>% filter(!is.na(NPmean)) %$% abline(h = max(NPmean), lty = 2)
abline(h = 0)
for(i in 1:1e3) curve(a.prior[i] + b.prior[i] * x, # loop
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = col.alpha("black", 0.1))
# slope and intercept are way too variable and predictions cover improbable space
# beyond the data margins

b.prior <- rnorm(n = 1e5, mean = b.mu, sd = 0.1)
# this is more probable since this prior is suspicious of positive slopes

ggplot() +
  geom_density(aes(x = b.prior)) +
  geom_vline(aes(xintercept = c(b.mu - 2 * 0.1,
                                b.mu + 2 * 0.1))) +
  theme_minimal()

a.prior <- rnorm(n = 1e5, mean = a.mu, sd = 3)
# this constrains the intercept somewhat while still allowing for discrepancies in
# Ecklonia photosynthesis between the present and previous studies

ggplot() +
  geom_density(aes(x = a.prior)) +
  geom_vline(aes(xintercept = c(a.mu - 2 * 3,
                                a.mu + 2 * 3))) +
  theme_minimal()

plot(NULL, xlim = c(0, 120), ylim = c(-30, 30), # base plot is better with loops
     xlab = "d", ylab = "µmol g^-1 h^-1") # empty plot
NP.summary.df %>% filter(!is.na(NPmean)) %$% abline(h = min(NPmean), lty = 2) # data range
NP.summary.df %>% filter(!is.na(NPmean)) %$% abline(h = max(NPmean), lty = 2)
abline(h = 0)
for(i in 1:1e3) curve(a.prior[i] + b.prior[i] * x, # loop
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = col.alpha("black", 0.1))
# still lots of improbable predictions but these priors allow for net heterotrophy
# after ~20 d (cf. Wright et al. 2022) while indicating the general negative trend
# evident from the literature

# due to the difficulty of measuring mass and photosynthesis of disintegrated samples,
# whenever a sample fully decomposed, photosynthesis was not observed (see Meta.csv)
# this means that unfortunately data are not missing at random since disintegration caused
# the inability to measure photosynthesis, making disintegration conditional upon detrital age
# and cessation of photosynthesis
# missing data conditional on the observed response variable are impossible to impute and
# it is not straightforward to model disintegration as a function of detrital age within
# the missing data model
# I have therefore decided to model the additional information contained in the disintegration
# variable separately as part of series of binomial models on photosynthetic death

# so the planned photosynthesis models require complete cases (which is already the case for photo.summary.df)
# and need to be in list format; it is also a good idea to standardise the variables for which the intercept
# is not meaningful (in this case this applies to all potential confounders)

NP.summary.l <- NP.summary.df %>%
  mutate(Tsd = Tsd / sd(Tmean),
         O2sd = O2sd / sd(O2mean),
         Tmean = standardize(Tmean),
         O2mean = standardize(O2mean),
         Pmean = standardize(Pmean),
         Smean = standardize(Smean),
         Mass = standardize(Mass)) %>%
  as.list()

NP.summary.l$N <- nrow(NP.summary.df)

# check for multicollinearity
pairs(~ Day + O2mean + Tmean + Pmean + Smean + Mass, data = NP.summary.l) # none evident

# run multilevel models with partial pooling for all categorical variables first (default)

NP.mod.pp <- ulam(
  alist(
    # response measurement error
    NPmean ~ dnorm(mean = NP, sd = NPsd), # here the observed measurement error is introduced
    vector[N]:NP ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP

    # linear model
    mu <- a + b[Treatment] * Day + # treatment effect
      br[Tank] * Day + # tank effect
      bO2 * O2[i] + bT * T_[i] + bP * Pmean + # confounders
      bS * Smean + bM * Mass,

    # predictor measurement error
    O2mean ~ dnorm(mean = O2, sd = O2sd),
    vector[N]:O2 ~ dnorm(mean = 0, sd = 1), # prior for unobserved O2
    Tmean ~ dnorm(mean = T_, sd = Tsd),
    vector[N]:T_ ~ dnorm(mean = 0, sd = 1), # prior for unobserved T

    # higher level priors for partial pooling
    b[Treatment] ~ dnorm(mean = b_, sd = sb),
    br[Tank] ~ dnorm(mean = 0, sd = sbr),

    # lower level priors
    a ~ dnorm(mean = 15.39553, sd = 3), # intercept defined by the literature
    b_ ~ dnorm(mean = -0.1374649, sd = 0.1), # slope defined by the literature
    c(bO2, bT, bP, bS, bM) ~ dnorm(mean = 0, sd = 1), # confounder slopes
    c(sigma, sb, sbr) ~ dexp(rate = 1) # standard deviations
  ),
  data = NP.summary.l, chains = 8, cores = parallel::detectCores(), iter = 1e4,
)

# note divergent transitions (~5%)
traceplot(NP.mod.pp) # healthy chains according to traceplot
trankplot(NP.mod.pp) # but trankplot highlights some issues
dev.off()
plot(precis(NP.mod.pp, depth = 1)) # including tanks explains little additional variability
options(max.print = 2000)
precis(NP.mod.pp, depth = 2) # note low variability indicated by br[1:16] and sbr
# treatments are also heavily regularised, which might exacerbate the issue of uneven missing data
# ok Rhat4 but some seriously low n_eff
# -> try non-centred priors

NP.mod.pp.nc <- ulam(
  alist(
    # response measurement error
    NPmean ~ dnorm(mean = NP, sd = NPsd), # here the observed measurement error is introduced
    vector[N]:NP ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP

    # linear model
    mu <- a + b[Treatment] * Day + # treatment effect
      br[Tank] * Day + # tank effect
      bO2 * O2[i] + bT * T_[i] + bP * Pmean + # confounders
      bS * Smean + bM * Mass,

    # predictor measurement error
    O2mean ~ dnorm(mean = O2, sd = O2sd),
    vector[N]:O2 ~ dnorm(mean = 0, sd = 1),
    Tmean ~ dnorm(mean = T_, sd = Tsd),
    vector[N]:T_ ~ dnorm(mean = 0, sd = 1),

    # define multilevel coefficients using z-scores that HMC can sample better
    save> vector[4]:b <<- b_ + zb * sb,
    save> vector[16]:br <<- 0 + zbr * sbr,

    # priors
    vector[4]:zb ~ dnorm(mean = 0, sd = 1), # z-scores
    vector[16]:zbr ~ dnorm(mean = 0, sd = 1),
    a ~ dnorm(mean = 15.39553, sd = 3), # intercept
    b_ ~ dnorm(mean = -0.1374649, sd = 0.1), # slope
    c(bO2, bT, bP, bS, bM) ~ dnorm(mean = 0, sd = 1), # confounder slopes
    c(sigma, sb, sbr) ~ dexp(rate = 1) # standard deviations
  ),
  data = NP.summary.l, chains = 8, cores = parallel::detectCores(), iter = 1e4,
)

# note lack of divergent transitions (~0%)
traceplot(NP.mod.pp.nc)
trankplot(NP.mod.pp.nc) # healthier chains
dev.off()
plot(precis(NP.mod.pp.nc, depth = 1))
options(max.print = 2000)
precis(NP.mod.pp.nc, depth = 2) # good Rhat4 and n_eff

# visual cross-validation
precis_c <- precis(NP.mod.pp, depth = 2)
precis_nc <- precis(NP.mod.pp.nc, depth = 2)

pars <- c(paste("NP[",1:57,"]",sep=""), paste("O2[",1:57,"]",sep=""), paste("T_[",1:57,"]",sep=""),
          paste("b[",1:4,"]",sep=""), paste("br[",1:16,"]",sep=""), "a" , "b_",
          "bO2", "bT", "bP", "bS", "bM", "sb", "sbr", "sigma")

neff_table <- cbind(precis_c[pars,"n_eff"] , precis_nc[pars,"n_eff"] )
dev.off()
plot( neff_table , xlim=range(neff_table) , ylim=range(neff_table) ,
      xlab="n_eff (centred)" , ylab="n_eff (non-centred)" , lwd=2 )
abline( a=0 , b=1 , lty=2 )

# non-centred model clearly wins

# treatment will be coded without partial pooling because regularisation would lead to
# even more bias towards the most complete data (treatment 2) due to non-random missing data
# keep tank partial pooling for now

NP.mod <- ulam(
  alist(
    # response measurement error
    NPmean ~ dnorm(mean = NP, sd = NPsd), # here the observed measurement error is introduced
    vector[N]:NP ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP

    # linear model
    mu <- a + b[Treatment] * Day + # treatment effect
      br[Tank] * Day + # tank effect
      bO2 * O2[i] + bT * T_[i] + bP * Pmean + # confounders
      bS * Smean + bM * Mass,

    # predictor measurement error
    O2mean ~ dnorm(mean = O2, sd = O2sd),
    vector[N]:O2 ~ dnorm(mean = 0, sd = 1),
    Tmean ~ dnorm(mean = T_, sd = Tsd),
    vector[N]:T_ ~ dnorm(mean = 0, sd = 1),

    # higher level priors for partial pooling
    br[Tank] ~ dnorm(mean = 0, sd = sbr),

    # priors
    a ~ dnorm(mean = 15.39553, sd = 3), # intercept
    b[Treatment] ~ dnorm(mean = -0.1374649, sd = 0.1), # slope
    c(bO2, bT, bP, bS, bM) ~ dnorm(mean = 0, sd = 1), # confounder slopes
    c(sigma, sbr) ~ dexp(rate = 1) # standard deviations
  ),
  data = NP.summary.l, chains = 8, cores = parallel::detectCores(), iter = 1e4,
)

# note divergent transitions (~1%)
traceplot(NP.mod) # unhealthy chains according to traceplot
trankplot(NP.mod) # and trankplot
dev.off()
plot(precis(NP.mod, depth = 1))
options(max.print = 2000)
precis(NP.mod, depth = 2) # ok Rhat4 but some low n_eff
# -> try non-centred priors

NP.mod.nc <- ulam(
  alist(
    # response measurement error
    NPmean ~ dnorm(mean = NP, sd = NPsd), # here the observed measurement error is introduced
    vector[N]:NP ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP

    # linear model
    mu <- a + b[Treatment] * Day + # treatment effect
      br[Tank] * Day + # tank effect
      bO2 * O2[i] + bT * T_[i] + bP * Pmean + # confounders
      bS * Smean + bM * Mass,

    # predictor measurement error
    O2mean ~ dnorm(mean = O2, sd = O2sd),
    vector[N]:O2 ~ dnorm(mean = 0, sd = 1),
    Tmean ~ dnorm(mean = T_, sd = Tsd),
    vector[N]:T_ ~ dnorm(mean = 0, sd = 1),

    # define multilevel coefficients using z-scores that HMC can sample better
    save> vector[16]:br <<- 0 + zbr * sbr,

    # priors
    vector[16]:zbr ~ dnorm(mean = 0, sd = 1), # z-score
    a ~ dnorm(mean = 15.39553, sd = 3), # intercept
    b[Treatment] ~ dnorm(mean = -0.1374649, sd = 0.1), # slopes
    c(bO2, bT, bP, bS, bM) ~ dnorm(mean = 0, sd = 1), # confounder slopes
    c(sigma, sbr) ~ dexp(rate = 1) # standard deviations
  ),
  data = NP.summary.l, chains = 8, cores = parallel::detectCores(), iter = 1e4,
)

# note lack of divergent transitions (~0%)
traceplot(NP.mod.nc)
trankplot(NP.mod.nc) # healthier chains
dev.off()
plot(precis(NP.mod.nc, depth = 1))
options(max.print = 2000)
precis(NP.mod.nc, depth = 2) # good Rhat4 and n_eff

# visual cross-validation
precis_c <- precis(NP.mod, depth = 2)
precis_nc <- precis(NP.mod.nc, depth = 2)

pars <- c(paste("NP[",1:57,"]",sep=""), paste("O2[",1:57,"]",sep=""), paste("T_[",1:57,"]",sep=""),
          paste("b[",1:4,"]",sep=""), paste("br[",1:16,"]",sep=""), "a",
          "bO2", "bT", "bP", "bS", "bM", "sbr", "sigma")

neff_table <- cbind(precis_c[pars,"n_eff"] , precis_nc[pars,"n_eff"] )
dev.off()
plot( neff_table , xlim=range(neff_table) , ylim=range(neff_table) ,
      xlab="n_eff (centred)" , ylab="n_eff (non-centred)" , lwd=2 )
abline( a=0 , b=1 , lty=2 )

# non-centred model clearly wins


# Prior-posterior comparison
NP.mod.nc.posterior <- extract.samples(NP.mod.nc)
NP.mod.nc.prior <- extract.prior(NP.mod.nc, n = 2e4) # for some reason prior simulation went wrong
# -> manually simulate for parameters of interest (a, b[1:4], bO2, bT, bP, bS, bM, O2[1:57], T_[1:57])

NP.mod.nc.prior.df <- data.frame(n = rep(1:1e4, 124),
                                 coefficient = factor(c(rep("a", 1e4), rep("b", 1e4*4),
                                                        rep(c("bO2", "bT", "bP", "bS", "bM"), each = 1e4),
                                                        rep(c("O2", "T_"), each = 1e4*57))),
                                 level = factor(c(rep(NA, 1e4), rep(1:4, each = 1e4), rep(NA, 5*1e4),
                                                  rep(rep(1:57, each = 1e4), 2))),
                                 sample = c(rnorm(n = 1e4, mean = 15.39553, sd = 3),
                                            rep(rnorm(n = 1e4, mean = -0.1374649, sd = 0.1), 4),
                                            rep(rnorm(n = 1e4, mean = 0, sd = 1), 5+57*2)))

NP.mod.nc.posterior.df <- data.frame(NP.mod.nc.posterior) %>%
  rownames_to_column(var = "n") %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "level")) %>%
  mutate(level = as.integer(level)) %>%
  arrange(coefficient, level, n) %>%
  mutate(level = factor(level), coefficient = factor(coefficient))


NP.mod.nc.prior.posterior.df <- data.frame(rbind(NP.mod.nc.prior.df,
                                                 NP.mod.nc.posterior.df %>%
                                                   filter(!coefficient %in% c("NP", "zbr", "sigma", "sbr", "br")) %>%
                                                   droplevels()),
                                           distribution = rep(c("prior", "posterior"), each = 1240000)) %>%
  mutate(distribution = factor(distribution))

require(ggh4x)
ggplot(data = NP.mod.nc.prior.posterior.df) +
  geom_density(aes(x = sample, colour = distribution)) +
  facet_nested_wrap(~ coefficient + level, scales = "free",
                    nest_line = element_line()) +
  theme_minimal() +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")


# Posterior prediction

# Contrasts
brsim <- NP.mod.nc.posterior %$% rnorm(n = length(sbr), mean = 0, sd = sbr)
NP.comp <- NP.mod.nc.posterior.df %>%
  filter(coefficient == "b") %>%
  mutate(br = rep(brsim, 4),
         b = sample + br) %>% # add variation across tanks to posterior slope
  select(-c(sample, br)) %>%
  pivot_wider(names_from = c(coefficient, level), values_from = b) %>%
  mutate(b1_2 = b_1 - b_2, b1_3 = b_1 - b_3, b1_4 = b_1 - b_4,
         b2_3 = b_2 - b_3, b2_4 = b_2 - b_4, b3_4 = b_3 - b_4)

NP.comp %>%
  filter(b1_2 < 0) %>% pull(b1_2) %>% length() /
  NP.comp %>% pull(b1_2) %>% length()
# 92% probability that b1 < b2

NP.comp %>%
  filter(b1_3 < 0) %>% pull(b1_3) %>% length() /
  NP.comp %>% pull(b1_3) %>% length()
# 58% probability that b1 < b3

NP.comp %>%
  filter(b1_4 > 0) %>% pull(b1_4) %>% length() /
  NP.comp %>% pull(b1_4) %>% length()
# 52% probability that b1 > b4

NP.comp %>%
  filter(b2_3 > 0) %>% pull(b2_3) %>% length() /
  NP.comp %>% pull(b2_3) %>% length()
# 86% probability that b2 > b3

NP.comp %>%
  filter(b2_4 > 0) %>% pull(b2_4) %>% length() /
  NP.comp %>% pull(b2_4) %>% length()
# 86% probability that b2 > b4

NP.comp %>%
  filter(b3_4 > 0) %>% pull(b3_4) %>% length() /
  NP.comp %>% pull(b3_4) %>% length()
# 58% probability that b3 > b4

# Probability mass below zero
NP.comp %>%
  filter(b_1 < 0) %>% pull(b_1) %>% length() /
  NP.comp %>% pull(b_1) %>% length()
# 98% probability that b1 < 0

NP.comp %>%
  filter(b_2 < 0) %>% pull(b_2) %>% length() /
  NP.comp %>% pull(b_2) %>% length()
# 91% probability that b2 < 0

NP.comp %>%
  filter(b_3 < 0) %>% pull(b_3) %>% length() /
  NP.comp %>% pull(b_3) %>% length()
# 96% probability that b3 < 0

NP.comp %>%
  filter(b_4 < 0) %>% pull(b_4) %>% length() /
  NP.comp %>% pull(b_4) %>% length()
# 95% probability that b4 < 0

# Coefficients
NP.b <- NP.mod.nc.posterior.df %>%
  filter(coefficient %in% c("b", "bT", "bO2", "bP", "bS", "bM")) %>%
  mutate(br = rep(brsim, 9),
         b = sample + br) %>%
  group_by(coefficient, level) %>%
  summarise(mean = mean(b),
            sd = sd(b))

NP.coef <- NP.mod.nc.posterior.df %>%
  group_by(coefficient, level) %>%
  summarise(mean = mean(sample),
            sd = sd(sample)) %>%
  left_join(NP.b %>%
              rename(meanr = "mean", sdr = "sd"), 
            by = c("coefficient", "level"))
# a = 6.4 ± 0.67
# b1 = -0.12 ± 0.06
# b2 = -0.05 ± 0.04
# b3 = -0.11 ± 0.06
# b4 = -0.12 ± 0.08

# unstandardised bO2 =
NP.coef %>% filter(coefficient == "bO2") %>%
  pull(mean) / sd(NP.summary.df$O2mean) # 0.08
# unstandardised bT =
NP.coef %>% filter(coefficient == "bT") %>%
  pull(mean) / sd(NP.summary.df$Tmean) # 0.36
# unstandardised bP =
NP.coef %>% filter(coefficient == "bP") %>%
  pull(mean) / sd(NP.summary.df$Pmean) # 0.1
# unstandardised bS =
NP.coef %>% filter(coefficient == "bS") %>%
  pull(mean) / sd(NP.summary.df$Smean) # -0.21
# unstandardised bM =
NP.coef %>% filter(coefficient == "bM") %>%
  pull(mean) / sd(NP.summary.df$Mass) # 0.6

# unstandardised aO2 =



# Equation

NP.annotation <- NP.summary.df %>%
  group_by(Treatment) %>%
  summarise(n = length(Treatment))
str(NP.annotation)

NP.annotation %<>%
  mutate(n = as.character(c(expression(italic("n ")*"= 15"),
                            expression(italic("n ")*"= 22"),
                            expression(italic("n ")*"= 13"),
                            expression(italic("n ")*"= 7"))),
         equation = as.character(c(expression(italic("y ")*"= –0.12"*italic("x ")*"+ 6.4"),
                                   expression(italic("y ")*"= –0.05"*italic("x ")*"+ 6.4"),
                                   expression(italic("y ")*"= –0.11"*italic("x ")*"+ 6.4"),
                                   expression(italic("y ")*"= –0.12"*italic("x ")*"+ 6.4"))))

# Lines and intervals
predictor.df <- data.frame(Day = c(seq(0, 37, length.out = 200),
                                   seq(0, 119, length.out = 200),
                                   seq(0, 37, length.out = 200),
                                   seq(0, 27, length.out = 200)),
                           Treatment = rep(1:4, each = 200),
                           O2mean = 0, Tmean = 0, Pmean = 0, Smean = 0, Mass = 0)
# note that all confounders are set to 0, which is their mean because they are standardised

NP.mod.nc.mu <- link(NP.mod.nc, data = predictor.df) # regular link function fails when
# Tank is ignored

# # calculate posterior prediction for the average Tank, i.e. ignoring Tank
# link_ <- function(Treatment, Day){
#   mu <- with(NP.mod.nc.posterior, a + b[,Treatment] * Day) # mean Tank effect is 0
#   return(mu)
# }
# 
# NP.mod.nc.mu <- predictor.df %$% mapply(link_, Treatment, Day) 

# calculate posterior prediction for new simulated Tanks, i.e. incorporating Tank s.d.
link_sim <- function(Treatment, Day){
  mu <- with(NP.mod.nc.posterior, a + b[,Treatment] * Day + brsim * Day)
  return(mu)
}

NP.mod.nc.mu <- predictor.df %$% mapply(link_sim, Treatment, Day) # use mapply instead of sapply
# because there are two explanatory variables

# summarise posterior probability distributions
predictor.df$NPmu.mean <- apply(NP.mod.nc.mu, 2, mean)
predictor.df$NPmu.pi.lwr <- t(apply(NP.mod.nc.mu, 2, PI, prob = 0.5))[,1]
predictor.df$NPmu.pi.upr <- t(apply(NP.mod.nc.mu, 2, PI, prob = 0.5))[,2]
predictor.df$NPmu.pi.lwr2 <- t(apply(NP.mod.nc.mu, 2, PI, prob = 0.8))[,1]
predictor.df$NPmu.pi.upr2 <- t(apply(NP.mod.nc.mu, 2, PI, prob = 0.8))[,2]
predictor.df$NPmu.pi.lwr3 <- t(apply(NP.mod.nc.mu, 2, PI, prob = 0.9))[,1]
predictor.df$NPmu.pi.upr3 <- t(apply(NP.mod.nc.mu, 2, PI, prob = 0.9))[,2]

# predictor.df$NPmu.median <- apply(NP.mod.nc.mu, 2, median)
# predictor.df$NPmu.mode <- apply(NP.mod.nc.mu, 2, chainmode)
# predictor.df$NPmu.hpdi.lwr <- t(apply(NP.mod.nc.mu, 2, HPDI, prob = 0.98))[,1]
# predictor.df$NPmu.hpdi.upr <- t(apply(NP.mod.nc.mu, 2, HPDI, prob = 0.98))[,2]


mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .5, .2, .2),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15, hjust = 0),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.key.size = unit(.3, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.spacing.x = unit(.1, "cm"),
                 legend.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title.align = 0,
                 legend.title = element_text(size = 12, face = "bold"),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 15, hjust = 0),
                 panel.spacing = unit(.5, "cm"),
                 text = element_text(family = "Helvetica Neue"))

Treatment_names <- c(
  '1' = "Dark 15°C",
  '2' = "Light 15°C",
  '3' = "Light 20°C",
  '4' = "Light 25°C"
)

# plot percentile intervals for unknown Tanks
NPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df, 
              aes(Day, NP, fill = factor(Treatment), colour = factor(Treatment), group = ID),
              alpha = 0.2, position = "identity", width = 8) +
  geom_line(data = predictor.df, aes(Day, NPmu.mean, colour = factor(Treatment))) +
  geom_ribbon(data = predictor.df, 
              aes(Day, ymin = NPmu.pi.lwr, ymax = NPmu.pi.upr, fill = factor(Treatment)),
              alpha = 0.5) +
  geom_ribbon(data = predictor.df, 
              aes(Day, ymin = NPmu.pi.lwr2, ymax = NPmu.pi.upr2, fill = factor(Treatment)),
              alpha = 0.4) +
  geom_ribbon(data = predictor.df, 
              aes(Day, ymin = NPmu.pi.lwr3, ymax = NPmu.pi.upr3, fill = factor(Treatment)),
              alpha = 0.3) +
  geom_text(data = NP.annotation, aes(0, -20, label = n),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = 0) +
  geom_text(data = NP.annotation, aes(0, -26, label = equation),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = 0) +
  scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  facet_grid(~Treatment, space = "free_x", scales = "free_x", labeller = labeller(Treatment = Treatment_names)) +
  labs(y = expression(italic(P[n])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = "Detrital age (d)") +
  scale_x_continuous(breaks = seq(0, 120, by = 10)) +
  coord_cartesian(ylim = c(-30, 20), expand = FALSE, clip = "off") +
  mytheme
NPp


# plot confounders
# calculate posterior prediction for new simulated Tanks, i.e. incorporating Tank s.d.
predictor.O2.df <- NP.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             O2mean = seq(min(O2mean), max(O2mean), length.out = 200)) %>%
  mutate(O2mean.us = O2mean * sd(NP.summary.df$O2mean) + mean(NP.summary.df$O2mean)) # unstandardised

link_O2 <- function(Day, O2mean){
  mu <- with(NP.mod.nc.posterior, a + brsim * Day + bO2 * O2mean)
  return(mu)
}

NP.mod.nc.mu.O2 <- predictor.O2.df %$% mapply(link_O2, Day, O2mean)


predictor.O2.df$mu.mean <- apply(NP.mod.nc.mu.O2, 2, mean)
predictor.O2.df$mu.pi.lwr <- t(apply(NP.mod.nc.mu.O2, 2, PI, prob = 0.5))[,1]
predictor.O2.df$mu.pi.upr <- t(apply(NP.mod.nc.mu.O2, 2, PI, prob = 0.5))[,2]
predictor.O2.df$mu.pi.lwr2 <- t(apply(NP.mod.nc.mu.O2, 2, PI, prob = 0.8))[,1]
predictor.O2.df$mu.pi.upr2 <- t(apply(NP.mod.nc.mu.O2, 2, PI, prob = 0.8))[,2]
predictor.O2.df$mu.pi.lwr3 <- t(apply(NP.mod.nc.mu.O2, 2, PI, prob = 0.9))[,1]
predictor.O2.df$mu.pi.upr3 <- t(apply(NP.mod.nc.mu.O2, 2, PI, prob = 0.9))[,2]

require(ggdensity)
O2NPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_hdr(data = photo.density.df, aes(O2NP, NP, group = ID), n = 500, method = "mvnorm", probs = 0.999,
           fill = "#363538", alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.O2.df, aes(O2mean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.O2.df, aes(O2mean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.O2.df, aes(O2mean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.O2.df, aes(O2mean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[n])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = expression("Initial O"[2]*" ("*mu*"M)")) +
  scale_x_continuous(breaks = seq(150, 350, by = 100)) +
  coord_cartesian(xlim = c(150, 350), ylim = c(-30, 20), clip = "off", expand = FALSE) +
  mytheme
O2NPp




predictor.T.df <- NP.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Tmean = seq(min(Tmean), max(Tmean), length.out = 200)) %>%
  mutate(Tmean.us = Tmean * sd(NP.summary.df$Tmean) + mean(NP.summary.df$Tmean)) # unstandardised

link_T <- function(Day, Tmean){
  mu <- with(NP.mod.nc.posterior, a + brsim * Day + bT * Tmean)
  return(mu)
}

NP.mod.nc.mu.T <- predictor.T.df %$% mapply(link_T, Day, Tmean)


predictor.T.df$mu.mean <- apply(NP.mod.nc.mu.T, 2, mean)
predictor.T.df$mu.pi.lwr <- t(apply(NP.mod.nc.mu.T, 2, PI, prob = 0.5))[,1]
predictor.T.df$mu.pi.upr <- t(apply(NP.mod.nc.mu.T, 2, PI, prob = 0.5))[,2]
predictor.T.df$mu.pi.lwr2 <- t(apply(NP.mod.nc.mu.T, 2, PI, prob = 0.8))[,1]
predictor.T.df$mu.pi.upr2 <- t(apply(NP.mod.nc.mu.T, 2, PI, prob = 0.8))[,2]
predictor.T.df$mu.pi.lwr3 <- t(apply(NP.mod.nc.mu.T, 2, PI, prob = 0.9))[,1]
predictor.T.df$mu.pi.upr3 <- t(apply(NP.mod.nc.mu.T, 2, PI, prob = 0.9))[,2]

TNPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_hdr(data = photo.density.df, aes(TNP, NP, group = ID), n = 500, method = "mvnorm", probs = 0.999,
           fill = "#363538", alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.T.df, aes(Tmean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.T.df, aes(Tmean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.T.df, aes(Tmean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.T.df, aes(Tmean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[n])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = "Temperature (°C)") +
  scale_x_continuous(breaks = seq(16.5, 18.5, by = 1)) +
  coord_cartesian(xlim = c(16.5, 18.5), ylim = c(-30, 20), clip = "off", expand = FALSE) +
  mytheme
TNPp



predictor.P.df <- NP.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Pmean = seq(min(Pmean), max(Pmean), length.out = 200)) %>%
  mutate(Pmean.us = Pmean * sd(NP.summary.df$Pmean) + mean(NP.summary.df$Pmean)) # unstandardised

link_P <- function(Day, Pmean){
  mu <- with(NP.mod.nc.posterior, a + brsim * Day + bP * Pmean)
  return(mu)
}

NP.mod.nc.mu.P <- predictor.P.df %$% mapply(link_P, Day, Pmean)


predictor.P.df$mu.mean <- apply(NP.mod.nc.mu.P, 2, mean)
predictor.P.df$mu.pi.lwr <- t(apply(NP.mod.nc.mu.P, 2, PI, prob = 0.5))[,1]
predictor.P.df$mu.pi.upr <- t(apply(NP.mod.nc.mu.P, 2, PI, prob = 0.5))[,2]
predictor.P.df$mu.pi.lwr2 <- t(apply(NP.mod.nc.mu.P, 2, PI, prob = 0.8))[,1]
predictor.P.df$mu.pi.upr2 <- t(apply(NP.mod.nc.mu.P, 2, PI, prob = 0.8))[,2]
predictor.P.df$mu.pi.lwr3 <- t(apply(NP.mod.nc.mu.P, 2, PI, prob = 0.9))[,1]
predictor.P.df$mu.pi.upr3 <- t(apply(NP.mod.nc.mu.P, 2, PI, prob = 0.9))[,2]

PNPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df %>%
                left_join(photo.summary.df %>%
                            select(ID, PNPmean), by = "ID"),
              aes(PNPmean, NP, group = ID), fill = "#363538", width = 1,
              alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.P.df, aes(Pmean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.P.df, aes(Pmean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.P.df, aes(Pmean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.P.df, aes(Pmean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[n])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = "Pressure (hPa)") +
  scale_x_continuous(breaks = seq(1005, 1035, by = 10)) +
  coord_cartesian(xlim = c(1005, 1035), ylim = c(-30, 20), clip = "off", expand = FALSE) +
  mytheme
PNPp



predictor.S.df <- NP.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Smean = seq(min(Smean), max(Smean), length.out = 200)) %>%
  mutate(Smean.us = Smean * sd(NP.summary.df$Smean) + mean(NP.summary.df$Smean)) # unstandardised

link_S <- function(Day, Smean){
  mu <- with(NP.mod.nc.posterior, a + brsim * Day + bS * Smean)
  return(mu)
}

NP.mod.nc.mu.S <- predictor.S.df %$% mapply(link_S, Day, Smean)


predictor.S.df$mu.mean <- apply(NP.mod.nc.mu.S, 2, mean)
predictor.S.df$mu.pi.lwr <- t(apply(NP.mod.nc.mu.S, 2, PI, prob = 0.5))[,1]
predictor.S.df$mu.pi.upr <- t(apply(NP.mod.nc.mu.S, 2, PI, prob = 0.5))[,2]
predictor.S.df$mu.pi.lwr2 <- t(apply(NP.mod.nc.mu.S, 2, PI, prob = 0.8))[,1]
predictor.S.df$mu.pi.upr2 <- t(apply(NP.mod.nc.mu.S, 2, PI, prob = 0.8))[,2]
predictor.S.df$mu.pi.lwr3 <- t(apply(NP.mod.nc.mu.S, 2, PI, prob = 0.9))[,1]
predictor.S.df$mu.pi.upr3 <- t(apply(NP.mod.nc.mu.S, 2, PI, prob = 0.9))[,2]

SNPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df %>%
                left_join(photo.summary.df %>%
                            select(ID, Smean), by = "ID"),
              aes(Smean, NP, group = ID), fill = "#363538", width = 0.2,
              alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.S.df, aes(Smean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.S.df, aes(Smean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.S.df, aes(Smean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.S.df, aes(Smean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[n])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = "Salinity (‰)") +
  scale_x_continuous(breaks = seq(33, 36, by = 1)) +
  coord_cartesian(xlim = c(33, 36), ylim = c(-30, 20), clip = "off", expand = FALSE) +
  mytheme
SNPp



predictor.M.df <- NP.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Mass = seq(min(Mass), max(Mass), length.out = 200)) %>%
  mutate(Mass.us = Mass * sd(NP.summary.df$Mass) + mean(NP.summary.df$Mass)) # unstandardised

link_M <- function(Day, Mass){
  mu <- with(NP.mod.nc.posterior, a + brsim * Day + bM * Mass)
  return(mu)
}

NP.mod.nc.mu.M <- predictor.M.df %$% mapply(link_M, Day, Mass)


predictor.M.df$mu.mean <- apply(NP.mod.nc.mu.M, 2, mean)
predictor.M.df$mu.pi.lwr <- t(apply(NP.mod.nc.mu.M, 2, PI, prob = 0.5))[,1]
predictor.M.df$mu.pi.upr <- t(apply(NP.mod.nc.mu.M, 2, PI, prob = 0.5))[,2]
predictor.M.df$mu.pi.lwr2 <- t(apply(NP.mod.nc.mu.M, 2, PI, prob = 0.8))[,1]
predictor.M.df$mu.pi.upr2 <- t(apply(NP.mod.nc.mu.M, 2, PI, prob = 0.8))[,2]
predictor.M.df$mu.pi.lwr3 <- t(apply(NP.mod.nc.mu.M, 2, PI, prob = 0.9))[,1]
predictor.M.df$mu.pi.upr3 <- t(apply(NP.mod.nc.mu.M, 2, PI, prob = 0.9))[,2]

MNPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df, aes(Mass, NP, group = ID), fill = "#363538", width = 0.2,
              alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.M.df, aes(Mass.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.M.df, aes(Mass.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.M.df, aes(Mass.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.M.df, aes(Mass.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[n])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = "Mass (g)") +
  scale_x_continuous(breaks = seq(0, 5, by = 1)) +
  coord_cartesian(xlim = c(0, 5), ylim = c(-30, 20), clip = "off", expand = FALSE) +
  mytheme
MNPp







# gross photosynthesis model based on blotted mass
GP.summary.df <- photo.summary.df %>%
  select(Tank, Treatment, Day, # main predictor variables
         GPmean, GPsd, # response variable
         O2NPRmean, O2NPRsd, TNPRmean, TNPRsd, # potential confounders with measurement error
         PNPRmean, Smean, Mass) %>% # potential confounders without measurement error
  rename(O2mean = "O2NPRmean", O2sd = "O2NPRsd", Tmean = "TNPRmean",
         Tsd = "TNPRsd", Pmean = "PNPRmean")

# Prior selection
# Wright et al. 2022 (doi: 10.1111/gcb.16299) provide detrital photosynthesis over detrital age
# slopes for three Laminaria species: 0.007, -0.016 and -0.035 mg C g^-1 dry mass h^-1 d^-1
# Frontier et al. 2021 (doi: 10.1016/j.marenvres.2021.105277) is the only similarly useful prior
# study but unfortunately does not provide any slopes due to categorical analysis
r <- 0.127003814 # Wright et al. 2022 also provide this mean dry to wet mass ratio
g.mol <- 12.0107 # g mol^-1 for carbon

# note that somewhat counter-intuitively their rates per dry mass need to be multiplied by r to get
# rates per wet mass because a given dry mass photosynthesises more than the same amount of wet mass
b.mu <- mean(c(0.007, -0.016, -0.035)) * r / g.mol * 1e3
b.mu # -0.1550886 µmol O2 g^-1 wet mass h^-1 d^-1

# Staehr & Wernberg 2009 and Wernberg et al. 2016 provide baseline dry mass photosynthesis rates
# but no mass conversion, but I can use our mean dry to blotted mass ratio:
r <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Meta.csv") %>%
  mutate(MR = Dry / Blotted) %>% filter(!is.na(MR)) %>% pull(MR) %>% mean()
g.mol <- 15.9994 # g mol^-1 for oxygen

a.mu <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Prior.csv") %>%
  mutate(GP = NP + R) %>%
  filter(Temperature >= 15, Temperature <= 20) %>%
  # these are closest to my experimental temperature
  pull(GP) %>% mean() * r / g.mol * 1e3
a.mu # 22.30601 µmol O2 g^-1 wet mass h^-1


b.prior <- rnorm(n = 1e5, mean = b.mu, sd = 0.1)

ggplot() +
  geom_density(aes(x = b.prior)) +
  geom_vline(aes(xintercept = c(b.mu - 2 * 0.1,
                                b.mu + 2 * 0.1))) +
  theme_minimal()

a.prior <- rnorm(n = 1e5, mean = a.mu, sd = 4)
# more uncertainty because more parameters (NP + R) and greater mean

ggplot() +
  geom_density(aes(x = a.prior)) +
  geom_vline(aes(xintercept = c(a.mu - 2 * 4,
                                a.mu + 2 * 4))) +
  theme_minimal()

plot(NULL, xlim = c(0, 120), ylim = c(-10, 30), # base plot is better with loops
     xlab = "d", ylab = "µmol g^-1 h^-1") # empty plot
GP.summary.df %>% filter(!is.na(GPmean)) %$% abline(h = min(GPmean), lty = 2) # data range
GP.summary.df %>% filter(!is.na(GPmean)) %$% abline(h = max(GPmean), lty = 2)
abline(h = 0)
for(i in 1:1e3) curve(a.prior[i] + b.prior[i] * x, # loop
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = col.alpha("black", 0.1))
# still lots of improbable predictions (GP theoretically cannot be < 0) but I do not want to
# go narrower considering the narrow range of GP compared to estimates from the literature


GP.summary.l <- GP.summary.df %>%
  mutate(Tsd = Tsd / sd(Tmean),
         O2sd = O2sd / sd(O2mean),
         Tmean = standardize(Tmean),
         O2mean = standardize(O2mean),
         Pmean = standardize(Pmean),
         Smean = standardize(Smean),
         Mass = standardize(Mass)) %>%
  as.list()

GP.summary.l$N <- nrow(GP.summary.df)

# check for multicollinearity
pairs(~ Day + O2mean + Tmean + Pmean + Smean + Mass, data = GP.summary.l) # none evident


GP.mod <- ulam(
  alist(
    # response measurement error
    GPmean ~ dnorm(mean = GP, sd = GPsd), # here the observed measurement error is introduced
    vector[N]:GP ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP

    # linear model
    mu <- a + b[Treatment] * Day + # treatment effect
      br[Tank] * Day + # tank effect
      bO2 * O2[i] + bT * T_[i] + bP * Pmean + # confounders
      bS * Smean + bM * Mass,

    # predictor measurement error
    O2mean ~ dnorm(mean = O2, sd = O2sd),
    vector[N]:O2 ~ dnorm(mean = 0, sd = 1),
    Tmean ~ dnorm(mean = T_, sd = Tsd),
    vector[N]:T_ ~ dnorm(mean = 0, sd = 1),

    # higher level priors for partial pooling
    br[Tank] ~ dnorm(mean = 0, sd = sbr),

    # priors
    a ~ dnorm(mean = 22.30601, sd = 4), # intercept
    b[Treatment] ~ dnorm(mean = -0.1550886, sd = 0.1), # slope
    c(bO2, bT, bP, bS, bM) ~ dnorm(mean = 0, sd = 1), # confounder slopes
    c(sigma, sbr) ~ dexp(rate = 1) # standard deviations
  ),
  data = GP.summary.l, chains = 8, cores = parallel::detectCores(), iter = 1e4,
)

# note divergent transitions (~3%)
traceplot(GP.mod) # unhealthy chains according to traceplot
trankplot(GP.mod) # and trankplot
dev.off()
plot(precis(GP.mod, depth = 1))
options(max.print = 2000)
precis(GP.mod, depth = 2) # ok Rhat4 but some low n_eff
# -> try non-centred priors

GP.mod.nc <- ulam(
  alist(
    # response measurement error
    GPmean ~ dnorm(mean = GP, sd = GPsd), # here the observed measurement error is introduced
    vector[N]:GP ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP

    # linear model
    mu <- a + b[Treatment] * Day + # treatment effect
      br[Tank] * Day + # tank effect
      bO2 * O2[i] + bT * T_[i] + bP * Pmean + # confounders
      bS * Smean + bM * Mass,

    # predictor measurement error
    O2mean ~ dnorm(mean = O2, sd = O2sd),
    vector[N]:O2 ~ dnorm(mean = 0, sd = 1),
    Tmean ~ dnorm(mean = T_, sd = Tsd),
    vector[N]:T_ ~ dnorm(mean = 0, sd = 1),

    # define multilevel coefficients using z-scores that HMC can sample better
    save> vector[16]:br <<- 0 + zbr * sbr,

    # priors
    vector[16]:zbr ~ dnorm(mean = 0, sd = 1), # z-score
    a ~ dnorm(mean = 22.30601, sd = 4), # intercept
    b[Treatment] ~ dnorm(mean = -0.1550886, sd = 0.1), # slopes
    c(bO2, bT, bP, bS, bM) ~ dnorm(mean = 0, sd = 1), # confounder slopes
    c(sigma, sbr) ~ dexp(rate = 1) # standard deviations
  ),
  data = GP.summary.l, chains = 8, cores = parallel::detectCores(), iter = 1e4,
)

# note lack of divergent transitions (~0%)
traceplot(GP.mod.nc)
trankplot(GP.mod.nc) # healthier chains
dev.off()
plot(precis(GP.mod.nc, depth = 1))
options(max.print = 2000)
precis(GP.mod.nc, depth = 2) # good Rhat4 and n_eff

# visual cross-validation
precis_c <- precis(GP.mod, depth = 2)
precis_nc <- precis(GP.mod.nc, depth = 2)

pars <- c(paste("GP[",1:57,"]",sep=""), paste("O2[",1:57,"]",sep=""), paste("T_[",1:57,"]",sep=""),
          paste("b[",1:4,"]",sep=""), paste("br[",1:16,"]",sep=""), "a",
          "bO2", "bT", "bP", "bS", "bM", "sbr", "sigma")

neff_table <- cbind(precis_c[pars,"n_eff"] , precis_nc[pars,"n_eff"] )
dev.off()
plot( neff_table , xlim=range(neff_table) , ylim=range(neff_table) ,
      xlab="n_eff (centred)" , ylab="n_eff (non-centred)" , lwd=2 )
abline( a=0 , b=1 , lty=2 )

# non-centred model clearly wins


# Prior-posterior comparison
GP.mod.nc.posterior <- extract.samples(GP.mod.nc)
# manually simulate priors for parameters of interest (a, b[1:4], bO2, bT, bP, bS, bM, O2[1:57], T_[1:57])

GP.mod.nc.prior.df <- data.frame(n = rep(1:1e4, 124),
                                 coefficient = factor(c(rep("a", 1e4), rep("b", 1e4*4),
                                                        rep(c("bO2", "bT", "bP", "bS", "bM"), each = 1e4),
                                                        rep(c("O2", "T_"), each = 1e4*57))),
                                 level = factor(c(rep(NA, 1e4), rep(1:4, each = 1e4), rep(NA, 5*1e4),
                                                  rep(rep(1:57, each = 1e4), 2))),
                                 sample = c(rnorm(n = 1e4, mean = 22.30601, sd = 4),
                                            rep(rnorm(n = 1e4, mean = -0.1550886, sd = 0.1), 4),
                                            rep(rnorm(n = 1e4, mean = 0, sd = 1), 5+57*2)))

GP.mod.nc.posterior.df <- data.frame(GP.mod.nc.posterior) %>%
  rownames_to_column(var = "n") %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "level")) %>%
  mutate(level = as.integer(level)) %>%
  arrange(coefficient, level, n) %>%
  mutate(level = factor(level), coefficient = factor(coefficient))


GP.mod.nc.prior.posterior.df <- data.frame(rbind(GP.mod.nc.prior.df,
                                                 GP.mod.nc.posterior.df %>%
                                                   filter(!coefficient %in% c("GP", "zbr", "sigma", "sbr", "br")) %>%
                                                   droplevels()),
                                           distribution = rep(c("prior", "posterior"), each = 1240000)) %>%
  mutate(distribution = factor(distribution))

ggplot(data = GP.mod.nc.prior.posterior.df) +
  geom_density(aes(x = sample, colour = distribution)) +
  facet_nested_wrap(~ coefficient + level, scales = "free",
                    nest_line = element_line()) +
  theme_minimal() +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")


# Posterior prediction
# Contrasts
brsim <- GP.mod.nc.posterior %$% rnorm(n = length(sbr), mean = 0, sd = sbr)
GP.comp <- GP.mod.nc.posterior.df %>%
  filter(coefficient == "b") %>%
  mutate(br = rep(brsim, 4),
         b = sample + br) %>% # add variation across tanks to posterior slope
  select(-c(sample, br)) %>%
  pivot_wider(names_from = c(coefficient, level), values_from = b) %>%
  mutate(b1_2 = b_1 - b_2, b1_3 = b_1 - b_3, b1_4 = b_1 - b_4,
         b2_3 = b_2 - b_3, b2_4 = b_2 - b_4, b3_4 = b_3 - b_4)


GP.comp %>%
  filter(b1_2 < 0) %>% pull(b1_2) %>% length() /
  GP.comp %>% pull(b1_2) %>% length()
# 88% probability that b1 < b2

GP.comp %>%
  filter(b1_3 > 0) %>% pull(b1_3) %>% length() /
  GP.comp %>% pull(b1_3) %>% length()
# 74% probability that b1 > b3

GP.comp %>%
  filter(b1_4 > 0) %>% pull(b1_4) %>% length() /
  GP.comp %>% pull(b1_4) %>% length()
# 69% probability that b1 > b4

GP.comp %>%
  filter(b2_3 > 0) %>% pull(b2_3) %>% length() /
  GP.comp %>% pull(b2_3) %>% length()
# 96% probability that b2 > b3

GP.comp %>%
  filter(b2_4 > 0) %>% pull(b2_4) %>% length() /
  GP.comp %>% pull(b2_4) %>% length()
# 91% probability that b2 > b4

GP.comp %>%
  filter(b3_4 < 0) %>% pull(b3_4) %>% length() /
  GP.comp %>% pull(b3_4) %>% length()
# 53% probability that b3 < b4

# Probability mass below zero
GP.comp %>%
  filter(b_1 < 0) %>% pull(b_1) %>% length() /
  GP.comp %>% pull(b_1) %>% length()
# 95% probability that b1 < 0

GP.comp %>%
  filter(b_2 < 0) %>% pull(b_2) %>% length() /
  GP.comp %>% pull(b_2) %>% length()
# 87% probability that b2 < 0

GP.comp %>%
  filter(b_3 < 0) %>% pull(b_3) %>% length() /
  GP.comp %>% pull(b_3) %>% length()
# 98% probability that b3 < 0

GP.comp %>%
  filter(b_4 < 0) %>% pull(b_4) %>% length() /
  GP.comp %>% pull(b_4) %>% length()
# 96% probability that b4 < 0





# Coefficients
GP.b <- GP.mod.nc.posterior.df %>%
  filter(coefficient %in% c("b", "bT", "bO2", "bP", "bS", "bM")) %>%
  mutate(br = rep(brsim, 9),
         b = sample + br) %>%
  group_by(coefficient, level) %>%
  summarise(mean = mean(b),
            sd = sd(b))

GP.coef <- GP.mod.nc.posterior.df %>%
  group_by(coefficient, level) %>%
  summarise(mean = mean(sample),
            sd = sd(sample)) %>%
  left_join(GP.b %>%
              rename(meanr = "mean", sdr = "sd"), 
            by = c("coefficient", "level"))
# a = 8.9 ± 0.58
# b1 = -0.08 ± 0.05
# b2 = -0.03 ± 0.04
# b3 = -0.12 ± 0.05
# b4 = -0.11 ± 0.07

# unstandardised bO2 =
GP.coef %>% filter(coefficient == "bO2") %>%
  pull(mean) / sd(GP.summary.df$O2mean) # 0.05
# unstandardised bT =
GP.coef %>% filter(coefficient == "bT") %>%
  pull(mean) / sd(GP.summary.df$Tmean) # 0.25
# unstandardised bP =
GP.coef %>% filter(coefficient == "bP") %>%
  pull(mean) / sd(GP.summary.df$Pmean) # 0.06
# unstandardised bS =
GP.coef %>% filter(coefficient == "bS") %>%
  pull(mean) / sd(GP.summary.df$Smean) # -0.22
# unstandardised bM =
GP.coef %>% filter(coefficient == "bM") %>%
  pull(mean) / sd(GP.summary.df$Mass) # 0.44

# unstandardised aO2 =





GP.annotation <- GP.summary.df %>%
  group_by(Treatment) %>%
  summarise(n = length(Treatment))
str(GP.annotation)

GP.annotation %<>%
  mutate(n = as.character(c(expression(italic("n ")*"= 15"),
                            expression(italic("n ")*"= 22"),
                            expression(italic("n ")*"= 13"),
                            expression(italic("n ")*"= 7"))),
         equation = as.character(c(expression(italic("y ")*"= –0.08"*italic("x ")*"+ 8.9"),
                                   expression(italic("y ")*"= –0.03"*italic("x ")*"+ 8.9"),
                                   expression(italic("y ")*"= –0.12"*italic("x ")*"+ 8.9"),
                                   expression(italic("y ")*"= –0.11"*italic("x ")*"+ 8.9"))))







# Lines and intervals
# calculate posterior prediction for new simulated Tanks, i.e. incorporating Tank s.d.
link_sim <- function(Treatment, Day){
  mu <- with(GP.mod.nc.posterior, a + b[,Treatment] * Day + brsim * Day)
  return(mu)
}

GP.mod.nc.mu <- predictor.df %$% mapply(link_sim, Treatment, Day)

# summarise posterior probability distributions
predictor.df$GPmu.mean <- apply(GP.mod.nc.mu, 2, mean)
predictor.df$GPmu.pi.lwr <- t(apply(GP.mod.nc.mu, 2, PI, prob = 0.5))[,1]
predictor.df$GPmu.pi.upr <- t(apply(GP.mod.nc.mu, 2, PI, prob = 0.5))[,2]
predictor.df$GPmu.pi.lwr2 <- t(apply(GP.mod.nc.mu, 2, PI, prob = 0.8))[,1]
predictor.df$GPmu.pi.upr2 <- t(apply(GP.mod.nc.mu, 2, PI, prob = 0.8))[,2]
predictor.df$GPmu.pi.lwr3 <- t(apply(GP.mod.nc.mu, 2, PI, prob = 0.9))[,1]
predictor.df$GPmu.pi.upr3 <- t(apply(GP.mod.nc.mu, 2, PI, prob = 0.9))[,2]


GPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df, 
              aes(Day, GP, colour = factor(Treatment), fill = factor(Treatment), group = ID),
              alpha = 0.2, position = "identity", width = 8) +
  geom_line(data = predictor.df, aes(Day, GPmu.mean, colour = factor(Treatment))) +
  geom_ribbon(data = predictor.df, 
              aes(Day, ymin = GPmu.pi.lwr, ymax = GPmu.pi.upr, fill = factor(Treatment)),
              alpha = 0.5) +
  geom_ribbon(data = predictor.df, 
              aes(Day, ymin = GPmu.pi.lwr2, ymax = GPmu.pi.upr2, fill = factor(Treatment)),
              alpha = 0.4) +
  geom_ribbon(data = predictor.df, 
              aes(Day, ymin = GPmu.pi.lwr3, ymax = GPmu.pi.upr3, fill = factor(Treatment)),
              alpha = 0.3) +
  geom_text(data = GP.annotation, aes(0, -12, label = n),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = 0) +
  geom_text(data = GP.annotation, aes(0, -16.5, label = equation),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = 0) +
  scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  facet_grid(~Treatment, space = "free_x", scales = "free_x", labeller = labeller(Treatment = Treatment_names)) +
  labs(y = expression(italic(P[g])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = "Detrital age (d)") +
  scale_x_continuous(breaks = seq(0, 120, by = 10)) +
  coord_cartesian(ylim = c(-20, 20), expand = FALSE, clip = "off") +
  mytheme
GPp


# plot confounders
# calculate posterior prediction for new simulated Tanks, i.e. incorporating Tank s.d.
predictor.O2.df <- GP.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             O2mean = seq(min(O2mean), max(O2mean), length.out = 200)) %>%
  mutate(O2mean.us = O2mean * sd(GP.summary.df$O2mean) + mean(GP.summary.df$O2mean)) # unstandardised

link_O2 <- function(Day, O2mean){
  mu <- with(GP.mod.nc.posterior, a + brsim * Day + bO2 * O2mean)
  return(mu)
}

GP.mod.nc.mu.O2 <- predictor.O2.df %$% mapply(link_O2, Day, O2mean)


predictor.O2.df$mu.mean <- apply(GP.mod.nc.mu.O2, 2, mean)
predictor.O2.df$mu.pi.lwr <- t(apply(GP.mod.nc.mu.O2, 2, PI, prob = 0.5))[,1]
predictor.O2.df$mu.pi.upr <- t(apply(GP.mod.nc.mu.O2, 2, PI, prob = 0.5))[,2]
predictor.O2.df$mu.pi.lwr2 <- t(apply(GP.mod.nc.mu.O2, 2, PI, prob = 0.8))[,1]
predictor.O2.df$mu.pi.upr2 <- t(apply(GP.mod.nc.mu.O2, 2, PI, prob = 0.8))[,2]
predictor.O2.df$mu.pi.lwr3 <- t(apply(GP.mod.nc.mu.O2, 2, PI, prob = 0.9))[,1]
predictor.O2.df$mu.pi.upr3 <- t(apply(GP.mod.nc.mu.O2, 2, PI, prob = 0.9))[,2]

O2GPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_hdr(data = photo.density.df, aes(O2NPR, GP, group = ID), n = 500, method = "mvnorm", probs = 0.999,
           fill = "#363538", alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.O2.df, aes(O2mean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.O2.df, aes(O2mean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.O2.df, aes(O2mean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.O2.df, aes(O2mean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[g])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = expression("Initial O"[2]*" ("*mu*"M)")) +
  scale_x_continuous(breaks = seq(150, 350, by = 100)) +
  coord_cartesian(xlim = c(150, 350), ylim = c(-20, 20), clip = "off", expand = FALSE) +
  mytheme
O2GPp




predictor.T.df <- GP.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Tmean = seq(min(Tmean), max(Tmean), length.out = 200)) %>%
  mutate(Tmean.us = Tmean * sd(GP.summary.df$Tmean) + mean(GP.summary.df$Tmean)) # unstandardised

link_T <- function(Day, Tmean){
  mu <- with(GP.mod.nc.posterior, a + brsim * Day + bT * Tmean)
  return(mu)
}

GP.mod.nc.mu.T <- predictor.T.df %$% mapply(link_T, Day, Tmean)


predictor.T.df$mu.mean <- apply(GP.mod.nc.mu.T, 2, mean)
predictor.T.df$mu.pi.lwr <- t(apply(GP.mod.nc.mu.T, 2, PI, prob = 0.5))[,1]
predictor.T.df$mu.pi.upr <- t(apply(GP.mod.nc.mu.T, 2, PI, prob = 0.5))[,2]
predictor.T.df$mu.pi.lwr2 <- t(apply(GP.mod.nc.mu.T, 2, PI, prob = 0.8))[,1]
predictor.T.df$mu.pi.upr2 <- t(apply(GP.mod.nc.mu.T, 2, PI, prob = 0.8))[,2]
predictor.T.df$mu.pi.lwr3 <- t(apply(GP.mod.nc.mu.T, 2, PI, prob = 0.9))[,1]
predictor.T.df$mu.pi.upr3 <- t(apply(GP.mod.nc.mu.T, 2, PI, prob = 0.9))[,2]

TGPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_hdr(data = photo.density.df, aes(TNPR, GP, group = ID), n = 500, method = "mvnorm", probs = 0.999,
           fill = "#363538", alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.T.df, aes(Tmean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.T.df, aes(Tmean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.T.df, aes(Tmean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.T.df, aes(Tmean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[g])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = "Temperature (°C)") +
  scale_x_continuous(breaks = seq(16.5, 18.5, by = 1)) +
  coord_cartesian(xlim = c(16.5, 18.5), ylim = c(-20, 20), clip = "off", expand = FALSE) +
  mytheme
TGPp



predictor.P.df <- GP.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Pmean = seq(min(Pmean), max(Pmean), length.out = 200)) %>%
  mutate(Pmean.us = Pmean * sd(GP.summary.df$Pmean) + mean(GP.summary.df$Pmean)) # unstandardised

link_P <- function(Day, Pmean){
  mu <- with(GP.mod.nc.posterior, a + brsim * Day + bP * Pmean)
  return(mu)
}

GP.mod.nc.mu.P <- predictor.P.df %$% mapply(link_P, Day, Pmean)


predictor.P.df$mu.mean <- apply(GP.mod.nc.mu.P, 2, mean)
predictor.P.df$mu.pi.lwr <- t(apply(GP.mod.nc.mu.P, 2, PI, prob = 0.5))[,1]
predictor.P.df$mu.pi.upr <- t(apply(GP.mod.nc.mu.P, 2, PI, prob = 0.5))[,2]
predictor.P.df$mu.pi.lwr2 <- t(apply(GP.mod.nc.mu.P, 2, PI, prob = 0.8))[,1]
predictor.P.df$mu.pi.upr2 <- t(apply(GP.mod.nc.mu.P, 2, PI, prob = 0.8))[,2]
predictor.P.df$mu.pi.lwr3 <- t(apply(GP.mod.nc.mu.P, 2, PI, prob = 0.9))[,1]
predictor.P.df$mu.pi.upr3 <- t(apply(GP.mod.nc.mu.P, 2, PI, prob = 0.9))[,2]

PGPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df %>%
                left_join(photo.summary.df %>%
                            select(ID, PNPRmean), by = "ID"),
              aes(PNPRmean, GP, group = ID), fill = "#363538", width = 1,
              alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.P.df, aes(Pmean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.P.df, aes(Pmean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.P.df, aes(Pmean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.P.df, aes(Pmean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[g])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = "Pressure (hPa)") +
  scale_x_continuous(breaks = seq(1005, 1035, by = 10)) +
  coord_cartesian(xlim = c(1005, 1035), ylim = c(-20, 20), clip = "off", expand = FALSE) +
  mytheme
PGPp



predictor.S.df <- GP.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Smean = seq(min(Smean), max(Smean), length.out = 200)) %>%
  mutate(Smean.us = Smean * sd(GP.summary.df$Smean) + mean(GP.summary.df$Smean)) # unstandardised

link_S <- function(Day, Smean){
  mu <- with(GP.mod.nc.posterior, a + brsim * Day + bS * Smean)
  return(mu)
}

GP.mod.nc.mu.S <- predictor.S.df %$% mapply(link_S, Day, Smean)


predictor.S.df$mu.mean <- apply(GP.mod.nc.mu.S, 2, mean)
predictor.S.df$mu.pi.lwr <- t(apply(GP.mod.nc.mu.S, 2, PI, prob = 0.5))[,1]
predictor.S.df$mu.pi.upr <- t(apply(GP.mod.nc.mu.S, 2, PI, prob = 0.5))[,2]
predictor.S.df$mu.pi.lwr2 <- t(apply(GP.mod.nc.mu.S, 2, PI, prob = 0.8))[,1]
predictor.S.df$mu.pi.upr2 <- t(apply(GP.mod.nc.mu.S, 2, PI, prob = 0.8))[,2]
predictor.S.df$mu.pi.lwr3 <- t(apply(GP.mod.nc.mu.S, 2, PI, prob = 0.9))[,1]
predictor.S.df$mu.pi.upr3 <- t(apply(GP.mod.nc.mu.S, 2, PI, prob = 0.9))[,2]

SGPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df %>%
                left_join(photo.summary.df %>%
                            select(ID, Smean), by = "ID"),
              aes(Smean, GP, group = ID), fill = "#363538", width = 0.2,
              alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.S.df, aes(Smean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.S.df, aes(Smean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.S.df, aes(Smean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.S.df, aes(Smean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[g])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = "Salinity (‰)") +
  scale_x_continuous(breaks = seq(33, 36, by = 1)) +
  coord_cartesian(xlim = c(33, 36), ylim = c(-20, 20), clip = "off", expand = FALSE) +
  mytheme
SGPp



predictor.M.df <- GP.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Mass = seq(min(Mass), max(Mass), length.out = 200)) %>%
  mutate(Mass.us = Mass * sd(GP.summary.df$Mass) + mean(GP.summary.df$Mass)) # unstandardised

link_M <- function(Day, Mass){
  mu <- with(GP.mod.nc.posterior, a + brsim * Day + bM * Mass)
  return(mu)
}

GP.mod.nc.mu.M <- predictor.M.df %$% mapply(link_M, Day, Mass)


predictor.M.df$mu.mean <- apply(GP.mod.nc.mu.M, 2, mean)
predictor.M.df$mu.pi.lwr <- t(apply(GP.mod.nc.mu.M, 2, PI, prob = 0.5))[,1]
predictor.M.df$mu.pi.upr <- t(apply(GP.mod.nc.mu.M, 2, PI, prob = 0.5))[,2]
predictor.M.df$mu.pi.lwr2 <- t(apply(GP.mod.nc.mu.M, 2, PI, prob = 0.8))[,1]
predictor.M.df$mu.pi.upr2 <- t(apply(GP.mod.nc.mu.M, 2, PI, prob = 0.8))[,2]
predictor.M.df$mu.pi.lwr3 <- t(apply(GP.mod.nc.mu.M, 2, PI, prob = 0.9))[,1]
predictor.M.df$mu.pi.upr3 <- t(apply(GP.mod.nc.mu.M, 2, PI, prob = 0.9))[,2]

MGPp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df, aes(Mass, GP, group = ID), fill = "#363538", width = 0.2,
              alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.M.df, aes(Mass.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.M.df, aes(Mass.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.M.df, aes(Mass.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.M.df, aes(Mass.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[g])*" ("*mu*"mol O"[2]*" g"^-1*" h"^-1*")"),
       x = "Mass (g)") +
  scale_x_continuous(breaks = seq(0, 5, by = 1)) +
  coord_cartesian(xlim = c(0, 5), ylim = c(-20, 20), clip = "off", expand = FALSE) +
  mytheme
MGPp







# daily net photosynthesis model based on blotted mass
NPD.summary.df <- photo.summary.df %>%
  select(Tank, Treatment, Day, # main predictor variables
         NPDmean, NPDsd, # response variable
         O2NPRmean, O2NPRsd, TNPRmean, TNPRsd, # potential confounders with measurement error
         PNPRmean, Smean, Mass) %>% # potential confounders without measurement error
  rename(O2mean = "O2NPRmean", O2sd = "O2NPRsd", Tmean = "TNPRmean",
         Tsd = "TNPRsd", Pmean = "PNPRmean")

# Prior selection
# Wright et al. 2022 (doi: 10.1111/gcb.16299) provide detrital photosynthesis over detrital age
# slopes for three Laminaria species for net (-0.001, -0.009, -0.029) and gross (0.007, -0.016 and -0.035 
# mg C g^-1 dry mass h^-1 d^-1) photosynthesis
# daily net photosynthesis can be estimated by subtracting net from gross photosynthesis and the subtracting
# the difference (respiration) from net photosynthesis before multiplying by 12
r <- 0.127003814 # Wright et al. 2022 also provide this mean dry to wet mass ratio
g.mol <- 12.0107 # g mol^-1 for carbon

# note that somewhat counter-intuitively their rates per dry mass need to be multiplied by r to get
# rates per wet mass because a given dry mass photosynthesises more than the same amount of wet mass
b.mu <- mean(c(-0.001, -0.009, -0.029) - (c(0.007, -0.016, -0.035) - c(-0.001, -0.009, -0.029))) * 
  12 * r / g.mol * 1e3
b.mu # -1.438094 µmol O2 g^-1 wet mass d^-1 d^-1

# Staehr & Wernberg 2009 and Wernberg et al. 2016 provide baseline dry mass photosynthesis rates
# but no mass conversion, but I can use our mean dry to blotted mass ratio:
r <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Meta.csv") %>%
  mutate(MR = Dry / Blotted) %>% filter(!is.na(MR)) %>% pull(MR) %>% mean()
g.mol <- 15.9994 # g mol^-1 for oxygen

a.mu <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Prior.csv") %>%
  mutate(NPD = (NP - R) * 12) %>%
  filter(Temperature >= 15, Temperature <= 20) %>%
  # these are closest to my experimental temperature
  pull(NPD) %>% mean() * r / g.mol * 1e3
a.mu # 101.8206 µmol O2 g^-1 wet mass d^-1


b.prior <- rnorm(n = 1e5, mean = b.mu, sd = 2)

ggplot() +
  geom_density(aes(x = b.prior)) +
  geom_vline(aes(xintercept = c(b.mu - 2 * 2,
                                b.mu + 2 * 2))) +
  theme_minimal()

a.prior <- rnorm(n = 1e5, mean = a.mu, sd = 30)
# more uncertainty because more parameters (NP + R) and greater mean

ggplot() +
  geom_density(aes(x = a.prior)) +
  geom_vline(aes(xintercept = c(a.mu - 2 * 30,
                                a.mu + 2 * 30))) +
  theme_minimal()

plot(NULL, xlim = c(0, 120), ylim = c(-400, 200), # base plot is better with loops
     xlab = "d", ylab = "µmol g^-1 h^-1") # empty plot
NPD.summary.df %>% filter(!is.na(NPDmean)) %$% abline(h = min(NPDmean), lty = 2) # data range
NPD.summary.df %>% filter(!is.na(NPDmean)) %$% abline(h = max(NPDmean), lty = 2)
abline(h = 0)
for(i in 1:1e3) curve(a.prior[i] + b.prior[i] * x, # loop
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = col.alpha("black", 0.1))
# still lots of improbable predictions but narrow enough


NPD.summary.l <- NPD.summary.df %>%
  mutate(Tsd = Tsd / sd(Tmean),
         O2sd = O2sd / sd(O2mean),
         Tmean = standardize(Tmean),
         O2mean = standardize(O2mean),
         Pmean = standardize(Pmean),
         Smean = standardize(Smean),
         Mass = standardize(Mass)) %>%
  as.list()

NPD.summary.l$N <- nrow(NPD.summary.df)

# check for multicollinearity
pairs(~ Day + O2mean + Tmean + Pmean + Smean + Mass, data = NPD.summary.l) # none evident


NPD.mod <- ulam(
  alist(
    # response measurement error
    NPDmean ~ dnorm(mean = NPD, sd = NPDsd), # here the observed measurement error is introduced
    vector[N]:NPD ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP

    # linear model
    mu <- a + b[Treatment] * Day + # treatment effect
      br[Tank] * Day + # tank effect
      bO2 * O2[i] + bT * T_[i] + bP * Pmean + # confounders
      bS * Smean + bM * Mass,

    # predictor measurement error
    O2mean ~ dnorm(mean = O2, sd = O2sd),
    vector[N]:O2 ~ dnorm(mean = 0, sd = 1),
    Tmean ~ dnorm(mean = T_, sd = Tsd),
    vector[N]:T_ ~ dnorm(mean = 0, sd = 1),

    # higher level priors for partial pooling
    br[Tank] ~ dnorm(mean = 0, sd = sbr),

    # priors
    a ~ dnorm(mean = 101.8206, sd = 30), # intercept
    b[Treatment] ~ dnorm(mean = -1.438094, sd = 2), # slope
    c(bO2, bT, bP, bS, bM) ~ dnorm(mean = 0, sd = 15), # confounder slopes
    c(sigma, sbr) ~ dexp(rate = 1) # standard deviations
  ),
  data = NPD.summary.l, chains = 8, cores = parallel::detectCores(), iter = 1e4,
)

# note lack of divergent transitions (~0%)
traceplot(NPD.mod) # fairly healthy chains according to traceplot
trankplot(NPD.mod) # and trankplot
dev.off()
plot(precis(NPD.mod, depth = 1))
options(max.print = 2000)
precis(NPD.mod, depth = 2) # good Rhat4 but some low n_eff
# -> try non-centred priors

NPD.mod.nc <- ulam(
  alist(
    # response measurement error
    NPDmean ~ dnorm(mean = NPD, sd = NPDsd), # here the observed measurement error is introduced
    vector[N]:NPD ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP

    # linear model
    mu <- a + b[Treatment] * Day + # treatment effect
      br[Tank] * Day + # tank effect
      bO2 * O2[i] + bT * T_[i] + bP * Pmean + # confounders
      bS * Smean + bM * Mass,

    # predictor measurement error
    O2mean ~ dnorm(mean = O2, sd = O2sd),
    vector[N]:O2 ~ dnorm(mean = 0, sd = 1),
    Tmean ~ dnorm(mean = T_, sd = Tsd),
    vector[N]:T_ ~ dnorm(mean = 0, sd = 1),

    # define multilevel coefficients using z-scores that HMC can sample better
    save> vector[16]:br <<- 0 + zbr * sbr,

    # priors
    vector[16]:zbr ~ dnorm(mean = 0, sd = 1), # z-score
    a ~ dnorm(mean = 101.8206, sd = 30), # intercept
    b[Treatment] ~ dnorm(mean = -1.438094, sd = 2), # slopes
    c(bO2, bT, bP, bS, bM) ~ dnorm(mean = 0, sd = 15), # confounder slopes
    c(sigma, sbr) ~ dexp(rate = 1) # standard deviations
  ),
  data = NPD.summary.l, chains = 8, cores = parallel::detectCores(), iter = 1e4,
)

# note complete lack of divergent transitions
traceplot(NPD.mod.nc)
trankplot(NPD.mod.nc) # healthier chains
dev.off()
plot(precis(NPD.mod.nc, depth = 1))
options(max.print = 2000)
precis(NPD.mod.nc, depth = 2) # good Rhat4 and better n_eff

# visual cross-validation
precis_c <- precis(NPD.mod, depth = 2)
precis_nc <- precis(NPD.mod.nc, depth = 2)

pars <- c(paste("NPD[",1:57,"]",sep=""), paste("O2[",1:57,"]",sep=""), paste("T_[",1:57,"]",sep=""),
          paste("b[",1:4,"]",sep=""), paste("br[",1:16,"]",sep=""), "a",
          "bO2", "bT", "bP", "bS", "bM", "sbr", "sigma")

neff_table <- cbind(precis_c[pars,"n_eff"] , precis_nc[pars,"n_eff"] )
dev.off()
plot( neff_table , xlim=range(neff_table) , ylim=range(neff_table) ,
      xlab="n_eff (centred)" , ylab="n_eff (non-centred)" , lwd=2 )
abline( a=0 , b=1 , lty=2 )

# non-centred model clearly wins!


# Prior-posterior comparison
NPD.mod.nc.posterior <- extract.samples(NPD.mod.nc)
# manually simulate priors for parameters of interest (a, b[1:4], bO2, bT, bP, bS, bM, O2[1:57], T_[1:57])

NPD.mod.nc.prior.df <- data.frame(n = rep(1:1e4, 124),
                                 coefficient = factor(c(rep("a", 1e4), rep("b", 1e4*4),
                                                        rep(c("bO2", "bT", "bP", "bS", "bM"), each = 1e4),
                                                        rep(c("O2", "T_"), each = 1e4*57))),
                                 level = factor(c(rep(NA, 1e4), rep(1:4, each = 1e4), rep(NA, 5*1e4),
                                                  rep(rep(1:57, each = 1e4), 2))),
                                 sample = c(rnorm(n = 1e4, mean = 101.8206, sd = 30),
                                            rep(rnorm(n = 1e4, mean = -1.438094, sd = 2), 4),
                                            rep(rnorm(n = 1e4, mean = 0, sd = 15), 5),
                                            rep(rnorm(n = 1e4, mean = 0, sd = 1), 57*2)))

NPD.mod.nc.posterior.df <- data.frame(NPD.mod.nc.posterior) %>%
  rownames_to_column(var = "n") %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "level")) %>%
  mutate(level = as.integer(level)) %>%
  arrange(coefficient, level, n) %>%
  mutate(level = factor(level), coefficient = factor(coefficient))


NPD.mod.nc.prior.posterior.df <- data.frame(rbind(NPD.mod.nc.prior.df,
                                                 NPD.mod.nc.posterior.df %>%
                                                   filter(!coefficient %in% c("NPD", "zbr", "sigma", "sbr", "br")) %>%
                                                   droplevels()),
                                           distribution = rep(c("prior", "posterior"), each = 1240000)) %>%
  mutate(distribution = factor(distribution))

ggplot(data = NPD.mod.nc.prior.posterior.df) +
  geom_density(aes(x = sample, colour = distribution)) +
  facet_nested_wrap(~ coefficient + level, scales = "free",
                    nest_line = element_line()) +
  theme_minimal() +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")


# Posterior prediction

# Contrasts
brsim <- NPD.mod.nc.posterior %$% rnorm(n = length(sbr), mean = 0, sd = sbr)
NPD.comp <- NPD.mod.nc.posterior.df %>%
  filter(coefficient == "b") %>%
  mutate(br = rep(brsim, 4),
         b = sample + br) %>% # add variation across tanks to posterior slope
  select(-c(sample, br)) %>%
  pivot_wider(names_from = c(coefficient, level), values_from = b) %>%
  mutate(b1_2 = b_1 - b_2, b1_3 = b_1 - b_3, b1_4 = b_1 - b_4,
         b2_3 = b_2 - b_3, b2_4 = b_2 - b_4, b3_4 = b_3 - b_4)


NPD.comp %>%
  filter(b1_2 < 0) %>% pull(b1_2) %>% length() /
  NPD.comp %>% pull(b1_2) %>% length()
# 96% probability that b1 < b2

NPD.comp %>%
  filter(b1_3 < 0) %>% pull(b1_3) %>% length() /
  NPD.comp %>% pull(b1_3) %>% length()
# 93% probability that b1 < b3

NPD.comp %>%
  filter(b1_4 < 0) %>% pull(b1_4) %>% length() /
  NPD.comp %>% pull(b1_4) %>% length()
# 70% probability that b1 < b4

NPD.comp %>%
  filter(b2_3 < 0) %>% pull(b2_3) %>% length() /
  NPD.comp %>% pull(b2_3) %>% length()
# 56% probability that b2 < b3

NPD.comp %>%
  filter(b2_4 > 0) %>% pull(b2_4) %>% length() /
  NPD.comp %>% pull(b2_4) %>% length()
# 73% probability that b2 > b4

NPD.comp %>%
  filter(b3_4 > 0) %>% pull(b3_4) %>% length() /
  NPD.comp %>% pull(b3_4) %>% length()
# 74% probability that b3 > b4

# Probability mass below zero
NPD.comp %>%
  filter(b_1 < 0) %>% pull(b_1) %>% length() /
  NPD.comp %>% pull(b_1) %>% length()
# 98% probability that b1 < 0

NPD.comp %>%
  filter(b_2 < 0) %>% pull(b_2) %>% length() /
  NPD.comp %>% pull(b_2) %>% length()
# 85% probability that b2 < 0

NPD.comp %>%
  filter(b_3 < 0) %>% pull(b_3) %>% length() /
  NPD.comp %>% pull(b_3) %>% length()
# 69% probability that b3 < 0

NPD.comp %>%
  filter(b_4 < 0) %>% pull(b_4) %>% length() /
  NPD.comp %>% pull(b_4) %>% length()
# 86% probability that b4 < 0



# Coefficients
NPD.b <- NPD.mod.nc.posterior.df %>%
  filter(coefficient %in% c("b", "bT", "bO2", "bP", "bS", "bM")) %>%
  mutate(br = rep(brsim, 9),
         b = sample + br) %>%
  group_by(coefficient, level) %>%
  summarise(mean = mean(b),
            sd = sd(b))

NPD.coef <- NPD.mod.nc.posterior.df %>%
  group_by(coefficient, level) %>%
  summarise(mean = mean(sample),
            sd = sd(sample)) %>%
  left_join(NPD.b %>%
              rename(meanr = "mean", sdr = "sd"), 
            by = c("coefficient", "level"))
# a = 38.34 ± 7.97
# b1 = -1.68 ± 0.87
# b2 = -0.54 ± 0.7
# b3 = -0.45 ± 0.92
# b4 = -1.12 ± 1.09

# unstandardised bO2 =
NPD.coef %>% filter(coefficient == "bO2") %>%
  pull(mean) / sd(NPD.summary.df$O2mean) # 1.55
# unstandardised bT =
NPD.coef %>% filter(coefficient == "bT") %>%
  pull(mean) / sd(NPD.summary.df$Tmean) # -2.84
# unstandardised bP =
NPD.coef %>% filter(coefficient == "bP") %>%
  pull(mean) / sd(NPD.summary.df$Pmean) # 0.33
# unstandardised bS =
NPD.coef %>% filter(coefficient == "bS") %>%
  pull(mean) / sd(NPD.summary.df$Smean) # -0.1
# unstandardised bM =
NPD.coef %>% filter(coefficient == "bM") %>%
  pull(mean) / sd(NPD.summary.df$Mass) # 16.41

# unstandardised aO2 =





NPD.annotation <- NPD.summary.df %>%
  group_by(Treatment) %>%
  summarise(n = length(Treatment))
str(NPD.annotation)

NPD.annotation %<>%
  mutate(n = as.character(c(expression(italic("n ")*"= 15"),
                            expression(italic("n ")*"= 22"),
                            expression(italic("n ")*"= 13"),
                            expression(italic("n ")*"= 7"))),
         equation = as.character(c(expression(italic("y ")*"= –1.68"*italic("x ")*"+ 38.34"),
                                   expression(italic("y ")*"= –0.54"*italic("x ")*"+ 38.34"),
                                   expression(italic("y ")*"= –0.45"*italic("x ")*"+ 38.34"),
                                   expression(italic("y ")*"= –1.12"*italic("x ")*"+ 38.34"))))











# Lines and intervals
# calculate posterior prediction for new simulated Tanks, i.e. incorporating Tank s.d.

link_sim <- function(Treatment, Day){
  mu <- with(NPD.mod.nc.posterior, a + b[,Treatment] * Day + brsim * Day)
  return(mu)
}

NPD.mod.nc.mu <- predictor.df %$% mapply(link_sim, Treatment, Day)

# summarise posterior probability distributions
predictor.df$NPDmu.mean <- apply(NPD.mod.nc.mu, 2, mean)
predictor.df$NPDmu.pi.lwr <- t(apply(NPD.mod.nc.mu, 2, PI, prob = 0.5))[,1]
predictor.df$NPDmu.pi.upr <- t(apply(NPD.mod.nc.mu, 2, PI, prob = 0.5))[,2]
predictor.df$NPDmu.pi.lwr2 <- t(apply(NPD.mod.nc.mu, 2, PI, prob = 0.8))[,1]
predictor.df$NPDmu.pi.upr2 <- t(apply(NPD.mod.nc.mu, 2, PI, prob = 0.8))[,2]
predictor.df$NPDmu.pi.lwr3 <- t(apply(NPD.mod.nc.mu, 2, PI, prob = 0.9))[,1]
predictor.df$NPDmu.pi.upr3 <- t(apply(NPD.mod.nc.mu, 2, PI, prob = 0.9))[,2]


NPDp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df, 
              aes(Day, NPD, colour = factor(Treatment), fill = factor(Treatment), group = ID),
              alpha = 0.2, position = "identity", width = 8) +
  geom_line(data = predictor.df, aes(Day, NPDmu.mean, colour = factor(Treatment))) +
  geom_ribbon(data = predictor.df, 
              aes(Day, ymin = NPDmu.pi.lwr, ymax = NPDmu.pi.upr, fill = factor(Treatment)),
              alpha = 0.5) +
  geom_ribbon(data = predictor.df, 
              aes(Day, ymin = NPDmu.pi.lwr2, ymax = NPDmu.pi.upr2, fill = factor(Treatment)),
              alpha = 0.4) +
  geom_ribbon(data = predictor.df, 
              aes(Day, ymin = NPDmu.pi.lwr3, ymax = NPDmu.pi.upr3, fill = factor(Treatment)),
              alpha = 0.3) +
  geom_text(data = NPD.annotation, aes(0, -440, label = n),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = 0) +
  geom_text(data = NPD.annotation, aes(0, -535, label = equation),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = 0) +
  scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  facet_grid(~Treatment, space = "free_x", scales = "free_x", labeller = labeller(Treatment = Treatment_names)) +
  labs(y = expression(italic(P[d])*" ("*mu*"mol O"[2]*" g"^-1*" d"^-1*")"),
       x = "Detrital age (d)") +
  scale_x_continuous(breaks = seq(0, 120, by = 10)) +
  coord_cartesian(ylim = c(-600, 200), expand = FALSE, clip = "off") +
  mytheme
NPDp


# plot confounders
# calculate posterior prediction for new simulated Tanks, i.e. incorporating Tank s.d.
predictor.O2.df <- NPD.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             O2mean = seq(min(O2mean), max(O2mean), length.out = 200)) %>%
  mutate(O2mean.us = O2mean * sd(NPD.summary.df$O2mean) + mean(NPD.summary.df$O2mean)) # unstandardised

link_O2 <- function(Day, O2mean){
  mu <- with(NPD.mod.nc.posterior, a + brsim * Day + bO2 * O2mean)
  return(mu)
}

NPD.mod.nc.mu.O2 <- predictor.O2.df %$% mapply(link_O2, Day, O2mean)


predictor.O2.df$mu.mean <- apply(NPD.mod.nc.mu.O2, 2, mean)
predictor.O2.df$mu.pi.lwr <- t(apply(NPD.mod.nc.mu.O2, 2, PI, prob = 0.5))[,1]
predictor.O2.df$mu.pi.upr <- t(apply(NPD.mod.nc.mu.O2, 2, PI, prob = 0.5))[,2]
predictor.O2.df$mu.pi.lwr2 <- t(apply(NPD.mod.nc.mu.O2, 2, PI, prob = 0.8))[,1]
predictor.O2.df$mu.pi.upr2 <- t(apply(NPD.mod.nc.mu.O2, 2, PI, prob = 0.8))[,2]
predictor.O2.df$mu.pi.lwr3 <- t(apply(NPD.mod.nc.mu.O2, 2, PI, prob = 0.9))[,1]
predictor.O2.df$mu.pi.upr3 <- t(apply(NPD.mod.nc.mu.O2, 2, PI, prob = 0.9))[,2]

O2NPDp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_hdr(data = photo.density.df, aes(O2NPR, NPD, group = ID), n = 500, method = "mvnorm", probs = 0.999,
           fill = "#363538", alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.O2.df, aes(O2mean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.O2.df, aes(O2mean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.O2.df, aes(O2mean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.O2.df, aes(O2mean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[d])*" ("*mu*"mol O"[2]*" g"^-1*" d"^-1*")"),
       x = expression("Initial O"[2]*" ("*mu*"M)")) +
  scale_x_continuous(breaks = seq(150, 350, by = 100)) +
  scale_y_continuous(breaks = seq(-600, 200, by = 200)) +
  coord_cartesian(xlim = c(150, 350), ylim = c(-600, 200), clip = "off", expand = FALSE) +
  mytheme
O2NPDp




predictor.T.df <- NPD.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Tmean = seq(min(Tmean), max(Tmean), length.out = 200)) %>%
  mutate(Tmean.us = Tmean * sd(NPD.summary.df$Tmean) + mean(NPD.summary.df$Tmean)) # unstandardised

link_T <- function(Day, Tmean){
  mu <- with(NPD.mod.nc.posterior, a + brsim * Day + bT * Tmean)
  return(mu)
}

NPD.mod.nc.mu.T <- predictor.T.df %$% mapply(link_T, Day, Tmean)


predictor.T.df$mu.mean <- apply(NPD.mod.nc.mu.T, 2, mean)
predictor.T.df$mu.pi.lwr <- t(apply(NPD.mod.nc.mu.T, 2, PI, prob = 0.5))[,1]
predictor.T.df$mu.pi.upr <- t(apply(NPD.mod.nc.mu.T, 2, PI, prob = 0.5))[,2]
predictor.T.df$mu.pi.lwr2 <- t(apply(NPD.mod.nc.mu.T, 2, PI, prob = 0.8))[,1]
predictor.T.df$mu.pi.upr2 <- t(apply(NPD.mod.nc.mu.T, 2, PI, prob = 0.8))[,2]
predictor.T.df$mu.pi.lwr3 <- t(apply(NPD.mod.nc.mu.T, 2, PI, prob = 0.9))[,1]
predictor.T.df$mu.pi.upr3 <- t(apply(NPD.mod.nc.mu.T, 2, PI, prob = 0.9))[,2]

TNPDp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_hdr(data = photo.density.df, aes(TNPR, NPD, group = ID), n = 500, method = "mvnorm", probs = 0.999,
           fill = "#363538", alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.T.df, aes(Tmean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.T.df, aes(Tmean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.T.df, aes(Tmean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.T.df, aes(Tmean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[d])*" ("*mu*"mol O"[2]*" g"^-1*" d"^-1*")"),
       x = "Temperature (°C)") +
  scale_x_continuous(breaks = seq(16.5, 18.5, by = 1)) +
  scale_y_continuous(breaks = seq(-600, 200, by = 200)) +
  coord_cartesian(xlim = c(16.5, 18.5), ylim = c(-600, 200), clip = "off", expand = FALSE) +
  mytheme
TNPDp



predictor.P.df <- NPD.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Pmean = seq(min(Pmean), max(Pmean), length.out = 200)) %>%
  mutate(Pmean.us = Pmean * sd(NPD.summary.df$Pmean) + mean(NPD.summary.df$Pmean)) # unstandardised

link_P <- function(Day, Pmean){
  mu <- with(NPD.mod.nc.posterior, a + brsim * Day + bP * Pmean)
  return(mu)
}

NPD.mod.nc.mu.P <- predictor.P.df %$% mapply(link_P, Day, Pmean)


predictor.P.df$mu.mean <- apply(NPD.mod.nc.mu.P, 2, mean)
predictor.P.df$mu.pi.lwr <- t(apply(NPD.mod.nc.mu.P, 2, PI, prob = 0.5))[,1]
predictor.P.df$mu.pi.upr <- t(apply(NPD.mod.nc.mu.P, 2, PI, prob = 0.5))[,2]
predictor.P.df$mu.pi.lwr2 <- t(apply(NPD.mod.nc.mu.P, 2, PI, prob = 0.8))[,1]
predictor.P.df$mu.pi.upr2 <- t(apply(NPD.mod.nc.mu.P, 2, PI, prob = 0.8))[,2]
predictor.P.df$mu.pi.lwr3 <- t(apply(NPD.mod.nc.mu.P, 2, PI, prob = 0.9))[,1]
predictor.P.df$mu.pi.upr3 <- t(apply(NPD.mod.nc.mu.P, 2, PI, prob = 0.9))[,2]

PNPDp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df %>%
                left_join(photo.summary.df %>%
                            select(ID, PNPRmean), by = "ID"),
              aes(PNPRmean, NPD, group = ID), fill = "#363538", width = 1,
              alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.P.df, aes(Pmean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.P.df, aes(Pmean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.P.df, aes(Pmean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.P.df, aes(Pmean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[d])*" ("*mu*"mol O"[2]*" g"^-1*" d"^-1*")"),
       x = "Pressure (hPa)") +
  scale_x_continuous(breaks = seq(1005, 1035, by = 10)) +
  scale_y_continuous(breaks = seq(-600, 200, by = 200)) +
  coord_cartesian(xlim = c(1005, 1035), ylim = c(-600, 200), clip = "off", expand = FALSE) +
  mytheme
PNPDp



predictor.S.df <- NPD.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Smean = seq(min(Smean), max(Smean), length.out = 200)) %>%
  mutate(Smean.us = Smean * sd(NPD.summary.df$Smean) + mean(NPD.summary.df$Smean)) # unstandardised

link_S <- function(Day, Smean){
  mu <- with(NPD.mod.nc.posterior, a + brsim * Day + bS * Smean)
  return(mu)
}

NPD.mod.nc.mu.S <- predictor.S.df %$% mapply(link_S, Day, Smean)


predictor.S.df$mu.mean <- apply(NPD.mod.nc.mu.S, 2, mean)
predictor.S.df$mu.pi.lwr <- t(apply(NPD.mod.nc.mu.S, 2, PI, prob = 0.5))[,1]
predictor.S.df$mu.pi.upr <- t(apply(NPD.mod.nc.mu.S, 2, PI, prob = 0.5))[,2]
predictor.S.df$mu.pi.lwr2 <- t(apply(NPD.mod.nc.mu.S, 2, PI, prob = 0.8))[,1]
predictor.S.df$mu.pi.upr2 <- t(apply(NPD.mod.nc.mu.S, 2, PI, prob = 0.8))[,2]
predictor.S.df$mu.pi.lwr3 <- t(apply(NPD.mod.nc.mu.S, 2, PI, prob = 0.9))[,1]
predictor.S.df$mu.pi.upr3 <- t(apply(NPD.mod.nc.mu.S, 2, PI, prob = 0.9))[,2]

SNPDp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df %>%
                left_join(photo.summary.df %>%
                            select(ID, Smean), by = "ID"),
              aes(Smean, NPD, group = ID), fill = "#363538", width = 0.2,
              alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.S.df, aes(Smean.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.S.df, aes(Smean.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.S.df, aes(Smean.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.S.df, aes(Smean.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[d])*" ("*mu*"mol O"[2]*" g"^-1*" d"^-1*")"),
       x = "Salinity (‰)") +
  scale_x_continuous(breaks = seq(33, 36, by = 1)) +
  scale_y_continuous(breaks = seq(-600, 200, by = 200)) +
  coord_cartesian(xlim = c(33, 36), ylim = c(-600, 200), clip = "off", expand = FALSE) +
  mytheme
SNPDp



predictor.M.df <- NPD.summary.l %$%
  data.frame(Day = mean(Day), # all other predictor variables have mean = 0
             Mass = seq(min(Mass), max(Mass), length.out = 200)) %>%
  mutate(Mass.us = Mass * sd(NPD.summary.df$Mass) + mean(NPD.summary.df$Mass)) # unstandardised

link_M <- function(Day, Mass){
  mu <- with(NPD.mod.nc.posterior, a + brsim * Day + bM * Mass)
  return(mu)
}

NPD.mod.nc.mu.M <- predictor.M.df %$% mapply(link_M, Day, Mass)


predictor.M.df$mu.mean <- apply(NPD.mod.nc.mu.M, 2, mean)
predictor.M.df$mu.pi.lwr <- t(apply(NPD.mod.nc.mu.M, 2, PI, prob = 0.5))[,1]
predictor.M.df$mu.pi.upr <- t(apply(NPD.mod.nc.mu.M, 2, PI, prob = 0.5))[,2]
predictor.M.df$mu.pi.lwr2 <- t(apply(NPD.mod.nc.mu.M, 2, PI, prob = 0.8))[,1]
predictor.M.df$mu.pi.upr2 <- t(apply(NPD.mod.nc.mu.M, 2, PI, prob = 0.8))[,2]
predictor.M.df$mu.pi.lwr3 <- t(apply(NPD.mod.nc.mu.M, 2, PI, prob = 0.9))[,1]
predictor.M.df$mu.pi.upr3 <- t(apply(NPD.mod.nc.mu.M, 2, PI, prob = 0.9))[,2]

MNPDp <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_violin(data = photo.density.df, aes(Mass, NPD, group = ID), fill = "#363538", width = 0.2,
              alpha = 0.2, colour = "#363538", position = "identity") +
  geom_line(data = predictor.M.df, aes(Mass.us, mu.mean), colour = "#363538") +
  geom_ribbon(data = predictor.M.df, aes(Mass.us, ymin = mu.pi.lwr, ymax = mu.pi.upr),
              alpha = 0.5, fill = "#363538") +
  geom_ribbon(data = predictor.M.df, aes(Mass.us, ymin = mu.pi.lwr2, ymax = mu.pi.upr2),
              alpha = 0.4, fill = "#363538") +
  geom_ribbon(data = predictor.M.df, aes(Mass.us, ymin = mu.pi.lwr3, ymax = mu.pi.upr3),
              alpha = 0.3, fill = "#363538") +
  labs(y = expression(italic(P[d])*" ("*mu*"mol O"[2]*" g"^-1*" d"^-1*")"),
       x = "Mass (g)") +
  scale_x_continuous(breaks = seq(0, 5, by = 1)) +
  scale_y_continuous(breaks = seq(-600, 200, by = 200)) +
  coord_cartesian(xlim = c(0, 5), ylim = c(-600, 200), clip = "off", expand = FALSE) +
  mytheme
MNPDp


require(patchwork)
Fig.2 <- (NPp + theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank())) /
         (GPp + theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      strip.text = element_blank())) /
         (NPDp + theme(strip.text = element_blank())) +
  plot_annotation(tag_levels = "a") & # the & symbol is important so that theme is applied to all tags
  theme(plot.tag = element_text(family = "Helvetica Neue",
                                size = 20, face = "bold"))
Fig.2 # 8 x 13 in

Fig.S2 <- (TNPp + theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank()) | 
           O2NPp + theme(axis.title.x = element_blank(),
                         axis.text.x = element_blank()) | 
           PNPp + theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank()) | 
           SNPp + theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank()) | 
           MNPp + theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank())) /
          (TGPp + theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank()) | 
           O2GPp + theme(axis.title.x = element_blank(),
                         axis.text.x = element_blank()) | 
           PGPp + theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank()) | 
           SGPp + theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank()) | 
           MGPp + theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank())) /
          (TNPDp | O2NPDp | PNPDp | SNPDp | MNPDp) +
  plot_annotation(tag_levels = list(c("a", rep("", 4), "b", rep("", 4), "c", rep("", 4)))) &
  theme(plot.tag = element_text(family = "Helvetica Neue",
                                size = 20, face = "bold"))
Fig.S2 # 8 x 13 in

# # try to fix with non-centred priors, i.e. priors that are not hierarchical
# 
# NP.mod.nc <- ulam(
#   alist(
#     # likelihood
#     NPmean ~ dnorm(mean = mu, sd = NPsd), # here the observed measurement error is introduced
#     vector[N]:NP ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP
#     
#     # linear model
#     mu <- a[Treatment] + b[Treatment] * Day + # treatment of interest (temp/PAR)
#       at[Tank] + bt[Tank] * Day + # tank effect
#       ac + bct * T_[i] + bco * O2[i] + bcp * Pmean + 
#       bcs * Salinity + bcm * Mass, # potential confounds
#     
#     # predictor measurement error
#     Tmean ~ dnorm(mean = T_, sd = Tsd),
#     vector[N]:T_ ~ dnorm(mean = 17.5, sd = 0.2),
#     O2mean ~ dnorm(mean = O2, sd = O2sd),
#     vector[N]:O2 ~ dnorm(mean = 232.5, sd = 30),
#     
#     # define multilevel coefficients using z-scores that HMC can sample better
#     save> vector[4]:a <<- a_ + za * sa,
#     save> vector[4]:b <<- b_ + zb * sb,
#     save> vector[16]:at <<- at_ + zat * sat,
#     save> vector[16]:bt <<- bt_ + zbt * sbt,
#     
#     # priors
#     vector[4]:za ~ dnorm(mean = 0, sd = 1), # z-score priors
#     vector[4]:zb ~ dnorm(mean = 0, sd = 1),
#     vector[16]:zat ~ dnorm(mean = 0, sd = 1),
#     vector[16]:zbt ~ dnorm(mean = 0, sd = 1),
#     
#     c(a_,at_,ac) ~ dnorm(mean = 15.39553, sd = 1), # intercept priors
#     c(b_,bt_) ~ dnorm(mean = -0.1374649, sd = 0.05), # slope priors
#     c(bct,bco,bcp,bcs,bcm) ~ dnorm(mean = 0, sd = 0.1), # potential confound slope priors
#     c(sigma,sa,sb,sat,sbt) ~ dexp(rate = 1) # standard deviation priors
#   ), 
#   data = meta.NP.cc.l, chains = 8, cores = 8, iter = 1e3,
#   control = list(adapt_delta = 0.99),
# )
# 
# traceplot(NP.mod.nc)
# trankplot(NP.mod.nc) # chains are even worse! who would have thought that possible?!
# precis(NP.mod.nc, depth = 2) # extraordinarily small n_eff and large Rhat4
# 
# 
# # the only remaining option is to model the effects of interest and 
# # the potential confounds separately
# 
# 
# #########################################
# NP.mod.main <- ulam(
#   alist(
#     # likelihood 
#     NPmean ~ normal(mean = NP, sd = NPsd), # here the observed measurement error is introduced
#     vector[N]:NP ~ normal(mean = mu, sd = sigma), # this describes the true, unobserved NP
#     
#     # linear model
#     mu <- atreat[Treatment] + btreat[Treatment] * Day + # treatment of interest (temp/PAR)
#       atank[Tank] + btank[Tank] * Day, # tank effect
#     
#     # adaptive priors
#     c(atreat, btreat)[Treatment] ~ multi_normal(c(a_treat, b_treat), rhotreat, sigmatreat),
#     c(atank, btank)[Tank] ~ multi_normal(Mu = c(a_tank, b_tank), rhotank, sigmatank),
#     
#     # hyper-priors
#     c(a_treat,a_tank) ~ normal(mean = 15.39553, sd = 1), # intercepts
#     c(b_treat,b_tank) ~ normal(mean = -0.1374649, sd = 0.05), # slopes across detrital age
#     c(sigma,sigmatreat,sigmatank) ~ exponential(rate = 1), # standard deviations
#     c(rhotreat, rhotank) ~ lkj_corr(eta = 2)
#   ), 
#   data = meta.NP.cc.l, chains = 8, cores = 8, iter = 1e3,
#   control = list(adapt_delta = 0.99),
# )
# 
# traceplot(NP.mod.main)
# trankplot(NP.mod.main) # chains are quite healthy
# precis(NP.mod.main, depth = 2) # ok n_eff and no extraordinarily large large Rhat4;
# # but note that the model has a hard time sampling the nested, centred priors
# # -> try to fix with non-centred priors
# 
# NP.mod.main.nc <- ulam(
#   alist(
#     # likelihood
#     NPmean ~ dnorm(mean = mu, sd = NPsd), # here the observed measurement error is introduced
#     vector[N]:NP ~ dnorm(mean = mu, sd = sigma), # this describes the true, unobserved NP
#     
#     # linear model
#     mu <- a[Treatment] + b[Treatment] * Day + # treatment of interest (temp/PAR)
#       at[Tank] + bt[Tank] * Day, # tank effect
#     
#     # define multilevel coefficients using z-scores that HMC can sample better
#     save> vector[4]:a <<- a_ + za * sa,
#     save> vector[4]:b <<- b_ + zb * sb,
#     save> vector[16]:at <<- at_ + zat * sat,
#     save> vector[16]:bt <<- bt_ + zbt * sbt,
#     
#     # priors
#     vector[4]:za ~ dnorm(mean = 0, sd = 1), # z-score priors
#     vector[4]:zb ~ dnorm(mean = 0, sd = 1),
#     vector[16]:zat ~ dnorm(mean = 0, sd = 1),
#     vector[16]:zbt ~ dnorm(mean = 0, sd = 1),
#     
#     c(a_,at_) ~ dnorm(mean = 15.39553, sd = 1), # intercept priors
#     c(b_,bt_) ~ dnorm(mean = -0.1374649, sd = 0.05), # slope priors
#     c(sigma,sa,sb,sat,sbt) ~ dexp(rate = 1) # standard deviation priors
#   ), 
#   data = meta.NP.cc.l, chains = 8, cores = 8, iter = 1e3,
#   control = list(adapt_delta = 0.99),
# )
# 
# traceplot(NP.mod.main.nc)
# trankplot(NP.mod.main.nc) # chains are ok
# precis(NP.mod.main.nc, depth = 2) # less good n_eff and some large Rhat4;
# # note that this model is better at some aspects but worse at others
# # hard to tell so compare n_eff between models visually
# 
# 
# # visual cross-validation
# precis_c <- precis(NP.mod.main, depth = 2)
# precis_nc <- precis(NP.mod.main.nc, depth = 2)
# 
# pars <- c(paste("NP[",1:52,"]",sep="") , paste("a[",1:4,"]",sep=""),
#           paste("b[",1:4,"]",sep="") , paste("at[",1:16,"]",sep=""),
#           paste("bt[",1:16,"]",sep=""), "at_" , "a_", "bt_", "b_",
#           "sbt", "sat", "sb", "sa", "sigma")
# 
# neff_table <- cbind(precis_c[pars,"n_eff"] , precis_nc[pars,"n_eff"] )
# dev.off() 
# plot( neff_table , xlim=range(neff_table) , ylim=range(neff_table) ,
# xlab="n_eff (centred)" , ylab="n_eff (non-centred)" , lwd=2 )
# abline( a=0 , b=1 , lty=2 )
# 
# # the centred prior model clearly wins!
#   




require(patchwork)





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

# data("chimpanzees")
# d <- chimpanzees
# d$treatment <- 1 + d$prosoc_left + 2*d$condition
# 
# dat_list <- list(
# pulled_left = d$pulled_left,
# actor = d$actor,
# block_id = d$block,
# treatment = as.integer(d$treatment) )
# 
# m13.4nc <- ulam(
# alist(
# pulled_left ~ dbinom( 1 , p ) ,
# logit(p) <- a_bar + z[actor]*sigma_a + # actor intercepts
#   x[block_id]*sigma_g +
#   b[treatment] ,
# b[treatment] ~ dnorm( 0 , 0.5 ),
# z[actor] ~ dnorm( 0 , 1 ),
# x[block_id] ~ dnorm( 0 , 1 ),
# a_bar ~ dnorm( 0 , 1.5 ),
# sigma_a ~ dexp(1),
# # block intercepts
# sigma_g ~ dexp(1),
# gq> vector[actor]:a <<- a_bar + z*sigma_a,
# gq> vector[block_id]:g <<- x*sigma_g
# ) , data=dat_list , chains=4 , cores=4 )
# 
# post <- extract.samples(m13.4nc)
# 
# p_link <- function( treatment , actor=2 , block_id=1 ) {
#   logodds <- with( post ,
#   a[,actor] + g[,block_id] + b[,treatment] )
#   return( inv_logit(logodds) )
# }
# 
# p_raw <- sapply( 1:4, p_link)
# p_mu <- apply( p_raw , 2 , mean )
# p_ci <- apply( p_raw , 2 , PI )
# 
# pd <- data.frame(t(p_ci), 
#                  mean = p_mu,
#                  treatment = 1:4)
# 
# ggplot() + geom_pointrange(data = pd, aes(x = treatment, y = mean, ymin = X5., ymax = X94.))
# 
# p_link_abar <- function( treatment ) {
#   logodds <- with( post , a_bar + b[,treatment] )
#   return( inv_logit(logodds) )
# }
# 
# p_raw <- sapply( 1:4 , p_link_abar)
# p_mu <- apply(p_raw , 2 , mean )
# p_ci <- apply(p_raw , 2 , PI )
# 
# pd <- data.frame(t(p_ci), 
#                  mean = p_mu,
#                  treatment = 1:4)
# 
# ggplot() + geom_pointrange(data = pd, aes(x = treatment, y = mean, ymin = X5., ymax = X94.))
# 
# 
# a_sim <- with( post , rnorm( length(post$a_bar) , a_bar , sigma_a ) )
# ggplot() + geom_density(aes(a_sim))
# 
# p_link_asim <- function( treatment ) {
#   logodds <- with( post , a_sim + b[,treatment] )
#   return( logodds )
# }
# p_raw_asim <- sapply( 1:4 , p_link_asim)
# 
# p_mu <- inv_logit(apply(p_raw_asim , 2 , mean ))
# p_ci <- inv_logit(apply(p_raw_asim , 2 , PI ))
# 
# pd <- data.frame(t(p_ci), 
#                  mean = p_mu,
#                  treatment = 1:4)
# 
# ggplot() + geom_pointrange(data = pd, aes(x = treatment, y = mean, ymin = X5., ymax = X94.))

