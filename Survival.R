### Packages ####
require(tidyverse) # this is my preferred data processing and visualisation framework 
require(magrittr) # this extends the pipe operator (%>%) of the tidyverse
require(rethinking) # this provides my Bayesian data analysis framework of choice

# load progress 
photo.summary.df <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/photo.density.df.csv") %>%
  select(-1) %>%
  group_by(ID) %>% 
  summarise(Tank = mean(Tank),
            Treatment = mean(Treatment),
            Day = mean(Day),
            NPmean = mean(NP),
            GPmean = mean(GP),
            NPDmean = mean(NPD),
            Samples = length(n))

meta <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Meta.csv") %>% 
  mutate(Treatment = ifelse(Temperature == 15 & PAR == 0, 1, # numbering follows increasing  
                            ifelse(Temperature == 15 & PAR == 8, 2, # PAR and temperature
                                   ifelse(Temperature == 20, 3, 4))))

meta$Treatment[is.na(meta$Treatment)] <- photo.summary.df$Treatment[photo.summary.df$Day == 0]
meta$Tank[is.na(meta$Tank)] <- photo.summary.df$Tank[photo.summary.df$Day == 0]

survival <- meta %>% 
  select(ID, Day, Tank, Treatment, Disintegration) %>%
  left_join(photo.summary.df %>% select(-Samples), by = c("ID", "Day", "Tank", "Treatment")) %>%
  mutate(AS = ifelse(Disintegration == 0 & NPmean > 0, 1, 0), # autotrophic survival
         PS = ifelse(Disintegration == 0 & GPmean > 0, 1, 0), # photosynthetic survival
         ADS = ifelse(Disintegration == 0 & NPDmean > 0, 1, 0)) # daily autotrophic survival
# disintegrated samples are assumed not to be alive so zeros are recorded for all survival variables

AS.df <- survival %>%
  select(Day, Tank, Treatment, AS)
AS.l <- as.list(AS.df)

# assign generic flat intercept and neutral slope priors first
a.prior <- rnorm(n = 1e4, mean = 0, sd = 1.5)
ggplot() +
  geom_density(aes(inv_logit(a.prior))) +
  theme_minimal()

b.prior <- rnorm(n = 1e4, mean = 0, sd = 0.05)
ggplot() +
  geom_density(aes(inv_logit(b.prior))) +
  theme_minimal()


predictor <- AS.df %$% seq(min(Day), max(Day), length.out = 100)
plot(NULL, xlim = c(0, 120), ylim = c(0, 1), # base plot is better with loops
     xlab = "d", ylab = "probability") # empty plot
for(i in 1:1e3) curve(inv_logit(a.prior[i] + b.prior[i] * x), # loop
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = col.alpha("black", 0.1))
# nearly all relationships are possible
# but we already know that we are modelling survival, which necessarily starts at 100%
# so the prior intercept should be around 1, which on the logit scale is approached 
# with near certainty at 5
# survival curves also pre-exempt positive slopes so the prior slope should be negative 
# or neutral

a.prior <- rnorm(n = 1e4, mean = 5, sd = 2.5) 
ggplot() +
  geom_density(aes(inv_logit(a.prior))) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() # intercepts below p = 0.75 very unlikely

b.prior <- rnorm(n = 1e4, mean = -0.1, sd = 0.2)
ggplot() +
  geom_density(aes(inv_logit(b.prior))) +
  geom_vline(xintercept = 0.5) +
  theme_minimal() # positive slopes less unlikely


predictor <- AS.df %$% seq(min(Day), max(Day), length.out = 100)
plot(NULL, xlim = c(0, 120), ylim = c(0, 1), # base plot is better with loops
     xlab = "d", ylab = "probability") # empty plot
for(i in 1:1e3) curve(inv_logit(a.prior[i] + b.prior[i] * x), # loop
                      from = min(predictor),
                      to = max(predictor),
                      add = T, col = col.alpha("black", 0.1))
# looks better while still allowing some improbable relationships


AS.mod <- ulam(
  alist(
    AS ~ dbinom(size = 1, prob = p), # likelihood same as Bernoulli
    logit(p) <- a + b[Treatment] * Day + br[Tank] * Day, # generalised linear model
    br[Tank] ~ dnorm(mean = 0, sd = sbr), # higher level prior for partial pooling
    a ~ dnorm(mean = 5, sd = 2.5), # intercept prior
    b[Treatment] ~ dnorm(mean = -0.1, sd = 0.2), # slope prior
    sbr ~ dexp(rate = 1) # standard deviation prior for mean Tank
  ), data = AS.l, chains = 8, cores = parallel::detectCores(), iter = 1e4
)
# note warnings: chains finished unsuccessfully and
# there are some divergent transitions
traceplot(AS.mod)
trankplot(AS.mod) # some unhealthy chains
precis(AS.mod, depth = 2) # some low n_eff
dev.off()
plot(precis(AS.mod, depth = 2))
# try to fix with non-centred prior

AS.mod.nc <- ulam(
  alist(
    AS ~ dbinom(size = 1, prob = p), # likelihood same as Bernoulli
    logit(p) <- a + b[Treatment] * Day + br[Tank] * Day, # generalised linear model
    save> vector[16]:br <<- 0 + zbr * sbr, # non-centred expression using z-score
    vector[16]:zbr ~ dnorm(mean = 0, sd = 1), # z-score prior
    a ~ dnorm(mean = 5, sd = 2.5), # intercept prior
    b[Treatment] ~ dnorm(mean = -0.1, sd = 0.2), # slope prior
    sbr ~ dexp(rate = 1) # standard deviation prior for mean Tank
  ), data = AS.l, chains = 8, cores = parallel::detectCores(), iter = 1e4
)

# no warnings (if warnings, run again)
traceplot(AS.mod.nc)
trankplot(AS.mod.nc) # healthier chains
precis(AS.mod.nc, depth = 2) # better n_eff
dev.off()
plot(precis(AS.mod.nc, depth = 2))

# visual cross-validation
precis_c <- precis(AS.mod, depth = 2)
precis_nc <- precis(AS.mod.nc, depth = 2)

pars <- c(paste("b[",1:4,"]",sep=""), paste("br[",1:16,"]",sep=""), "a", "sbr")

neff_table <- cbind(precis_c[pars,"n_eff"] , precis_nc[pars,"n_eff"] )
dev.off()
plot( neff_table , xlim=range(neff_table) , ylim=range(neff_table) ,
      xlab="n_eff (centred)" , ylab="n_eff (non-centred)" , lwd=2 )
abline( a=0 , b=1 , lty=2 )

# non-centred model clearly wins!


AS.mod.nc.prior <- extract.prior(AS.mod.nc, iter = 2e4)
AS.mod.nc.posterior <- extract.samples(AS.mod.nc)

AS.mod.nc.prior.df <- data.frame(AS.mod.nc.prior) %>%
  rownames_to_column(var = "n") %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "level")) %>%
  mutate(level = as.integer(level)) %>%
  arrange(coefficient, level, n) %>%
  mutate(level = factor(level), coefficient = factor(coefficient))

AS.mod.nc.posterior.df <- data.frame(AS.mod.nc.posterior) %>%
  rownames_to_column(var = "n") %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "level")) %>%
  mutate(level = as.integer(level)) %>%
  arrange(coefficient, level, n) %>%
  mutate(level = factor(level), coefficient = factor(coefficient))

AS.mod.nc.prior.posterior.df <- data.frame(rbind(AS.mod.nc.prior.df,
                                                 AS.mod.nc.posterior.df),
                                        distribution = rep(c("prior", "posterior"), each = 380000)) %>%
  mutate(distribution = factor(distribution))


ggplot(data = AS.mod.nc.prior.posterior.df %>%
         filter(coefficient %in% c("a", "b"))) +
  geom_density(aes(x = sample, colour = distribution)) +
  facet_wrap(~ coefficient + level, scales = "free") +
  theme_minimal() +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")


# Posterior prediction
# Parameters
# not only can log-odds be interpreted directly, making conversion to the
# probability scale unnecessary at this stage, log-odds parameters are also
# the coefficients in the logistic equation and meaningful as such

# only two parameters are really of interest: the slope (b) which is the logistic
# decay rate k and the intercept (a) from which the inflection point of the sigmoid
# (µ, x at 50% y) can be derived as a / b

brsim <- AS.mod.nc.posterior %$% rnorm(n = length(sbr), mean = 0, sd = sbr)

AS.para <- AS.mod.nc.posterior.df %>%
  filter(coefficient == "b") %>%
  mutate(br = rep(brsim, 4),
         b = sample + br) %>% # add variation across tanks to posterior slope
  select(-c(coefficient, sample)) %>%
  left_join(AS.mod.nc.posterior.df %>%
              filter(coefficient == "a") %>%
              select(-c(coefficient, level)),
            by = "n") %>%
  rename(a = "sample") %>%
  mutate(µ = -(a / b))


# Contrasts
AS.comp <- AS.para %>%
  select(-c(a, br)) %>%
  pivot_wider(names_from = level, values_from = c(b, µ)) %>%
  mutate(b1_2 = b_1 - b_2, b1_3 = b_1 - b_3, b1_4 = b_1 - b_4,
         b2_3 = b_2 - b_3, b2_4 = b_2 - b_4, b3_4 = b_3 - b_4,
         µ1_2 = µ_1 - µ_2, µ1_3 = µ_1 - µ_3, µ1_4 = µ_1 - µ_4,
         µ2_3 = µ_2 - µ_3, µ2_4 = µ_2 - µ_4, µ3_4 = µ_3 - µ_4)

# b
AS.comp %>%
  filter(b1_2 < 0) %>% pull(b1_2) %>% length() /
  AS.comp %>% pull(b1_2) %>% length()
# 97% probability that b1 < b2

AS.comp %>%
  filter(b1_3 > 0) %>% pull(b1_3) %>% length() /
  AS.comp %>% pull(b1_3) %>% length()
# 61% probability that b1 > b3

AS.comp %>%
  filter(b1_4 > 0) %>% pull(b1_4) %>% length() /
  AS.comp %>% pull(b1_4) %>% length()
# 94% probability that b1 > b4

AS.comp %>%
  filter(b2_3 > 0) %>% pull(b2_3) %>% length() /
  AS.comp %>% pull(b2_3) %>% length()
# 98% probability that b2 > b3

AS.comp %>%
  filter(b2_4 > 0) %>% pull(b2_4) %>% length() /
  AS.comp %>% pull(b2_4) %>% length()
# 100% probability that b2 > b4

AS.comp %>%
  filter(b3_4 > 0) %>% pull(b3_4) %>% length() /
  AS.comp %>% pull(b3_4) %>% length()
# 90% probability that b3 > b4

# µ
AS.comp %>%
  filter(µ1_2 < 0) %>% pull(µ1_2) %>% length() /
  AS.comp %>% pull(µ1_2) %>% length()
# 79% probability that µ1 < µ2

AS.comp %>%
  filter(µ1_3 > 0) %>% pull(µ1_3) %>% length() /
  AS.comp %>% pull(µ1_3) %>% length()
# 60% probability that µ1 > µ3

AS.comp %>%
  filter(µ1_4 > 0) %>% pull(µ1_4) %>% length() /
  AS.comp %>% pull(µ1_4) %>% length()
# 91% probability that µ1 > µ4

AS.comp %>%
  filter(µ2_3 > 0) %>% pull(µ2_3) %>% length() /
  AS.comp %>% pull(µ2_3) %>% length()
# 79% probability that µ2 > µ3

AS.comp %>%
  filter(µ2_4 > 0) %>% pull(µ2_4) %>% length() /
  AS.comp %>% pull(µ2_4) %>% length()
# 79% probability that µ2 > µ4

AS.comp %>%
  filter(µ3_4 > 0) %>% pull(µ3_4) %>% length() /
  AS.comp %>% pull(µ3_4) %>% length()
# 88% probability that µ3 > µ4

# Probability mass of b < 0
AS.comp %>%
  filter(b_1 < 0) %>% pull(b_1) %>% length() /
  AS.comp %>% pull(b_1) %>% length() # 97%

AS.comp %>%
  filter(b_2 < 0) %>% pull(b_2) %>% length() /
  AS.comp %>% pull(b_2) %>% length() # 78%

AS.comp %>%
  filter(b_3 < 0) %>% pull(b_3) %>% length() /
  AS.comp %>% pull(b_3) %>% length() # 97%

AS.comp %>%
  filter(b_4 < 0) %>% pull(b_4) %>% length() /
  AS.comp %>% pull(b_4) %>% length() # 99%

# Equations
AS.coef <- AS.para %>%
  group_by(level) %>%
  summarise(a.mean = mean(a),
            a.sd = sd(a),
            b.mean = mean(b),
            b.sd = sd(b),
            b.lwr = PI(b, prob = 0.9)[1],
            b.upr = PI(b, prob = 0.9)[2], 
            µ.median = median(µ),
            µ.mean = mean(µ)) %>% # quotient distributions are near impossible to 
  mutate(µ.mean.derived = -(a.mean / b.mean)) # describe with central tendencies
# a = 3.12 ± 0.63
# b1 = -0.12 ± 0.07
# b2 = -0.04 ± 0.06
# b3 = -0.13 ± 0.07
# b4 = -0.2 ± 0.08
# µ1 = 26
# µ2 = 85
# µ3 = 24
# µ4 = 15

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

ggplot() +
  geom_vline(data = AS.para, aes(xintercept = mode(µ))) +
  geom_vline(data = AS.para, aes(xintercept = mean(µ)),
             colour = "orange") +
  geom_vline(data = AS.coef, aes(xintercept = µ.mean.derived),
             colour = "red") +
  geom_density(data = AS.para, aes(µ)) +
  facet_grid(~level) +
  scale_x_continuous(limits = c(0, 120), 
                     oob = scales::oob_keep) 
# impossible to describe µ for treatment 2


AS.annotation <- AS.df %>%
  group_by(Treatment) %>%
  summarise(n = length(Treatment))
str(AS.annotation)

AS.annotation %<>%
  mutate(n = as.character(c(expression(italic("n ")*"= 26"),
                            expression(italic("n ")*"= 25"),
                            expression(italic("n ")*"= 25"),
                            expression(italic("n ")*"= 25"))),
         equation = as.character(c(expression(italic("y ")*"= "*frac(1, 1+e^{0.12*"("*italic("x ")*"– 26)"})),
                                   expression(italic("y ")*"= "*frac(1, 1+e^{0.04*"("*italic("x ")*"– 85)"})),
                                   expression(italic("y ")*"= "*frac(1, 1+e^{0.13*"("*italic("x ")*"– 24)"})),
                                   expression(italic("y ")*"= "*frac(1, 1+e^{0.2*"("*italic("x ")*"– 15)"})))))



# Lines and intervals
predictor.df <- data.frame(Day = rep(seq(0, 119, length.out = 200), 4),
                           Treatment = rep(1:4, each = 200))

# calculate posterior prediction for new simulated Tanks, i.e. incorporating Tank s.d.

link_sim <- function(Treatment, Day){
  logit_p <- with(AS.mod.nc.posterior, a + b[,Treatment] * Day + brsim * Day)
  return(logit_p)
}

AS.mod.nc.logit_p <- predictor.df %$% mapply(link_sim, Treatment, Day)

# summarise posterior probability distributions
predictor.df$AS.p.mean <- inv_logit(apply(AS.mod.nc.logit_p, 2, mean))
predictor.df$AS.p.pi.lwr <- inv_logit(t(apply(AS.mod.nc.logit_p, 2, PI, prob = 0.5))[,1])
predictor.df$AS.p.pi.upr <- inv_logit(t(apply(AS.mod.nc.logit_p, 2, PI, prob = 0.5))[,2])
predictor.df$AS.p.pi.lwr2 <- inv_logit(t(apply(AS.mod.nc.logit_p, 2, PI, prob = 0.8))[,1])
predictor.df$AS.p.pi.upr2 <- inv_logit(t(apply(AS.mod.nc.logit_p, 2, PI, prob = 0.8))[,2])
predictor.df$AS.p.pi.lwr3 <- inv_logit(t(apply(AS.mod.nc.logit_p, 2, PI, prob = 0.9))[,1])
predictor.df$AS.p.pi.upr3 <- inv_logit(t(apply(AS.mod.nc.logit_p, 2, PI, prob = 0.9))[,2])


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
                 panel.spacing = unit(.8, "cm"),
                 text = element_text(family = "Helvetica Neue"))

Treatment_names <- c(
  '1' = "Dark 15°C",
  '2' = "Light 15°C",
  '3' = "Light 20°C",
  '4' = "Light 25°C"
)

require(ggdist)
ASp <- ggplot() +
  # geom_point(data = AS.df %>%
  #              group_by(Treatment, Day) %>%
  #              summarise(Day = mean(Day),
  #                        AS = mean(AS)),
  #            aes(Day, AS, colour = factor(Treatment))) +
  geom_dots(data = AS.df, aes(Day, AS, side = factor(AS), fill = factor(Treatment)),
            scale = 0.25, colour = NA) +
  geom_line(data = predictor.df, aes(Day, AS.p.mean, colour = factor(Treatment))) +
  geom_ribbon(data = predictor.df, aes(Day, ymin = AS.p.pi.lwr, ymax = AS.p.pi.upr,
                                       fill = factor(Treatment)), alpha = 0.5) +
  geom_ribbon(data = predictor.df, aes(Day, ymin = AS.p.pi.lwr2, ymax = AS.p.pi.upr2,
                                       fill = factor(Treatment)),  alpha = 0.4) +
  geom_ribbon(data = predictor.df, aes(Day, ymin = AS.p.pi.lwr3, ymax = AS.p.pi.upr3,
                                       fill = factor(Treatment)),  alpha = 0.3) +
  geom_point(data = AS.coef %>% rename(Treatment = "level"),
             aes(µ.mean.derived, 0.5, fill = factor(Treatment)),
             shape = 21, size = 3.5, colour = "#ffffff") +
  geom_text(data = AS.annotation, aes(120, 1.075, label = n),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = 1) +
  geom_text(data = AS.annotation, aes(c(120, 0, 120, 120), c(0.88, 0.25, 0.88, 0.88), label = equation),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = c(1, 0, 1, 1)) +
  scale_side_mirrored(guide = "none") +
  scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  facet_grid(~Treatment, labeller = labeller(Treatment = Treatment_names)) +
  labs(y = expression("Autotrophic ("*italic(p)*")"),
       x = "Detrital age (d)") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 120),
                  clip = "off", expand = FALSE) +
  mytheme
ASp


PS.df <- survival %>%
  select(Day, Tank, Treatment, PS)
PS.l <- as.list(PS.df)

PS.mod <- ulam(
  alist(
    PS ~ dbinom(size = 1, prob = p), # likelihood same as Bernoulli
    logit(p) <- a + b[Treatment] * Day + br[Tank] * Day, # generalised linear model
    br[Tank] ~ dnorm(mean = 0, sd = sbr), # higher level prior for partial pooling
    a ~ dnorm(mean = 5, sd = 2.5), # intercept prior
    b[Treatment] ~ dnorm(mean = -0.1, sd = 0.2), # slope prior
    sbr ~ dexp(rate = 1) # standard deviation prior for mean Tank
  ), data = PS.l, chains = 8, cores = parallel::detectCores(), iter = 1e4
)
# note warnings: chains finished unsuccessfully and
# there are some divergent transitions
traceplot(PS.mod)
trankplot(PS.mod) # some unhealthy chains
precis(PS.mod, depth = 2) # some low n_eff
dev.off()
plot(precis(PS.mod, depth = 2))
# try to fix with non-centred prior

PS.mod.nc <- ulam(
  alist(
    PS ~ dbinom(size = 1, prob = p), # likelihood same as Bernoulli
    logit(p) <- a + b[Treatment] * Day + br[Tank] * Day, # generalised linear model
    save> vector[16]:br <<- 0 + zbr * sbr, # non-centred expression using z-score
    vector[16]:zbr ~ dnorm(mean = 0, sd = 1), # z-score prior
    a ~ dnorm(mean = 5, sd = 2.5), # intercept prior
    b[Treatment] ~ dnorm(mean = -0.1, sd = 0.2), # slope prior
    sbr ~ dexp(rate = 1) # standard deviation prior for mean Tank
  ), data = PS.l, chains = 8, cores = parallel::detectCores(), iter = 1e4
)

# no warnings (if warnings, run again)
traceplot(PS.mod.nc)
trankplot(PS.mod.nc) # healthier chains
precis(PS.mod.nc, depth = 2) # better n_eff
dev.off()
plot(precis(PS.mod.nc, depth = 2))

# visual cross-validation
precis_c <- precis(PS.mod, depth = 2)
precis_nc <- precis(PS.mod.nc, depth = 2)

pars <- c(paste("b[",1:4,"]",sep=""), paste("br[",1:16,"]",sep=""), "a", "sbr")

neff_table <- cbind(precis_c[pars,"n_eff"] , precis_nc[pars,"n_eff"] )
dev.off()
plot( neff_table , xlim=range(neff_table) , ylim=range(neff_table) ,
      xlab="n_eff (centred)" , ylab="n_eff (non-centred)" , lwd=2 )
abline( a=0 , b=1 , lty=2 )

# non-centred model clearly wins!


PS.mod.nc.prior <- extract.prior(PS.mod.nc, iter = 2e4)
PS.mod.nc.posterior <- extract.samples(PS.mod.nc)

PS.mod.nc.prior.df <- data.frame(PS.mod.nc.prior) %>%
  rownames_to_column(var = "n") %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "level")) %>%
  mutate(level = as.integer(level)) %>%
  arrange(coefficient, level, n) %>%
  mutate(level = factor(level), coefficient = factor(coefficient))

PS.mod.nc.posterior.df <- data.frame(PS.mod.nc.posterior) %>%
  rownames_to_column(var = "n") %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "level")) %>%
  mutate(level = as.integer(level)) %>%
  arrange(coefficient, level, n) %>%
  mutate(level = factor(level), coefficient = factor(coefficient))

PS.mod.nc.prior.posterior.df <- data.frame(rbind(PS.mod.nc.prior.df,
                                                 PS.mod.nc.posterior.df),
                                        distribution = rep(c("prior", "posterior"), each = 380000)) %>%
  mutate(distribution = factor(distribution))


ggplot(data = PS.mod.nc.prior.posterior.df %>%
         filter(coefficient %in% c("a", "b"))) +
  geom_density(aes(x = sample, colour = distribution)) +
  facet_wrap(~ coefficient + level, scales = "free") +
  theme_minimal() +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")


# Posterior prediction
# Parameters
# not only can log-odds be interpreted directly, making conversion to the
# probability scale unnecessary at this stage, log-odds parameters are also
# the coefficients in the logistic equation and meaningful as such

# only two parameters are really of interest: the slope (b) which is the logistic
# decay rate k and the intercept (a) from which the inflection point of the sigmoid
# (µ, x at 50% y) can be derived as a / k

brsim <- PS.mod.nc.posterior %$% rnorm(n = length(sbr), mean = 0, sd = sbr)
PS.para <- PS.mod.nc.posterior.df %>%
  filter(coefficient == "b") %>%
  mutate(br = rep(brsim, 4),
         b = sample + br) %>% # add variation across tanks to posterior slope
  select(-c(coefficient, sample)) %>%
  left_join(PS.mod.nc.posterior.df %>%
              filter(coefficient == "a") %>%
              select(-c(coefficient, level)),
            by = "n") %>%
  rename(a = "sample") %>%
  mutate(µ = -(a / b))


# Contrasts
PS.comp <- PS.para %>%
  select(-c(a, br)) %>%
  pivot_wider(names_from = level, values_from = c(b, µ)) %>%
  mutate(b1_2 = b_1 - b_2, b1_3 = b_1 - b_3, b1_4 = b_1 - b_4,
         b2_3 = b_2 - b_3, b2_4 = b_2 - b_4, b3_4 = b_3 - b_4,
         µ1_2 = µ_1 - µ_2, µ1_3 = µ_1 - µ_3, µ1_4 = µ_1 - µ_4,
         µ2_3 = µ_2 - µ_3, µ2_4 = µ_2 - µ_4, µ3_4 = µ_3 - µ_4)

# b
PS.comp %>%
  filter(b1_2 < 0) %>% pull(b1_2) %>% length() /
  PS.comp %>% pull(b1_2) %>% length()
# 99% probability that b1 < b2

PS.comp %>%
  filter(b1_3 > 0) %>% pull(b1_3) %>% length() /
  PS.comp %>% pull(b1_3) %>% length()
# 70% probability that b1 > b3

PS.comp %>%
  filter(b1_4 > 0) %>% pull(b1_4) %>% length() /
  PS.comp %>% pull(b1_4) %>% length()
# 99% probability that b1 > b4

PS.comp %>%
  filter(b2_3 > 0) %>% pull(b2_3) %>% length() /
  PS.comp %>% pull(b2_3) %>% length()
# 100% probability that b2 > b3

PS.comp %>%
  filter(b2_4 > 0) %>% pull(b2_4) %>% length() /
  PS.comp %>% pull(b2_4) %>% length()
# 100% probability that b2 > b4

PS.comp %>%
  filter(b3_4 > 0) %>% pull(b3_4) %>% length() /
  PS.comp %>% pull(b3_4) %>% length()
# 98% probability that b3 > b4

# µ
PS.comp %>%
  filter(µ1_2 < 0) %>% pull(µ1_2) %>% length() /
  PS.comp %>% pull(µ1_2) %>% length()
# 90% probability that µ1 < µ2

PS.comp %>%
  filter(µ1_3 > 0) %>% pull(µ1_3) %>% length() /
  PS.comp %>% pull(µ1_3) %>% length()
# 70% probability that µ1 > µ3

PS.comp %>%
  filter(µ1_4 > 0) %>% pull(µ1_4) %>% length() /
  PS.comp %>% pull(µ1_4) %>% length()
# 99% probability that µ1 > µ4

PS.comp %>%
  filter(µ2_3 > 0) %>% pull(µ2_3) %>% length() /
  PS.comp %>% pull(µ2_3) %>% length()
# 91% probability that µ2 > µ3

PS.comp %>%
  filter(µ2_4 > 0) %>% pull(µ2_4) %>% length() /
  PS.comp %>% pull(µ2_4) %>% length()
# 91% probability that µ2 > µ4

PS.comp %>%
  filter(µ3_4 > 0) %>% pull(µ3_4) %>% length() /
  PS.comp %>% pull(µ3_4) %>% length()
# 98% probability that µ3 > µ4

# Probability mass of b < 0
PS.comp %>%
  filter(b_1 < 0) %>% pull(b_1) %>% length() /
  PS.comp %>% pull(b_1) %>% length() # 99%

PS.comp %>%
  filter(b_2 < 0) %>% pull(b_2) %>% length() /
  PS.comp %>% pull(b_2) %>% length() # 91%

PS.comp %>%
  filter(b_3 < 0) %>% pull(b_3) %>% length() /
  PS.comp %>% pull(b_3) %>% length() # 99%

PS.comp %>%
  filter(b_4 < 0) %>% pull(b_4) %>% length() /
  PS.comp %>% pull(b_4) %>% length() # 100%

# Equations
PS.coef <- PS.para %>%
  group_by(level) %>%
  summarise(a.mean = mean(a),
            a.sd = sd(a),
            b.mean = mean(b),
            b.sd = sd(b),
            b.lwr = PI(b, prob = 0.9)[1],
            b.upr = PI(b, prob = 0.9)[2], 
            µ.median = median(µ),
            µ.mean = mean(µ)) %>% # quotient distributions are near impossible to 
  mutate(µ.mean.derived = -(a.mean / b.mean)) # describe with central tendencies
# a = 4.82 ± 0.92
# b1 = -0.13 ± 0.06
# b2 = -0.05 ± 0.05
# b3 = -0.16 ± 0.06
# b4 = -0.27 ± 0.08
# µ1 = 36
# µ2 = 99
# µ3 = 31
# µ4 = 18



PS.annotation <- PS.df %>%
  group_by(Treatment) %>%
  summarise(n = length(Treatment))
str(PS.annotation)

PS.annotation %<>%
  mutate(n = as.character(c(expression(italic("n ")*"= 26"),
                            expression(italic("n ")*"= 25"),
                            expression(italic("n ")*"= 25"),
                            expression(italic("n ")*"= 25"))),
         equation = as.character(c(expression(italic("y ")*"= "*frac(1, 1+e^{0.13*"("*italic("x ")*"– 36)"})),
                                   expression(italic("y ")*"= "*frac(1, 1+e^{0.05*"("*italic("x ")*"– 99)"})),
                                   expression(italic("y ")*"= "*frac(1, 1+e^{0.16*"("*italic("x ")*"– 31)"})),
                                   expression(italic("y ")*"= "*frac(1, 1+e^{0.27*"("*italic("x ")*"– 18)"})))))



# Lines and intervals

# calculate posterior prediction for new simulated Tanks, i.e. incorporating Tank s.d.

link_sim <- function(Treatment, Day){
  logit_p <- with(PS.mod.nc.posterior, a + b[,Treatment] * Day + brsim * Day)
  return(logit_p)
}

PS.mod.nc.logit_p <- predictor.df %$% mapply(link_sim, Treatment, Day)

# summarise posterior probability distributions
predictor.df$PS.p.mean <- inv_logit(apply(PS.mod.nc.logit_p, 2, mean))
predictor.df$PS.p.pi.lwr <- inv_logit(t(apply(PS.mod.nc.logit_p, 2, PI, prob = 0.5))[,1])
predictor.df$PS.p.pi.upr <- inv_logit(t(apply(PS.mod.nc.logit_p, 2, PI, prob = 0.5))[,2])
predictor.df$PS.p.pi.lwr2 <- inv_logit(t(apply(PS.mod.nc.logit_p, 2, PI, prob = 0.8))[,1])
predictor.df$PS.p.pi.upr2 <- inv_logit(t(apply(PS.mod.nc.logit_p, 2, PI, prob = 0.8))[,2])
predictor.df$PS.p.pi.lwr3 <- inv_logit(t(apply(PS.mod.nc.logit_p, 2, PI, prob = 0.9))[,1])
predictor.df$PS.p.pi.upr3 <- inv_logit(t(apply(PS.mod.nc.logit_p, 2, PI, prob = 0.9))[,2])


PSp <- ggplot() +
  geom_dots(data = PS.df, aes(Day, PS, side = factor(PS), fill = factor(Treatment)),
            scale = 0.25, colour = NA) +
  geom_line(data = predictor.df, aes(Day, PS.p.mean, colour = factor(Treatment))) +
  geom_ribbon(data = predictor.df, aes(Day, ymin = PS.p.pi.lwr, ymax = PS.p.pi.upr,
                                       fill = factor(Treatment)), alpha = 0.5) +
  geom_ribbon(data = predictor.df, aes(Day, ymin = PS.p.pi.lwr2, ymax = PS.p.pi.upr2,
                                       fill = factor(Treatment)),  alpha = 0.4) +
  geom_ribbon(data = predictor.df, aes(Day, ymin = PS.p.pi.lwr3, ymax = PS.p.pi.upr3,
                                       fill = factor(Treatment)),  alpha = 0.3) +
  geom_point(data = PS.coef %>% rename(Treatment = "level"),
             aes(µ.mean.derived, 0.5, fill = factor(Treatment)),
             shape = 21, size = 3.5, colour = "#ffffff") +
  geom_text(data = PS.annotation, aes(120, 1.075, label = n),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = 1) +
  geom_text(data = PS.annotation, aes(c(120, 0, 120, 120), c(0.88, 0.15, 0.88, 0.88), label = equation),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = c(1, 0, 1, 1)) +
  scale_side_mirrored(guide = "none") +
  scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  facet_grid(~Treatment, labeller = labeller(Treatment = Treatment_names)) +
  labs(y = expression("Photosynthetic ("*italic(p)*")"),
       x = "Detrital age (d)") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 120),
                  clip = "off", expand = FALSE) +
  mytheme
PSp



ADS.df <- survival %>%
  select(Day, Tank, Treatment, ADS)
ADS.l <- as.list(ADS.df)

ADS.mod <- ulam(
  alist(
    ADS ~ dbinom(size = 1, prob = p), # likelihood same as Bernoulli
    logit(p) <- a + b[Treatment] * Day + br[Tank] * Day, # generalised linear model
    br[Tank] ~ dnorm(mean = 0, sd = sbr), # higher level prior for partial pooling
    a ~ dnorm(mean = 5, sd = 2.5), # intercept prior
    b[Treatment] ~ dnorm(mean = -0.1, sd = 0.2), # slope prior
    sbr ~ dexp(rate = 1) # standard deviation prior for mean Tank
  ), data = ADS.l, chains = 8, cores = parallel::detectCores(), iter = 1e4
)
# note warnings: chains finished unsuccessfully and
# there are some divergent transitions
traceplot(ADS.mod)
trankplot(ADS.mod) # some unhealthy chains
precis(ADS.mod, depth = 2) # some low n_eff
dev.off()
plot(precis(ADS.mod, depth = 2))
# try to fix with non-centred prior

ADS.mod.nc <- ulam(
  alist(
    ADS ~ dbinom(size = 1, prob = p), # likelihood same as Bernoulli
    logit(p) <- a + b[Treatment] * Day + br[Tank] * Day, # generalised linear model
    save> vector[16]:br <<- 0 + zbr * sbr, # non-centred expression using z-score
    vector[16]:zbr ~ dnorm(mean = 0, sd = 1), # z-score prior
    a ~ dnorm(mean = 5, sd = 2.5), # intercept prior
    b[Treatment] ~ dnorm(mean = -0.1, sd = 0.2), # slope prior
    sbr ~ dexp(rate = 1) # standard deviation prior for mean Tank
  ), data =ADS.l, chains = 8, cores = parallel::detectCores(), iter = 1e4
)

# no warnings (if warnings, run again)
traceplot(ADS.mod.nc)
trankplot(ADS.mod.nc) # healthier chains
precis(ADS.mod.nc, depth = 2) # better n_eff
dev.off()
plot(precis(ADS.mod.nc, depth = 2))

# visual cross-validation
precis_c <- precis(ADS.mod, depth = 2)
precis_nc <- precis(ADS.mod.nc, depth = 2)

pars <- c(paste("b[",1:4,"]",sep=""), paste("br[",1:16,"]",sep=""), "a", "sbr")

neff_table <- cbind(precis_c[pars,"n_eff"] , precis_nc[pars,"n_eff"] )
dev.off()
plot( neff_table , xlim=range(neff_table) , ylim=range(neff_table) ,
      xlab="n_eff (centred)" , ylab="n_eff (non-centred)" , lwd=2 )
abline( a=0 , b=1 , lty=2 )

# non-centred model clearly wins!


ADS.mod.nc.prior <- extract.prior(ADS.mod.nc, iter = 2e4)
ADS.mod.nc.posterior <- extract.samples(ADS.mod.nc)

ADS.mod.nc.prior.df <- data.frame(ADS.mod.nc.prior) %>%
  rownames_to_column(var = "n") %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "level")) %>%
  mutate(level = as.integer(level)) %>%
  arrange(coefficient, level, n) %>%
  mutate(level = factor(level), coefficient = factor(coefficient))

ADS.mod.nc.posterior.df <- data.frame(ADS.mod.nc.posterior) %>%
  rownames_to_column(var = "n") %>%
  pivot_longer(cols = !n, names_to = "coefficient", values_to = "sample") %>%
  mutate(n = as.integer(n)) %>%
  separate_wider_delim(col = coefficient, delim = ".", too_few = "align_start",
                       names = c("coefficient", "level")) %>%
  mutate(level = as.integer(level)) %>%
  arrange(coefficient, level, n) %>%
  mutate(level = factor(level), coefficient = factor(coefficient))

ADS.mod.nc.prior.posterior.df <- data.frame(rbind(ADS.mod.nc.prior.df,
                                                  ADS.mod.nc.posterior.df),
                                        distribution = rep(c("prior", "posterior"), each = 380000)) %>%
  mutate(distribution = factor(distribution))


ggplot(data = ADS.mod.nc.prior.posterior.df %>%
         filter(coefficient %in% c("a", "b"))) +
  geom_density(aes(x = sample, colour = distribution)) +
  facet_wrap(~ coefficient + level, scales = "free") +
  theme_minimal() +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")


# Posterior prediction
# Parameters
# not only can log-odds be interpreted directly, making conversion to the
# probability scale unnecessary at this stage, log-odds parameters are also
# the coefficients in the logistic equation and meaningful as such

# only two parameters are really of interest: the slope (b) which is the logistic
# decay rate k and the intercept (a) from which the inflection point of the sigmoid
# (µ, x at 50% y) can be derived as a / k

brsim <- ADS.mod.nc.posterior %$% rnorm(n = length(sbr), mean = 0, sd = sbr)
ADS.para <- ADS.mod.nc.posterior.df %>%
  filter(coefficient == "b") %>%
  mutate(br = rep(brsim, 4),
         b = sample + br) %>% # add variation across tanks to posterior slope
  select(-c(coefficient, sample)) %>%
  left_join(ADS.mod.nc.posterior.df %>%
              filter(coefficient == "a") %>%
              select(-c(coefficient, level)),
            by = "n") %>%
  rename(a = "sample") %>%
  mutate(µ = -(a / b))


# Contrasts
ADS.comp <- ADS.para %>%
  select(-c(a, br)) %>%
  pivot_wider(names_from = level, values_from = c(b, µ)) %>%
  mutate(b1_2 = b_1 - b_2, b1_3 = b_1 - b_3, b1_4 = b_1 - b_4,
         b2_3 = b_2 - b_3, b2_4 = b_2 - b_4, b3_4 = b_3 - b_4,
         µ1_2 = µ_1 - µ_2, µ1_3 = µ_1 - µ_3, µ1_4 = µ_1 - µ_4,
         µ2_3 = µ_2 - µ_3, µ2_4 = µ_2 - µ_4, µ3_4 = µ_3 - µ_4)


# b
ADS.comp %>%
  filter(b1_2 < 0) %>% pull(b1_2) %>% length() /
  ADS.comp %>% pull(b1_2) %>% length()
# 99% probability that b1 < b2

ADS.comp %>%
  filter(b1_3 < 0) %>% pull(b1_3) %>% length() /
  ADS.comp %>% pull(b1_3) %>% length()
# 85% probability that b1 < b3

ADS.comp %>%
  filter(b1_4 > 0) %>% pull(b1_4) %>% length() /
  ADS.comp %>% pull(b1_4) %>% length()
# 56% probability that b1 > b4

ADS.comp %>%
  filter(b2_3 > 0) %>% pull(b2_3) %>% length() /
  ADS.comp %>% pull(b2_3) %>% length()
# 92% probability that b2 > b3

ADS.comp %>%
  filter(b2_4 > 0) %>% pull(b2_4) %>% length() /
  ADS.comp %>% pull(b2_4) %>% length()
# 99% probability that b2 > b4

ADS.comp %>%
  filter(b3_4 > 0) %>% pull(b3_4) %>% length() /
  ADS.comp %>% pull(b3_4) %>% length()
# 88% probability that b3 > b4

# µ
ADS.comp %>%
  filter(µ1_2 < 0) %>% pull(µ1_2) %>% length() /
  ADS.comp %>% pull(µ1_2) %>% length()
# 64% probability that µ1 < µ2

ADS.comp %>%
  filter(µ1_3 < 0) %>% pull(µ1_3) %>% length() /
  ADS.comp %>% pull(µ1_3) %>% length()
# 77% probability that µ1 > µ3

ADS.comp %>%
  filter(µ1_4 > 0) %>% pull(µ1_4) %>% length() /
  ADS.comp %>% pull(µ1_4) %>% length()
# 55% probability that µ1 > µ4

ADS.comp %>%
  filter(µ2_3 > 0) %>% pull(µ2_3) %>% length() /
  ADS.comp %>% pull(µ2_3) %>% length()
# 65% probability that µ2 > µ3

ADS.comp %>%
  filter(µ2_4 > 0) %>% pull(µ2_4) %>% length() /
  ADS.comp %>% pull(µ2_4) %>% length()
# 63% probability that µ2 > µ4

ADS.comp %>%
  filter(µ3_4 > 0) %>% pull(µ3_4) %>% length() /
  ADS.comp %>% pull(µ3_4) %>% length()
# 79% probability that µ3 > µ4

# Probability mass of b < 0
ADS.comp %>%
  filter(b_1 < 0) %>% pull(b_1) %>% length() /
  ADS.comp %>% pull(b_1) %>% length() # 96%

ADS.comp %>%
  filter(b_2 < 0) %>% pull(b_2) %>% length() /
  ADS.comp %>% pull(b_2) %>% length() # 61%

ADS.comp %>%
  filter(b_3 < 0) %>% pull(b_3) %>% length() /
  ADS.comp %>% pull(b_3) %>% length() # 88%

ADS.comp %>%
  filter(b_4 < 0) %>% pull(b_4) %>% length() /
  ADS.comp %>% pull(b_4) %>% length() # 97%

# Equations
ADS.coef <- ADS.para %>%
  group_by(level) %>%
  summarise(a.mean = mean(a),
            a.sd = sd(a),
            b.mean = mean(b),
            b.sd = sd(b),
            b.lwr = PI(b, prob = 0.9)[1],
            b.upr = PI(b, prob = 0.9)[2], 
            µ.median = median(µ),
            µ.mean = mean(µ)) %>% # quotient distributions are near impossible to 
  mutate(µ.mean.derived = -(a.mean / b.mean)) # describe with central tendencies
# a = 1.97 ± 0.52
# b1 = -0.18 ± 0.12
# b2 = -0.02 ± 0.11
# b3 = -0.11 ± 0.11
# b4 = -0.19 ± 0.12
# µ1 = 11
# µ2 = 94
# µ3 = 19
# µ4 = 10


ADS.annotation <- ADS.df %>%
  group_by(Treatment) %>%
  summarise(n = length(Treatment))
str(ADS.annotation)

ADS.annotation %<>%
  mutate(n = as.character(c(expression(italic("n ")*"= 26"),
                            expression(italic("n ")*"= 25"),
                            expression(italic("n ")*"= 25"),
                            expression(italic("n ")*"= 25"))),
         equation = as.character(c(expression(italic("y ")*"= "*frac(1, 1+e^{0.18*"("*italic("x ")*"– 11)"})),
                                   expression(italic("y ")*"= "*frac(1, 1+e^{0.02*"("*italic("x ")*"– 94)"})),
                                   expression(italic("y ")*"= "*frac(1, 1+e^{0.11*"("*italic("x ")*"– 19)"})),
                                   expression(italic("y ")*"= "*frac(1, 1+e^{0.19*"("*italic("x ")*"– 10)"})))))



# Lines and intervals

# calculate posterior prediction for new simulated Tanks, i.e. incorporating Tank s.d.

link_sim <- function(Treatment, Day){
  logit_p <- with(ADS.mod.nc.posterior, a + b[,Treatment] * Day + brsim * Day)
  return(logit_p)
}

ADS.mod.nc.logit_p <- predictor.df %$% mapply(link_sim, Treatment, Day)

# summarise posterior probability distributions
predictor.df$ADS.p.mean <- inv_logit(apply(ADS.mod.nc.logit_p, 2, mean))
predictor.df$ADS.p.pi.lwr <- inv_logit(t(apply(ADS.mod.nc.logit_p, 2, PI, prob = 0.5))[,1])
predictor.df$ADS.p.pi.upr <- inv_logit(t(apply(ADS.mod.nc.logit_p, 2, PI, prob = 0.5))[,2])
predictor.df$ADS.p.pi.lwr2 <- inv_logit(t(apply(ADS.mod.nc.logit_p, 2, PI, prob = 0.8))[,1])
predictor.df$ADS.p.pi.upr2 <- inv_logit(t(apply(ADS.mod.nc.logit_p, 2, PI, prob = 0.8))[,2])
predictor.df$ADS.p.pi.lwr3 <- inv_logit(t(apply(ADS.mod.nc.logit_p, 2, PI, prob = 0.9))[,1])
predictor.df$ADS.p.pi.upr3 <- inv_logit(t(apply(ADS.mod.nc.logit_p, 2, PI, prob = 0.9))[,2])


ADSp <- ggplot() +
  geom_dots(data = ADS.df, aes(Day, ADS, side = factor(ADS), fill = factor(Treatment)),
            scale = 0.25, colour = NA) +
  geom_line(data = predictor.df, aes(Day, ADS.p.mean, colour = factor(Treatment))) +
  geom_ribbon(data = predictor.df, aes(Day, ymin = ADS.p.pi.lwr, ymax = ADS.p.pi.upr,
                                       fill = factor(Treatment)), alpha = 0.5) +
  geom_ribbon(data = predictor.df, aes(Day, ymin = ADS.p.pi.lwr2, ymax = ADS.p.pi.upr2,
                                       fill = factor(Treatment)),  alpha = 0.4) +
  geom_ribbon(data = predictor.df, aes(Day, ymin = ADS.p.pi.lwr3, ymax = ADS.p.pi.upr3,
                                       fill = factor(Treatment)),  alpha = 0.3) +
  geom_point(data = ADS.coef %>% rename(Treatment = "level"),
             aes(µ.mean.derived, 0.5, fill = factor(Treatment)),
             shape = 21, size = 3.5, colour = "#ffffff") +
  geom_text(data = ADS.annotation, aes(120, 1.075, label = n),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = 1) +
  geom_text(data = ADS.annotation, aes(c(120, 0, 120, 120), c(0.88, 0.25, 0.88, 0.88), label = equation),
            family = "Helvetica Neue", parse = TRUE, size = 4.2, hjust = c(1, 0, 1, 1)) +
  scale_side_mirrored(guide = "none") +
  scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
  facet_grid(~Treatment, labeller = labeller(Treatment = Treatment_names)) +
  labs(y = expression("Daily autotrophic ("*italic(p)*")"),
       x = "Detrital age (d)") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 120),
                  clip = "off", expand = FALSE) +
  mytheme
ADSp


require(patchwork)
Fig.3 <- (ASp + theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank())) /
         (PSp + theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      strip.text = element_blank())) /
         (ADSp + theme(strip.text = element_blank())) +
  plot_annotation(tag_levels = "a") & # the & symbol is important so that theme is applied to all tags
  theme(plot.tag = element_text(family = "Helvetica Neue",
                                size = 20, face = "bold"))
Fig.3 # 8 x 13 in



# Contrasts between responses
S.comp <- AS.para %>%
  left_join(PS.para, by = c("n", "level")) %>%
  left_join(ADS.para, by = c("n", "level")) %>%
  select(level, a.x, a.y, a) %>%
  rename(a_AS = "a.x", a_PS = "a.y", a_ADS = "a") %>%
  filter(level == 1) %>% # all levels have the same intercept so select one
  mutate(aAS_PS = a_AS - a_PS, aAS_ADS = a_AS - a_ADS, aPS_ADS = a_PS - a_ADS)


S.comp %>%
  filter(aAS_PS < 0) %>% pull(aAS_PS) %>% length() /
  S.comp %>% pull(aAS_PS) %>% length()
# 94% probability that aAS < aPS

S.comp %>%
  filter(aAS_ADS > 0) %>% pull(aAS_ADS) %>% length() /
  S.comp %>% pull(aAS_ADS) %>% length()
# 92% probability that aAS > aADS

S.comp %>%
  filter(aPS_ADS > 0) %>% pull(aPS_ADS) %>% length() /
  S.comp %>% pull(aPS_ADS) %>% length()
# 100% probability that aPS > aADS















