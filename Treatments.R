require(tidyverse)
require(magrittr)
gis.df <- map_data(map = "world")

require(rgbif)
key <- name_backbone("Ecklonia Hornemann, 1828")$usageKey
occ_download(pred("taxonKey", key),
             pred_not(pred("verbatimScientificName", "Ecklonia maxima")),
             pred("occurrenceStatus", "PRESENT"),
             pred("hasCoordinate", TRUE),
             pred("geometry", "POLYGON((20 -60,180 -60,180 -20,20 -20,20 -60))"),
             format = "SIMPLE_CSV",
             user = "lukaseamus",
             email = "luka@wright.it",
             pwd = "wrighty4862")
# Citation: GBIF Occurrence Download https://doi.org/10.15468/dl.peg33d 
# Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-05-10
occ_download_wait("0265881-230224095556074")
gbif.df <- occ_download_get("0265881-230224095556074") %>%
  occ_download_import()
gbif.df %>% filter(verbatimScientificName == "Ecklonia maxima") 
gbif.df %>% filter(verbatimScientificName == "Ecklonia maxima (Osbeck) Papenf.")
# double-checked that there isn't any Ecklonia maxima

require(robis)
obis.df <- occurrence(scientificname = "Ecklonia radiata")

NPP.df <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/NPP Dataset.csv") %>%
  filter(Species == "Ecklonia radiata")
# NPP.df <- NPP.df %>% mutate(stdev_NPP_kg_C_m2_y = replace_na(stdev_NPP_kg_C_m2_y, 1e-10))
# NPP.l <- NPP.df %$% list(Reference = Reference,
#                          Latitude = Latitude_decimal_degrees - mean(Latitude_decimal_degrees),
#                          NPPmean = Avg_NPP_kg_C_m2_y - mean(Avg_NPP_kg_C_m2_y),
#                          NPPsem = stdev_NPP_kg_C_m2_y)
# NPP.l$N <- nrow(NPP.df)
# 
# knots <- 15
# knot_list <- quantile(NPP.l$Latitude , probs = seq(0, 1, length.out = knots))
# require(splines)
# B <- bs(NPP.l$Latitude,
#         knots = knot_list[-c(1, knots)],
#         degree = 3, intercept = TRUE)
# NPP.l$B <- B

# require(rethinking)
# NPP.mod <- ulam(
#   alist(
#     # NPPmean ~ dnorm(mean = NPP, sd = NPPsem),
#     # vector[N]:NPP ~ dnorm(mean = mu, sd = sigma),
#     NPPmean ~ dnorm(mean = mu, sd = sigma),
#     mu <- a + B %*% w,
#     a ~ dnorm(mean = 0, sd = 5),
#     w ~ dnorm(mean = 0, sd = 1),
#     sigma ~ dexp(rate = 1)
#   ), data = NPP.l, start = list(w = rep(0, ncol(B))), iter = 1e4, 
#      chains = 8, cores = parallel::detectCores()
# )




T.df <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Tolerance clean.csv") %>%
  group_by(Unit) %>%
  mutate(mean.s = (mean - mean(mean)) / sd(mean))


E.df <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Light.csv") %>%
  pivot_longer(cols = c(Ec.mean, Ek.mean, Ec.sem, Ek.sem), 
               names_to = c("Parameter", ".value"), # .value enables creation of multiple columns
               names_sep = 2) %>%
  rename(mean = ".mean", sem = ".sem")


P.df <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Prior.csv") %>%
  filter(Site == "Marmion") %>%
  mutate(GPmean = NPmean + Rmean,
         GPsem = sqrt(NPsem^2 + Rsem^2)) %>% # apply variance sum law
  pivot_longer(cols = c(NPmean, GPmean, NPsem, GPsem), 
               names_to = c("Parameter", ".value"), # .value enables creation of multiple columns
               names_sep = 2)





mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .5, .2, .2), "cm"),
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
                 legend.title = element_blank(),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 15, hjust = 0),
                 panel.spacing = unit(.5, "cm"),
                 text = element_text(family = "Helvetica Neue"))


map <- ggplot() +
  geom_point(data = gbif.df, aes(decimalLongitude,
                                 decimalLatitude),
             shape = 16, size = 4, colour = "#C3B300", alpha = 0.2) +
  geom_point(data = obis.df, aes(decimalLongitude, 
                                 decimalLatitude),
             shape = 16, size = 4, colour = "#C3B300", alpha = 0.2) +
  geom_map(data = gis.df, map = gis.df,
           colour = "#ffffff", fill = "#363538", size = 0.1,
           aes(map_id = region)) +
  annotate("point", x = 115.679017, y = -31.789367, shape = 1, size = 4) +
  annotate("segment", x = 110, xend = 114.5, y = -31.789367, 
           yend = -31.789367) +
  annotate("text", x = 97, y = -31.789367, label = "study site", 
           size = 4.2, family = "Helvetica Neue", hjust = 0) +
  scale_y_continuous(breaks = seq(-50, -10, by = 20), expand = c(0, 0),
                     labels = c("50°S", "30°S", "10°S")) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_fixed(xlim = c(10, 180), ylim = c(-50, -10)) +
  mytheme +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 1),
        plot.margin = unit(c(.2, 0, 0, .2), "cm"))

map

NPP.mod <- loess(Avg_NPP_kg_C_m2_y ~ Latitude_decimal_degrees, 
                 span = 1.5, data = NPP.df)
NPP.opt <- data.frame(Latitude_decimal_degrees = seq(-45, -32, 0.01)) %>% # loess does not extrapolate so stay within range
  mutate(Avg_NPP_kg_C_m2_y = predict(NPP.mod, newdata = Latitude_decimal_degrees))

NPP.opt %>% 
  filter(Avg_NPP_kg_C_m2_y == max(Avg_NPP_kg_C_m2_y)) # NPP is maximal at 34.72°S

NPP <- ggplot() +
  geom_point(data = NPP.df, aes(x = Latitude_decimal_degrees, 
                                y = Avg_NPP_kg_C_m2_y),
             shape = 16, size = 2, colour = "#C3B300", alpha = 0.2) +
  geom_smooth(data = NPP.df, aes(x = Latitude_decimal_degrees, 
                                 y = Avg_NPP_kg_C_m2_y),
              method = "loess", se = FALSE, span = 1.5, colour = "#C3B300") +
  geom_rug(data = NPP.opt %>% filter(Avg_NPP_kg_C_m2_y == max(Avg_NPP_kg_C_m2_y)), 
           aes(Latitude_decimal_degrees), colour = "#C3B300", length = unit(.25, "cm")) +
  ylab(expression(italic("NPP ")*"(kg C m"^-2*" yr"^-1*")")) +
  scale_x_continuous(breaks = seq(-50, -10, by = 20)) +
  scale_y_continuous(breaks = seq(0, 2, by = 1)) +
  coord_flip(xlim = c(-50, -10), ylim = c(0, 2), 
             expand = FALSE, clip = "off") +
  mytheme +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.length.y = unit(-0.25, "cm"),
        plot.margin = unit(c(.2, .5, 0, 0), "cm"))
NPP


T.mod <- loess(mean.s ~ Temperature, 
               span = 0.8, data = T.df)
T.opt <- data.frame(Temperature = seq(10, 30, 0.01)) %>%
  mutate(mean.s = predict(T.mod, newdata = Temperature))

T.opt %>% 
  filter(mean.s == max(mean.s)) # performance is maximal at 20.45°C

Ts <- ggplot() +
  geom_vline(xintercept = c(15, 20, 25), colour = "#DADADC") +
  geom_point(data = T.df, aes(Temperature, mean.s),
             colour = "#C3B300", shape = 16, size = 2, alpha = 0.2) +
  geom_smooth(data = T.df, aes(Temperature, mean.s), method = "loess", 
              colour = "#C3B300", se = FALSE, span = 0.8) +
  geom_rug(data = T.opt %>% filter(mean.s == max(mean.s)), 
           aes(Temperature), colour = "#C3B300", length = unit(.25, "cm")) +
  labs(x = expression("Temperature (°C)"^-1), # "^-1" added for spacing
       y = expression("Performance ("*italic(z)*"-score)")) +
  scale_x_continuous(labels = c("0°C", "10°C", "20°C", "30°C")) +
  scale_y_continuous(breaks = seq(-4, 16, 4)) +
  coord_cartesian(xlim = c(0, 30), ylim = c(-4, 16), expand = FALSE, clip = "off") +
  mytheme +
  theme(plot.margin = unit(c(0, .5, 0, .2), "cm"),
        axis.title.x = element_text(colour = NA)) # colour = NA rather than element_blank for spacing
Ts



options(pillar.sigfig = 4) # 4 significant figures
E.df %>%
  filter(!is.na(mean)) %>%
  group_by(Parameter) %>%
  summarise(µ = mean(mean),
            σ = sd(mean))
# Parameter     µ      σ
# Ec        13.82  9.240
# Ek        78.91  44.03 

require(ggridges)
require(ggnewscale)

L <- ggplot() +
  geom_vline(xintercept = 8, colour = "#DADADC") +
  geom_density_ridges_gradient(data = E.df %>% filter(Parameter == "Ek"),
                               aes(x = mean, y = 2, fill = factor(stat(quantile)), colour = Parameter),
                               scale = 15, bandwidth = 8, calc_ecdf = TRUE,
                               jittered_points = TRUE, point_alpha = 0.2, point_shape = 16, point_size = 2,
                               position = "raincloud", quantile_lines = TRUE, quantiles = c(0.05, 0.5, 0.95)) +
  scale_fill_manual(values = alpha("#C3B300", c(0.2, 0.5, 0.5, 0.2)), guide = "none") +
  new_scale_fill() +
  geom_density_ridges_gradient(data = E.df %>% filter(Parameter == "Ec"),
                               aes(x = mean, y = 1, fill = factor(stat(quantile)), colour = Parameter),
                               scale = 15, bandwidth = 8, calc_ecdf = TRUE,
                               jittered_points = TRUE, point_alpha = 0.2, point_shape = 16, point_size = 2,
                               position = "raincloud", quantile_lines = TRUE, quantiles = c(0.05, 0.5, 0.95)) +
  scale_fill_manual(values = alpha("#686000", c(0.2, 0.5, 0.5, 0.2)), guide = "none") +
  geom_segment(data = E.df %>%
                 filter(!is.na(mean)) %>%
                 group_by(Parameter) %>%
                 summarise(value = c(mean(mean), quantile(mean, probs = c(0.05, 0.5, 0.95))),
                           stat = c("mean", rep("quantile", 3))),
               aes(x = value, xend = value, y = c(1.545, rep(1, 3), 2.115, rep(2, 3)),
                   yend = rep(c(0.54, 1.54), each = 4), colour = Parameter, lty = stat)) +
  scale_linetype_manual(values = c(5, 1), guide = "none") +
  scale_colour_manual(values = c("#686000", "#C3B300"),
                      labels = c(expression(italic("E"[c])), expression(italic("E"[k])))) +
  xlab(expression(italic("E ")*"("*mu*"mol photons m"^-2*" s"^-1*")")) +
  coord_flip(xlim = c(0, 300), expand = FALSE, clip = "off") +
  mytheme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = c(0.1, 0.93),
        plot.margin = unit(c(0, .5, 0, .2), "cm"))
L



NP.mod <- loess(mean ~ Temperature, 
                span = 1, data = P.df %>% filter(Parameter == "NP"))
GP.mod <- loess(mean ~ Temperature, 
                span = 1, data = P.df %>% filter(Parameter == "GP"))

P.pred <- data.frame(Temperature = seq(10, 30, 0.01)) %>%
  mutate(NP = predict(NP.mod, newdata = Temperature),
         GP = predict(GP.mod, newdata = Temperature)) %>%
  pivot_longer(cols = c(NP, GP), names_to = "Parameter",
               values_to = "mean")

P.pred %>% 
  filter(Parameter == "NP") %>%
  filter(mean == max(mean)) %>%
  pull(Temperature) # performance is maximal at 23.51°C

P.pred %>% 
  filter(Parameter == "GP") %>%
  filter(mean == max(mean)) %>%
  pull(Temperature) # performance is maximal at 26.6°C


P.opt <- P.pred %>%
  group_by(Parameter) %>%
  filter(mean == max(mean)) %>%
  summarise(Temperature = Temperature)

P <- ggplot() +
  geom_vline(xintercept = c(15, 20, 25), colour = "#DADADC") +
  geom_pointrange(data = P.df, aes(Temperature, mean, ymin = mean - sem,
                                   ymax = mean + sem, shape = Reference, colour = Parameter),
                  size = 0.5, position = position_dodge(width = 2)) +
  geom_smooth(data = P.df, aes(Temperature, mean, colour = Parameter), method = "loess",
              se = FALSE, span = 1) +
  geom_rug(data = P.opt, aes(Temperature, colour = Parameter), length = unit(.25, "cm")) +
  scale_colour_manual(values = c("#686000", "#C3B300"),
                      labels = c(expression(italic(P[g])), expression(italic(P[n])))) +
  scale_shape_manual(values = c(15, 16), guide = "none") +
  scale_x_continuous(labels = c("10°C", "15°C", "20°C", "25°C", "30°C")) +
  labs(x = "Temperature (°C)",
       y = expression(italic(P[max]*" ")*"(mg O"[2]*" g"^-1*" h"^-1*")")) +
  coord_cartesian(xlim = c(8, 30), ylim = c(0, 5), expand = FALSE, clip = "off") +
  mytheme +
  theme(legend.position = c(0.1, 0.93),
        plot.margin = unit(c(0, .5, 0, .2), "cm"),
        axis.title.x = element_blank())
P

E <- ggplot() +
  geom_vline(xintercept = c(15, 20, 25), colour = "#DADADC") +
  geom_hline(yintercept = 8, colour = "#DADADC") +
  geom_pointrange(data = E.df %>% 
                    filter(Reference == "Staehr & Wernberg 2009",
                           Site == "Marmion Lagoon, Western Australia, Australia"), 
                  aes(Temperature, mean, ymin = mean - sem,
                      ymax = mean + sem, colour = fct_rev(Parameter)), 
                  shape = 15, size = 0.5) +
  geom_smooth(data = E.df %>% 
                filter(Reference == "Staehr & Wernberg 2009",
                       Site == "Marmion Lagoon, Western Australia, Australia"), 
              aes(Temperature, mean, colour = fct_rev(Parameter)), method = "loess",
              se = FALSE, span = 1.2) +
  scale_colour_manual(values = c("#C3B300", "#686000"),
                      labels = c(expression(italic(E[k])), expression(italic(E[c])))) +
  scale_x_continuous(labels = c("10°C", "15°C", "20°C", "25°C", "30°C")) +
  labs(x = "Temperature (°C)",
       y = expression(italic("E ")*"("*mu*"mol photons m"^-2*" s"^-1*")")) +
  coord_cartesian(xlim = c(8, 30), ylim = c(0, 150), expand = FALSE, clip = "off") +
  mytheme +
  theme(plot.margin = unit(c(0, .5, 0, .2), "cm"),
        axis.title.x = element_blank(),
        legend.position = c(0.1, 0.93))
E













# data(cherry_blossoms)
# d <- cherry_blossoms
# precis(d)
# d2 <- d[complete.cases(d$temp),]
# 
# knots <- 15
# knot_list <- quantile(d2$year , probs = seq(0, 1, length.out = knots))
# require(splines)
# B <- bs(d2$year,
#         knots = knot_list[-c(1, knots)],
#         degree = 3, intercept = TRUE)
# 
# 
# m4.7 <- quap(
#   alist(
#     Temp ~ dnorm(mu , sigma),
#     mu <- a + B %*% w ,
#     a ~ dnorm(6,10),
#     w ~ dnorm(0,1),
#     sigma ~ dexp(1)
#   ),
#   data=list( Temp=d2$temp , B=B ) ,
#   start=list( w=rep( 0 , ncol(B))))
# precis(m4.7)
# 
# 
# 
# d3 <- d2 %>% mutate(B = B)
# 
# require(brms)
# b4.8 <- 
#   brm(data = d3,
#       family = gaussian,
#       temp ~ 1 + B,
#       prior = c(prior(normal(6, 10), class = Intercept),
#                 prior(normal(0, 10), class = b),
#                 prior(exponential(1), class = sigma)),
#       iter = 2000, warmup = 1000, chains = 4, cores = 4,
#       seed = 4)
# stancode(b4.8)
# 
# 
# d <- data.frame(x = rep(seq(0, 5 * pi, length.out = 100), 2),
#                 group = rep(c("a", "b"), each = 100))
# d$y[1:100] <- 2 * sin(2 * d$x[1:100]) + rnorm(100, 0, 1)
# d$y[101:200] <- 5 * sin(3 * d$x[101:200]) + rnorm(100, 0, 0.5)
# 
# knots <- 15
# knot_list <- quantile(d$x , probs = seq(0, 1, length.out = knots))
# require(splines)
# B <- bs(d$x,
#         knots = knot_list[-c(1, knots)],
#         degree = 3, intercept = TRUE)
# 
# require(rethinking)
# r.mod <- quap(
#   alist(
#     y ~ dnorm(mu, sigma),
#     mu <- a + B %*% w,
#     a ~ dnorm(0, 1),
#     w ~ dnorm(0, 1),
#     sigma ~ dexp(1)
#   ), data = list(y = d$y, B = B),
#   start = list(w = rep(0, ncol(B)))
# )
# precis(r.mod) # model produces parameter estimates
# 
# 
# r.mod.stan <- ulam( # note ulam instead of quap
#   alist(
#     y ~ dnorm(mu, sigma),
#     mu <- a + B %*% w,
#     a ~ dnorm(0, 1),
#     w ~ dnorm(0, 1),
#     sigma ~ dexp(1)
#   ), data = list(y = d$y, B = B),
#   start = list(w = rep(0, ncol(B))),
#   chains = 8, cores = 8, iter = 1e4
# )
# 
# # model fails due to apparently semantic error:
# # Ill-typed arguments supplied to assignment operator =: lhs has type real and rhs has type row_vector
# 
# require(brms)
# d$B <- B
# b.mod.stan <- brm(data = d,
#                   family = gaussian,
#                   y ~ 1 + B,
#                   prior = c(prior(normal(0, 1), class = Intercept),
#                             prior(normal(0, 1), class = b),
#                             prior(exponential(1), class = sigma)),
#                   iter = 1e4, warmup = 5e3, chains = 4, cores = 4,
#                   seed = 4)
# summary(b.mod.stan) # same estimates as quap model




### Realised treatments ###

Tr.df <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Temperature.csv")

# Make numeric time sequence
start <- mdy_hms("06/29/22 12:00:00", tz = "Australia/Perth")
end <- Tr.df %>% filter(row_number() == n()) %>% 
  mutate(Datetime = paste(Date, Time)) %$% 
  mdy_hms(Datetime, tz = "Australia/Perth")
experiment <- start %--% end
experiment / ddays() # experiment ran for 119 d
experiment / dminutes(15) + 1 # temperature was measured in 15-min intervals,
# so the experimental interval would be split into 11425 measurements (1 is added
# because it is the difference between the number of intervals and the number of cases)

# but temperature recording was started a few hours late
startT <- Tr.df %>% filter(row_number() == 1) %>% 
  mutate(Datetime = paste(Date, Time)) %$% 
  mdy_hms(Datetime, tz = "Australia/Perth")
record <- startT %--% end
record / ddays() # record ran for 117.8333 d
record / dminutes(15) + 1 # = 11313 measurements in reality

# this can be cross-validated
Tr.df %>% filter(Tank == 2) %$% length(Temperature) # = 11313

# while all recordings started at the same time, they did not finish at the same time
# while we could calculate the number of days of each recording using lubridate, 
# the easiest way forward is to summarise the lengths of each Tank vector
Tr.df %>% group_by(Tank) %>%
  summarise(length = length(Temperature))
# only Tank 12 has a shorter measurement series: 3459
# again, cross-validate
end12 <- Tr.df %>% filter(Tank == 12) %>% filter(row_number() == n()) %>%
  mutate(Datetime = paste(Date, Time)) %$% 
  mdy_hms(Datetime, tz = "Australia/Perth")
record12 <- startT %--% end12
record12 / dminutes(15) + 1 # = 3459 measurments

# now we create a sequence for the duration of the 119-d experiment
d <- seq(0, 119, length.out = 11425)
d[2] - d[1] # the interval between measurements = 0.01041667 d
# this can again be cross-validated
dminutes(15) / ddays() # 15 min = 0.01041667 d


# now we trim the first part to match the temperature record
d <- d[-c(1:(11425 - 11313))] 
d[1] # the first temperature measurement was taken 
# 1.166667 d after the start of the experiment
# again, cross-validate
experiment / ddays() - record / ddays() # = 1.166667 d
length(d) - (record / dminutes(15) + 1) # d now has the correct length
d[11313] # = 119 d
# cross-validate 
record / ddays() + d[1] # = 119 d

# now we trim the last part to match the temperature record in Tank 12
d12 <- d[1:3459]
d12[3459] # = 37.1875 d
# cross-validate
record12 / ddays() + d12[1] # = 37.1875 d

# add d and d12 to dataframe as new variable
# double-check dimensions
length(d) * 7 + length(d12) - T.df %$% length(Temperature)

Tr.df <- Tr.df %>%
  mutate(Day = c(rep(d, 5), d12, rep(d, 2)))

# plot
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .5, .2, .2), "cm"),
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
                 legend.title = element_blank(),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 15, hjust = 0),
                 panel.spacing = unit(.5, "cm"),
                 text = element_text(family = "Helvetica Neue"))

T. <- ggplot() +
  geom_line(data = Tr.df, aes(x = Day, y = Temperature, colour = factor(Tank))) +
  mytheme
T.
# the initial spike in Tank 6 and 8 is due to the loggers being left outside until
# 2nd July at 7:45 and measuring air rather than seawater temperature

start6.8 <- mdy_hms("07/02/22 07:45:00", tz = "Australia/Perth")
record6.8 <- start6.8 %--% end
record6.8 / ddays() # = 116.1771 d
record6.8 / dminutes(15) + 1 # = 11154 measurements
  
# filter out the erroneous measurements
Tr.df <- Tr.df %>% 
  filter(Tank %in% c(6, 8)) %>%
  group_by(Tank) %>%
  slice_tail(n = record6.8 / dminutes(15) + 1) %>%
  bind_rows(Tr.df %>% filter(!Tank %in% c(6, 8))) %>%
  arrange(Tank)

Tr.df <- Tr.df %>%
  mutate(Treatment = ifelse(Tank %in% c(2, 4), "Light 15°C",
                            ifelse(Tank %in% c(6, 8), "Dark 15°C",
                                   ifelse(Tank %in% c(10, 12), "Light 20°C", "Light 25°C"))),
         Order = ifelse(Tank %in% c(2, 6, 10, 14), "Second tank", "Fourth tank"))

Tr <- ggplot() +
  geom_hline(yintercept = c(15, 20, 25), colour = "#DADADC") +
  geom_line(data = Tr.df, aes(x = Day, y = Temperature, 
                             colour = Treatment, alpha = fct_rev(Order))) +
  scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"),
                      guide = guide_legend(nrow = 2)) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_y_continuous(labels = c("10°C", "15°C", "20°C", "25°C", "30°C")) +
  scale_x_continuous(breaks = seq(0, 120, by = 15)) + 
  xlab("Time since start of experiment (d)") +
  coord_cartesian(ylim = c(10, 30), xlim = c(0, 120), expand = FALSE) +
  mytheme +
  theme(legend.position = c(0.85, 0.93),
        legend.box = "horizontal",
        plot.margin = unit(c(0, 0, .2, .2), "cm"),
        axis.title.y = element_blank())
Tr

require(ggridges)
require(ggnewscale)
Tdens <- ggplot() +
  geom_vline(xintercept = c(15, 20, 25), colour = "#DADADC") +
  geom_density_ridges_gradient(data = Tr.df %>% filter(Treatment == "Light 25°C"),
                               aes(x = Temperature, y = 1.1, fill = factor(stat(quantile)), colour = Treatment),
                               scale = 1.3, bandwidth = 0.2, calc_ecdf = TRUE, quantile_lines = TRUE, 
                               quantiles = c(0.05, 0.5, 0.95)) +
  scale_fill_manual(values = alpha("#d1750c", c(0.2, 0.5, 0.5, 0.2)), guide = "none") +
  new_scale_fill() +
  geom_density_ridges_gradient(data = Tr.df %>% filter(Treatment == "Light 20°C"),
                               aes(x = Temperature, y = 1.05, fill = factor(stat(quantile)), colour = Treatment),
                               scale = 1.3, bandwidth = 0.2, calc_ecdf = TRUE, quantile_lines = TRUE, 
                               quantiles = c(0.05, 0.5, 0.95)) +
  scale_fill_manual(values = alpha("#f5a54a", c(0.2, 0.5, 0.5, 0.2)), guide = "none") +
  new_scale_fill() +
  geom_density_ridges_gradient(data = Tr.df %>% filter(Treatment == "Light 15°C"),
                               aes(x = Temperature, y = 1, fill = factor(stat(quantile)), colour = Treatment),
                               scale = 1.3, bandwidth = 0.2, calc_ecdf = TRUE, quantile_lines = TRUE, 
                               quantiles = c(0.05, 0.5, 0.95)) +
  scale_fill_manual(values = alpha("#6a98b4", c(0.2, 0.5, 0.5, 0.2)), guide = "none") +
  new_scale_fill() +
  geom_density_ridges_gradient(data = Tr.df %>% filter(Treatment == "Dark 15°C"),
                               aes(x = Temperature, y = 0, fill = factor(stat(quantile)), colour = Treatment),
                               scale = 1.3, bandwidth = 0.2, calc_ecdf = TRUE, quantile_lines = TRUE, 
                               quantiles = c(0.05, 0.5, 0.95)) +
  scale_fill_manual(values = alpha("#2e4a5b", c(0.2, 0.5, 0.5, 0.2)), guide = "none") +
  geom_segment(data = Tr.df %>%
                 group_by(Treatment) %>%
                 summarise(mean = mean(Temperature)),
               aes(x = mean, xend = mean, y = c(0, 1, 1.05, 1.1),
                   yend = c(0.98, 1.97, 2.64, 1.975), colour = Treatment), lty = 5) +
  scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"),
                      guide = "none") +
  coord_flip(xlim = c(10, 30), expand = FALSE, clip = "off") +
  mytheme +
  theme(plot.margin = unit(c(0, .5, .2, 0), "cm"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.length.y = unit(0, "cm"), # axis.ticks.y = element_blank() does not remove the space for some reason
        axis.line.y = element_blank(),
        plot.background = element_blank())
Tdens




# calculate descriptive stats
options(pillar.sigfig = 4) # 4 significant figures
Tr.df %>% 
  group_by(Treatment) %>%
  summarise(mean = mean(Temperature),
            sd = sd(Temperature),
            n = length(Temperature))
# Treatment   mean     sd     n
# Dark 15°C  15.23 0.7255 22308
# Light 15°C 15.39 0.5981 22626
# Light 20°C 19.87 0.3054 14772
# Light 25°C 24.90 1.174  22626

Er.df <- read.csv("~/Desktop/PhD/Papers/Detrital photosynthesis/Data/Emax.csv")

start <- dmy("29.06.22", tz = "Australia/Perth")

Er.df <- Er.df %>% 
  mutate(Day = start %--% dmy(Date, tz = "Australia/Perth") / ddays(),
         Treatment = ifelse(Tank %in% 1:4, "Light 15°C",
                            ifelse(Tank %in% 5:8, "Dark 15°C",
                                   ifelse(Tank %in% 9:12, "Light 20°C", "Light 25°C"))),
         Order = ifelse(Tank %in% c(1, 5, 9, 13), "First tank", 
                        ifelse(Tank %in% c(2, 6, 10, 14), "Second tank", 
                               ifelse(Tank %in% c(3, 7, 11, 15), "Third tank", "Fourth tank"))))

ggplot() +
  geom_point(data = Er.df %>% filter(Treatment != "Dark 15°C"), 
             aes(x = Day, y = PAR, colour = Treatment, shape = Order)) +
  mytheme
# too sparse to plot by day
# note dips in PAR which are due to measurements not being timed correctly 
# in relation to the peak of the parabola

# remove days 6, 13, 82 and 92 where PAR was not measured at noon
ggplot() +
  geom_density(data = Er.df %>% filter(Treatment != "Dark 15°C",
                                      !Day %in% c(6, 13, 82, 92)), 
               aes(x = PAR, colour = Treatment), bw = 1) +
  mytheme
# not worth a plot

# just calculate descriptive stats
options(pillar.sigfig = 3) # 3 significant figures
Er.df %>% 
  filter(Treatment != "Dark 15°C",
         !Day %in% c(6, 13, 82, 92)) %>% 
  group_by(Treatment) %>%
  summarise(mean = mean(PAR),
            sd = sd(PAR),
            n = length(PAR))
# Treatment   mean    sd     n
# Light 15°C  7.88 0.719    16
# Light 20°C  8    0.918    20
# Light 25°C  7.75 0.967    20

Er.df %>% 
  filter(Treatment != "Dark 15°C",
         !Day %in% c(6, 13, 82, 92)) %>% 
  summarise(mean = mean(PAR),
            sd = sd(PAR),
            n = length(PAR))
# mean         sd  n
# 7.875 0.8751623 56




require(patchwork)
( (map | NPP) + plot_layout(widths = c(1, 0.18)) ) / (Ts | L | P | E) / ( (Tr | Tdens) + plot_layout(widths = c(1, 0.1)) ) + 
  plot_annotation(tag_levels = list(c("a", "", "b", "", "c", "", "d"))) & # the & symbol is important so that theme is applied to all tags
  theme(plot.tag = element_text(family = "Helvetica Neue",
                                size = 20, face = "bold")) # 9.78 x 13 in

