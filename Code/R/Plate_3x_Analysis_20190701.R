# Well Plate Full Data Analysis
# Gabriel Odom
# 2019-07-01


######  Overview  #############################################################
# Grace called me on 29 April asking why her data wasn't normal and what
#   methods can she use on this data. I asked her for the raw data to explore.
#   I sent described to her the exploration and analysis in
#   Code/Plate_Analysis_20190430.R on 2019-06-06.
# The data cleaning was handled in Code/Plate_3x_Cleaning_20190606.R. This
#   script analyse the ROS, LDH, and MTT data, as they all share similar forms.

library(tidyverse)
.dataDir <- "Google Drive/R Projects/Grace Data/Data/clean2019/"
dataDir <- .dataDir

#


######  ROS  ##################################################################
###  Read Clean Data  ###
rosPlate_df <- read_csv(paste0(dataDir, "ROS_x_DEP_exposure_20190606.csv"))
str(rosPlate_df)


###  Inspect Fluorescence  ###
summary(rosPlate_df$fluorescence)
hist(rosPlate_df$fluorescence)
hist(log(rosPlate_df$fluorescence))
# This is clearly bi-modal, but we'll see if the residuals are as well. The
#   treatment may very well be the cause.

plot(
  2 * rosPlate_df$ENDO + 3 * rosPlate_df$MICRO + 4 * rosPlate_df$HAPI
)
# There are 5 total groups: ENDO, ENDO + HAPI, HAPI, ENDO + MICRO, and MICRO.
#   Let's make ENDO the baseline.
rosPlate2_df <- 
  rosPlate_df %>% 
  mutate(only_HAPI  = (HAPI & !ENDO)) %>% 
  mutate(only_MICRO = (MICRO & !ENDO)) %>% 
  mutate(HAPI_ENDO  = (HAPI & ENDO)) %>% 
  mutate(MICRO_ENDO = (MICRO & ENDO)) %>% 
  select(-ENDO, -HAPI, -MICRO)

rosPlate2_df %>% 
  select(-fluorescence, -dose_ug) %>% 
  rowSums() %>% 
  plot
# There are rows with all 0 values. These are the ENDO cells alone (no CTRLs).


###  Linear Model  ###
ros_mod <- lm(log(fluorescence) ~ ., data = rosPlate2_df)
par(mfrow = c(2, 2))
plot(ros_mod)
par(mfrow = c(1, 1))

summary(ros_mod)
# The intercept has the ENDO effect. The control effects make sense, and we see
#   a positive (but small) dose effect regardless of the cell type. Further, it
#   appears that changing from ENDO to any other cell type--or even adding
#   another cell type--decreases fluorescence. The worst offender is MICRO.



###  Plot Dose Effect within Cell Type  ###
rosPlate3_df <- 
  rosPlate_df %>%  
  mutate(only_HAPI  = (HAPI & !ENDO)) %>% 
  mutate(only_MICRO = (MICRO & !ENDO)) %>%
  mutate(only_ENDO  = (ENDO & (!MICRO & !HAPI))) %>%  
  mutate(HAPI_ENDO  = (HAPI & ENDO)) %>% 
  mutate(MICRO_ENDO = (MICRO & ENDO)) %>% 
  select(-MICRO, -HAPI, -ENDO) %>% 
  gather(Type, toKeep, -negCTRL, -posCTRL, -dose_ug, -fluorescence) %>% 
  filter(toKeep) %>%
  select(-toKeep)

# Treatment vs Control
rosPlate3_df %>% 
  filter(!posCTRL) %>% 
  filter(!negCTRL) %>%
  mutate(
    dose_ug = as.factor(dose_ug)
  ) %>% 
  mutate(Type = factor(Type)) %>% 
  ggplot() +
    aes(x = dose_ug, y = fluorescence) +
    scale_y_log10() +
    geom_boxplot(aes(fill = Type))
# we see a small dose effect (as expected)


# Dose Effect within Cell Type
rosPlate3_df %>% 
  filter(!negCTRL) %>% 
  filter(!posCTRL) %>% 
  mutate(Type = factor(Type)) %>% 
  ggplot() +
    aes(x = dose_ug, y = fluorescence, group = Type) +
    scale_y_log10() +
    geom_jitter(aes(color = Type), width = 0.1, alpha = 0.5) +
    stat_smooth(aes(color = Type), method = "lm", se = FALSE) +
    ggtitle("Log Fluorescence by Dose over Cell Type")

# As mentioned previously, the dose effect is small, positive, significant, and
#   is expressed across groups. However, this effect is dwarfed by the cell type
#   effect.



###  Repeated Analysis After Cell-type Normalization  ###
# Grace called me this morning (2019-07-02) and asked me to "standardize" the 
#   response treatment values within cell type to the "blank control". In math,
#   this means that
#   for each cell type {
#     calculate the mean response at dose = 0 (ignoring al CTRLs); save as mu
#     for each dose > 0 {
#       divide each response by mu 
#     }
#   }
# This will transform all response values to percent changes from the dose = 0
#   "blank control".
rosBlankCTRL_df <- 
  rosPlate3_df %>% 
  filter(!posCTRL) %>% 
  filter(!negCTRL) %>% 
  filter(dose_ug == 0) %>% 
  group_by(Type) %>% 
  summarise(
    meanFLUO   = mean(fluorescence),
    stdDevFLUO = sd(fluorescence)
  ) %>% 
  mutate(Type = as.character(Type))

rosActiveDFs_ls <- 
  rosPlate3_df %>% 
  filter(!posCTRL) %>% 
  filter(!negCTRL) %>% 
  select(-posCTRL, -negCTRL) %>% 
  split(., .$Type)

rosTreat_df <- 
  lapply(names(rosActiveDFs_ls), function(name){
    # browser()
    
    mu_num <- 
      rosBlankCTRL_df %>% 
      filter(Type == name) %>% 
      select(meanFLUO) %>% 
      pull
    sigma_num <- 
      rosBlankCTRL_df %>% 
      filter(Type == name) %>% 
      select(stdDevFLUO) %>% 
      pull
    rosActiveDFs_ls[[name]] %>% 
      mutate(scaledFLUO = fluorescence / mu_num) %>% 
      mutate(centredFLUO = (fluorescence - mu_num) / sigma_num) %>% 
      select(-fluorescence)
    
  }) %>% 
  bind_rows()


###  Linear Model for "Scaled" Fluorescence  ###
# FLUO / ave(FLUO) within Cell Type
rosScale_mod <- lm(log(scaledFLUO) ~ . -centredFLUO, data = rosTreat_df)
par(mfrow = c(2, 2))
plot(rosScale_mod)
par(mfrow = c(1, 1))
# Well that's an absolute nightmare. EDIT: REMOVE fluorescence YOU IDIOT.

summary(rosScale_mod)
# Great. Small dose effect. Small ENDO and ENDO + MICRO effect. Really low R^2
#   now (after removing FLUO).

# normalized FLUO within Cell Type
rosCentred_mod <- lm(centredFLUO ~ . -scaledFLUO, data = rosTreat_df)
par(mfrow = c(2, 2))
plot(rosCentred_mod)
par(mfrow = c(1, 1))
# Worse residuals

hist(residuals(rosCentred_mod))

summary(rosCentred_mod)
# But better fit.


###  Plots for "Scaled" Fluorescence  ###
# Treatment vs Control
rosTreat_df %>% 
  mutate(
    dose_ug = as.factor(dose_ug)
  ) %>% 
  mutate(Type = factor(Type)) %>% 
  ggplot() +
    aes(x = dose_ug, y = scaledFLUO) +
    geom_boxplot(aes(fill = Type)) +
    ggtitle("Fluorescence / Average Fluorescence within Cell Type by Dose")
# we see a small dose effect (as expected)
rosTreat_df %>% 
  mutate(
    dose_ug = as.factor(dose_ug)
  ) %>% 
  mutate(Type = factor(Type)) %>% 
  ggplot() +
    aes(x = dose_ug, y = centredFLUO) +
    geom_boxplot(aes(fill = Type)) +
    ggtitle("Within Cell-Type Standardized Fluorescence by Dose")
# When we standardize to the cell type, we see a strange effect from the
#   MICRO + ENDO cell combination. This combination, while middle-of-the-pack
#   for raw changes in fluorescence, has marked increase in standardized changes


# Dose Effect within Cell Type
rosTreat_df %>%
  mutate(Type = factor(Type)) %>% 
  ggplot() +
    aes(x = dose_ug, y = scaledFLUO, group = Type) +
    geom_jitter(aes(color = Type), width = 0.1, alpha = 0.5) +
    stat_smooth(aes(color = Type), method = "lm", se = FALSE) +
    ggtitle("Fluorescence / Average Fluorescence within Cell Type by Dose")
# Sharpest slopes are ENDO alone and ENDO + MICRO
rosTreat_df %>%
  mutate(Type = factor(Type)) %>% 
  ggplot() +
    aes(x = dose_ug, y = centredFLUO, group = Type) +
    geom_jitter(aes(color = Type), width = 0.1, alpha = 0.5) +
    stat_smooth(aes(color = Type), method = "lm", se = FALSE) +
    ggtitle("Within Cell-Type Standardized Fluorescence by Dose")


# rm(list = ls())
# dataDir <- .dataDir



######  LDH  ##################################################################
###  Read Clean Data  ###
ldhPlate_df <- read_csv(paste0(dataDir, "LDH_x_DEP_exposure_20190606.csv"))
str(ldhPlate_df)

###  Inspect absorbance  ###
summary(ldhPlate_df$absorbance)
hist(ldhPlate_df$absorbance)
hist(log(ldhPlate_df$absorbance))
# This is clearly bi-modal, but we'll see if the residuals are as well. The
#   treatment may very well be the cause.

plot(
  2 * ldhPlate_df$ENDO + 3 * ldhPlate_df$MICRO + 4 * ldhPlate_df$HAPI
)
# There are 5 total groups: ENDO, ENDO + HAPI, HAPI, ENDO + MICRO, and MICRO.
#   Let's make ENDO the baseline.
ldhPlate2_df <- 
  ldhPlate_df %>% 
  mutate(only_HAPI  = (HAPI & !ENDO)) %>% 
  mutate(only_MICRO = (MICRO & !ENDO)) %>% 
  mutate(HAPI_ENDO  = (HAPI & ENDO)) %>% 
  mutate(MICRO_ENDO = (MICRO & ENDO)) %>% 
  select(-ENDO, -HAPI, -MICRO)

ldhPlate2_df %>% 
  select(-absorbance, -dose_ug) %>% 
  rowSums() %>% 
  plot
# There are rows with all 0 values. These are the ENDO cells alone (no CTRLs).


###  Linear Model  ###
ldh_mod <- lm(
  log(absorbance) ~ .,
  data = ldhPlate2_df
)
par(mfrow = c(2, 2))
plot(ldh_mod)
par(mfrow = c(1, 1))
# The tails of those residuals are problematic...

# So, the residuals are not bi-modal, but they are lighter-tailed than they
#   should be.
par(mfrow = c(2, 1))
hist(
  rnorm(nrow(ldhPlate2_df), sd = sd(residuals(ldh_mod))),
  xlim = c(-1, 1),
  main = "Random Normal Density to Match Variance of LDH Model Residuals"
)
hist(residuals(ldh_mod), xlim = c(-1, 1), main = "LDH Model Residuals")
par(mfrow = c(1, 1))

summary(ldh_mod)
# The intercept has the ENDO effect. The control effects make sense, and we see
#   a positive (but small) dose effect regardless of the cell type. Further, it
#   appears that changing from ENDO to MICRO decreases response.



###  Plot Dose Effect within Cell Type  ###
ldhPlate3_df <- 
  ldhPlate_df %>%  
  mutate(only_HAPI  = (HAPI & !ENDO)) %>% 
  mutate(only_MICRO = (MICRO & !ENDO)) %>%
  mutate(only_ENDO  = (ENDO & (!MICRO & !HAPI))) %>%  
  mutate(HAPI_ENDO  = (HAPI & ENDO)) %>% 
  mutate(MICRO_ENDO = (MICRO & ENDO)) %>% 
  select(-MICRO, -HAPI, -ENDO) %>% 
  gather(Type, toKeep, -LPS, -posCTRL, -dose_ug, -absorbance) %>% 
  filter(toKeep) %>%
  select(-toKeep) %>% 
  mutate(Type = factor(Type))

# Treatment vs Control
ldhPlate3_df %>% 
  filter(!posCTRL) %>% 
  filter(!LPS) %>%
  mutate(
    dose_ug = as.factor(dose_ug)
  ) %>% 
  mutate(Type = factor(Type)) %>% 
  ggplot() +
  aes(x = dose_ug, y = absorbance) +
  scale_y_log10() +
  geom_boxplot(aes(fill = Type))
# we see a small dose effect for all three ENDO cell types


# Dose Effect within Cell Type
ldhPlate3_df %>% 
  filter(!LPS) %>% 
  filter(!posCTRL) %>% 
  ggplot() +
  aes(x = dose_ug, y = absorbance, group = Type) +
  scale_y_log10() +
  geom_jitter(aes(color = Type), width = 0.1, alpha = 0.5) +
  stat_smooth(aes(color = Type), method = "lm", se = FALSE) +
  ggtitle("Log absorbance by Dose over Cell Type")

# As mentioned previously, the dose effect is small, positive, significant, and
#   is expressed across the three cell types that included ENDO


###  Repeated Analysis After Cell-type Normalization  ###
ldhBlankCTRL_df <- 
  ldhPlate3_df %>% 
  filter(!posCTRL) %>% 
  filter(!LPS) %>% 
  filter(dose_ug == 0) %>% 
  group_by(Type) %>% 
  summarise(
    meanABSORB   = mean(absorbance),
    stdDevABSORB = sd(absorbance)
  ) %>% 
  mutate(Type = as.character(Type))

ldhActiveDFs_ls <- 
  ldhPlate3_df %>% 
  filter(!posCTRL) %>% 
  filter(!LPS) %>% 
  select(-posCTRL, -LPS) %>% 
  split(., .$Type)

ldhTreat_df <- 
  lapply(names(ldhActiveDFs_ls), function(name){
    # browser()
    
    mu_num <- 
      ldhBlankCTRL_df %>% 
      filter(Type == name) %>% 
      select(meanABSORB) %>% 
      pull
    sigma_num <- 
      ldhBlankCTRL_df %>% 
      filter(Type == name) %>% 
      select(stdDevABSORB) %>% 
      pull
    ldhActiveDFs_ls[[name]] %>% 
      mutate(scaledABSORB = absorbance / mu_num) %>% 
      mutate(centredABSORB = (absorbance - mu_num) / sigma_num) %>% 
      select(-absorbance)
    
  }) %>% 
  bind_rows()


###  Linear Model for "Scaled" absorbance  ###
# ABSORB / ave(ABSORB) within Cell Type
ldhScale_mod <- lm(log(scaledABSORB) ~ . -centredABSORB, data = ldhTreat_df)
par(mfrow = c(2, 2))
plot(ldhScale_mod)
par(mfrow = c(1, 1))
# Well that's an absolute nightmare

summary(ldhScale_mod)
# I don't want to interpret this at all.

# normalized ABSORB within Cell Type
ldhCentred_mod <- lm(centredABSORB ~ . -scaledABSORB, data = ldhTreat_df)
par(mfrow = c(2, 2))
plot(ldhCentred_mod)
par(mfrow = c(1, 1))
# Well that's also an absolute nightmare

summary(ldhCentred_mod)
# I don't want to interpret this one either.


###  Plots for "Scaled" absorbance  ###
# Treatment vs Control
ldhTreat_df %>% 
  mutate(
    dose_ug = as.factor(dose_ug)
  ) %>% 
  mutate(Type = factor(Type)) %>% 
  ggplot() +
  aes(x = dose_ug, y = scaledABSORB) +
  geom_boxplot(aes(fill = Type)) +
  ggtitle("absorbance / Average absorbance within Cell Type by Dose")
# we see a small dose effect for all the ENDO-type cells
ldhTreat_df %>% 
  mutate(
    dose_ug = as.factor(dose_ug)
  ) %>% 
  mutate(Type = factor(Type)) %>% 
  ggplot() +
  aes(x = dose_ug, y = centredABSORB) +
  geom_boxplot(aes(fill = Type)) +
  ggtitle("Within Cell-Type Standardized absorbance by Dose")
# Once again, we see a small dose effect for all the ENDO-type cells, but MICRO
#   may also be slightly sensitive to larger doses of treatment


# Dose Effect within Cell Type
ldhTreat_df %>%
  mutate(Type = factor(Type)) %>% 
  ggplot() +
  aes(x = dose_ug, y = scaledABSORB, group = Type) +
  geom_jitter(aes(color = Type), width = 0.1, alpha = 0.5) +
  stat_smooth(aes(color = Type), method = "lm", se = FALSE) +
  ggtitle("absorbance / Average absorbance within Cell Type by Dose")
# MICRO might be flat overall; HAPI decreases; all ENDO cell types have a tight
#   similar effect to treatment

# # DELETE FIGURE
# ldhTreat_df %>%
#   mutate(Type = factor(Type)) %>% 
#   ggplot() +
#   aes(x = dose_ug, y = centredABSORB, group = Type) +
#   geom_jitter(aes(color = Type), width = 0.1, alpha = 0.5) +
#   stat_smooth(aes(color = Type), method = "lm", se = FALSE) +
#   ggtitle("Within Cell-Type Standardized absorbance by Dose")



# rm(list = ls())
# dataDir <- .dataDir



######  MTT  ##################################################################
###  Read Clean Data  ###
mttPlate_df <- read_csv(paste0(dataDir, "MTT_x_DEP_exposure_20190606.csv"))
str(mttPlate_df)

###  Inspect absorbance  ###
summary(mttPlate_df$absorbance)
hist(mttPlate_df$absorbance)
# This is clearly bi-modal, but we'll see if the residuals are as well. The
#   treatment may very well be the cause.

plot(
  2 * mttPlate_df$ENDO + 3 * mttPlate_df$MICRO + 4 * mttPlate_df$HAPI
)
# There are 5 total groups: ENDO, ENDO + HAPI, HAPI, ENDO + MICRO, and MICRO.
#   Let's make ENDO the baseline.
mttPlate2_df <- 
  mttPlate_df %>% 
  mutate(only_HAPI  = (HAPI & !ENDO)) %>% 
  mutate(only_MICRO = (MICRO & !ENDO)) %>% 
  mutate(HAPI_ENDO  = (HAPI & ENDO)) %>% 
  mutate(MICRO_ENDO = (MICRO & ENDO)) %>% 
  select(-ENDO, -HAPI, -MICRO)

mttPlate2_df %>% 
  select(-absorbance, -dose_ug) %>% 
  rowSums() %>% 
  plot
# There are rows with all 0 values. These are the ENDO cells alone (no CTRLs).


###  Linear Model  ###
mtt_mod <- lm(absorbance ~ ., data = mttPlate2_df)
par(mfrow = c(2, 2))
plot(mtt_mod)
par(mfrow = c(1, 1))
# The "residuals v fitted" plot may show a problem.

# So, the residuals are not bi-modal, but they are lighter-tailed than they
#   should be.
par(mfrow = c(2, 1))
hist(
  rnorm(nrow(mttPlate2_df), sd = sd(residuals(mtt_mod))),
  xlim = c(-0.5, 0.5),
  main = "Random Normal Density to Match Variance of MTT Model Residuals"
)
hist(residuals(mtt_mod), xlim = c(-0.5, 0.5), main = "MTT Model Residuals")
par(mfrow = c(1, 1))
# Once again, the kurtosis is strong, but I don't know if that will yield
#   problems with the inference.

summary(mtt_mod)
# Is this a problem? The LPS and positive controls have negative signs compared
#   to ENDO.



###  Plot Dose Effect within Cell Type  ###
mttPlate3_df <- 
  mttPlate_df %>%  
  mutate(only_HAPI  = (HAPI & !ENDO)) %>% 
  mutate(only_MICRO = (MICRO & !ENDO)) %>%
  mutate(only_ENDO  = (ENDO & (!MICRO & !HAPI))) %>%  
  mutate(HAPI_ENDO  = (HAPI & ENDO)) %>% 
  mutate(MICRO_ENDO = (MICRO & ENDO)) %>% 
  select(-MICRO, -HAPI, -ENDO) %>% 
  gather(Type, toKeep, -LPS, -posCTRL, -dose_ug, -absorbance) %>% 
  filter(toKeep) %>%
  select(-toKeep) %>% 
  mutate(Type = factor(Type))

# Treatment vs Control
mttPlate3_df %>% 
  filter(!posCTRL) %>% 
  filter(!LPS) %>%
  mutate(
    dose_ug = as.factor(dose_ug)
  ) %>% 
  mutate(Type = factor(Type)) %>% 
  ggplot() +
  aes(x = dose_ug, y = absorbance) +
  geom_boxplot(aes(fill = Type))
# We see a small negative dose effect for ENDO + HAPI (and ENDO). Technically, 
#   the dose effect isn't significant at 0.05, so this isn't real. Moreover,
#   there may be confounding effects due to cell type. For this metric in
#   particular, scaling or centering within cell type will be valuable.


# Dose Effect within Cell Type
mttPlate3_df %>% 
  filter(!LPS) %>% 
  filter(!posCTRL) %>% 
  ggplot() +
  aes(x = dose_ug, y = absorbance, group = Type) +
  geom_jitter(aes(color = Type), width = 0.1, alpha = 0.5) +
  stat_smooth(aes(color = Type), method = "lm", se = FALSE) +
  ggtitle("Log absorbance by Dose over Cell Type")

# It looks like the potential negative effect to dose is driven entirely by the
#   ENDO + HAPI group.



###  Repeated Analysis After Cell-type Normalization  ###
mttBlankCTRL_df <- 
  mttPlate3_df %>% 
  filter(!posCTRL) %>% 
  filter(!LPS) %>% 
  filter(dose_ug == 0) %>% 
  group_by(Type) %>% 
  summarise(
    meanABSORB   = mean(absorbance),
    stdDevABSORB = sd(absorbance)
  ) %>% 
  mutate(Type = as.character(Type))

mttActiveDFs_ls <- 
  mttPlate3_df %>% 
  filter(!posCTRL) %>% 
  filter(!LPS) %>% 
  select(-posCTRL, -LPS) %>% 
  split(., .$Type)

mttTreat_df <- 
  lapply(names(mttActiveDFs_ls), function(name){
    # browser()
    
    mu_num <- 
      mttBlankCTRL_df %>% 
      filter(Type == name) %>% 
      select(meanABSORB) %>% 
      pull
    sigma_num <- 
      mttBlankCTRL_df %>% 
      filter(Type == name) %>% 
      select(stdDevABSORB) %>% 
      pull
    mttActiveDFs_ls[[name]] %>% 
      mutate(scaledABSORB = absorbance / mu_num) %>% 
      mutate(centredABSORB = (absorbance - mu_num) / sigma_num) %>% 
      select(-absorbance)
    
  }) %>% 
  bind_rows()


###  Linear Model for "Scaled" absorbance  ###
# ABSORB / ave(ABSORB) within Cell Type
mttScale_mod <- lm(log(scaledABSORB) ~ . -centredABSORB, data = mttTreat_df)
par(mfrow = c(2, 2))
plot(mttScale_mod)
par(mfrow = c(1, 1))
# Well that's an absolute nightmare

summary(mttScale_mod)
# I don't want to interpret this at all. The dose effect completely vanishes

# normalized ABSORB within Cell Type
mttCentred_mod <- lm(centredABSORB ~ . -scaledABSORB, data = mttTreat_df)
par(mfrow = c(2, 2))
plot(mttCentred_mod)
par(mfrow = c(1, 1))
# Well that's also an absolute nightmare

summary(mttCentred_mod)
# I don't want to interpret this one either. The dose effect is back?


###  Plots for "Scaled" absorbance  ###
# Treatment vs Control
mttTreat_df %>% 
  mutate(
    dose_ug = as.factor(dose_ug)
  ) %>% 
  mutate(Type = factor(Type)) %>% 
  ggplot() +
  aes(x = dose_ug, y = scaledABSORB) +
  geom_boxplot(aes(fill = Type)) +
  ggtitle("absorbance / Average absorbance within Cell Type by Dose")
# Variance kills...
mttTreat_df %>% 
  mutate(
    dose_ug = as.factor(dose_ug)
  ) %>% 
  mutate(Type = factor(Type)) %>% 
  ggplot() +
  aes(x = dose_ug, y = centredABSORB) +
  geom_boxplot(aes(fill = Type)) +
  ggtitle("Within Cell-Type Standardized absorbance by Dose")
# To me, this plot is the most telling. There may be some overall negative
#   effect of the dose, but all normalized values are within 2.5 standard
#   deviations. It appears that dose has no effect on MTT within these cell
#   types.


# Dose Effect within Cell Type
mttTreat_df %>%
  mutate(Type = factor(Type)) %>% 
  ggplot() +
  aes(x = dose_ug, y = scaledABSORB, group = Type) +
  geom_jitter(aes(color = Type), width = 0.1, alpha = 0.5) +
  stat_smooth(aes(color = Type), method = "lm", se = FALSE) +
  ggtitle("absorbance / Average absorbance within Cell Type by Dose")
# Negative dose effect for ENDO and HAPI + ENDO. Positive dose effect for
#   MICRO + ENDO (but no effect for MICRO). Once again, recall the p-value of
#   the dose effet for this model: 0.966.

# DELETE FIGURE
mttTreat_df %>%
  mutate(Type = factor(Type)) %>%
  ggplot() +
  aes(x = dose_ug, y = centredABSORB, group = Type) +
  geom_jitter(aes(color = Type), width = 0.1, alpha = 0.5) +
  stat_smooth(aes(color = Type), method = "lm", se = FALSE) +
  ggtitle("Within Cell-Type Standardized absorbance by Dose")
# Negative dose effect for ENDO and HAPI + ENDO. It's super small though.

# rm(list = ls())

