ldhPlate_df <- read_csv("/Users/amjad_dabi/Google Drive/R Projects/Grace Data/LDH_x_DEP_exposure_20190606.csv")
str(ldhPlate_df)

ldhPlate2_df <- 
  ldhPlate_df %>% 
  mutate(only_HAPI  = (HAPI & !ENDO)) %>% 
  mutate(only_MICRO = (MICRO & !ENDO)) %>% 
  mutate(HAPI_ENDO  = (HAPI & ENDO)) %>% 
  mutate(MICRO_ENDO = (MICRO & ENDO)) %>% 
  select(-ENDO, -HAPI, -MICRO)

# sampleSizes_int <- 
#   ldhPlate2_df %>% 
#   filter(!LPS, !posCTRL, only_MICRO) %>% 
#   count(dose_ug) %>% 
#   pull(n)
# sampleSizes_int
# # 48 samples
# 
# # Create sample indices
# whichSamps <- sample(1:48, replace = TRUE)
# whichSamps


######  LDH MICRO  ############################################################
# ldhPlateMICRO_nCTRL_df <- 
#   ldhPlate2_df %>% 
#   filter(!LPS, !posCTRL, only_MICRO) %>% 
#   select(dose_ug, absorbance)
# 
# # Guts of loop:
# samples_idx <- unlist(
#   lapply(seq_along(sampleSizes_int), function(i){
#     
#     (sampleSizes_int[i] * (i - 1)) + 
#       sample(1:sampleSizes_int[i], replace = TRUE)
#     
#   })
# )
# 
# ldhMICROBoot_df <- 
#   ldhPlateMICRO_nCTRL_df[samples_idx, ]
# 
# # lm(log(absorbance) ~ dose_ug, data = ldhMICROBoot_df)
# coef(lm(log(absorbance) ~ dose_ug, data = ldhMICROBoot_df))[2]
# 
# 
# slopes_ls <- replicate(1000, expr = {
#   
#   samples_idx <- unlist(
#     lapply(seq_along(sampleSizes_int), function(i){
#       
#       (sampleSizes_int[i] * (i - 1)) + 
#         sample(1:sampleSizes_int[i], replace = TRUE)
#       
#     })
#   )
#   
#   ldhMICROBoot_df <- 
#     ldhPlateMICRO_nCTRL_df[samples_idx, ]
#   
#   # lm(log(absorbance) ~ dose_ug, data = ldhMICROBoot_df)
#   coef(lm(log(absorbance) ~ dose_ug, data = ldhMICROBoot_df))[2]
#   
# }, simplify = FALSE)
# slopesMicro_num <- unlist(slopes_ls)
# hist(slopesMicro_num)


######  LDH HAPI  #############################################################
sampleSizesHAPI_int <-
  ldhPlate2_df %>%
  filter(!LPS, !posCTRL, only_HAPI) %>%
  count(dose_ug) %>%
  pull(n)
sampleSizesHAPI_int

ldhPlateHAPI_nCTRL_df <-
  ldhPlate2_df %>%
  filter(!LPS, !posCTRL, only_HAPI) %>%
  select(dose_ug, absorbance)

slopesHAPI_ls <- replicate(1000, expr = {

  samples_idx <- unlist(
    lapply(seq_along(sampleSizesHAPI_int), function(i){

      (sampleSizesHAPI_int[i] * (i - 1)) +
        sample(1:sampleSizesHAPI_int[i], replace = TRUE)

    })
  )

  ldhHAPIBoot_df <-
    ldhPlateHAPI_nCTRL_df[samples_idx, ]

  coef(lm(log(absorbance) ~ dose_ug, data = ldhHAPIBoot_df))[2]

}, simplify = FALSE)
slopesHapi_num <- unlist(slopesHAPI_ls)
hist(slopesHapi_num)

t.test(slopesHapi_num, slopesMicro_num)
