########################### preparation ###########################
library(tidyverse)
# if default install gives wrong version, verify using direct v0.2.2 install
# pak::pak("spCF@0.2.2")
library(spCF)
library(stringr)
#library(fastDummies)
library(ggcorrplot)
library(viridis)
library(sf)
library(ggpattern)

# set the working directory to the main R project directory
setwd(here::here())

########################### correlation plot ###########################
seFin <- read_csv("seFin.csv")

corrPlot <- ggcorrplot(cor(seFin %>% 
                 select(c("samplingYr",
                          "samplingM",
                          "samplingSiteX",
                          "samplingSiteY",
                          "maxRich",
                          "murphyRich",
                          "tmiRich",
                          "tmi",
                          "observationCount",
                          "waterDepthMax")) %>% 
                   rename(Year = "samplingYr",
                          Month = "samplingM",
                          X = "samplingSiteX",
                          Y = "samplingSiteY",
                          `Full Richness` = "maxRich",
                          `Murphy Richness` = "murphyRich",
                          `TMI Richness` = "tmiRich",
                          TMI = "tmi",
                          Count = "observationCount",
                          `Survey Depth` = "waterDepthMax")),
           type = "upper",
           lab = TRUE,
           digits = 2) +
  scale_fill_gradientn(colors = viridis(256, option = 'D'))

#### now for water quality too
# no need to use both max and murphy richness --> better to just use just murphy
# TMI strongly related to alkalinity, conductivity, and pH
# TMI richness weakly related to just pH
# otherwise, no major correlations between chemistry and biodiv (pH is the highest)
# Alk, pH, and conductivity are very close
# TP and TN are also close (and somewhat linked to Chl-a)
# TOC somewhat linked to TN
corrPlotFull <- ggcorrplot(cor(seFin %>% 
                 select(c("samplingYr",
                          "samplingM",
                          "samplingSiteX",
                          "samplingSiteY",
                          "maxRich",
                          "murphyRich",
                          "tmiRich",
                          "tmi",
                          "observationCount",
                          "waterDepthMax",
                          "Alk",
                          "Abs_F420",
                          "Kfyll",
                          "Kond_25",
                          "TOC",
                          "Tot_P",
                          "Temp",
                          "pH",
                          "TN")) %>% 
                   relocate(c(samplingM, samplingYr,
                              samplingSiteX, samplingSiteY,
                              maxRich, murphyRich, tmiRich, tmi,
                              observationCount, waterDepthMax,
                              Abs_F420, Alk, Kfyll, Kond_25, TOC, pH,
                              Tot_P, Temp, TN)) %>% 
                          rename(Year = "samplingYr",
                                 Month = "samplingM",
                                 X = "samplingSiteX",
                                 Y = "samplingSiteY",
                                 `Full Richness` = "maxRich",
                                 `Murphy Richness` = "murphyRich",
                                 `TMI Richness` = "tmiRich",
                                 TMI = "tmi",
                                 `# Obs` = "observationCount",
                                 `Survey Depth` = "waterDepthMax",
                                 Absorbance = "Abs_F420",
                                 Alkalinity = "Alk",
                                 `Chlorophyl-a` = "Kfyll",
                                 Conductivity = "Kond_25",
                                 TOC = "TOC",
                                 TP = "Tot_P",
                                 `T` = "Temp",
                                 pH = "pH",
                                 TN = "TN"),
               # keep NA/missing values (i.e., get full correlations for richnesses even without chemistry)
               use = "pairwise.complete.obs"),
           type = "upper",
           lab = TRUE,
           digits = 2) +
  scale_fill_viridis_c(begin = 0.3) +
  theme(legend.position = "none")

ggsave("figS1.jpeg",
       corrPlotFull,
       dpi = 300,
       scale = 2.7,
       width = 8,
       height = 8.3,
       units = "cm")

########################### null model on full data ###########################

# basic information from example at https://dmuraka.r-universe.dev/articles/spCF/spCF_glm.html
# fit the actual CF-GLMM
# specify model elements separately
coords <- seFin %>% select(c(samplingSiteX, samplingSiteY)) %>% as.matrix()
x <- seFin %>% select(c(samplingYr)) %>% 
  # code below to make dummy variables, unused
  #dummy_cols(select_columns = "samplingYr",
  # otherwise perfect multicollinearity, 2007 as baseline
  #           remove_first_dummy = TRUE) %>% 
  #select(-samplingYr) %>% 
  as.matrix()
y <- seFin %>% pull(maxRich)
offset <- seFin %>% pull(observationCount)

# fit the basic, full model with holdout validation
mod_hv <- cf_glm_hv(y = y, #x = x, 
                    #offset = offset, # looks like model doesn't find suitable coefficients with offset
                    coords = coords, 
                    family = negbin())

# then train the full model using the cf_glm function
mod <- cf_glm(y = y, #x = x, 
              #offset = offset,
              coords = coords, mod_hv = mod_hv)

# no summary function; inspect results
mod

############### residual plots ############### 

# basic residuals check
nullResults <- seFin %>% 
  bind_cols(# Predictive mean
    mod$pred$pred, 
    # Predictive SD
    mod$pred$pred_sd) %>% 
  rename(predRich = `...39`,
         predSd = `...40`) %>% 
  mutate(rawResid = maxRich - predRich)

# residuals vs fitted
residFit <- ggplot(nullResults, aes(x = predRich, y = rawResid)) + 
  geom_point() +
  theme_bw() +
  labs(x = "Error",
       y = "Predicted richness")

# temporal sensitivity
residTemp <- ggplot(nullResults, aes(x = samplingYr, y = rawResid)) + 
  geom_point() +
  theme_bw() +
  scale_x_continuous(minor_breaks = seq(2007,2024,
                                        by = 1)) +
  labs(x = "Year",
       y = "Error")

# spatial sensitivity
# make data sf
nullSpatial <- st_as_sf(nullResults, 
                        coords = c("samplingSiteX", "samplingSiteY"), 
                        crs = 3006)

# get basemap as terra SpatVector polygon
seMap <- geodata::gadm(country = "SWE", level = 0, path = tempdir())
# make SWEREF 99
seMap <- st_as_sf(seMap, crs = 4326) %>% 
  st_transform(crs = 3009)

# spatially plot residuals
residSpat <- ggplot() + 
  geom_sf(data = seMap) +
  geom_sf(data = nullSpatial, aes(colour = rawResid)) +
  theme_bw() +
  scale_colour_viridis_c() +
  labs(colour = "Error")

# format plots in replicable cowplot format
resids <- cowplot::plot_grid(residFit, 
                             cowplot::plot_grid(residTemp, residSpat,
                                                labels = c('b', 'c'), 
                                                label_size = 12,
                                                rel_widths = c(1.5, 1),
                                                align = "h"
                                                ),
                             labels = c("a", ""),
                             label_size = 12,
                             nrow = 2,
                             align = "h") 
cowplot::save_plot("fig4.jpeg", resids, 
                   nrow = 2,
                   dpi = 300,
                   base_width = 7)

############### spatial feature extraction ############### 
# we can (and should) change the bandwidth on these to reflect the model results
# can then extract the specific preditive means for just those bandwidth features of a certain scale
mod_s1 <- sp_scalewise(mod,bw_range=c(100000,Inf)) # Large scale (100+ km), 15 scales - 1
mod_s2 <- sp_scalewise(mod,bw_range=c(10000,100000)) # medium scale (10-100 km), 21 scales - 2
mod_s3 <- sp_scalewise(mod,bw_range=c(0,10000)) # medium scale (under 10 km), 14 scales - 2

nullPlot <- nullSpatial %>% 
  bind_cols(mod_s1$pred$pred,
            mod_s2$pred$pred,
            mod_s3$pred$pred) %>% 
  rename(large = `...41`,
         medium = `...42`,
         small = `...43`) %>% 
  pivot_longer(cols = c("small", "medium", "large"),
               names_to = "scale",
               values_to = "predSpat") %>% 
  mutate(scale = factor(scale,
                        levels = c("small", "medium", "large")))

spatDecompPlot <- ggplot() + 
  geom_sf(data = seMap) +
  geom_sf(data = nullPlot, aes(colour = predSpat)) +
  theme_bw() +
  scale_colour_viridis_c() +
  facet_wrap(~scale) +
  labs(colour = "Adjustment \nfrom baseline")

ggsave("fig3.jpeg", 
       spatDecompPlot,
       width = 24,
       height = 17,
       units = "cm",
       #scale = 1,
       dpi = 300)

############### spatial process checks ############### 
# the goal here is to find the covariates that distort the spatial process the least
# we evaluate these by 4 metrics: (1) min, (2) max bandwidth
# (3) average distance between bandwidth
# (4) correlation between current model predictions and those of the null model

# save the model output
nullResRef <- nullResults %>% 
  select(c("samplingYr", "samplingSiteId", "predRich", "predSd")) %>% 
  # rename to simplify things later
  rename(predRichF = "predRich",
         predSdF = "predSd")

# first, get these for the null model
subsetSum <- data.frame(
  # number of observations
  n = nrow(nullResRef),
  # min bandwidth
  minBand = min(mod$bands), 
  # max bandwidth
  maxBand = max(mod$bands), 
  # number of unique bands
  nBand = length(mod$bands),
  # geometric average difference between bands
  # geometric mean because of the non-linear scaling of band increases
  bandDiff = mod$bands %>% sort() %>% diff() %>% log() %>% mean() %>% exp(),
  # correlation between predicted values
  modCorr = cor(nullResRef$predRichF, nullResRef$predRichF),
  # predictive score metrics
  r2 = mod$e_summary %>% filter(str_detect(stat, "R2")) %>% pull(value),
  rmse = mod$e_summary %>% filter(str_detect(stat, "RMSE")) %>% pull(value),
  mae = mod$e_summary %>% filter(str_detect(stat, "MAE")) %>% pull(value),
  # model name
  mod = "Full"
)

# get a list of all unique covariate combinations
varList <- c("secchi", "waterDepthMax", "Abs_F420",
             "Alk", "Kfyll", "Kond_25", "TOC", "Tot_P",
             "Temp", "pH", "TN") %>% sort()

varComb <- combinations <- unlist(
  lapply(1:length(varList), function(i) {
    combn(varList, i, simplify = FALSE)
  }), 
  recursive = FALSE
)

# build helper function to get all the output
# variables must be vector of strings
getModSum <- function(variables){
  seFilt <- seFin %>% select(all_of(c("samplingSiteX", "samplingSiteY", "maxRich",
                                      # ID info for later
                                      "samplingYr", "samplingSiteId",
                                      variables))) %>% 
    drop_na()
  # specify model elements separately
  coordsTest <- seFilt %>% select(c(samplingSiteX, samplingSiteY)) %>% as.matrix()
  yTest <- seFilt %>% pull(maxRich)
  
  # fit the basic, full model with holdout validation
  test_hv <- cf_glm_hv(y = yTest,
                       coords = coordsTest, 
                       family=negbin()) %>% 
    # quiet the output
    suppressMessages()
  
  # then train the full model using the cf_glm function
  testMod <- cf_glm(y = yTest,
                    coords = coordsTest, mod_hv = test_hv) %>% 
    # quiet the output
    suppressMessages()
  
  testRes <- seFilt %>% 
    # select bare minimum for consistent column numbers
    select("samplingSiteX", "samplingSiteY", "maxRich",
           "samplingYr", "samplingSiteId") %>% 
    bind_cols(
      testMod$pred$pred, # Predictive mean
      testMod$pred$pred_sd # Predictive SD
    ) %>% 
    rename(predRich = `...6`,
           predSd = `...7`) %>% 
    select(c("samplingYr", "samplingSiteId", "predRich", "predSd")) %>% 
    left_join(nullResRef, by = c("samplingYr", "samplingSiteId")) %>% 
    # quiet the output
    suppressMessages()
  
  # do the same summary as before
  return(
    data.frame(
      n = nrow(testRes),
      minBand = min(testMod$bands),
      maxBand = max(testMod$bands),
      nBand = length(testMod$bands),
      bandDiff = testMod$bands %>% sort() %>% diff() %>% log() %>% mean() %>% exp(),
      modCorr = cor(testRes$predRich, testRes$predRichF),
      r2 = testMod$e_summary %>% filter(str_detect(stat, "R2")) %>% pull(value),
      rmse = testMod$e_summary %>% filter(str_detect(stat, "RMSE")) %>% pull(value),
      mae = testMod$e_summary %>% filter(str_detect(stat, "MAE")) %>% pull(value),
      mod = paste(variables, collapse = ";")
    ) %>% 
      # quiet the output
      suppressMessages()
  )
}

for(i in 1:length(varComb)){
  invisible({capture.output({
    tempSum <- getModSum(varComb[[i]])
  })})
  # combine with total summary
  subsetSum <- subsetSum %>% bind_rows(tempSum)
}

cleanSubsetSum <- subsetSum %>% 
  group_by(across(all_of(colnames(subsetSum %>% select(-mod))))) %>% 
  slice_max(nchar(mod), n = 1)

# export for writing
write.csv(cleanSubsetSum,
          "tableS2.csv")

########################### filter w/ chemistry data ###########################
# re-run with just water sites with water chemistry data
# method is sensitive to NAs which must be removed
# filter only for values not NA in selected variables
# done manually but these can be found in entry 3 of cleanSubsetSum
# pH;Temp;TN;TOC;Tot_P;waterDepthMax
seFinFilt <- seFin %>% drop_na(any_of(c(
  #"secchi",
  #"Abs_F420",
  #"Alk",
  #"Kfyll",
  #"Kond_25",
  "pH",
  "Temp",
  "TN",
  "TOC",
  "Tot_P",
  "waterDepthMax"))) %>% 
  select(-c("secchi",
            "Abs_F420",
            "Alk",
            "Kfyll",
            "Kond_25",
            #"TOC",
            #"Tot_P",
            #"Temp",
            #"pH",
            #"TN",
            #"waterDepthMax"
  )) %>% 
  mutate(
    pH_sq = pH^2,
    Temp_sq = Temp^2,
    TN_sq = TN^2,
    TOC_sq = TOC^2,
    Tot_P_sq = Tot_P^2,
    waterDepthMax_sq = waterDepthMax^2
  )

# prepare data for package
coordFilt <- seFinFilt %>% select(c(samplingSiteX, samplingSiteY)) %>% as.matrix()
yFilt <- seFinFilt %>% pull(maxRich)

# fit the null model again
nullFiltHV <- cf_glm_hv(y = yFilt,
                        coords = coordFilt, 
                        family = negbin())
nullFilt <- cf_glm(y = yFilt,
                   coords = coordFilt, 
                   mod_hv = nullFiltHV)

# general inspection
#nullFilt
# spatial trends are in
#nullFilt$bands
# covariates are in nullFilt$beta
#nullFilt$beta
# summary statistics in
#nullFilt$e_summary

# basic residuals check
nullFiltResults <- seFinFilt %>% 
  bind_cols(# Predictive mean
    nullFilt$pred$pred, 
    # Predictive SD
    nullFilt$pred$pred_sd) %>% 
  rename(predRich = `...40`,
         predSd = `...41`) %>% 
  mutate(rawResid = maxRich - predRich)


# residuals vs fitted
# more systematic overprediction than underprediction but nothing else really?
ggplot(nullFiltResults, aes(x = predRich, y = rawResid)) + geom_point()
# absolutely lethal misfit where residuals  increase as richness increases
# this is pretty typical for a misfit model and not uniquely horrible
ggplot(nullFiltResults, aes(x = maxRich, y = rawResid)) + geom_point()

# temporal sensitivity
# middle years are especially undersampled but nothing structural
ggplot(nullFiltResults, aes(x = samplingYr, y = rawResid)) + geom_point()

# spatial sensitivity
nullFiltSpatial <- st_as_sf(nullFiltResults, 
                            coords = c("samplingSiteX", "samplingSiteY"), 
                            crs = 3006)
# nothing hugely concerning spatially
# slightly higher errors in the mountains and around cities but nothing huge
ggplot() + 
  geom_sf(data = seMap) +
  geom_sf(data = nullFiltSpatial, aes(colour = rawResid)) +
  theme_bw() +
  scale_colour_viridis_c() +
  labs(colour = "Raw residuals")

########################### manual chem model dredge ###########################

# explore the variety of combined models
# helper function to run things faster
runQuickCF <- function(covList) {
  coordFilt <- seFinFilt %>% select(c(samplingSiteX, samplingSiteY)) %>% as.matrix()
  yFilt <- seFinFilt %>% pull(maxRich)
  xFilt <- seFinFilt %>% select(all_of(covList)) %>% as.matrix()
  
  modHV <- cf_glm_hv(y = yFilt, x = xFilt, 
                     coords = coordFilt, 
                     family = negbin()) %>% 
    # silence output
    suppressMessages()
  
  modFit <- cf_glm(y = yFilt, x = xFilt, 
                   coords = coordFilt, 
                   mod_hv = modHV) %>% 
    suppressMessages()
  
  return(modFit)
}

# make heatmap for visualisation
# empty df with each column as some bandwidth
modSums <- data.frame(band = nullFiltHV$other$bands_all)
# add null model
modSums <- modSums %>% 
  mutate(Null = case_when(
    band %in% nullFilt$band ~ nullFilt$e_summary %>% 
      filter(str_detect(stat, "R2")) %>% 
      pull(value),
    .default = 0
  ))

# more summary statistics
# all relevant to the validation set
# just get it directly from the summary object
modScores <- as.data.frame(nullFilt$e_summary) %>%
  pivot_wider(names_from = stat,
              values_from = value) %>% 
  mutate(mod = "Null")

# finally the fixed effects
# left empty because we don't care about the null model
modFE = data.frame(#var = covListFull,
  coef = numeric(),
  isSig = numeric(),
  mod = character()
)

#get full list of covariates
finVars <- c("pH",
             "Temp",
             "TN",
             "TOC",
             "Tot_P",
             "waterDepthMax")

# generate a matrix of states to simplify computation, where
# 0 = omit; 1 = incl linear term only; 2 = incl linear & sqare terms
varsGrid <- expand.grid(rep(list(0:2), length(finVars)))

chemVarsList <- apply(varsGrid, # apply over the grid
                      1, # for rows
                      function(i) {
                        unlist(Map(function(v, # for some variable
                                            s) { # and for some state
                          if (s == 0) {
                            # include nothing
                            return(NULL) 
                          }
                          else if (s == 1) {
                            # include variable linearly
                            return(v) 
                          }
                          else {
                            # include variable linearly and as square
                            c(v, paste0(v, "_sq")) 
                          }
                        }, finVars, i)) %>% 
                          # remove element names
                          unname() 
                      })

# loop through all our candidate sets
# start at 2 because 1 is our null model
for(i in 2:length(chemVarsList)){
  # get the set of covariates we're working with now
  covVect <- chemVarsList[[i]]
  # run the model
  modTemp <- runQuickCF(covVect)
  # add to our temp spatial results vector
  modSums <- modSums %>% 
    mutate(tempBand = case_when(
      band %in% modTemp$band ~ modTemp$e_summary$value[1],
      .default = 0
    ))
  # rename temp column name with model name
  colnames(modSums)[length(colnames(modSums))] <- paste(covVect,
                                                        collapse = ";")
  
  # get FE coefficients and significance
  # 1-variable case is a little different
  if(length(covVect) == 1){
    feTemp <- modTemp$beta %>% 
      # model output no longer provides any variable names
      # assume that first row of the beta matrix is always the intercept
      slice(-1) %>% 
      mutate(var = covVect[1],
             isSig = as.numeric(sign(lower_95CI) == sign(upper_95CI))) %>% 
      select(var, coef, isSig)
  } else {
    # otherwise for multiple coef models
    feTemp <- modTemp$beta %>% rownames_to_column() %>% rename(var = "rowname") %>% 
      filter(str_length(var) > 1) %>% 
      mutate(var = str_remove_all(var, "x"),
             isSig = as.numeric(sign(lower_95CI) == sign(upper_95CI))) %>% 
      select(var, coef, isSig)
  }
  
  # add column indicating the model name (just in case)
  feTemp <- feTemp %>% mutate(mod = paste(covVect,
                                          collapse = ";"))
  
  # combine with temp counter
  #modFE <- left_join(modFE, feTemp, by = "var")
  modFE <- bind_rows(modFE, feTemp)
  
  # finally, add the model scores
  modScoreTemp <- as.data.frame(modTemp$e_summary) %>%
    pivot_wider(names_from = stat,
                values_from = value) %>% 
    mutate(mod = paste(covVect,
                       collapse = ";"))
  # append to the dataset
  modScores <- bind_rows(modScores, modScoreTemp)
}

########################### clean & plot results ###########################

# clean up model scores
modScoresExp <- modScores %>% 
  # regex magic to clean things up
  mutate(modClean = str_replace_all(mod,
                                    # capture the complete second occurrence of a duplicate pair
                                    # leaving the underscore suffix
                                    "(^|;)([^;]+);(\\2(?:_\\d+)?)",
                                    # using backreferences to 1st and 3rd groups (though 3rd never exists)
                                    "\\1\\3") %>% 
           str_replace_all(";", " + ") %>% 
           str_replace_all("_sq", "^2") %>% 
           str_replace_all("Temp", "T") %>% 
           str_replace_all("Tot_P", "TP") %>% 
           str_replace_all("waterDepthMax", "Depth")
         ) %>% 
  select(-mod) %>% 
  relocate(modClean) %>% 
  mutate(across(where(is.numeric), ~round(., 7)))

# save for manuscript
# some weird floating point issues with long training 00000000003 but it's fine
write_csv(modScoresExp, "tableS3.csv")


# check the full FE numbers
modFE %>% 
  filter(var != "Intercept") %>%
  mutate(var = var %>% str_replace_all("waterDepthMa$", "waterDepthMax") %>% 
           str_replace_all("waterDepthMa_", "waterDepthMax_"),
         containsSq = str_detect(mod, paste0(var, "_sq")),
         isSigWith = as.numeric(isSig & containsSq)) %>%
  summarise(.by = var,
            n = n(),
            nSig = sum(isSig),
            nSigWith = sum(isSigWith))

# get only the models without square terms for reporting
modSumsVis <- modSums %>% select(!contains("_sq"))
modFeVis <- modFE %>% filter(!str_detect(mod, "_sq"))

# get the minimum distance for all columns
minDist <- modSumsVis %>%
  filter(if_any(-band, ~ . != 0)) %>%
  summarise(minDist = min(band)) %>%
  pull(minDist)

# set clean model order for readability (first step)
modelOrderClean <- colnames(modSumsVis)[-1] %>% 
  str_replace_all("Temp", "T") %>% 
  str_replace_all("Tot_P", "TP") %>% 
  str_replace_all("waterDepthMax", "Depth")

# format for plot
modOrderLab <- modelOrderClean %>% 
  # another regex to remove all the ; and duplicates with +
  # kind of irrelevant now that there are no square terms but it's fine
  str_replace_all("_sq", "(^2") %>% 
  str_remove_all("\\;[:alpha:]{1,2}\\(") %>%
  str_replace_all("\\;", " + ")

modelOrderFinal <- c(modelOrderClean %>% 
                       as_tibble() %>% 
                       filter(!str_detect(value, "pH")) %>% 
                       pull(value),
                     modelOrderClean %>% 
                       as_tibble() %>% 
                       filter(str_detect(value, "pH")) %>% 
                       pull(value))
modelLabelFinal <- c(modOrderLab %>% 
                       as_tibble() %>% 
                       filter(!str_detect(value, "pH")) %>% 
                       pull(value),
                     modOrderLab %>% 
                       as_tibble() %>% 
                       filter(str_detect(value, "pH")) %>% 
                       pull(value))

modPlot <- modSumsVis %>% 
  # remove all bandwidths below the minimum detection for all models
  filter(band >= minDist) %>%
  # convert from m to km
  mutate(band = band/1000) %>% 
  pivot_longer(cols = colnames(modSumsVis)[-1],
               names_to = "modName",
               values_to = "R2") %>% 
  # make names more consistent & readable
  mutate(modName = modName %>% 
           str_replace_all("Temp", "T") %>% 
           str_replace_all("Tot_P", "TP") %>% 
           str_replace_all("waterDepthMax", "Depth")) %>% 
  # format factor levels manually to make plot more readable
  mutate(distance = as.factor(format(round(band, 1), nsmall = 1)),
         modName = factor(modName,
                          levels = modelOrderFinal,
                          labels = modelLabelFinal))

# plot
spatPlot <- ggplot(modPlot, aes(modName, distance, fill = R2)) +
  geom_tile() +
  theme_bw() + 
  labs(x = "",
       y = "Bandwitdh (km)",
       fill = bquote(R^2)) +
  # make limits & labels more readable
  scale_fill_viridis_c(limits = c(0, 
                                  # get slightly above the rounded value just in case
                                  round(max(modPlot$R2),2) + 0.01),
                       # round sequence values to avoid trailing decimals
                       breaks = round(seq(0, round(max(modPlot$R2),2) + 0.01,
                                          length.out = 4), 2)) +
  scale_x_discrete(labels = scales::parse_format(),
                   sec.axis = dup_axis(labels = NULL)) +
  #scale_y_discrete(limits = rev) + 
  # make x axis readable
  theme(
    axis.text.x = element_text(angle = 55,
                               hjust = 1),
    plot.margin = margin(t = 0.1, b = 5.5, l = 5.5, r = 5.5)
  )

# add another heatmap for FE coefficients showing
modReFE <- modFeVis %>% 
  # make names more consistent & readable
  mutate(mod = mod %>% 
           str_replace_all("Temp", "T") %>% 
           str_replace_all("Tot_P", "TP") %>% 
           str_replace_all("waterDepthMax", "Depth"),
         # and again for the variables
         var = var %>% 
           str_replace_all("Temp", "T") %>% 
           str_replace_all("Tot_P", "TP") %>% 
           str_replace_all("waterDepthMax", "Depth") %>% 
           str_replace_all("waterDepthMa", "Depth"),
         isNeg = as.numeric(coef < 0)
  ) %>%
  # remove intercept
  filter(var != "Intercept")

modReFE <- modReFE %>%
  # add null model
  bind_rows(data.frame(coef = NA,
                       isSig = 0,
                       mod = "Null",
                       var = "TN",
                       meanEff = NA,
                       coefResc = NA,
                       isNeg = NA
                     )) %>% 
  mutate(# set factor levels to improve readability
    mod = factor(mod, 
                 levels = modelOrderFinal,
                 labels = modelLabelFinal),
    var = factor(var)) 
  

# plot
fePlot <- ggplot(data = modReFE, aes(mod, var, 
                                     pattern = as.character(isNeg), 
                                     fill = coef)) + 
  geom_tile() +
  # fill for the fixed effect estimates
  scale_fill_viridis_c(na.value = "transparent") +
  labs(fill = bquote(hat(beta))) +
  scale_y_discrete(labels = scales::parse_format()) +
  geom_tile_pattern(pattern_colour = NA,
                    pattern_fill = "grey",
                    pattern_angle = 45,
                    pattern_density = 0.5,
                    pattern_spacing = 0.025,
                    pattern_key_scale_factor = 1) +
  scale_pattern_manual(values = c("1" = "circle", "0" = "none"),
                       guide = "none") +
  geom_tile(data = modReFE %>% filter(isSig == 1), 
            aes(mod, var,
                # d is 0 = insignificant; 1 = significant
                colour = as.factor(isSig)), 
            linewidth = 2, fill=NA) +
  # add outline to boxes corresponding to significant coefficients
  scale_colour_manual(values = c("grey"),
                      guide = "none") +
  scale_x_discrete(sec.axis = dup_axis(labels = NULL)) +
  theme_bw() +
  labs(x = "",
       y = "") + 
  theme(axis.text.x = element_blank(),
        plot.margin = margin(t = 5.5, b = 0.1, l = 5.5, r = 5.5)
  ) 

modSum <- cowplot::plot_grid(fePlot, spatPlot, 
                   labels = c("a", "b"),
                   ncol = 1,
                   rel_heights = c(0.4, 1.4),
                   align = "v") 
cowplot::save_plot("fig5.jpeg", modSum, 
          nrow = 2,
          dpi = 300,
          base_height = 5,
          base_asp = 2.2)







########################### old stuff (not used) ###########################
### unused diagnostics
# raw vs fitted; quasi-QQ plot
ggplot(nullResults, aes(x = maxRich, y = predRich)) + geom_point() +
  geom_abline(intercept = 0, slope = 1)
# slightly less brute QQplot
qqplot(nullResults$maxRich, nullResults$predRich)
# residual histogram
ggplot(nullResults, aes(x = rawResid)) + geom_histogram()

## LLM output automatic function for calculating dominance
# roughly a measure of predictive variable importance
# unused; verification incomplete; included only for interested readers

# both general (equally averaged) and conditional (weighted based on model size)
# non-bootstrapped (unlike in literature) and perpetuates bias in model space
# evaluates squared effects: "incremental benefits of adding non-linearity"
# for more, see https://psycnet.apa.org/doi/10.1037/1082-989X.8.2.129
dominance_contributions <- function(data,
                                    model_col = 1,
                                    metric_col = 2,
                                    higher_is_better = FALSE,
                                    sep = ";",
                                    sq_suffix = "_sq") {
  
  ###########  1. Standardise input  ########### 
  dat <- data %>%
    transmute(model = as.character(data[[model_col]]),
              metric = as.numeric(data[[metric_col]])) %>%
    mutate(model = str_trim(model),
           model = if_else(is.na(model), "", model))
  
  ###########  2. Parse model strings into sets of terms  ########### 
  model_terms <- dat %>%
    mutate(
      terms = str_split(model, fixed(sep)),
      terms = map(terms, ~ .x[.x != ""]),
      terms = map(terms, ~ sort(unique(.x))),
      key = map_chr(terms, paste, collapse = sep),
      # Number of terms in the model; model size used for conditional dominance
      model_size = map_int(terms, length)
    )
  
  # Check for duplicate specifications
  duplicated_models <- model_terms %>%
    count(key) %>%
    filter(n > 1)
  
  if (nrow(duplicated_models) > 0) {
    stop(
      "Duplicate model specifications detected. ",
      "Each unique model must occur only once."
    )
  }
  
  # Lookup table for metric values
  lookup <- model_terms %>%
    select(key, metric)
  
  ###########  3. Expand models into model/term combinations  ########### 
  expanded <- model_terms %>%
    select(key, terms, model_size) %>%
    unnest_longer(terms, values_to = "term")
  
  ########### 4. Construct nested model by removing each term  ########### 
  comparisons <- expanded %>%
    mutate(
      is_squared = str_ends(term, fixed(sq_suffix)),
      
      base_term = if_else(
        is_squared,
        str_remove(term, fixed(sq_suffix)),
        term
      ),
      
      reduced_terms = map2(
        key,
        term,
        ~ {
          current <- str_split(.x, fixed(sep))[[1]]
          current <- current[current != ""]
          sort(setdiff(current, .y))
        }
      ),
      
      reduced_key = map_chr(
        reduced_terms,
        paste,
        collapse = sep
      ),
      
      # Size of the simpler/nested model.
      reduced_model_size = model_size - 1L
    )
  
  ########### 5. Match each comparison to its nested model  ########### 
  comparisons <- comparisons %>%
    left_join(
      lookup %>%
        rename(
          reduced_key = key,
          metric_without = metric
        ),
      by = "reduced_key"
    ) %>%
    left_join(
      lookup %>%
        rename(
          key = key,
          metric_with = metric
        ),
      by = "key"
    ) %>%
    filter(
      !is.na(metric_without),
      !is.na(metric_with)
    )
  
  ###########  6. Enforce hierarchical treatment of squared terms ########### 
  #
  # For var1_sq:
  #
  #   var1 + var1_sq + A
  #             vs
  #   var1 + A
  #
  # rather than:
  #
  #   var1 + var1_sq + A
  #             vs
  #   A
  comparisons <- comparisons %>%
    filter(
      !is_squared |
        map2_lgl(
          reduced_terms,
          base_term,
          ~ .y %in% .x
        )
    )
  
  ###########  7. Calculate marginal improvement  ########### 
  #
  # For RMSE-like metrics:
  #
  #   improvement = metric_without - metric_with
  #
  # Therefore:
  #
  #   positive = improvement
  #   negative = deterioration
  comparisons <- comparisons %>%
    mutate(
      delta = if (higher_is_better) {
        metric_with - metric_without
      } else {
        metric_without - metric_with
      }
    )
  
  ###########  8. CONDITIONAL DOMINANCE  ########### 
  #
  # Conditional dominance is calculated separately for every
  # model size.
  #
  # model_size refers to the larger model:
  #
  #   A + B + C  -> model_size = 3
  #   A + B      -> reduced size = 2
  #
  # Thus the contribution is "adding the variable to a model
  # of size model_size - 1".
  conditional <- comparisons %>%
    group_by(term, base_term, is_squared, model_size) %>%
    summarise(
      n_comparisons = n(),
      conditional_mean = mean(delta, na.rm = TRUE),
      
      conditional_sd = if (n() > 1) {
        sd(delta, na.rm = TRUE)
      } else {
        NA_real_
      },
      # This is the uncertainty of the mean within the
      # available model contexts at this model size.
      conditional_se = if (n() > 1) {
        sd(delta, na.rm = TRUE) / sqrt(n())
      } else {
        NA_real_
      },
      
      .groups = "drop"
    ) %>%
    mutate(
      term_type = if_else(
        is_squared,
        "squared/additional",
        "linear"
      )
    ) %>%
    arrange(term, model_size)
  
  ###########  9. GENERAL DOMINANCE  ########### 
  #
  # Each model size receives equal weight.
  #
  # G_j = mean(C_jk)
  #
  # where C_jk is the conditional contribution at model size k.
  #
  # General uncertainty is the SD of C_jk across model sizes.
  # This is explicitly MODEL-SPACE variability, not sampling SE.
  general <- conditional %>%
    group_by(term, base_term, is_squared, term_type) %>%
    summarise(
      n_model_sizes = n(),
      general_mean = mean(conditional_mean, na.rm = TRUE),
      general_sd = if (n() > 1) {
        sd(conditional_mean, na.rm = TRUE)
      } else {
        NA_real_
      },
      total_comparisons = sum(n_comparisons),
      min_model_size = min(model_size),
      max_model_size = max(model_size),
      .groups = "drop") %>%
    arrange(desc(general_mean))
  
  ########### 10. Return all results  ########### 
  list(
    # primary, equal-weighted comparison
    general = general,
    # slightly awkward summary where each row is a different # parameters too
    # needs post-processing to be re-weighted into a single metric
    conditional = conditional,
    # more just an internal summary of the comparisons
    comparisons = comparisons
  )
}

# run to test with RMSE (but can re-run with MAE, for instance)
domTest <- dominance_contributions(data = modScores %>% 
                                     # get columns in the right order
                                     select(c(mod, validation_RMSE)) %>% 
                                     # format the null model approriately
                                     mutate(mod = mod %>% str_remove_all("Null")))

# general dominance is probably the only usable one in this state
#domTest$general