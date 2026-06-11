########################### preparation ###########################
library(tidyverse)
library(spCF)
library(stringr)
#library(fastDummies)
library(ggcorrplot)
library(viridis)
library(sf)
library(ggpattern)

setwd(here::here())
setwd("dataTS/lake/Sweden/MacrophyteCFGLMM")

########################### correlation plot ###########################
corrPlot <- ggcorrplot(cor(seFin %>% 
                 select(c("samplingYr",
                          "samplingM",
                          "samplingSiteX",
                          "samplingSiteY",
                          "rich",
                          "tmiRich",
                          "tmi",
                          "murphyRich",
                          "maxRich",
                          "observationCount")) %>% 
                   rename(Year = "samplingYr",
                          Month = "samplingM",
                          X = "samplingSiteX",
                          Y = "samplingSiteY",
                          Richness = "rich",
                          `TMI Richness` = "tmiRich",
                          TMI = "tmi",
                          `Murphy Richness` = "murphyRich",
                          `Full Richness` = "maxRich",
                          Count = "observationCount")),
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
                          "rich",
                          "tmiRich",
                          "tmi",
                          "murphyRich",
                          "maxRich",
                          "observationCount",
                          "Alk",
                          "Kfyll",
                          "Kond_25",
                          "TOC",
                          "Tot_P",
                          "Temp",
                          "pH",
                          "TN")) %>% 
                   relocate(c(samplingM, samplingYr,
                              samplingSiteX, samplingSiteY,
                              rich, murphyRich, maxRich, tmiRich, tmi,
                              observationCount,
                              Alk, Kfyll, Kond_25, TOC, pH,
                              Tot_P, Temp, TN)) %>% 
                          rename(Year = "samplingYr",
                                 Month = "samplingM",
                                 X = "samplingSiteX",
                                 Y = "samplingSiteY",
                                 `WFD Richness` = "rich",
                                 `TMI Richness` = "tmiRich",
                                 TMI = "tmi",
                                 `Murphy Richness` = "murphyRich",
                                 `Full Richness` = "maxRich",
                                 `# Obs` = "observationCount",
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
y <- seFin %>% pull(rich)
offset <- seFin %>% pull(observationCount)

# fit the basic, full model with holdout validation
mod_hv <- cf_glm_hv(y = y, #x = x, 
                    #offset = offset, # looks like model doesn't find suitable coefficients with offset
                    coords = coords, 
                    family=poisson()) # poisson and quasipoisson are identical

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
  mutate(rawResid = rich - predRich)

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
cowplot::save_plot("fig6.jpeg", resids, 
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

ggsave("fig5.jpeg", 
       spatDecompPlot,
       width = 24,
       height = 17,
       units = "cm",
       #scale = 1,
       dpi = 300)


########################### filter w/ chemistry data ###########################

# re-run with just water sites with water chemistry data
# method is sensitive to NAs which must be removed

# summarise all NA values again
View(
  seFin %>% ungroup() %>% 
    summarise(across(everything(), ~ sum(is.na(.))))
)
# filter only for values not NA in a few values
seFinFilt <- seFin %>% drop_na(any_of(c(#"secchi", # missing values here are different from the chem
  #"Abs_F420", # drops R2 to 74%
  #"Alk",
  #"Kfyll",
  #"Kond_25",
  #"TOC", # drops R2 to 63%
  "Tot_P",
  "Temp",
  #"pH",
  "TN"))) %>% 
  select(-c("secchi",
            "Abs_F420",
            "Alk",
            "Kfyll",
            "Kond_25",
            "TOC",
            #"Tot_P",
            #"Temp",
            "pH",
            #"TN"
  )) %>% 
  mutate(TPsq = Tot_P^2,
         Tsq = Temp^2,
         TNsq = TN^2)

# prepare data for package
coordFilt <- seFinFilt %>% select(c(samplingSiteX, samplingSiteY)) %>% as.matrix()
yFilt <- seFinFilt %>% pull(rich)
# for later test
#offsetFilt <- seFinFilt %>% pull(observationCount)


nullFiltHV <- cf_glm_hv(y = yFilt,
                        coords = coordFilt, 
                        family = poisson()) # poisson and quasipoisson are identical
nullFilt <- cf_glm(y = yFilt,
                   coords = coordFilt, 
                   mod_hv = nullFiltHV)
# model when filtering only for all rows with some chem value is horrible
# best R2 at 85,6% when we have TN, Temp, and TP
nullFilt

# spatial trends are in
nullFilt$bands
# covariates are in nullFilt$beta
nullFilt$beta
# summary statistics in
nullFilt$e_summary


# basic residuals check
nullFiltResults <- seFinFilt %>% 
  bind_cols(# Predictive mean
    nullFilt$pred$pred, 
    # Predictive SD
    nullFilt$pred$pred_sd) %>% 
  rename(predRich = `...35`,
         predSd = `...36`) %>% 
  mutate(rawResid = rich - predRich)


# residuals vs fitted
# clearly some misfit with more extreme underpredictions as predicted richness increases
ggplot(nullFiltResults, aes(x = predRich, y = rawResid)) + geom_point()
# absolutely lethal misfit for total and just TN data
# true richness vs residual
ggplot(nullFiltResults, aes(x = rich, y = rawResid)) + geom_point()

# temporal sensitivity
# missing a few years later
ggplot(nullFiltResults, aes(x = samplingYr, y = rawResid)) + geom_point()

# spatial sensitivity
nullFiltSpatial <- st_as_sf(nullFiltResults, 
                            coords = c("samplingSiteX", "samplingSiteY"), 
                            crs = 3006)
# all majorly concerning residuals situated in eastern Svealand, southern coast, or extreme northwest
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
  yFilt <- seFinFilt %>% pull(rich)
  xFilt <- seFinFilt %>% select(all_of(covList)) %>% as.matrix()
  
  modHV <- cf_glm_hv(y = yFilt, x = xFilt, 
                     coords = coordFilt, 
                     family = poisson()) # poisson and quasipoisson are identical
  modFit <- cf_glm(y = yFilt, x = xFilt, 
                   coords = coordFilt, 
                   mod_hv = modHV)
  return(modFit)
}

# make heatmap for visualisation
# empty df with each column as some bandwidth
modSums <- data.frame(band = nullFiltHV$other$bands_all)
# add null model
modSums <- modSums %>% 
  mutate(Null = case_when(
    band %in% nullFilt$band ~ nullFilt$e_summary$value[1],
    .default = 0
  ))

# more summary statistics
# all relevant to the validation set
# just get it directly from the summary object
modScores <- as.data.frame(nullFilt$e_summary) %>%
  pivot_wider(names_from = stat,
              values_from = value) %>% 
  mutate(mod = "Null")

#full list of covariates
covListFull <- c("TN",
                 "Temp",
                 "Tot_P",
                 "TNsq",
                 "Tsq",
                 "TPsq")

modFE = data.frame(#var = covListFull,
  coef = numeric(),
  isSig = numeric(),
  mod = character()
)

# fully enumerate all candidate wq models
cov1 <- c("TN")
cov2 <- c("Temp")
cov3 <- c("Tot_P")
cov4 <- c("TN", "TNsq")
cov5 <- c("Temp", "Tsq")
cov6 <- c("Tot_P", "TPsq")
cov7 <- c("TN", "Temp")
cov8 <- c("TN", "TNsq", "Temp")
cov9 <- c("TN", "Temp", "Tsq")
cov10 <- c("TN", "TNsq", "Temp", "Tsq")

cov11 <- c("Tot_P", "Temp")
cov12 <- c("Tot_P", "TPsq", "Temp")
cov13 <- c("Tot_P", "Temp", "Tsq")
cov14 <-  c("Tot_P", "TPsq", "Temp", "Tsq")

cov15 <- c("TN", "Tot_P")
cov16 <- c("TN", "TNsq", "Tot_P")
cov17 <- c("TN", "Tot_P", "TPsq")
cov18 <- c("TN", "TNsq", "Tot_P", "TPsq")

cov19 <- c("TN", "Temp", "Tot_P")
cov20 <- c("TN", "TNsq", "Temp", "Tot_P")
cov21 <- c("TN", "Temp", "Tsq", "Tot_P") 
cov22 <- c("TN", "Temp", "Tot_P", "TPsq") 

cov23 <- c("TN", "TNsq", "Temp", "Tsq", "Tot_P")
cov24 <- c("TN", "TNsq", "Temp", "Tot_P", "TPsq")
cov25 <- c("TN", "Temp", "Tsq", "Tot_P", "TPsq")
cov26 <- c("TN", "TNsq", "Temp", "Tsq", "Tot_P", "TPsq")

# loop through all our candidate sets
for(i in 1:26){
  # get the set of covariates we're working with now
  # eval to parse for function use 
  covVect <- eval(parse(text = paste0("cov",i)))
  # silently run its model
  invisible({capture.output({
    modTemp <- runQuickCF(covVect)
  })})
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
    feTemp <- modTemp$beta %>% rownames_to_column() %>% rename(var = "rowname") %>% 
      filter(var == "x2") %>% 
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

# get the minimum distance for all columns
minDist <- modSums %>%
  filter(if_any(-band, ~ . != 0)) %>%
  summarise(minDist = min(band)) %>%
  pull(minDist)

# set clean model order for readability
modelOrderClean <- str_replace_all(
  str_replace_all(colnames(modSums)[-1],
                  "Temp","T"),
  "Tot_P","TP")

# format for plot
modOrderLab <- str_replace_all(modelOrderClean, "sq", "(^2") %>% 
  str_remove_all("\\;[:alpha:]{1,2}\\(") %>%
  str_replace_all("\\;", " + ")
  # probably a better way to do this, but this is fine
  #str_replace_all("\\^", "\\(^")

modPlot <- modSums %>% 
  # remove all bandwidths below the minimum detection for all models
  filter(band >= minDist) %>%
  # convert from m to km
  mutate(band = band/1000) %>% 
  pivot_longer(cols = colnames(modSums)[-1],
               names_to = "modName",
               values_to = "R2") %>% 
  # make names more consistent & readable
  mutate(modName = str_replace_all(modName,
                                   "Tot_P",
                                   "TP"),
         # do this twice bc combining too much work
         modName = str_replace_all(modName,
                                   "Temp",
                                   "T")) %>% 
  # format factor levels manually to make plot more readable
  mutate(distance = as.factor(format(round(band, 1), nsmall = 1)),
         modName = factor(modName,
                          levels = modelOrderClean,
                          labels = modOrderLab))

# plot
spatPlot <- ggplot(modPlot, aes(modName, distance, fill = R2)) +
  geom_tile() +
  theme_bw() + 
  labs(x = "",
       y = "Bandwitdh (km)",
       fill = bquote(R^2)) +
  # make limits & labels more readable
  scale_fill_viridis_c(limits = c(0, round(max(modPlot$R2),2)),
                       breaks = seq(0, round(max(modPlot$R2),2), 
                                    length.out = 4)) +
  scale_x_discrete(labels = scales::parse_format()) +
  # make x axis readable
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# add another heatmap for FE coefficients showing
# rescale by the fixed effect mean to get results on the same scale
# first get the means from each covariate
meanFE <- as.data.frame(t(seFinFilt %>% select(c("TPsq", "Tot_P", "Tsq", "Temp", "TNsq", "TN")) %>% 
                            summarise(across(everything(), ~ mean(., na.rm = TRUE))))) %>% rename(meanEff = "V1") %>% 
  rownames_to_column(var = "var")
# combine with the original data
modReFE <- left_join(modFE, meanFE, by = "var") %>% 
  mutate(coefResc = exp(coef * meanEff))

# modReFE the factors to our desited level
modReFE <- modReFE %>% 
  # make names more consistent & readable
  mutate(mod = str_replace_all(mod, "Tot_P", "TP"),
         # do this twice bc combining too much work
         mod = str_replace_all(mod, "Temp", "T"),
         # and again for the variables
         var = var %>% replace_when(
           var == "Tot_P" ~ "TP",
           var == "Temp" ~ "T"
         ),
         isNeg = as.numeric(coef < 0)
  )

modReFE <- bind_rows(modReFE,
                     data.frame(
                       coef = NA,
                       isSig = 0,
                       mod = "Null",
                       var = "TN",
                       meanEff = NA,
                       coefResc = NA,
                       isNeg = NA
                     )) %>% 
  mutate(# set factor levels to improve readability
    mod = factor(mod, levels = modelOrderClean),
    var = factor(var,
                 levels = c("TPsq", "TP", "Tsq", "T", "TNsq", "TN"),
                 labels = c("TP^2", "TP", "T^2", "T", "TN^2", "TN")))

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
  theme_bw() +
  labs(x = "",
       y = "") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
  ) 

modSum <- cowplot::plot_grid(fePlot, spatPlot, 
                   labels = c("a", "b"),
                   ncol = 1,
                   rel_heights = c(0.4, 1.4),
                   align = "v") 
cowplot::save_plot("fig1.jpeg", modSum, 
          nrow = 2,
          dpi = 300,
          base_height = 5)

# possible extensions:
# explore the survey effort fixed effects
#       -our goal is to find fixed effects that give us similar bandwidths before and after site filtering




















########################### old stuff (not used) ###########################

# spCF does not provide a dispersion parameter for quasipoisson
# we can either do tests with base glm() function, or
# in basis-type additive models (e.g., GAMs) dispersion is best estimated from residuals
# Pearson's goodness-of-fit for Poisson
pearsonGOF <- sum((nullResults$rawResid)^2 / nullResults$predRich)
# dispersion by dividing by degrees of freedom --> not available
# can get DHARMa to work if there's a way to get simulated residuals

### unused diagnostics
# raw vs fitted; quasi-QQ plot
ggplot(nullResults, aes(x = rich, y = predRich)) + geom_point() +
  geom_abline(intercept = 0, slope = 1)
# slightly less brute QQplot
qqplot(nullResults$rich, nullResults$predRich)
# residual histogram
ggplot(nullResults, aes(x = rawResid)) + geom_histogram()