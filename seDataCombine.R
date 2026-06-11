########################### preparation ###########################
library(tidyverse)
library(readxl)

setwd(here::here())
setwd("dataTS/lake/Sweden")

# import from API downloads
seSpecImp <- read_csv("seCommSpat.csv")
seChemImp <- read_csv("seChemClean.csv")

########################### tidy chemistry data ###########################
# keep bare minimum columns
seChemImp <- seChemImp %>% 
  select(c(samplingSiteId, samplingDate, propertyCode, val)) %>% 
  mutate(samplingDateD = as.Date(samplingDate),
         samplingYr = year(samplingDateD))

# get only unique site-year combinations and
obsList <- seSpecImp %>% select(samplingSiteId, samplingYr)
# filter the chemistry data by only sites we want
seChemFilt <- left_join(obsList, seChemImp)
# almost half the sites have no chemistry at all
nrow(unique(seChemFilt %>% filter(is.na(val)) %>% 
              select(samplingSiteId, samplingYr)))
# get the empty NAs out for now
seChemFilt <- seChemFilt %>% filter(!is.na(val))

# available values
table(seChemImp$propertyCode)
# remove some variables we don't need
seChemMin <- seChemFilt %>% filter(!(propertyCode %in% c("NH4_N",
                                                         "NO2_NO3_N",
                                                         "PO4_P",
                                                         "Kjeldahl_N",
                                                         "NO2_N",
                                                         "NO3_N")))
# how many measurements per year for observations
View(seChemMin %>% 
  mutate(sampDate = as.Date(samplingDate),
         sampYr = year(sampDate)) %>% 
  summarise(.by = c("samplingSiteId", "sampYr", "propertyCode"),
            dateS = min(sampDate),
            dateE = max(sampDate),
            nSamp = n()))



seST <- seChemMin %>% mutate(monthChem = month(samplingDate),
                        inSzn = monthChem %in% 6:9)

# look out how many are out-of-season
seErr <- seST %>% summarise(.by = c("samplingSiteId", "samplingYr", "propertyCode"),
                            nSzn = sum(inSzn, na.rm = TRUE)) %>% 
  filter(nSzn == 0)



# we want data as average of measurements during summer
# a weighted average where June-Aug data is weighted 3
# and May-Sep data is weighted 2
# that way if we only have 2 observations (e.g., 1 july and 1 may) then the july data will be more important
# this is important due to the "true" growing season conditions in northern parts of sweden
seComp <- seST %>% mutate(inSummer = monthChem %in% 6:8,
                          inSample = monthChem %in% c(5,9)) %>% 
  summarise(.by = c("samplingSiteId", "samplingYr", 
                    "propertyCode"), # add some more columns later?
            val = case_when(
              sum(inSummer, na.rm = TRUE) > 0 & sum(inSample, na.rm = TRUE) > 0 ~ 
                (3 * sum(val * inSummer, na.rm = TRUE)/sum(inSummer, na.rm = TRUE) + 
                   2 * sum(val * inSample, na.rm = TRUE)/sum(inSample, na.rm = TRUE))/5,
              sum(inSummer, na.rm = TRUE) > 0 ~ sum(val * inSummer, na.rm = TRUE)/sum(inSummer, na.rm = TRUE),
              sum(inSample, na.rm = TRUE) > 0 ~ sum(val * inSample, na.rm = TRUE)/sum(inSample, na.rm = TRUE),
              .default = NA))

# convert from long to wide using the old data
seCompChem <- seComp %>% 
  pivot_wider(names_from = propertyCode,
              values_from = val)


# compare this with "conservative" approach
# consistently a few more missing values, usually 2-10, per chem variable
# might be more accurate we ignore comparing winter chemistry to summer chemistry
View(
  seCompChem %>% summarise(across(everything(), ~ sum(is.na(.))))
  )

# remove columns with too many missing values (> 53 or >25% missing)
seCompChem <- seCompChem %>% select(-c("Abs_OF420",
                                       "Siktdjup",
                                       "Turb_FNU",
                                       "O2",
                                       "SaO2",
                                       "DOC",
                                       "Abs_F254",
                                       "Abs_F365",
                                       "Abs_436m",
                                       "Siktdjup_utan_kikare",
                                       "Siktdjup_kikare",
                                       "Farg",
                                       "P_ICP", # not available for any missing TP rows
                                       "Kond", # removed because unclear at which temperature
                                       "Salinitet", # limited salinity values are all > 0.5 PSU (ppt) which makes the rows that have them brackish??
                                       "Tot_P_F" # can't compare filtered to unfiltered
                                       ))


# infill TN using multiple columns
# results ARE method-dependent from pers comm w/ Lars Sonesten
# but maybe that level of precision doesn't matter
# prefer in the order of Tot_N_TNb > Tot_N_ps > Tot_N_summa > Tot_N
seChemFin <- seCompChem %>% 
  mutate(TN = case_when(
    !is.na(Tot_N_TNb) ~ Tot_N_TNb,
    !is.na(Tot_N_ps) ~ Tot_N_ps,
    !is.na(Tot_N_summa) ~ Tot_N_summa,
    !is.na(Tot_N) ~ Tot_N,
    .default = NA
  )) %>% 
  select(-c(Tot_N_TNb, Tot_N_ps, Tot_N_summa, Tot_N))





########################### biodiversity metrics ###########################
# clean up species to get biodiversity
# we can go back and get Shannon diversity and etc later but for now:
# 1) total richness
# 2) Murphy list richness
# also impute 0s for NAs in TMI --> assume it means none of the 105 species were found

# import Murphy list from https://ars.els-cdn.com/content/image/1-s2.0-S0304377019300300-mmc2.xlsx
# citation from https://doi.org/10.1016/j.aquabot.2019.06.006
murphyListRaw <- read_excel("MacrophyteCFGLMM/1-s2.0-S0304377019300300-mmc2.xlsx", 
                            sheet = "world", range = "A2:K3460")
murphyListHelper <- murphyListRaw %>% 
  # get just species-level, not subspecies or variant
  mutate(origClean = str_squish(str_to_sentence(word(Species, 1, 2))),
         synonymClean = str_squish(str_to_sentence(word(Synonym, 1, 2))),
         # we have a few identifications only at the genus level
         origGenus = word(origClean, 1),
         synonymGenus = word(synonymClean, 1))

# get all unique values from those 4 columns to make a macrophyte master list
# there's a few genera in the synonyms going "M." or "I."; hope it's fine
macroMasterList <- unique(c(murphyListHelper$origClean,
                            murphyListHelper$synonymClean,
                            murphyListHelper$origGenus,
                            murphyListHelper$origGenus))

# get all our species names to make life easier
seMacroNames <- colnames(seSpecImp %>% select(`Acorus calamus`:`Zygnemataceae`))
# now get our unique species names in the Murphy list
seMacroMurphy <- intersect(macroMasterList, seMacroNames)

# first, get "true" richness by converting everything to 0/1
seSpecProc <- seSpecImp %>% 
  mutate(across(all_of(seMacroNames), ~ as.numeric(.x > 0))) %>% 
  rowwise() %>% 
  # then, get richnesses
  mutate(maxRich = sum(c_across(all_of(seMacroNames)), na.rm = TRUE),
         murphyRich = sum(c_across(all_of(seMacroMurphy)), na.rm = TRUE)) %>% 
  # clean up the data
  select(-all_of(seMacroNames))

# we could use the difference between start and finish coordinates as another survey effort metric
# too much work to get working in R because of multiple CRS
# best way to do would be to remove RT90 rows, convert 4 coordinate pairs with sf to SWEREF99, and then bind_rows back
# but anyways, we remove these
seSpecFin <- seSpecProc %>% 
  select(-c("siteCoordStartXMin",
            "siteCoordStartYMin",
            "siteCoordEndXMin",
            "siteCoordEndYMin",
            "siteCoordStartXMax",
            "siteCoordStartYMax",
            "siteCoordEndXMax",
            "siteCoordEndYMax",
            "siteCoordCRS"))

# our final dataset:
seFin <- full_join(seSpecFin, seChemFin,
                   by = c("samplingSiteId", "samplingYr"))

########################### descriptive plotting ###########################
library(sf)

# edit the model for descriptive plotting
seFinPlot <- seFin %>% 
  mutate(inFilt = case_when(
    is.na(Tot_P) + is.na(TN) + is.na(Temp) == 0 ~ 1,
    .default = 0
  ))

# make spatial
seGeo <- st_as_sf(seFinPlot,
                  coords = c("samplingSiteX", "samplingSiteY"), 
                  crs = 3006)

#make WGS84 for reporting
seWGS <- seGeo %>% st_transform(crs = 4326) %>% 
  mutate(lon = unlist(map(geometry,1)),
         lat = unlist(map(geometry,2)))

# get basemap as vector polygon
seMap <- geodata::gadm(country = "SWE", level = 0, path = tempdir())
# make SWEREF 99
seMap <- st_as_sf(seMap, crs = 4326) %>% st_transform(crs = 3009)



# spatially plot richness
richPlot <- ggplot() + 
  geom_sf(data = seMap) +
  geom_sf(data = seGeo, aes(colour = rich,
                            shape = factor(inFilt,
                                           levels = c(1, 0),
                                           labels = c("Yes", "No")))) +
  theme_bw() +
  scale_colour_viridis_c() +
  labs(shape = "Chemistry",
       colour = "Richness") #+
  #theme(legend.position = c(.9, .1))


# spatiotemporal plot
stPlot <- ggplot() + 
  geom_sf(data = seMap) +
  geom_sf(data = seGeo, aes(colour = samplingYr,
                            shape = factor(inFilt,
                                           levels = c(1, 0),
                                           labels = c("Yes", "No")))) +
  theme_bw() +
  scale_colour_viridis_c() +
  labs(shape = "Chemistry",
       colour = "Year") #+
  #theme(legend.position = c(.9, .1))

# TN-T chloropleth map
library(biscale)
# classify and prepare data
seBiPlot <- bi_class(seFinPlot %>% filter(inFilt == 1), 
                     x = TN, 
                     y = Temp, 
                     style = "quantile", 
                     dim = 3) %>% 
  st_as_sf(coords = c("samplingSiteX", "samplingSiteY"), 
           crs = 3006)

# make custom palette
viridisPal <- c(
  "1-1" = "#E8F4F3", # low x, low y
  "2-1" = "#C2BDD6",
  "3-1" = "#9874A1", # high x, low y
  "1-2" = "#DDF2C0",
  "2-2" = "#85C2C0", # medium x, medium y
  "3-2" = "#6380A6",
  "1-3" = "#FEF287", # low x, high y
  "2-3" = "#72CF8E",
  "3-3" = "#21908D" # high x, high y
)

# make base bivariate plot
biPlot <- ggplot() + 
  geom_sf(data = seMap) +
  geom_sf(data = seBiPlot, aes(fill = bi_class), 
          colour = "black", 
          size = 2,
          pch = 21, 
          show.legend = FALSE
          ) +
  bi_scale_fill(pal = viridisPal, dim = 3) +
  theme_bw()
  
# make legend separately
biLeg <- bi_legend(pal = viridisPal,
                    dim = 3,
                    xlab = "High TN",
                    ylab = "High T",
                    size = 8) + 
  theme(plot.margin=unit(c(-0.05,0,0,0), "null")) # remove margin around plot

# combine plot together
#bivFinPlot <- cowplot::ggdraw(biPlot) #+
  #cowplot::draw_plot(biLeg, 0.5, .2, 0.2, 0.2)

# organise into final figure
sumGrid <- cowplot::plot_grid(richPlot, 
                              stPlot, 
                              biPlot,
                              labels = c("a", "b", "c"),
                              align = "h",
                              axis = "bt",
                              #rel_heights = c(1.2, 1, 1),
                              nrow = 1) +
  #cowplot::draw_plot(biLeg, 0.805, .0538, 0.18, 0.18) # old version w/o margin changes
  cowplot::draw_plot(biLeg, 0.855, .0545, 0.1, 0.12)

# save plot
cowplot::save_plot("MacrophyteCFGLMM/fig4.jpeg", sumGrid, 
                   nrow = 1,
                   base_height = 5,
                   base_asp = 2.2,
                   #base_width = 2,
                   dpi = 300)







########################### extra reporting stuff ###########################
# mostly just for reporting about programme names for the manuscript
# docs at https://www.slu.se/om-slu/organisation/institutioner/vatten-miljo/datavardskap/registersida/
table(seFinPlot$inFilt, seFinPlot$subProgramName)
table(seFinPlot$subProgramName, seFinPlot$county)