########################### preparation ###########################
library(tidyverse)
library(httr)
library(jsonlite)
library(stringr)

# set the working directory to the main R project directory
setwd(here::here())

################################ Chem API Call ################################
# docs https://miljodata.slu.se/api/docs/index.html

# replace the string with the API key
apiKey <- "api_key_here"

# first, get a list of all stations/sites with macrophyte observations
# build full API call data URL
urlMacroMeta <- paste0("https://miljodata.slu.se/api/observations-service/v2/stations/query?token=",
                       apiKey,
                       "&product=Makrofyter")
# Send http request
repMacroMeta <- httr::GET(urlMacroMeta)
# Read binary file
jsonMacroMeta <- readBin(repMacroMeta$content, "text")
# Parse JSON file
listMacroMeta <- jsonlite::fromJSON(jsonMacroMeta)
# process data to get list of names
dbMacroMeta <- listMacroMeta$stations #%>% filter(!is.na())

# now, get chemistry data only for stations with macrophyte observations
# URL has limits in the number of stations, so easier to do 1-by-1 via loop

# create placeholder dataframe
# set data types in advance to avoid type mismatch issues
fullChemRaw <- data.frame(sampleId = integer(),
                          sourceSampleId = character(),
                          isSensitive = logical(),
                          productCode = character(),
                          productName = character(),
                          productType = character(),
                          deliveryMethodType = character(),
                          deliveryTitle = character(),
                          samplingDate = character(),
                          sampleStatus = character(),
                          samplingMethod = character(),
                          samplingMethodStandard = logical(),
                          samplingMethodType = character(),
                          mediumTypeName = character(),
                          comment = character(),
                          sampleAccredited = character(),
                          analysisLabs = character(),
                          analysisMethods = character(),
                          flagCodes = logical(),
                          nationalStationId = character(),
                          nationalSiteId = character(),
                          samplingSiteId = integer(),
                          samplingSiteName = character(),
                          samplingSiteX = character(),
                          samplingSiteY = character(),
                          samplingSiteWkt = character(),
                          samplingSiteCoordinateSystem = character(),
                          relativeLocationType = character(),
                          minDepth = numeric(),
                          maxDepth = numeric(),
                          stationId = integer(),
                          stationName = character(),
                          stationType = character(),
                          stationEUID = character(),
                          stationCoordinateX = integer(),
                          stationCoordinateY = integer(),
                          stationCoordinateSystem = character(),
                          stationIsBlurred = logical(),
                          stationMetadata = list(),
                          county = character(),
                          municipality = character(),
                          studyName = character(),
                          subProgramName = character(),
                          surveyType = character(),
                          waterZone = logical(),
                          observationCount = integer(),
                          observations = list(),
                          groupCalculations = list(),
                          sampleIndices = logical(),
                          sampleMetadata = list(),
                          insertDate = character(),
                          editDate = character(),
                          lastEdit = character())

# pre-define the rest of the URL to call, station ID goes in the middle
# automatically add in the API key as above
urlPreamble <- paste0("https://miljodata.slu.se/api/observations-service/v2/full-samples/query?token=",
                      apiKey,
                      "&stationIds=")
urlPostamble <- "&productType=Chemistry&fromYear=2005&toYear2025"

# for each station with macrophytes
for(i in dbMacroMeta$stationId){ # use stationID bc nationalStationID has NAs
  # API call & basic conversion as above
  urlTemp <- paste0(urlPreamble, i, urlPostamble)
  responseTemp <- httr::GET(urlTemp)
  jsonTemp <- readBin(responseTemp$content, "text")
  dtaTemp <- jsonlite::fromJSON(jsonTemp)
  # get data out of the API output in dataframe format
  outputTemp <- dtaTemp$samples
  # append to placeholder db and repeat
  fullChemRaw <- bind_rows(fullChemRaw, outputTemp)
}

# save to avoid re-downloading later
write_rds(fullChemRaw, "seChemRawAPI.rds")

################################ Chem Cleaning ################################
# read-in again for restarting the script here
fullChemRaw <- read_rds("seChemRawAPI.rds")
View(fullChemRaw)

# isSensitive, stationIsBlurred is purely false
# samplingMethodStandard, flagCodes, waterZone, sampleIndices is purely NA
fullChem <- fullChemRaw %>% 
  select(-c(isSensitive, stationIsBlurred, samplingMethodStandard, flagCodes, 
            waterZone, sampleIndices)) %>% 
  rowwise() %>% 
  # what about these lists?
  mutate(statEle = length(stationMetadata),
         obsEle = length(observations),
         # no elements in any row
         groupEle = length(groupCalculations),
         sampleEle = length(sampleMetadata))

# messy columns are stationMetadata, observations, groupCalculations, sampleMetadata
# remove groupCalculations bc always empty
# station metadata is an enigma fullChem[[44]][[57]]
# sample metadata has meaning, but mostly duplicate/unneeded info fullChem[[47]][[8659]]
# observations is the real one; empty in 69 rows
# mostly due to bad ice/remoteness
View(fullChem %>% filter(obsEle == 0))
fullChem <- fullChem %>% select(-c(groupCalculations, groupEle))

fullChemTest <- unnest(fullChem, col = observations, names_repair = "unique")

# sample IDs are the same
all.equal(fullChemTest$sampleId...1,fullChemTest$sampleId...48)
# remove all columns with only one value
# SWEREF99 TM is the CRS (both site and station)
fullChemTest <- fullChemTest %>% select(-c(sampleId...48, 
                                           productCode,
                                           productName,
                                           productType,
                                           sampleStatus,
                                           samplingSiteCoordinateSystem,
                                           stationCoordinateSystem,
                                           propertyRepresentationType,
                                           accessRestriction,
                                           isQualitative,
                                           observationMetadata,
                                           obsEle)) %>% 
  rename(sampleId = sampleId...1) 

fullChemTest2 <- unnest(fullChemTest, col = observationValues, names_repair = "unique")


# observation IDs are the same
all.equal(fullChemTest2$observationId...41,fullChemTest2$observationId...50)
all.equal(fullChemTest2$name,fullChemTest2$propertyName)
# remove all columns with only one value or fully NA
# SWEREF99 TM is the CRS
fullChemTest2 <- fullChemTest2 %>% select(-c(observationId...50, 
                                           rank,
                                           valueType,
                                           isFromSubSample, 
                                           group,
                                           subSampleCd,
                                           name)) %>% 
  rename(observationId = observationId...41)

# get rid of the columns I can't understand for now
fullChemW <- fullChemTest2 %>% select(-c(stationMetadata,
                                         sampleMetadata,
                                         statEle, 
                                         sampleEle))

# 79 cases of nonstandard units, usually µg and mg doublets
wqDict <- unique(fullChemTest2 %>% select(propertyCode, propertyAbbrevName, propertyName, unit))
nrow(wqDict) - nrow(unique(fullChemTest2 %>% select(propertyCode, propertyAbbrevName, propertyName)))

wqDups <- wqDict %>% group_by(propertyCode) %>% filter(row_number() > 1)
wqDuplicates <- wqDict %>% filter(propertyCode %in% wqDups$propertyCode) %>% arrange(propertyCode)
#write_csv(wqDuplicates, "SwedenChemistryDuplicates.csv")




# harmonise units 
# https://miljodata.slu.se/mvm/DataContents/UnitConversion
# download dictionary as csv from ^
# make manual changes to give an acutal dictionary with no duplicated units
unitDict <- read_csv("Miljödata MVM - Administrera kontrollhalter.csv")
# match unit names to make life easier
colnames(unitDict) <- c("id", "propertyAbbrevName", "unitC", "unitAlt", 
                        "convFactRaw", "unitC_alt", "unitAlt_alt", "sigfig")
# replace commas with periods for numeric values
unitDict <- unitDict %>% mutate(convFact = as.numeric(str_replace_all(convFactRaw,
                                       ",",
                                       ".")))

# too many columns will be in the way
# how many distinct values are there: which rows can we throw out?
t(fullChemW %>% summarise(across(everything(), ~ n_distinct(.))))
t(fullChemW %>% summarise(across(everything(), ~ sum(is.na(.)))))

# doing this manually is impossible
# how many duplicates do we have if grouping by just ID and date
chemGroupT <- fullChemW %>% group_by(samplingSiteId, samplingDate) %>% 
  summarise(across(everything(), ~ n_distinct(.)))
View(chemGroupT %>% ungroup() %>% 
       summarise(across(everything(), ~ max(.))))

# keep only the bare minimum for visibility
chemCompact <- fullChemW %>% 
  select(-c("sourceSampleId", "nationalStationId", "nationalSiteId", "stationEUID", "observationId", # who needs another ID
            "deliveryMethodType", # good for uncertainty but we're not that picky for chemistry
            "deliveryTitle", "samplingMethod", "samplingMethodType", 
            "sampleAccredited", "analysisLabs", "analysisMethods",
            "studyName", "subProgramName", "surveyType", "observationCount",
            "insertDate", "editDate", #"lastEdit", 
            "isAccredited", "qualityCode", "analysisLab", "analysisMethod", "analysisDate", # ^^
            "mediumTypeName", "stationType", # might matter for spatial; all resurveys are lakes
            "comment", # too much work
            "samplingSiteWkt", # i think this is just geometry object; recalculate later
            "relativeLocationType", # no idea how to interpret so ignore
            "stationCoordinateX", "stationCoordinateY", # we're using sites
            "stationCoordinateN", "stationCoordinateE"))

# get only the variables we're interested in
inclVars <- c("Vattentemperatur", "Syrgashalt", "Totalfosfor",
              "Nitrit+Nitratkväve", "Fosfatfosfor",
              "Ammoniumkväve",
              "Konduktivitet 25°", "pH",
              "Syrgasmättnad",
              "Absorbans filtrerat (420 nm)",
              "Tot. org. kol TOC",
              "Turbiditet FNU",
              "Tot-N(TNb)",
              "Klorofyll a",
              "Alkalinitet/Acid.",
              "Tot-N(persulfat)",
              "Absorbans filtrerat (436 nm)",
              "Absorbans filtrerat (254 nm)",
              "Absorbans filtrerat (365 nm)",
              "Alkalinitet",
              "Absorbans ofiltrerat (420 nm)",
              "Siktdjup",
              "Siktdjup utan vattenkikare",
              "Salinitet",
              "Siktdjup med vattenkikare",
              "Tot-N",
              "Löst organiskt kol",
              "Färgtal",
              "Turbiditet NTU",
              "Kjeldahlkväve",
              "Nitratkväve", "Nitritkväve", "Kond",
              "Totalfosfor (ICP)",
              "Konduktivitet 20°",
              "Filtrerad totalfosfor",
              "Tot-N(summa)",
              "Totalfosfor filtrerad (ICP)",
              "Turb")

# filter the rest out
chemCompactF <- chemCompact %>% filter(propertyName %in% inclVars)

# replace commas with periods in value column
# no ranges to worry about
# anything with > and < thrown out for now
# table(str_replace_all(chemCompactF$value, "[:digit:]", ""))
chemCompactF <- chemCompactF %>% mutate(valueNew = as.numeric(str_replace_all(value,
                                                                              ",",
                                                                              ".")))
# when does our cleaning fail? usually with "<" and ">" signs
chemValsExp <- chemCompactF %>% filter(is.na(valueNew) & !is.na(value))
#View(chemValsExp)

# fix missing values
chemCompactF <- chemCompactF %>% rowwise() %>% 
  # no missing raw values (as far as we know)
  mutate(valueNew = ifelse(!is.na(valueNew),
                           valueNew,
                           # first, deal with below minimum detection limit
                           ifelse(str_detect(value, "<"),
                                  # get only the number
                                  as.numeric(str_replace_all(str_extract_all(value,
                                                                  "\\(?[0-9,.]+\\)?")[[1]],
                                                  ",",
                                                  # set to half the value as an assumption
                                                  "."))*0.5,
                                  as.numeric(str_replace_all(str_extract_all(value,
                                                                             "\\(?[0-9,.]+\\)?")[[1]],
                                                             ",",
                                                             # set to slightly more than the value as an assumption
                                                             "."))*1.1)))

# cleaning exposed missing values; purge them
chemCompactF <- chemCompactF %>% filter(!is.na(valueNew)) %>%
  # for date manipulation
  mutate(editDate = as.Date(lastEdit))

# check if we're going to have any duplicate problems
# we do, but way better than going by-station
View(chemCompactF %>% group_by(stationId, samplingDate, propertyAbbrevName, 
                               minDepth, maxDepth) %>% 
       summarise(n = n(),
                 nUnits = n_distinct(unit)) %>% filter(n > 1))

# so, let's harmonise units first
# we picked relatively easy variables
#View(unitDict %>% filter(propertyAbbrevName %in% unique(chemCompactF$propertyAbbrevName)))

# join with dictionary
fullChemUnits <- left_join(chemCompactF %>% rename(unitC = "unit"), 
                           unitDict %>% select(-c("unitC_alt", "unitAlt_alt",
                                                  "sigfig", "convFactRaw")), 
                           by = c("propertyAbbrevName", "unitC"))



# harmonise units
fullChemUnits2 <- fullChemUnits %>% 
  mutate(unitF = unitAlt,
         valF = valueNew * convFact) %>%
  # clean up unneeded columns
  # keep value for Alk/Acid cleaning
  select(-c("valueNew", "unitC", "unitAlt", "convFact", "id"))

# before averaging, fix the Alk_Acid column
# values with a "-" are pH
# values range from 0 to 0.1 which feels unlikely in freshwaters
# so we omit them for now (only 746 values)
chemFix <- fullChemUnits2 %>% 
  filter(propertyAbbrevName == "Alk/Acid" & str_detect(value, "-")) %>% 
  mutate(propertyCode = "pH",
         propertyAbbrevName = "pH",
         propertyName = "pH",
         unitF = NA)

# but any column with pH tested has 0 alkalinity (from pers comm)
fullChemUnits2 <- fullChemUnits2 %>%
  mutate(valF = replace_when(valF, 
                             propertyAbbrevName == "Alk/Acid" & str_detect(value, "-") ~ 0),
         # fix the names so that depth merging goes more smoothly
         propertyCode = replace_when(propertyCode, 
                                     propertyAbbrevName == "Alk/Acid" ~ "Alk"),
         propertyAbbrevName = replace_when(propertyAbbrevName, 
                                           propertyAbbrevName == "Alk/Acid" ~ "Alk."),
         propertyName = replace_when(propertyName, 
                                     propertyAbbrevName == "Alk/Acid" ~ "Alkalinitet"))

# now average values at the same depths
chemDepth <- fullChemUnits2 %>% 
  group_by(samplingSiteId, samplingDate, samplingSiteName,
           samplingSiteX, samplingSiteY, minDepth, maxDepth,
           stationId, stationName, county, municipality,
           propertyCode, propertyAbbrevName, propertyName,
           unitF) %>% 
  summarise(lastEditDate = max(editDate),
            val = mean(valF))

# what happens when the depths are unequal?
# lots of 0s involved but not only 0s
View(chemDepth %>% filter(minDepth != maxDepth))

# duplicate values indicating depth gradients
chemDups <- chemDepth %>% 
  group_by(samplingSiteId, samplingDate, propertyAbbrevName) %>% 
  summarise(n = n()) %>% filter(n > 1)
# half of our rows have depth duplicates
(sum(chemDups$n) - nrow(chemDups))/nrow(chemDepth)

# choose only the shallowest measure for a sampling occasion
# this matches other Nordic methodology
# use the shallowest maxDepth as tiebreaker
chemShallow <- chemDepth %>% 
  # throw out the remaining variables
  # they should be retained in the species data anyways
  group_by(samplingSiteId, samplingDate, propertyAbbrevName) %>%
  slice_min(tibble(minDepth, maxDepth), n = 1) 

# leaving in long form will make fusing with macrophyte data easier
write_csv(chemShallow, "seChemClean.csv")