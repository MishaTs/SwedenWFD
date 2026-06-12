########################### preparation ###########################
library(tidyverse)
library(httr)
library(jsonlite)
library(stringr)

# set the working directory to the main R project directory
setwd(here::here())

########################### API download ###########################
# replace the string with the API key
apiKey <- "api_key_here"

# build full macrophyte data URL
urlMacro <- paste0("https://miljodata.slu.se/api/observations-service/v2/full-samples/query?token=",
                   apiKey,
                   "&product=Makrofyter")
# Send http request
repMacro <- httr::GET(urlMacro)
# Read binary file
jsonMacro <- readBin(repMacro$content, "text")
# Parse JSON file
listMacro <- jsonlite::fromJSON(jsonMacro)
# save just the samples
# like the last API call, this has nested dfs and needs some loops to make usable
macroRaw <- listMacro$samples

# save temporarily
#write_rds(macroRaw, "seMacroRawAPI.rds")

########################### unpack and trim down data ##########################
#macroRaw <- read_rds("seMacroRawAPI.rds")

# remove fully NA columns
macroRaw <- macroRaw %>% select(where(~!all(is.na(.))))
# remove columns with only one unique value
# save the leftover columns in case they matter later
# for now, all that matters is that they're lakes (Sjö) and stations are unblurred
macroLeftover <- macroRaw %>% select(where(~n_distinct(.) <= 1))
macroRaw <- macroRaw %>% select(where(~n_distinct(.) > 1))

# transpose the long data to wide in stationMetadata
for(i in 1:nrow(macroRaw)){
  metaRow <- macroRaw$stationMetadata[[i]]
  if(dim(metaRow)[2] != 0){
    macroRaw$stationMetadata[[i]] <- metaRow %>% 
      pivot_wider(names_from = label, values_from = value)
  }
}

# stationMetadata corresponds to MS_CD columns; extract it
macroRaw2 <- unnest(macroRaw, col = stationMetadata, 
                    keep_empty = TRUE,
                    names_repair = "unique")

# remove unneeded columns
macroRaw2 <- macroRaw2 %>% select(-c(# only 1 non-NA min/max depth value (both 0.5)
                                     "minDepth", 
                                     "maxDepth",
                                     # just sf format of the existing coordinates
                                     "samplingSiteWkt",
                                     # coordinates either SWEREF99 TM or NA
                                     "samplingSiteCoordinateSystem",
                                     # we don't care about the name of the original Excel file
                                     "deliveryTitle",
                                     # duplicates of Y and X coordinates (respectively)
                                     "stationCoordinateN",
                                     "stationCoordinateE",
                                     # CD codes are gibberish
                                     "MS_CD C2",
                                     "MS_CD C3",
                                     "MS_CD C4"))

# do it all again for sampleMetadata
# transpose
for(i in 1:nrow(macroRaw2)){
  metaRow <- macroRaw2$sampleMetadata[[i]]
  if(dim(metaRow)[2] != 0){
    macroRaw2$sampleMetadata[[i]] <- metaRow %>% 
      pivot_wider(names_from = label, values_from = value)
  }
}
# extract
macroRaw3 <- unnest(macroRaw2, col = sampleMetadata, 
                    keep_empty = TRUE,
                    names_repair = "unique")
# fix names to make understandable
macroRaw3 <- macroRaw3 %>%
  rename(waterLevel = "Vattennivå",
         stdSamplingMethod = "Provtagningsmetodens standard (metodstandard för biologisk provtagning)",
         lakeDrainage = "Sjösänkning",
         hydromorphImpact = "Hydromorfologisk påverkan",
         lakeWholeOrPart = "Helsjö/Del av sjö",
         inventoryTransect = "Inventerade transekter",
         # technically counts for individual rake drags too??
         inventoryQuadrat = "Inventerade krattdrag/rutor",
         secchi = "Siktdjup (m)") %>% 
  # worse version of surveyType column
  select(-stdSamplingMethod)


# now the hard part begins; sampleIndices
# divided into two parts, indexQualityFlag and indexList
# first split them up
macroRaw4 <- unnest(macroRaw3, col = sampleIndices, 
                    keep_empty = TRUE,
                    names_repair = "unique") %>%
  # exactly the same as flagCodes column
  select(-indexQualityFlag)
# getting out the index list it makes sense to extract and then pivot_longer
macroRaw5 <- unnest(macroRaw4, col = indexList, 
                    keep_empty = TRUE,
                    names_repair = "unique")
# 3 indices: Macrophyte number of species, Macrophytes TMI-number of species, and TMI (Trophic Macrophyte Index)
# so, richness, TMI richness, and TMI
unique(macroRaw5$indexShortName)
macroRaw5 <- macroRaw5 %>% select(-c("indexName",
                                     # either v1 or NA; ignore
                                     "indexVersion",
                                     # always NA
                                     "unitOfMeasurement",
                                     # just another useless ID
                                     "indexValueId"))


# first, get inputValuesText as a value in IndexValue (name to be "detailTMI")
# then, pivot_wider with names_from indexShortName and values from IndexValue
# hacky way to do it but it does work
macroRaw6 <- macroRaw5 %>% pivot_longer(
  cols = c("indexValue", "inputValuesText"),
  values_to = "indVal") %>% 
  mutate(indexShortName = ifelse((!is.na(indVal) & name == "inputValuesText"), 
                                 "compTMI",
                                 indexShortName))

macroRaw7 <- macroRaw6 %>% select(-name) %>% pivot_wider(
  names_from = indexShortName,
  values_from = indVal,
  values_fn = ~ na_if(paste0(na.omit(.), collapse = '|'), ''))

# sanity check; if any of these aren't 0 then the approach needs to be refined
check1 <- sum(str_detect("|", na.omit(macroRaw7$`Mfr arter`)))
check2 <- sum(str_detect("|", na.omit(macroRaw7$`TMI arter`)))
check3 <- sum(str_detect("|", na.omit(macroRaw7$TMI)))
check4 <- sum(str_detect("|", na.omit(macroRaw7$compTMI)))
ifelse(check1 + check2 + check3 + check4 == 0, "SAFE", "TRY AGAIN")

# inputValuesText is two values like "Summa produkt=55.6;Summa viktfaktor=7.3"
# our goal is to get these into two separate columns
macroRaw7 <- macroRaw7 %>% 
  mutate(productSum = str_match(compTMI, 
                                "Summa produkt=\\s*(.*?)\\s*;Summa viktfaktor")[,2],
         weightFactorSum = str_replace(compTMI, 
                                       "(.*?)Summa viktfaktor=(.*?)", "")) %>% 
  select(-c("compTMI", "NA"))

# fix the numeric values of the 3 main vectors
macroRaw7 <- macroRaw7 %>% mutate(rich = str_extract(`Mfr arter`, "\\d+\\.*\\d*"),
                                  tmiRich = str_extract(`TMI arter`, "\\d+\\.*\\d*"),
                                  tmi = str_extract(TMI, "\\d+\\.*\\d*")) %>% 
  # remove the originals
  select(-c("Mfr arter", "TMI arter", "TMI"))
# save for later
write_rds(macroRaw7, "seMacroRawAPIOld.rds")



# now we get everything out of the observations column
macroRaw8 <- unnest(macroRaw7, col = observations, 
                    keep_empty = TRUE,
                    names_repair = "unique")

# remove NA and useless columns
macroRaw8 <- macroRaw8 %>% select(where(~!all(is.na(.))))
macroRaw8 <- macroRaw8 %>% select(-c("accessRestriction",
                                     "propertyRepresentationType",
                                     "isQualitative",
                                     "sampleId...24")) %>% 
  rename(sampleId = "sampleId...1")

# same deal
macroRaw9 <- unnest(macroRaw8, col = observationValues, 
                    keep_empty = TRUE,
                    names_repair = "unique")
# there's one observation from sample ID 355628 (Gransjön, Mitt) in 2012
# if you ignore its NAs, these rows are also all fully NA
macroRaw9 <- macroRaw9 %>% 
  select(-c("subSampleCd", "unit", "valueType",
            "name", "rank", "group", "isFromSubSample"))
# clean up IDs
macroRaw9 <- macroRaw9 %>% 
  select(-c("observationId...29", "observationValueId")) %>% 
  rename(observationId = "observationId...23")

# same deal, but unnest first and clean up after here
macroRaw10 <- unnest(macroRaw9, col = observationMetadata, 
                         keep_empty = TRUE,
                         names_repair = "unique")
# looping is inefficient and doing it after works fine too
macroRaw10 <- macroRaw10 %>% 
  pivot_wider(names_from = label, values_from = `value...31`) %>% 
  rename(value = "value...29")

#summary options of NAs and unique values
#View(macroRaw10 %>% summarise(across(everything(), ~ n_distinct(.))))
#View(macroRaw10 %>% summarise(across(everything(), ~ sum(is.na(.)))))

macroRaw10 <- macroRaw10 %>% select(-c("Transektnummer",
                                       "Krattdrag/ruta nummer",
                                       # duplicate of value but with -1, NA coding
                                       "TaxonId",
                                       # just another "nothing found" report for absences 
                                       #View(macroRaw10 %>% filter(!is.na(Analysresultat)))
                                       "Analysresultat",
                                       # some artefact from empty observations, i guess
                                       "NA")) %>%
  rename(waterDepth = "Vattendjup (m)",
         samplingMethod = "Provtagningsmetod (biologisk provtagning)",
         siteCoordStartX = "Provplatskoordinat START N/X",
         siteCoordStartY = "Provplatskoordinat START E/Y",
         siteCoordEndX = "Provplatskoordinat SLUT N/X",
         siteCoordEndY = "Provplatskoordinat SLUT E/Y",
         siteCoordCRS = "Provplatskoordinaternas koordinatsystem",
         shorelineDist = "Avstånd från strandlinjen (m)",
         dominantOrgSubstrate = "Bottensubstrat, organiskt, dominerande typ. D1",
         waterLevelDeviation = "Avvikelse från normal vattennivå (m)",
         dominantSubstrateFract = "Dominerande sedimentkaraktär angivet i huvudfraktioner",
         dominantSubstrate = "Dominerande substrattyp",
         dominantInorgSubstrateV1 = "Bottensubstrat, dominerande oorganiskt. D1: SS EN ISO 14688-1",
         dominantInorgSubstrateV3 = "Dominerande oorganiskt bottensubstrat enligt Undersökningstyp v3")

# save now that we're done with basic cleaning
#write_rds(macroRaw10, "seMacroRawAPI.rds")
macroRaw10 <- read_rds("seMacroRawAPI.rds")



# first, these are the new columns from observation data --> pivot_wider
# 828 is the ideal number of rows
#"observationId"            "qualityCode"             
#[25] "analysisDate"             "propertyCode"             "propertyAbbrevName"       "propertyName"            
#[29] "value"                    "waterDepth"               "samplingMethod"           "siteCoordStartX"         
#[33] "siteCoordStartY"          "siteCoordEndX"            "siteCoordEndY"            "siteCoordCRS"            
#[37] "shorelineDist"            "dominantOrgSubstrate"     "waterLevelDeviation"      "dominantSubstrateFract"  
#[41] "dominantSubstrate"        "dominantInorgSubstrateV1" "dominantInorgSubstrateV3"

# ignore: observationId, propertyCode, propertyAbbrevName


# useful information, but too much time to clean: waterLevelDeviation
macroRaw10 <- macroRaw10 %>% select(-c("observationId", 
                                       "propertyCode",
                                       "propertyAbbrevName",
                                       "waterLevelDeviation"))

# convert to numeric first
# value, waterDepth (commas to periods), coordinates, shorelineDist
# convert to fully lowercase
# dominantOrgSubstrate
# dominantSubstrate 0/1 for dominant substrate is organic (any value) 1; otherwise NA = 0
macroRaw11 <- macroRaw10 %>% 
  mutate(value = as.numeric(value),
         waterDepth = as.numeric(str_replace_all(waterDepth,
                                                 ",",
                                                 ".")),
         siteCoordStartX = as.numeric(siteCoordStartX),
         siteCoordStartY = as.numeric(siteCoordStartY),
         siteCoordEndX = as.numeric(siteCoordEndX),
         siteCoordEndY = as.numeric(siteCoordEndY),
         dominantOrgSubstrate = tolower(dominantOrgSubstrate),
         # 1 if the dominant substrate is organic; 0 otherwise
         dominantSubstrate = ifelse(is.na(dominantSubstrate), 0, 1),
         # 1 if analysis happens at a later date
         analysisDate = ifelse(analysisDate != samplingDate, 1, 0))

# maybe would be good to replace "N/A" with NA in qualityCode
# samplingMethod tolower --> convert to binary columns
# inorganic substrate do not overlap, so let's combine them together
# V1 is coarser than V3, so harmonise to V3 categories
# silt; ≤0,063 mm
# sand/fine gravel; >0,063-6,3 mm
# medium-coarse gravel; >6,3-63 mm" 
# stone; >63-200 mm"
# boulder; >200-630 mm"
# large boulder; >630 mm"
# dominantInorgSubstrateV1, dominantInorgSubstrateV3 --> convert to categories from small-to-large
macroRaw11 <- macroRaw11 %>% 
  mutate(qualityCode = ifelse(qualityCode == "N/A", NA, qualityCode),
         samplingMethod = tolower(samplingMethod),
         sampleRake = ifelse(samplingMethod %in% c("kratta", "räfsa", "krattning"), 1, 0),
         sampleSnorkel = ifelse(samplingMethod == "snorkling", 1, 0),
         sampleDive = ifelse(samplingMethod == "dykning", 1, 0),
         sampleBino = ifelse(samplingMethod == "vattenkikare", 1, 0),
         sampleGrid = ifelse(samplingMethod == "inventeringsruta", 1, 0),
         sampleHand = ifelse(samplingMethod == "utan hjälpmedel", 1, 0),
         inorgSubstrate1 = factor(tolower(dominantInorgSubstrateV1),
                                  levels = c("ler ≤ 0,002 mm",
                                             "silt 0,002-0,063 mm",
                                             
                                             "sand 0,063-2,0 mm",
                                             "fingrus 2,0-6,3 mm",
                                             
                                             "mellangrus 6,3-20 mm",
                                             # the 2,0 - 63 can be classed as either 2 or 3 but much more likely for it to be 3
                                             "grovgrus 20-63 mm", "grus 2,0-63 mm",
                                             
                                             "sten 63-200 mm",
                                             
                                             "block 200-630 mm",

                                             "stora block 630-2000 mm", "stora block >630 mm" ,
                                             "mycket stora block >2000 mm"),
                                  labels = c(1, 1,
                                             2, 2,
                                             3, 3, 3,
                                             4,
                                             5,
                                             6, 6, 6)),
         inorgSubstrate3 = factor(tolower(dominantInorgSubstrateV3),
                                  levels = c("ler-silt ≤0,063 mm",
                                             "sand-fingrus >0,063-6,3 mm",
                                             "mellan-grovgrus >6,3-63 mm",
                                             "sten >63-200 mm",
                                             "block >200-630 mm",
                                             "stora block >630 mm"),
                                  labels = c(1, 2, 3, 4, 5, 6))) %>% 
  # remove the ones we're done with
  select(-c(samplingMethod, dominantInorgSubstrateV1, dominantInorgSubstrateV3))

# this needs to be done rowwise
macroRaw11 <- macroRaw11 %>% rowwise() %>% 
  # also replace 0s with NAs because it creates a mess later on
  mutate(dominantInorgSubstrate = na_if(sum(as.numeric(inorgSubstrate1),
                                        as.numeric(inorgSubstrate3), 
                                        na.rm = TRUE), 0))
macroRaw11 <- macroRaw11 %>% select(-c(inorgSubstrate1, inorgSubstrate3))


########################### taxa cleaning ###########################
# helper import if starting the script from halfway, uncomment if needed
#macroRaw7 <- read_rds("seMacroRawAPIOld.rds")

extraCol <- colnames(macroRaw11)[!(colnames(macroRaw11) %in% colnames(macroRaw7))]
baseCol <- colnames(macroRaw11)[(colnames(macroRaw11) %in% colnames(macroRaw7))]

# hard decisions to get this into wide format with species
# first, we need to clean the species though
#n_distinct(macroRaw12$propertyName)
macroSpec <- macroRaw11 %>% group_by(across(all_of(c(baseCol, "propertyName")))) %>% 
  summarise(temp = n()) %>% group_by(propertyName) %>% summarise(n = n())
#View(macroSpec)

# also, go back to use the older data with 1 row per site/station pair
# get list of resampled stations and sites for reference
# first for sites
resampSeSite <- macroRaw7 %>% mutate(samplingDateD = as.Date(samplingDate),
                                     samplingYr = year(samplingDateD)) %>% 
  group_by(samplingSiteId) %>% 
  summarise(nYr = n_distinct(samplingYr)) %>% 
  filter(nYr > 1)

# test that there are no multiple coordinates for stations (deleted) and sites
# no problems, but some missing coordinates for sites
#summarise(xCoord = paste(na.omit(unique(samplingSiteX)), collapse = ";"),
#          yCoord = paste(na.omit(unique(samplingSiteX)), collapse = ";"))


# then for stations
resampSeStation <- macroRaw7 %>% mutate(samplingDateD = as.Date(samplingDate),
                                        samplingYr = year(samplingDateD)) %>% 
  group_by(stationId) %>% 
  summarise(nYr = n_distinct(samplingYr)) %>% 
  filter(nYr > 1)

# commented out because this can only be done later (once macroRaw12 is created)
# get a filtered data for only those stations
#macroRaw12Filt <- macroRaw12 %>% 
#  filter((stationId %in% resampSeStation$stationId) | (samplingSiteId %in% resampSeSite$samplingSiteId))
# and find its unique species for tiebreakers
#macroSpecFilt <- macroRaw12Filt %>% group_by(propertyName) %>% summarise(n = n())

# convert first letter to caps; remove leading/trailing whitespace first
# also turn that one NA value into Biota
# to NA: Algae, Biota, Cladophoraceae, Cladophora, Cladophora glomerata, 
#       Cyanobacteria, Dichothrix orsiniana, Dreissena polymorpha, Batrachospermum,
#       Batrachospermum gelatinosum, Chaetophora lobata, Chlorophyceae, Chlorophyta
#       Ephydatia fluviatilis, Hapalosiphon, Hildenbrandia, Hildenbrandia rubra,
#       Inget fynd, Mesogloia vermiculata, Microcystis, Mougeotia, Nostoc
#       Nostoc pruniforme, Nostoc zetterstedtii, Oedogonium, Paludicola keratophyta
#       Paludicola turfosa, Porifera, Rhizoclonium hieroglyphicum, Rhodophyta,
#       Spongilla, Spongilla lacustris, Spongillidae, Stigonema informe, 
#       Stigonema ocellatum var. globosum, Tolypothrix, Tolypothrix saviczii
#       trådalg, Ulothrix, Ulva, Ulva intestinalis, Ulva pilifera, Vaucheria, 
#       Vaucheria dichotoma
#       
# remove s.lat.  s. lat.
# anything with "×" " x " (all hybrids) to genus sp (keep only first word bc "Carex × saamica" is weird)
# remove agg. (assume somewhat certainty)
# "var" "subsp" remove everything after
# Chara globularis/virgata to Chara
# Chiloscyphus pallescens/polyanthos to Chiloscyphus
# knoppslinga --> Myriophyllum sibiricum
# möja --> Ranunculus
# Nitella flexilis/opaca --> Nitella
# trådstarr --> Carex lasiocarpa

# Barbilophozia lycopodioides only spec in Barbilophozia
# Cardamine pratensis only spec in Cardamine
# Deschampsia cespitosa only spec in Deschampsia
# Iris pseudacorus --> iris
# Jungermannia eucordifolia --> jungermannia
# Oncophorus integerrimus --> Oncophorus
# Philonotis fontana --> philonotis
# Pohlia wahlenbergii --> pohlia
# these next two make me rethink the rule (50+ observations at species level; 1 at genus)
# riccia fluitans --> Riccia
# Scutellaria galericulata --> Scutellaria
# Spirogyra groenlandica --> Spirogyra

# unclear about:
# Charophyceae = Chara??
# trådalg
# what to do about Ranunculus Batrachium agg.
# also Sphagnum subgenera e.g., Sphagnum subg. Subsecunda, Cuspidata, Sphagnum, Acutifolia

# checking with other datasets
# Denmark has Pohlia sphagnicola, Riccia fluitans
# no Sphagnum subgenera, Ranunculus Batrachium
# Finland has Iris pseudacorus, Philonotis fontana, Philonotis seriata, Pohlia wahlenbergii, Riccia fluitans 
# also no subgenera for Sphagnum or Ranunculus

macroSpec <- macroSpec %>% mutate(taxaClean = str_trim(str_to_sentence(propertyName)))
excludeTaxa <- c(NA, "Algae", "Biota", "Cladophoraceae", "Cladophora", "Cladophora glomerata",
                 "Cyanobacteria", "Dichothrix orsiniana", "Dreissena polymorpha", "Batrachospermum",
                 "Batrachospermum gelatinosum", "Chaetophora lobata", "Chlorophyceae", "Chlorophyta",
                 "Ephydatia fluviatilis", "Hapalosiphon", "Hildenbrandia", "Hildenbrandia rubra",
                 "Inget fynd", "Mesogloia vermiculata", "Microcystis", "Mougeotia", "Nostoc",
                 "Nostoc pruniforme", "Nostoc zetterstedtii", "Oedogonium", "Paludicola keratophyta",
                 "Paludicola turfosa", "Porifera", "Rhizoclonium hieroglyphicum", "Rhodophyta",
                 "Spongilla", "Spongilla lacustris", "Spongillidae", "Stigonema informe",
                 "Stigonema ocellatum var. globosum", "Tolypothrix", "Tolypothrix saviczii",
                 "Trådalg", "Ulothrix", "Ulva", "Ulva intestinalis", "Ulva pilifera", "Vaucheria",
                 "Vaucheria dichotoma", "Aegagropila linnaei")
macroSpec <- macroSpec %>% 
  # replace with the "Biota" placeholder if in  manually screened excludeTaxa list
  mutate(taxaClean2 = ifelse(taxaClean %in% excludeTaxa,
                             "Biota",
                             taxaClean),
         # limit to only first word if any of these characters exist in the strong
         taxaClean3 = ifelse(str_detect(taxaClean2, "×") | str_detect(taxaClean2, " x ") | str_detect(taxaClean2, "/"),
                             word(taxaClean2, 1),
                             taxaClean2),
         # remove these patterns bc they're always at the end of a string
         taxaClean4 = str_trim(str_replace_all(taxaClean3, c("s.lat."="",
                                                             "s. lat."="",
                                                             "agg."=""))),
         # remove everything after "var."
         taxaClean5 = str_trim(str_remove(taxaClean4, "(var\\.).*")),
         # probably can do this together with above but lazy
         taxaClean6 = str_trim(str_remove(taxaClean5, "(subsp\\.).*")),
         # replace some known issues with proper coding
         taxaClean7 = str_replace_all(taxaClean6, # $ is the string end delimiter
                                      c("Möja$"="Ranunculus",
                                        "Knoppslinga$"="Myriophyllum sibiricum",
                                        "Trådstarr$"="Carex lasiocarpa",
                                        # if resampled data has only species, then prefer species
                                        "Cardamine$" = "Cardamine pratensis",
                                        # when >50% of obs are in a single species, set genus to species (if only 1)
                                        "Iris$" = "Iris pseudacorus",
                                        "Riccia$" = "Riccia fluitans",
                                        "Scutellaria$" = "Scutellaria galericulata",
                                        # if only one of the two is found in another dataset, prefer that one
                                        "Philonotis$" = "Philonotis fontana",
                                        "Pohlia$" = "Pohlia wahlenbergii",
                                        # otherwise setting genus is more conservative 
                                        "Barbilophozia lycopodioides" = "Barbilophozia",
                                        "Deschampsia cespitosa" = "Deschampsia",
                                        "Jungermannia eucordifolia" = "Jungermannia",
                                        "Oncophorus integerrimus" = "Oncophorus",
                                        "Spirogyra groenlandica" = "Spirogyra",
                                        # get rid of subgenera manually for now
                                        "Sphagnum subgen. Acutifolia sect. Acutifolia" = "Sphagnum",
                                        "Sphagnum subg. Cuspidata" = "Sphagnum",
                                        "Sphagnum subg. Sphagnum" = "Sphagnum",
                                        "Sphagnum subg. Subsecunda" = "Sphagnum",
                                        # conflicting thoughts about Ranunculus bc many observations
                                        # do it to be consistent though
                                        "Ranunculus batrachium" = "Ranunculus")))

# see how much we've cut down
n_distinct(macroSpec$taxaClean7)
# prepare for re-merging
macroSpecClean <- macroSpec %>% select(propertyName, taxaClean7) %>% 
  rename(taxaClean = "taxaClean7")



# then merge with the original data
macroRaw11 <- left_join(macroRaw11, macroSpecClean, by = "propertyName")

########################### make data wide & clean covariates ##################
# concatenate: qualityCode
# sum: value, analysisDate
# min/max: waterDepth, coordinates (hopefully CRS is consistent)
# samplingMethod and dominant inorganic substrate were pre-processed
macroRaw12 <- macroRaw11 %>% 
  group_by(across(all_of(c(baseCol, "taxaClean")))) %>% 
  summarise(
    # value is always 1 for a recorded taxa, test yourself if in doubt
    #macroRaw11 %>% filter(value != 1)
    #macroRaw11 %>% filter(is.na(value))
    value = sum(value, na.rm = TRUE),
    analysisDate = sum(analysisDate, na.rm = TRUE),
    waterDepthMin = min(waterDepth, na.rm = TRUE),
    waterDepthMax = max(waterDepth, na.rm = TRUE),
    siteCoordStartXMin = min(siteCoordStartX, na.rm = TRUE),
    siteCoordStartYMin = min(siteCoordStartY, na.rm = TRUE),
    siteCoordEndXMin = min(siteCoordEndX, na.rm = TRUE),
    siteCoordEndYMin = min(siteCoordEndY, na.rm = TRUE),
    siteCoordStartXMax = max(siteCoordStartX, na.rm = TRUE),
    siteCoordStartYMax = max(siteCoordStartY, na.rm = TRUE),
    siteCoordEndXMax = max(siteCoordEndX, na.rm = TRUE),
    siteCoordEndYMax = max(siteCoordEndY, na.rm = TRUE),
    siteCoordCRS = paste(na.omit(unique(siteCoordCRS)), collapse = ";"),
    qualityCode = paste(na.omit(unique(qualityCode)), collapse = ";"),
    # updated due to variable recoding
    inorgSubstrateMin = min(dominantInorgSubstrate, na.rm = TRUE),
    inorgSubstrateMax = max(dominantInorgSubstrate, na.rm = TRUE),
    inorgSubstrateMed = median(dominantInorgSubstrate, na.rm = TRUE),
    inorgSubstrateMean = mean(dominantInorgSubstrate, na.rm = TRUE),
    sampleRake = sum(sampleRake, na.rm = TRUE),
    sampleSnorkel = sum(sampleSnorkel, na.rm = TRUE),
    sampleDive = sum(sampleDive, na.rm = TRUE),
    sampleBino = sum(sampleBino, na.rm = TRUE),
    sampleGrid = sum(sampleGrid, na.rm = TRUE),
    sampleHand = sum(sampleHand, na.rm = TRUE)) %>% 
  ungroup()

# replace Inf, -Inf, NaN, empty strings with NA
macroRaw12 <- macroRaw12 %>% mutate(across(where(is.double), ~na_if(., Inf)),
                                    across(where(is.double), ~na_if(., -Inf)),
                                    across(where(is.double), ~na_if(., NaN)),
                                    across(where(is.character), ~na_if(., "")))
# remove the 0.5 from medians by rounding in the direction of the mean
macroRaw12 <- macroRaw12 %>% rowwise() %>% 
  mutate(inorgSubstrateMed = ifelse(all.equal(inorgSubstrateMed, round(inorgSubstrateMed)) == TRUE, 
                                    inorgSubstrateMed,
                                    ifelse(inorgSubstrateMean > inorgSubstrateMed,
                                           ceiling(inorgSubstrateMed),
                                           floor(inorgSubstrateMed))))

# remove qualityCode
# sampling strategies code as binary rather than sum
# everything else keep as-is
# NOTE: macroRaw7 saved on line 155
extraCol2 <- colnames(macroRaw12)[!(colnames(macroRaw12) %in% colnames(macroRaw7))]
baseCol2 <- colnames(macroRaw12)[(colnames(macroRaw12) %in% colnames(macroRaw7))]

# you can't both pivot_wider and summarise, so my solution is to do separately and remerge after
# first, species
macroRaw13Spec <- macroRaw12 %>% select(c(baseCol2, "taxaClean", "value")) %>% 
  pivot_wider(names_from = taxaClean, values_from = value, 
              # assume all missing values are pure absences
              values_fill = 0) %>% 
  # no longer need this 0 observation placeholder
  # removed here and not earlier to retail 0 observation rows (of which there are none)
  select(-Biota) %>% 
  # order columns alphabetically
  select(baseCol2, sort(colnames(.)))

# then do the covariates, mostly same as before
macroRaw13Cov <- macroRaw12 %>% 
  group_by(across(all_of(c(baseCol2)))) %>% 
  summarise(analysisDate = sum(analysisDate, na.rm = TRUE),
            waterDepthMin = min(waterDepthMin, na.rm = TRUE),
            waterDepthMax = max(waterDepthMax, na.rm = TRUE),
            siteCoordStartXMin = min(siteCoordStartXMin, na.rm = TRUE),
            siteCoordStartYMin = min(siteCoordStartYMin, na.rm = TRUE),
            siteCoordEndXMin = min(siteCoordEndXMin, na.rm = TRUE),
            siteCoordEndYMin = min(siteCoordEndYMin, na.rm = TRUE),
            siteCoordStartXMax = max(siteCoordStartXMax, na.rm = TRUE),
            siteCoordStartYMax = max(siteCoordStartYMax, na.rm = TRUE),
            siteCoordEndXMax = max(siteCoordEndXMax, na.rm = TRUE),
            siteCoordEndYMax = max(siteCoordEndYMax, na.rm = TRUE),
            siteCoordCRS = paste(na.omit(unique(siteCoordCRS)), collapse = ";"),
            # the median of median and mean of mean are a bit weird, but whatever
            inorgSubstrateMin = min(inorgSubstrateMin, na.rm = TRUE),
            inorgSubstrateMax = max(inorgSubstrateMax, na.rm = TRUE),
            inorgSubstrateMed = median(inorgSubstrateMed, na.rm = TRUE),
            inorgSubstrateMean = mean(inorgSubstrateMean, na.rm = TRUE),
            sampleRake = sum(sampleRake, na.rm = TRUE),
            sampleSnorkel = sum(sampleSnorkel, na.rm = TRUE),
            sampleDive = sum(sampleDive, na.rm = TRUE),
            sampleBino = sum(sampleBino, na.rm = TRUE),
            sampleGrid = sum(sampleGrid, na.rm = TRUE),
            sampleHand = sum(sampleHand, na.rm = TRUE)) %>% 
  # always remember to ungroup
  ungroup()

# clean up covariates first
# get rid of Inf, -Inf, NaN, and empty strings again
macroRaw13Cov <- macroRaw13Cov %>% mutate(across(where(is.double), ~na_if(., Inf)),
                                          across(where(is.double), ~na_if(., -Inf)),
                                          across(where(is.double), ~na_if(., NaN)),
                                          across(where(is.character), ~na_if(., "")))
# again do the median rounding of inorganic substrates
macroRaw13Cov <- macroRaw13Cov %>% rowwise() %>% 
  mutate(inorgSubstrateMed = ifelse(all.equal(inorgSubstrateMed, round(inorgSubstrateMed)) == TRUE, 
                                    inorgSubstrateMed,
                                    ifelse(inorgSubstrateMean > inorgSubstrateMed,
                                           ceiling(inorgSubstrateMed),
                                           floor(inorgSubstrateMed)))) %>% 
  select(-inorgSubstrateMean)

# then merge both datasets together
macroRaw13 <- full_join(macroRaw13Cov, macroRaw13Spec, by = baseCol2)

# fix some outsanding covariate issues
macroRaw14 <- macroRaw13 %>% 
  mutate(hydromorphImpact = hydromorphImpact %>% 
           replace_when(str_detect(tolower(.),
                                   "saknas") ~ NA),
   lakeDrainage = str_to_sentence(lakeDrainage) %>% 
     replace_when(str_detect(., 
                             "saknas") ~ NA),
   # do this in 2 steps to be safe, first deal with non-numbers
   # done manually so please check if re-running with new API call
   secchi = str_trim(secchi) %>% 
     replace_values("2,36-2,76" ~ "2.56", # mean of upper and lower bound
                    "1,8-4,5" ~ "3.15",
                    "1,1-1,4" ~ "1.25",
                    "1,1-1,3" ~ "1.2",
                    "1-1,6" ~ "1.3",
                    "1,8-2,1" ~ "1.95",
                    "1 -1,25" ~ "1.125",
                    "1,7 meter" ~ "1.7",
                    "1,6-2,5" ~ "2.05",
                    "2,63-3,8" ~ "3.215"),
   secchi = case_when(
     # special case for greater than signs
     str_detect(secchi, ">") | str_detect(secchi, "≥") ~
       as.numeric(str_replace_all(str_extract_all(secchi, 
                                                  "\\(?[0-9,.]+\\)?"),
                                  # set to slightly more than the limit
                                  ",", "."))*1.1,
     # then, let NA by coercion fix 99% of the mess
     .default = as.numeric(str_replace_all(secchi, ",", "."))),
   rich = as.integer(rich),
   tmiRich = as.integer(tmiRich),
   # for temporal processing
   samplingDateD = as.Date(samplingDate),
   samplingYr = year(samplingDateD),
   samplingM = month(samplingDateD),
   siteIdUnique = paste0(samplingSiteX, samplingSiteY)) %>%
  # duplicate IDs
  select(-c("nationalStationId", "nationalSiteId")) %>% 
  ungroup()

########################### handle spatial uniqueness ##########################
# understand the process of duplicates 
tempSum <- macroRaw14 %>% 
  summarise(.by = siteIdUnique,
            nTot = n(),
            nYr = n_distinct(samplingYr),
            nDate = n_distinct(samplingDateD),
            nID = n_distinct(samplingSiteId))

# remove resurveys for analysis (issues with AR error)
# also errata
macroSpat <- macroRaw14 %>% 
  # only sites w/ coordinates
  # site 6679 is spatial duplicate and inaccurate (between basins)
  filter(siteIdUnique != "NANA" & samplingSiteId != 6679) %>% 
  # then, we choose surveys with
  group_by(samplingSiteId) %>% 
  # 1) earliest date for a site (to minimise experience impacts)
  slice_min(samplingDateD, n = 1) %>% 
  # 2) higher observation count (more survey effort) if multiple same-day surveys
  slice_max(observationCount, n = 1) %>% 
  # remove unneeded columns
  select(-c(sampleId,
            flagCodes,
            relativeLocationType, # only Gavel-Långsjön with Depth
            stationId,
            stationName,
            stationType, # only lakes
            stationEUID,
            stationCoordinateX,
            stationCoordinateY,
            waterLevel, # too many NAs
            lakeDrainage, # too many NAs
            hydromorphImpact, #too many NAs
            inventoryTransect,
            inventoryQuadrat,
            insertDate,
            editDate,
            lastEdit,
            productSum,
            weightFactorSum,
            analysisDate,
            inorgSubstrateMin, # too many NAs            
            inorgSubstrateMax,
            inorgSubstrateMed,
            siteIdUnique))
  
# export for analysis
# to inspect coordinates manually, check out https://rl.se/rt90
write_csv(macroSpat, "seCommSpat.csv")