#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(readxl))

# Travel Data from 2019 for overseas vists to the UK.

#This data is outdated (pre lockdowns ect).
## It might be worth looking at approved travel locations?

# Ensuring the UK_scheme_data dir is present, and if not creating it
if(length(list.files(pattern = "UK_scheme_data")) > 0){
  cat("UK_scheme_data directory found \n")
}else{
  dir.create("UK_scheme_data")
}

# DownLoading files
travel_file <- "overseas_travel_data_2019q1_to_2020q1.xlsx"

if(length(list.files("UK_scheme_data/",
                     pattern = travel_file)) > 0){
  cat("Already Downloaded Travel Data \n")
}else{
  download.file("https://www.ons.gov.uk/file?uri=%2fpeoplepopulationandcommunity%2fleisureandtourism%2fdatasets%2foverseastravelandtourism%2fcurrent/osq12020.xlsx",
                            destfile = paste0("UK_scheme_data/", travel_file))}

country_file <- "country_data.csv"
if(length(list.files("UK_scheme_data/",
                     pattern = country_file)) > 0){
  cat("Already Downloaded Country Data \n")
}else{
  download.file("https://raw.githubusercontent.com/lukes/ISO-3166-Countries-with-Regional-Codes/master/all/all.csv",
                destfile = paste0("UK_scheme_data/", country_file))
}

# Reading in the travel data.
## Some of these values are hard-coded. And are still working at 2021/10/19
suppressMessages(raw_data <- read_xlsx(paste0("UK_scheme_data/", travel_file), sheet = "Table 11",
                          col_names = FALSE))
colnames(raw_data) <- 1:ncol(raw_data)

# Selecting data for all of 2019 and q1 2020
data <- select(raw_data, 1,26,83)
colnames(data) <- c("name", "total_2019", "q1_2020")

# Removing some unwanted lines, and creating a total count
data2 <- data %>% filter(is.na(name) == FALSE) %>%
  .[-(1:2),] %>%
  mutate(across(.cols = 2:3, as.numeric)) %>%
  mutate(total = total_2019+q1_2020)

# To work out proportion, the total count is extacted
total_world <- data2 %>% 
  filter(name == "Total World") %>%
  select(total) %>%
  as.numeric()
  
data2<-mutate(data2, prop = total / total_world) %>%
  arrange(desc(prop))

suppressMessages(country_data <- read_csv(paste0("UK_scheme_data/", country_file)))

# Ensuring names are the same.
country_data$name[country_data$name == "United States of America"] <- "USA"
country_data$name[country_data$name == "Russian Federation"] <- "Russia"
country_data$name[country_data$name == "Ireland"] <- "Republic of Ireland"
country_data$name[country_data$name == "Czechia"] <- "Czech Republic"
country_data$name[country_data$name == "Hong Kong"] <- "Hong Kong (China)"
country_data$name[country_data$name == "Hong Kong"] <- "Hong Kong (China)"

# Joining the country data and the travel data
region_data <- inner_join(country_data, data2, by = "name") %>% 
  select(name, region, `sub-region`, total, prop)

# Grouping by sub-region
region_data2 <- region_data %>%
  group_by(`sub-region`) %>%
  summarise(prop = sum(prop),
            count = sum(total)) %>%
  arrange(desc(count))

percent <- sum(region_data2$prop)*100

cat(paste0(signif(percent,3), "% of travel data used \n"))
cat("Writing UK_scheme_data/region_data.csv")
write.csv(region_data2, "UK_scheme_data/region_data.csv")
