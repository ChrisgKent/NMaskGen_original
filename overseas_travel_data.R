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

# DownLoading files for overseas travel data.
overseas_vist_file <- "overseas_travel_data_2019q1_to_2020q1.xlsx"
if(length(list.files("UK_scheme_data/",
                     pattern = overseas_vist_file)) > 0){
  cat("Already Downloaded Travel Data \n")
}else{
  download.file("https://www.ons.gov.uk/file?uri=%2fpeoplepopulationandcommunity%2fleisureandtourism%2fdatasets%2foverseastravelandtourism%2fcurrent/osq12020.xlsx",
                            destfile = paste0("UK_scheme_data/", overseas_vist_file))}

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
suppressMessages(ov_raw_data <- read_xlsx(paste0("UK_scheme_data/", overseas_vist_file), sheet = "Table 11",
                          col_names = FALSE))
colnames(ov_raw_data) <- 1:ncol(ov_raw_data)

# Selecting data for all of 2019 and q1 2020
ov_data <- select(ov_raw_data, 1,26,83)
colnames(ov_data) <- c("name", "total_2019", "q1_2020")

# Removing some unwanted lines, and creating a total count
ov_data2 <- ov_data %>% filter(is.na(name) == FALSE) %>%
  .[-(1:2),] %>%
  mutate(across(.cols = 2:3, as.numeric)) %>%
  mutate(total = total_2019+q1_2020)

# To work out proportion, the total count is extacted
ov_total_world <- ov_data2$total[ov_data2$name == "Total World"]
  
ov_data2<-mutate(ov_data2, prop = total / ov_total_world) %>%
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
ov_region_data <- inner_join(country_data, ov_data2, by = "name") %>% 
  mutate(sub_region = `sub-region`) %>%
  select(name, region, sub_region, total, prop)

# Grouping by sub-region
ov_region_data2 <- ov_region_data %>%
  group_by(sub_region) %>%
  summarise(prop_overseas = sum(prop),
            count_overseas = sum(total)) %>%
  arrange(desc(count_overseas))

percent <- sum(ov_region_data2$prop_overseas)*100

cat(paste0(signif(percent,3), "% of overseas travel data used \n"))
cat("Writing UK_scheme_data/region_data.csv")
write.csv(ov_region_data2, "UK_scheme_data/region_data.csv")

# Starting on the tourism data
## Downloading the file if not present
tourism_file <- "ukq12020.xlsx"
if(length(list.files("UK_scheme_data/",
                     pattern = tourism_file)) > 0){
  cat("Already Downloaded Travel Data \n")
}else{
  download.file("https://www.ons.gov.uk/file?uri=%2fpeoplepopulationandcommunity%2fleisureandtourism%2fdatasets%2fukresidentsvisitsoverseasquarterly%2fquarter1jantomar2020/ukq12020.xlsx",
                destfile = paste0("UK_scheme_data/", tourism_file))}

tourism_data <- readxl::read_xlsx(paste0("UK_scheme_data/", tourism_file), 
                          sheet = "Table 9",
                          col_names = FALSE) %>%
  filter(is.na(.[,1]) == FALSE) %>%
  .[-(1:2),]

tourism_data_2019 <- tourism_data %>% select(1,26, 83)
colnames(tourism_data_2019) <- c("name", "q1_4_2019","q1_2020")
tourism_data_2019 <- tourism_data_2019 %>% 
  mutate(across(.cols = c(2,3), as.numeric)) %>%
  mutate(total = q1_4_2019+q1_2020)

total_tourism <- tourism_data_2019$total[tourism_data_2019$name == "Total World"]

tourism_data_2019 <- tourism_data_2019 %>% mutate(prop = total / total_tourism)

# Joining the country data and the travel data
tourism_region_data <- inner_join(country_data, tourism_data_2019, by = "name") %>% 
  mutate(sub_region = `sub-region`) %>% 
  select(name, region, sub_region, total, prop)

# Grouping by sub-region
tourism_region_data2 <- tourism_region_data %>%
  group_by(sub_region) %>%
  summarise(prop_tourism = sum(prop),
            count_tourism = sum(total)) %>%
  arrange(desc(count_tourism))


percent <- sum(tourism_region_data2$prop_tourism)*100

cat(paste0(signif(percent,3), "% of travel data used \n"))
cat("Writing UK_scheme_data/tourism_data.csv")
write.csv(tourism_region_data2, "UK_scheme_data/tourism_data.csv")


# Combining the data


data <- full_join(tourism_region_data2,ov_region_data2, by = "sub_region") %>%
  mutate(total_movement = count_overseas+count_tourism,
         total_prop = total_movement/sum(total_movement)) %>%
  arrange(desc(total_prop))

write_csv(data, "UK_scheme_data/total_travel_data.csv")



