library(rvest)

# link to NARS data
webpage <- read_html("https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys")

# get links from data column
data_links <- webpage %>% html_nodes("td:nth-child(3) a") %>% html_attr("href")
# get links from metadata column
metadata_links <- webpage %>% html_nodes("td:nth-child(4) a") %>% html_attr("href")

# download a specific file using url and link
url <- "https://www.epa.gov"
data_links_full <- paste0(url, data_links)
metadata_links_full <- paste0(url, metadata_links)

# this gets the info from the table
tbls <- html_nodes(webpage, "table")
tbls_ls <- webpage %>%
  html_nodes("table") %>%
  .[1] %>%
  html_table(fill = TRUE)
nars_data <- tbls_ls[[1]]

nars_data$data_link <- data_links_full
# nars_data$metadata_link <- metadata_links_full

write.csv(nars_data, "nars_data_table.csv", row.names = FALSE)
