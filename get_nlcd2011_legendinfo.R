library(rvest)
# this code gets the table from the NLCD webpage
# and rearranges it for columns with the class ID, name, description

# link to NLCD legend page
webpage <- read_html("https://www.mrlc.gov/nlcd11_leg.php")

legend_table <-
  webpage %>% 
  html_nodes("td table") %>% 
  .[1] %>% html_table(fill = TRUE)

legend_table <- legend_table[[1]] %>% 
  dplyr::slice(-1) %>% 
  dplyr::rename(Class = X1, Description = X2)

descriptions <- legend_table %>% 
  dplyr::pull(Description) %>%
  stringr::str_split(pattern = "-", 
  simplify = TRUE) %>% as.data.frame() %>%
  dplyr::filter(!V1 == "") %>%
  dplyr::mutate(description = str_c(V2, V3, V4)) %>%
  dplyr::rename(class = V1) %>% 
  dplyr::select(class, description)

legend_table$class_name <- legend_table %>% 
  dplyr::pull(Description) %>%
  stringr::str_split(pattern = "-", 
                     simplify = TRUE) %>% 
  as.data.frame() %>% dplyr::pull(V1)

legend_table <- legend_table %>%
  dplyr::select(Class, class_name) %>%
  dplyr::filter(class_name != "") %>%
  dplyr::left_join(descriptions, by = c("class_name" = "class"))

write.csv(legend_table, "nlcd_legend_2011.csv", row.names = FALSE)
