## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)
# Load raw data from .csv file
example_data <- read.csv("data-raw/example_data.csv")
# Apply preprocessing...
example_data <- dplyr::select(example_data, -X)

colnames(example_data) = c("MA1_1", "MA1_2", "MA2_1", "MA2_2")
rownames(example_data) = c("ES", "SE")
# Save the cleaned data in the required R package location
usethis::use_data(example_data)