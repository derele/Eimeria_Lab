## transfer the xml files from the nanodrop to your local storage
## open the console and use the following code to convert from xml to csv
## (if you don't have the package that includes the function ssconvert -->
# use the code 'apt install gnumeric')
## ssconvert sourcefile.xml destinationfile.csv

# Change your working directory to where your documents are located
#setwd("~/folder_with_csv_files")

# install libraries
library(dplyr)

# import files
## first we list all the file names and we save it in a list
list_of_names<- list.files(pattern = "*.csv")

# creade a function to read the csv files and remove the unecessary column X.
read_each_df <- function(x) {
  # read the file
  df1 <- read.csv(x)
  # delete column X.
  df1 <- df1 %>%
    dplyr::select(-X.)
}

## create a list with read data frames
list_of_df <- lapply(list_of_names, read_each_df)

# merge the data frames together
csv_product <- Reduce(rbind, list_of_df)


# if you need to remove a row (double measurements of a sample) then use
# csv_product <- csv_product(-c(1), )

# save output
# Change your working directory to where your documents are located
setwd(".....")

write.csv(x = csv_product, "~/Nano_drop_merged_files.csv", row.names = FALSE)
