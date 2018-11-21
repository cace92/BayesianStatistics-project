# -------------------------------
# Collection of various functions
# -------------------------------

# Function to clear the console
clc <-function() {
  cat("\f");
}

# Function that read the dataset
read_data <-function(directory) {
  setwd(directory);
  colClass <-c("character", "integer", "character", "integer", "character", "character", "character", "integer", 
                 rep("numeric", 12), "integer", rep("numeric", 20));
  data <-read.xlsx("data.xlsx", sheetName = "DATI", rowIndex = 3:1212, colIndex = 2:42, as.data.frame = TRUE,
                            header = FALSE, colClasses = colClass);
  colnames(data) <-c("lev", "sem", "camp", "cfu", "ssd", "cds", "sex", "age", "pm_i", "pf_i", "pi_i", "pf_i",
                     "pm_cs", "pf_cs", "pi_cs", "pf_cs", "pm_ac", "pf_ac", "pi_ac", "pf_ac", "ns",
                     "mq1", "mq2", "mq3", "mq4", "mq5", "mq6", "mq7", "mq8", "mq9", "mq10",
                     "mq11", "mq12", "mq13", "mq14", "mq15", "mq16", "mq17", "mq18", "mq19", "mq20");
  data$cds <-gsub("\n", "", data$cds, fixed = TRUE);
  save(data, file = "data.rda");
}