#### Trualling imputation without droppping sentinal lymph node

## breast
library(haven)
library(dplyr)

breast <- read_sas("/Users/Sam/Desktop/Imperial/Term_3/Week_7/NCT00041119_Breast/AllProvidedFiles_194/eval46_3_finala.sas7bdat")

dim(breast)

# Drop these columns TREAT_ASSIGNED,GROUP_ID, OH012,OH013,OH014, OH027, cod,
# event, Preamend, Agent, Receptor, survmos, OH032, OH037, OH016 (becasue everyone had had prior chemo) 
# Stra2, scase, Survstat, Sstat, ETHNIC_ID, Length [done]

# engineer these to have NAs for "unknown" and impute - OH006 [done]


# categories requireing one hot encoding - indrx, RACE_ID, tsize, agecat, OH002,
#OH005, OH006, OH028
#check NAs. 


sum(breast$dfsmos < 1)
sum(breast$dfsmos == 0)


tail(sort(breast$dfsmos), 5)

hist(breast$dfsmos)

dim(breast)
## REMOVE PEOPLE WHOS SURVIVAL EQUALLED EXACLT 0 I.E. THAT THEY SUPPOSEDLY HAD AN EVENT THE MOMENT THE TRIAL STARTED

breast <- breast %>%
  filter(dfsmos != 0)

summary(breast$dfsmos)

### Remove all ineligble patients
breast <- breast%>%
  filter(elig == 2)

dim(breast)

## Dropping columns
data <- breast%>%
  select(-c(TREAT_ASSIGNED,GROUP_ID, preamend, agent, receptor,OH012,OH013, 
            OH014, OH027,OH032, OH037,cod, event, survmos, OH016, stra2,
            scase, survstat, SSTAT, ETHNIC_ID, length))
dim(data)


## Removing ineligble participants 

data <- data%>%
  filter(elig == 2)

# remove elig
data <- data%>%
  select(-c(elig))

### Feature enginnering 
##HER status
table(data$OH006)


## CREATING A BINARY VARIABLE OF WHETHER HER2 IS POSITIVE OR NEGATIVE (OR MISSING)
# Assuming data$OH006 is the column containing the categories

# Create a new column to store the converted values
data$converted_OH006 <- NA
# Convert categories 1-4 to "1"
data$converted_OH006[data$OH006 %in% c(1, 2, 3, 4)] <- "1"
# Convert category 6 to "0"
data$converted_OH006[data$OH006 == 6] <- "0"
# Leave categories 5 and 7 as NA
data$converted_OH006[data$OH006 %in% c(5, 7)] <- NA
# Verify the updated values
table(data$converted_OH006)
# Replace the original column with the updated column
data$OH006 <- data$converted_OH006

table(data$OH006)
sum(is.na(data$OH006))

# Convert the OH006 column to a numeric
data$OH006 <- as.numeric(data$OH006)

# Verify the updated data type
table(data$OH006)

# remove converted_OH006
data <- data%>%
  select(-c(converted_OH006))

# keep for now but looks like it will haveto be reomved becuase of high missingness. 
table(data$OH036)
sum(is.na(data$OH036))


# Feature engineering of num_pos_nodes
# Feature enginnering  num_pos_nodes

table(data$num_pos_nodes)

data$new_num_pos_nodes <- ifelse(data$num_pos_nodes == 0, "0",
                                 ifelse(data$num_pos_nodes == 1, "1", "2+"))

table(data$new_num_pos_nodes)

data$num_pos_nodes <- data$new_num_pos_nodes
table(data$num_pos_nodes)

# Convert the num_pos_nodes column to a factor
data$num_pos_nodes <- as.factor(data$num_pos_nodes)


# remove new_num_pos_nodes
data <- data%>%
  select(-c(new_num_pos_nodes))


###plot the number
library(ggplot2)

ggplot(data, aes(x = dfsstat)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black") +
  xlab("Disease-Free Survival Status") +
  ylab("Frequency") +
  ggtitle("Histogram of Disease-Free Survival Status")

table(data$dfsstat)

# Rename columns:

colnames(data)

# Load dplyr library
library(dplyr)

# Rename columns
data <- rename(data,
               mask_id = mask_id,
               race_id = RACE_ID,
               treatment_assigned = indrx,
               menopause_status = stra1,
               tumor_side = OH002,
               receptor_status_er = OH003,
               receptor_status_pgrn = OH004,
               histologic_grade = OH005,
               her2_status = OH006,
               prior_hormonal_therapy = OH011,
               most_extensive_primary_surgery = OH028,
               tumor_size = tsize,
               age_category = agecat,
               disease_free_survival_status = dfsstat,
               disease_free_survival_months = dfsmos,
               num_pos_nodes = num_pos_nodes,
               sentinal_node_biopsy_positive =OH036)


# Convert the column to factor if it's not
data$race_id <- as.factor(data$race_id)

# Rename the levels
levels(data$race_id) <- c("White", "Black", "Asian", 
                          "Native_Hawaiian_American_Alaskan_or_Pacific_Islander", 
                          NA)

table(data$race_id)
sum(is.na(data$race_id))

### Plotting percentage missingess
# Load required packages
library(ggplot2)

# Calculate percentage of missing values for each variable
missing_values <- sapply(data, function(x) sum(is.na(x)) / length(x)) * 100

# Convert to a dataframe for plotting
missing_values_df <- data.frame(variable = names(missing_values), 
                                missing_percentage = as.numeric(missing_values))

# Order by descending missingness
missing_values_df <- missing_values_df[order(-missing_values_df$missing_percentage), ]

# Create barplot
ggplot(missing_values_df, aes(x=reorder(variable, missing_percentage), y=missing_percentage)) +
  geom_bar(stat='identity', fill='steelblue') +
  coord_flip() +  # to flip the axes
  theme_minimal() +
  labs(x = "Variable", y = "Missingness (%)", title = "Missingness Percentage of Each Variable")

## This time we;re not dropping sentinal_node_biopsy_positive due to missingness of 20%, but going to try and impute it. 
# data <- data%>%
#   select(-c(sentinal_node_biopsy_positive))

str(data)

##### Plotting the Kaplan Meier (KM) curves ####
library(survival)
library(survminer)
# Create  survival object
surv_obj <- Surv(data$disease_free_survival_months, 
                 data$disease_free_survival_status == 1)

# Fit survival curve
fit <- survfit(surv_obj ~ 1) 

# Plot train survival curve
ggsurvplot(fit, data = data,
           xlab = "Time in Months",
           ylab = "Survival Probability",
           title = "Kaplan-Meier Curve for Overall Disease-Free Survival",
           legend.labs = c("All Treatment Arms"),
           risk.table = T)


# Fit survival curves stratified by treatment
fit_stratified <- survfit(surv_obj ~ data$treatment_assigned)

# Plot  survival curves
ggsurvplot(fit_stratified, data = data,
           xlab = "Time in Months",
           ylab = "Survival Probability",
           title = "Kaplan-Meier Curves Stratified by Treatment",
           legend.labs = c("CC4","CC6", "T4", "T6"),
           risk.table = TRUE,
           palette = "jco",
           pval = TRUE) # included result of log rank test in the plot


##logrank tests
# Perform the log-rank test
logrank_train <- survdiff(surv_obj ~ data$treatment_assigned)

# Print the results
print(logrank_train)

## p = 0.2 N0 significant differences in the treatment


####### Split data into training & testing set #######


#remove.packages("caret")
#install.packages("caret")
library(caret)
set.seed(42)

# Split the data_encoded into training (80%) and testing (20%) sets, stratified by the outcome variable. 
trainIndex <- createDataPartition(data$disease_free_survival_status, p = 0.8, list = FALSE, times = 1)

trainSet <- data[trainIndex, ]

testSet <- data[-trainIndex, ]

##MICEranger imputation
library(miceRanger)
# Impute missing data in the train set
imputed_data_train <- miceRanger(
  trainSet,
  num.trees = 100,
  num.iterations = 5,
  verbose = FALSE
)

completed_train_data <- miceRanger::completeData(imputed_data_train)

train_data_imputed1 <- completed_train_data$Dataset_1

# Calculate the number of missing values in each column
missing_counts <- apply(trainSet, 2, function(x) sum(is.na(x)))

# Print the missing value counts
print(missing_counts)
dim(train_data_imputed1)

# Impute missing data in the test set using the predictorMatrix and varTypes from the train set
imputed_data_test <- miceRanger(
  testSet,
  predictorMatrix = imputed_data_train$predictorMatrix,
  varTypes = imputed_data_train$varTypes,
  num.trees = 100,
  verbose = FALSE
)

completed_test_data <- miceRanger::completeData(imputed_data_test)

test_data_imputed1 <- completed_test_data$Dataset_1
dim(test_data_imputed1)

colnames(train_data_imputed1)

### Make the "before and after imputation" plots
# sentinal_node_biopsy_positive
# sentinal_node_biopsy_positive in Training Set
# Recode sentinal_node_biopsy_positive in training and testing sets
trainSet$sentinal_node_biopsy_positive <- recode(trainSet$sentinal_node_biopsy_positive, `1` = "No", `2` = "Yes")
testSet$sentinal_node_biopsy_positive <- recode(testSet$sentinal_node_biopsy_positive, `1` = "No", `2` = "Yes")
train_data_imputed1$sentinal_node_biopsy_positive <- recode(train_data_imputed1$sentinal_node_biopsy_positive, `1` = "No", `2` = "Yes")
test_data_imputed1$sentinal_node_biopsy_positive <- recode(test_data_imputed1$sentinal_node_biopsy_positive, `1` = "No", `2` = "Yes")

# sentinal_node_biopsy_positive in Training Set
before_imputation_train <- data.frame(
  sentinal_node_biopsy_positive = trainSet$sentinal_node_biopsy_positive,
  Imputation = "Before",
  DataSet = "Train"
)

after_imputation_train <- data.frame(
  sentinal_node_biopsy_positive = train_data_imputed1$sentinal_node_biopsy_positive,
  Imputation = "After",
  DataSet = "Train"
)

# sentinal_node_biopsy_positive in Test Set
before_imputation_test <- data.frame(
  sentinal_node_biopsy_positive = testSet$sentinal_node_biopsy_positive,
  Imputation = "Before",
  DataSet = "Test"
)

after_imputation_test <- data.frame(
  sentinal_node_biopsy_positive = test_data_imputed1$sentinal_node_biopsy_positive,
  Imputation = "After",
  DataSet = "Test"
)

combined_data <- bind_rows(before_imputation_train, after_imputation_train, before_imputation_test, after_imputation_test)

# Create the combined bar plot with legend
ggplot(combined_data, aes(x = sentinal_node_biopsy_positive, fill = Imputation)) +
  geom_bar(position = "dodge") +
  facet_wrap(~DataSet, scales = "free", ncol = 2) +
  labs(title = "Sentinal Node Biopsy Positive (Before and After Imputation)") +
  xlab("Sentinal Node Biopsy Positive") +
  ylab("Count") +
  scale_fill_manual(values = c("Before" = "steelblue", "After" = "darkorange"), name = "Imputation") +
  scale_x_discrete(limits = c("No", "Yes"))


# Her2 Status in Training Set
before_imputation_train <- data.frame(
  her2_status = as.character(trainSet$her2_status),
  Imputation = "Before",
  DataSet = "Train"
)

after_imputation_train <- data.frame(
  her2_status = as.character(train_data_imputed1$her2_status),
  Imputation = "After",
  DataSet = "Train"
)

# Her2 Status in Test Set
before_imputation_test <- data.frame(
  her2_status = as.character(testSet$her2_status),
  Imputation = "Before",
  DataSet = "Test"
)

after_imputation_test <- data.frame(
  her2_status = as.character(test_data_imputed1$her2_status),
  Imputation = "After",
  DataSet = "Test"
)

combined_data <- bind_rows(before_imputation_train, after_imputation_train, before_imputation_test, after_imputation_test)

# Create the combined bar plot with legend
ggplot(combined_data, aes(x = her2_status, fill = Imputation)) +
  geom_bar(position = "dodge") +
  facet_wrap(~DataSet, scales = "free", ncol = 2) +
  labs(title = "Her2 Status (Before and After Imputation)") +
  xlab("Her2 Status") +
  ylab("Count") +
  scale_fill_manual(values = c("Before" = "steelblue", "After" = "darkorange"), name = "Imputation") +
  scale_x_discrete(limits = c("0", "1"))

#### Similar pre- & post- for race_id

# race_id in Training Set
before_imputation_train <- data.frame(
  race_id = as.character(trainSet$race_id),
  Imputation = "Before",
  DataSet = "Train"
)

after_imputation_train <- data.frame(
  race_id = as.character(train_data_imputed1$race_id),
  Imputation = "After",
  DataSet = "Train"
)

# race_id in Test Set
before_imputation_test <- data.frame(
  race_id = as.character(testSet$race_id),
  Imputation = "Before",
  DataSet = "Test"
)

after_imputation_test <- data.frame(
  race_id = as.character(test_data_imputed1$race_id),
  Imputation = "After",
  DataSet = "Test"
)

combined_data <- bind_rows(before_imputation_train, after_imputation_train, before_imputation_test, after_imputation_test)
combined_data_filtered <- combined_data[!is.na(combined_data$race_id), ] # required fro removing the NAs from the barplot. 

# Create the combined bar plot with legend
ggplot(combined_data_filtered, aes(x = race_id, fill = Imputation)) +
  geom_bar(position = "dodge") +
  facet_wrap(~DataSet, scales = "free", ncol = 2) +
  labs(title = "Race ID (Before and After Imputation)") +
  xlab("Race ID") +
  ylab("Count") +
  scale_fill_manual(values = c("Before" = "steelblue", "After" = "darkorange"), name = "Imputation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Now for histological grade too 

# Create the table for histologic_grade
hist_grade_table <- table(trainSet$histologic_grade)

# histologic_grade in Training Set
before_imputation_train <- data.frame(
  histologic_grade = as.character(trainSet$histologic_grade),
  Imputation = "Before",
  DataSet = "Train"
)

after_imputation_train <- data.frame(
  histologic_grade = as.character(train_data_imputed1$histologic_grade),
  Imputation = "After",
  DataSet = "Train"
)

# histologic_grade in Test Set
before_imputation_test <- data.frame(
  histologic_grade = as.character(testSet$histologic_grade),
  Imputation = "Before",
  DataSet = "Test"
)

after_imputation_test <- data.frame(
  histologic_grade = as.character(test_data_imputed1$histologic_grade),
  Imputation = "After",
  DataSet = "Test"
)

combined_data <- bind_rows(before_imputation_train, after_imputation_train, before_imputation_test, after_imputation_test)
combined_data_filtered <- combined_data[!is.na(combined_data$histologic_grade), ]

# Create the combined bar plot with legend
ggplot(combined_data_filtered, aes(x = histologic_grade, fill = Imputation)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ DataSet, scales = "free", ncol = 2) +
  labs(title = "Histologic Grade (Before and After Imputation)") +
  xlab("Histologic Grade") +
  ylab("Count") +
  scale_fill_manual(values = c("Before" = "steelblue", "After" = "darkorange"), name = "Imputation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





train_data_imputed1 <- as.data.frame(train_data_imputed1)
test_data_imputed1 <- as.data.frame(test_data_imputed1)

##### Descriptive statistics
library(tableone)
table_one <- CreateTableOne(data =train_data_imputed1)

print(table_one)


train_data_imputed1%>%
  select(-mask_id)%>%
  sumtable()



# install.packages("survival")
# install.packages("survminer")

library(survival)
library(survminer)


# there is a significant difference between thh treatments in the test set
# But maybe I shouldn't be 

### ONE_HOT_ENCODING 
library(fastDummies)
## one hot encoding  of training data

train_data_imputed1 <- as.data.frame(train_data_imputed1)
str(train_data_imputed1)

# Store columns to be excluded
cols_to_exclude <- c("mask_id", "disease_free_survival_status", "disease_free_survival_months")

# Create a vector of the columns to encode
cols_to_encode <- setdiff(names(train_data_imputed1), cols_to_exclude)

# Use the dummy_cols function to create the dummy variables
data_encoded_train <- dummy_cols(train_data_imputed1, select_columns = cols_to_encode, remove_first_dummy = TRUE)



# Identify the original columns that have been one-hot encoded
original_cols_to_drop <- setdiff(cols_to_encode, cols_to_exclude)

# Drop the original columns from the dataframe
data_encoded_train <- data_encoded_train[ , !(names(data_encoded_train) %in% original_cols_to_drop)]

colnames(data_encoded_train)


## one hot encoding  of test data

test_data_imputed1 <- as.data.frame(test_data_imputed1)
str(test_data_imputed1)

# Store columns to be excluded
cols_to_exclude <- c("mask_id", "disease_free_survival_status", "disease_free_survival_months")

# Create a vector of the columns to encode
cols_to_encode <- setdiff(names(test_data_imputed1), cols_to_exclude)

# Use the dummy_cols function to create the dummy variables
data_encoded_test <- dummy_cols(test_data_imputed1, select_columns = cols_to_encode, remove_first_dummy = TRUE)


# Identify the original columns that have been one-hot encoded
original_cols_to_drop <- setdiff(cols_to_encode, cols_to_exclude)

# Drop the original columns from the dataframe
data_encoded_test <- data_encoded_test[ , !(names(data_encoded_test) %in% original_cols_to_drop)]

colnames(data_encoded_test)
dim(data_encoded_test)
dim(data_encoded_train)
#######

## Renaming imputed, encoded datasets for saving and uploading back to python 

breast_train_data_imputed1_with_sentinal_node_data <- data_encoded_train
breast_test_data_imputed1_with_sentinal_node_data <- data_encoded_test

# Save train_data_imputed1
write.csv(breast_train_data_imputed1_with_sentinal_node_data, file = "/Users/Sam/Desktop/Imperial/Term_3/Week_7/Breast_imputed_data/breast_train_data_imputed1_with_sentinal_node_data.csv", row.names = FALSE)

# Save test_data_imputed1
write.csv(breast_test_data_imputed1_with_sentinal_node_data, file = "/Users/Sam/Desktop/Imperial/Term_3/Week_7/Breast_imputed_data/breast_test_data_imputed1_with_sentinal_node_data.csv", row.names = FALSE)

### Making a correlation heatmap of the variables. 

# Exclude columns from the correlation analysis
columns_to_exclude <- c("mask_id", "disease_free_survival_status", "disease_free_survival_months")
subset_data <- breast_train_data_imputed1_with_sentinal_node_data[, !colnames(breast_train_data_imputed1_with_sentinal_node_data) %in% columns_to_exclude]

# Calculate the correlation matrix
correlation_matrix <- cor(subset_data)

# Create the correlation heatmap
heatmap(correlation_matrix, 
        col = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Correlation Heatmap Predictor Variables")



library(gplots)

# install the gplots package if not already installed
if (!require(gplots)) {
  install.packages("gplots")
}

library(gplots)

# generate the heatmap with gplots' heatmap.2 function
heatmap.2(correlation_matrix,
          col = colorRampPalette(c("blue", "white", "red"))(100),
          xlab = "",  # remove x-axis labels
          dendrogram = "none",  # remove dendrograms
          main = "Correlation Heatmap Predictor Variables",
          trace = "none",  # remove trace lines
          key = TRUE,  # include a color scale
          keysize = 1.5,  # adjust the size of the color key
          density.info = "none",  # turn off density plot inside color key
          margins = c(5, 5)  # adjust margins
)

library(pheatmap)
pheatmap(correlation_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Correlation Heatmap Predictor Variables",
         cluster_rows = FALSE,  # remove row dendrogram
         cluster_cols = FALSE,  # remove column dendrogram
         display_numbers = FALSE,  # don't display numbers inside cells
         legend = TRUE  # include a color scale
)












