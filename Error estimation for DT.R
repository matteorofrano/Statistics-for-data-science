rm(list=ls()) 



library(rpart)
library(caret)
library(partykit)
library(reshape2)
library(ggplot2)
library(gridExtra)
set.seed(4)


#import diabets
df_diabets = read.csv("diabetes.csv", sep = ",", header = TRUE)

colnames(df_diabets)[colnames(df_diabets) == "Outcome"] ="Class"

df_diabets$Class = as.factor(df_diabets$Class)


Gen_error=function(tree, train){
  #Quantized Bayes Error
  leaf_assignments = predict(tree, train, type = "class")
  
  empirical_loss = sum(leaf_assignments != train$Class) / nrow(train)
  
  #SAMPLING ERROR 
  n_s = 0.5
  
  M = length(unique(train$Class))
  tree_party = as.party(tree)
  
  # Get the terminal node indices for each observation in the training set
  terminal_nodes_train = as.numeric(predict(tree_party, newdata = train, type = "node"))
  
  # Create a data frame combining terminal nodes and actual classes for the training set
  results_train = data.frame(Terminal_Node = terminal_nodes_train, Actual_Class = train$Class)
  
  # Count the occurrences of each class in each terminal node using table()
  #The with function is used to specify that the table should be created using
  #columns from the results_train data frame.
  class_counts_train = with(results_train, table(Terminal_Node, Actual_Class))
  
  # Convert class_counts_train to a data frame
  class_counts_df = as.data.frame(class_counts_train)
  
  # Pivot the data using dcast to get the desired matrix format
  k_i_y = dcast(data = class_counts_df, Terminal_Node ~ Actual_Class, value.var = "Freq")
  
  # Replace NA values with 0
  k_i_y[is.na(k_i_y)] = 0
  
  # Set the Terminal_Node column as row names
  rownames(k_i_y) = k_i_y$Terminal_Node
  k_i_y$Terminal_Node = NULL
  
  n_i =rowSums(k_i_y)
  
  n_i_y = (k_i_y+n_s)/(n_i+(M*n_s))
  
  VAR_n_i_y=(n_i_y*(1-n_i_y))/(1+n_i)
  
  BIAS_n_i_y=(n_i_y-(k_i_y/n_i))^2
  
  pmf=n_i/nrow(train)
  
  #generalization error estimation
  
  second_part=((VAR_n_i_y+BIAS_n_i_y)^0.5)*pmf
  
  sum=sum(second_part)
  
  GEN_ERROR= empirical_loss+sum
  return(GEN_ERROR)
}

CV<- function(dats, n.folds, depth_p){
  # This line creates an empty list called
  #folds to store the indices of the observations in each fold.
  folds = list() # flexible object for storing folds
  
  # This line calculates the size of each fold by dividing the 
  #number of rows in the dats data frame by the number of folds 
  #specified by the n.folds argument.
  fold.size = nrow(dats)/n.folds
  remain = 1:nrow(dats) # all obs are in
  
  # This line creates a matrix called 
  #results with n.folds rows and 2 columns, initialized with zeros.
  results = matrix(0, n.folds, 2)
  
  for (i in 1:n.folds){
    select = sample(remain, fold.size, replace = FALSE)
    #randomly sample “fold_size” from the ‘remaining observations’
    
    folds[[i]] = select # store indices
    
    
    if (i == n.folds){
      folds[[i]] = remain
    }
    
    #update remaining indices to reflect what was taken out
    remain = setdiff(remain, select)
    remain
  }
  
  for (i in 1:n.folds){
    # fold i
    indis = folds[[i]] #unpack into a vector
    train = dats[-indis, ] #split into train and test sets
    validation = dats[indis, ]
    
    tree = rpart(Class ~ ., data = train, method = "class", parms=list(split="gini"),
                 control = rpart.control(minsplit=2, minbucket=1, cp=0, maxdepth = depth_p)) 
    
    predicted = predict(tree, validation, type = "class")
    
    error = sum(predicted != validation$Class) / nrow(validation)
    depth = max(rpart:::tree.depth(as.numeric(rownames(tree$frame))))
    
    results[i, 1] = error
    results[i, 2] = depth
  }
  
  
  valori_buoni = results[, 1][which(results[, 2] == depth_p)]
  if (length(valori_buoni != 0)){
    mean_error = mean(valori_buoni)
  } else {
    mean_error = NA
  }
  
  return(mean_error)
}


##########################################################################
#DIABETES
dataframe = df_diabets
num_splits = 50
#THESE USED FOR THE DEPTHS PLOTS
jeff_results = list()
CV2_results=list()
CV5_results=list()
CV10_results=list()
test_results=list()


depths= seq.int(1, 25, by = 1)


r = c(0.1, 0.3, 0.5, 0.9)
plots_list = list()
for (prop in r) {
  for (depth in depths) {
    splits_jeff = list()
    splits_test=list()
    splits_CV2=list()
    splits_CV5=list()
    splits_CV10=list()
    for (i in 1:num_splits){
      # Generate random indices for the current split proportion
      num_samples = floor(prop * nrow(dataframe))
      indices = sample(1:nrow(dataframe), num_samples)
          
      # Create the training dataset for the current split
      train_data = dataframe[indices, ]
      test_data = dataframe[-indices, ]
          
     
      cv2_estimation = CV(train_data, 2, depth)
      cv5_estimation = CV(train_data, 5, depth)
      cv10_estimation = CV(train_data, 10, depth)
          
      splits_CV2[[i]]=list(estimation=cv2_estimation)
      splits_CV5[[i]]=list(estimation=cv5_estimation)
      splits_CV10[[i]]=list(estimation=cv10_estimation)
          
      # Perform classification using decision tree
      model = rpart(Class ~ ., data = train_data, method = "class", parms=list(split="gini"),
                    control = rpart.control(minsplit=2, minbucket=1, cp=0, maxdepth = depth))
          
    
          
    #?
      #predictions = vector("character", length(nrow(test_data)))    
      
      depth_model = max(rpart:::tree.depth(as.numeric(rownames(model$frame))))
      
      if (depth_model == depth){
        predictions=predict(model, newdata=test_data, type = "class")
        estimation=Gen_error(model, train_data)
      } else {
        predictions = NA
        estimation = NA
      }
      test_estimation = mean(predictions != test_data$Class)
      
      splits_jeff[[i]] = list(estimation=estimation)    
      splits_test[[i]] = list(estimation=test_estimation)
      }
    jeff_results[[as.character(depth)]] = splits_jeff
    test_results[[as.character(depth)]] = splits_test
    CV2_results[[as.character(depth)]]= splits_CV2
    CV5_results[[as.character(depth)]]= splits_CV5
    CV10_results[[as.character(depth)]]= splits_CV10
      
  }
  
  
  mean_performance_jeff = sapply(jeff_results, function(metrics) mean(sapply(metrics, function(split) split$estimation), na.rm = TRUE))
  mean_performance_test = sapply(test_results, function(metrics) mean(sapply(metrics, function(split) split$estimation), na.rm = TRUE))
  mean_performance_CV2 = sapply(CV2_results, function(metrics) mean(sapply(metrics, function(split) split$estimation), na.rm = TRUE))
  mean_performance_CV5 = sapply(CV5_results, function(metrics) mean(sapply(metrics, function(split) split$estimation), na.rm = TRUE))
  mean_performance_CV10 = sapply(CV10_results, function(metrics) mean(sapply(metrics, function(split) split$estimation), na.rm = TRUE))
  
  
  
  data <- data.frame(depths = depths,
                     Jeff = mean_performance_jeff,
                     Test = mean_performance_test,
                     CV2 = mean_performance_CV2,
                     CV5 = mean_performance_CV5,
                     CV10 = mean_performance_CV10)
  #remove NaN
  
  max_depth = max(max(data$depths[complete.cases(data$Jeff)]),
                   max(data$depths[complete.cases(data$Test)]),
                   max(data$depths[complete.cases(data$CV2)]),
                   max(data$depths[complete.cases(data$CV5)]),
                   max(data$depths[complete.cases(data$CV10)]))
  
  
  
  # Melt the data frame for easier plotting
  melted_data = melt(data, id.vars = "depths", variable.name = "Method")
  
  # Create the plot
  current_plot = ggplot(melted_data, aes(x = depths, y = value, color = Method)) +
    geom_line() +
    geom_point() +
    labs(x = "Tree depth", y = "Generalization Error", color = "Method") +
    ggtitle(paste("diabet (r:", prop, ")")) +
    scale_x_continuous(limits = c(1, max_depth)) +  # Set x-axis limits
    theme_minimal()
  
  plots_list[[as.character(prop)]] = current_plot
}

# Arrange the plots in a grid
grid.arrange(grobs = plots_list, ncol = 2)  # Adjust ncol as needed




#THESE USED FOR THE MSE/PEARSON PLOTS
r = c(0.1, 0.3, 0.5, 0.7, 0.9)
jeff_results2 = list()
CV2_results2=list()
CV5_results2=list()
CV10_results2=list()
test_results2=list()

pearson_jeff=c()
pearson_CV2=c()
pearson_CV5=c()
pearson_CV10=c()

for (prop in r) {
  splits_jeff_ = list()
  splits_CV2_=list()
  splits_CV5_=list()
  splits_CV10_=list()
  
  test_for_pearson=c()
  p_pearson_jeff=c()
  p_pearson_CV2=c()
  p_pearson_CV5=c()
  p_pearson_CV10=c()
  
  
  
  for (i in 1:num_splits) {
    # Generate random indices for the current split proportion
    num_samples = floor(prop * nrow(dataframe))
    indices = sample(1:nrow(dataframe), num_samples)
    
    # Create the training dataset for the current split
    train_data = dataframe[indices, ]
    test_data = dataframe[-indices, ]
    
    #model for ground truth
    model = rpart(Class ~ ., data = train_data, 
                  method = "class", parms=list(split="gini"))
    
    predictions=predict(model, newdata=test_data, type = "class")
    test_estimation = mean(predictions != test_data$Class)
    test_for_pearson=append(test_for_pearson, test_estimation)
    
    # Estimate cross-validation error 
    for (n in c(2,5,10)){
      train_ctrl = trainControl(method = "cv", number = n)  
      start_time_cv=Sys.time()
      cv_model = train(Class ~ ., data = train_data, method = "rpart",
                        trControl = train_ctrl, parms = list(split = "gini"))
      cv_estimation = 1-mean(cv_model$results$Accuracy)
      end_time_cv=Sys.time()
      cv_mse_obs=(test_estimation-cv_estimation)^2
      
      time_cv=end_time_cv-start_time_cv
      if (n==2){
        splits_CV2_[[i]]=list(estimation=cv_mse_obs, time=time_cv)
        p_pearson_CV2=append(p_pearson_CV2, cv_estimation)
      }else if(n==5){
        splits_CV5_[[i]]=list(estimation=cv_mse_obs, time=time_cv)
        p_pearson_CV5=append(p_pearson_CV5, cv_estimation)
      }else {
        splits_CV10_[[i]]=list(estimation=cv_mse_obs, time=time_cv)
        p_pearson_CV10=append(p_pearson_CV10, cv_estimation)}
    }
    
    
    #estimate generalization error
    start_time_jeff=Sys.time()
    estimation_jeff=Gen_error(model, train_data)
    end_time_jeff=Sys.time()
    
    p_pearson_jeff=append(p_pearson_jeff, estimation_jeff)
    mse_jeff_obs=(test_estimation-estimation_jeff)^2
    time_jeff= end_time_jeff - start_time_jeff
    
    splits_jeff_[[i]] = list(estimation=mse_jeff_obs, time=time_jeff)
    
  }
  jeff_results2[[as.character(prop)]] = splits_jeff_
  CV2_results2[[as.character(prop)]]= splits_CV2_
  CV5_results2[[as.character(prop)]]= splits_CV5_
  CV10_results2[[as.character(prop)]]= splits_CV10_
  
  
  cor_jeff=cor(test_for_pearson, p_pearson_jeff, method = "pearson")
  pearson_jeff=append(pearson_jeff, cor_jeff)
  cor_CV2=cor(test_for_pearson, p_pearson_CV2, method="pearson")
  pearson_CV2=append(pearson_CV2, cor_CV2)
  cor_CV5=cor(test_for_pearson, p_pearson_CV5, method="pearson")
  pearson_CV5=append(pearson_CV5, cor_CV5)
  cor_CV10=cor(test_for_pearson, p_pearson_CV10, method="pearson")
  pearson_CV10=append(pearson_CV10, cor_CV10)
  
}


# Function to extract MSE estimations for each proportion
calculate_average = function(nested_list) {
  estimations = sapply(nested_list, function(x) x$estimation)
  average = mean(estimations)
  return(average)
}

calculate_std = function(nested_list){
  estimations = sapply(nested_list, function(x) x$estimation)
  std=sd(estimations)
  return (std)
}

# Apply the function to each nested list within the main list
jeff = lapply(jeff_results2, calculate_average)
mse_jeff= unlist(jeff)

CV2 = lapply(CV2_results2, calculate_average)
mse_CV2=unlist(CV2)

CV5 = lapply(CV5_results2, calculate_average)
mse_CV5=unlist(CV5)

CV10 = lapply(CV10_results2, calculate_average)
mse_CV10 = unlist(CV10)

#MSE standard deviation
std_jeff=lapply(jeff_results2, calculate_std)
std_dev_mse_jeff=unlist(std_jeff)

std_CV2 = lapply(CV2_results2, calculate_std)
std_dev_mse_CV2=unlist(std_CV2)

std_CV5 = lapply(CV5_results2, calculate_std)
std_dev_mse_CV5=unlist(std_CV5)

std_CV10 = lapply(CV10_results2, calculate_std)
std_dev_mse_CV10=unlist(std_CV10)

# Create a data frame for the MSE plot
plot_data = data.frame(
  Proportion = r,
  MSE_Jeff = mse_jeff,
  SD_Jeff = std_dev_mse_jeff,
  MSE_CV2 = mse_CV2,
  SD_CV2 = std_dev_mse_CV2,
  MSE_CV5 = mse_CV5,
  SD_CV5 = std_dev_mse_CV5,
  MSE_CV10 = mse_CV10,
  SD_CV10 = std_dev_mse_CV10
)

# create the plot using ggplot2
plot0 = ggplot(plot_data, aes(x = Proportion)) +
  geom_line(aes(y = MSE_Jeff, color="MSE_Jeff"), linetype = "solid", size = 1) +
  geom_errorbar(aes(ymin = MSE_Jeff - SD_Jeff, ymax = MSE_Jeff + SD_Jeff, color="MSE_Jeff"),
                width = 0.1, size = 0.5) +
  geom_line(aes(y = MSE_CV2, color="MSE_CV2"), linetype = "solid", size = 1) +
  geom_errorbar(aes(ymin = MSE_CV2 - SD_CV2, ymax = MSE_CV2 + SD_CV2, color="MSE_CV2"),
                width = 0.1, size = 0.5) +
  geom_line(aes(y = MSE_CV5, color="MSE_CV5"), linetype = "solid", size = 1) +
  geom_errorbar(aes(ymin = MSE_CV5 - SD_CV5, ymax = MSE_CV5 + SD_CV5, color="MSE_CV5"),
                width = 0.1, size = 0.5) +
  geom_line(aes(y = MSE_CV10, color="MSE_CV10"), linetype = "solid", size = 1) +
  geom_errorbar(aes(ymin = MSE_CV10 - SD_CV10, ymax = MSE_CV10 + SD_CV10, color="MSE_CV10"),
                width = 0.1, size = 0.5) +
  scale_color_manual(values=c("MSE_Jeff"="blue", "MSE_CV2"="red", "MSE_CV5"="yellow", "MSE_CV10"="green")) +
  labs(title = "diabet mse",
       x = "#train/#data ratio",
       y = "MSE vs measured error") +
  scale_x_continuous(breaks = r) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title


print(plot0)

# Create a data frame for the pearson plot
plot_pearson = data.frame(
  Proportion = r,
  Correlation_Jeff = pearson_jeff,
  Correlation_CV2 = pearson_CV2,
  Correlation_CV5 = pearson_CV5,
  Correlation_CV10 = pearson_CV10
)

# create the plot 
plot2 = ggplot(plot_pearson, aes(x = Proportion)) +
  geom_line(aes(y = Correlation_Jeff, color="Correlation_Jeff"), linetype = "solid", size = 1) +
  geom_line(aes(y = Correlation_CV2, color="Correlation_CV2"), linetype = "solid", size = 1) +
  geom_line(aes(y = Correlation_CV5, color="Correlation_CV5"), linetype = "solid", size = 1) +
  geom_line(aes(y = Correlation_CV10, color="Correlation_CV10"), linetype = "solid", size = 1) +
  scale_color_manual(values=c("Correlation_Jeff"="blue", "Correlation_CV2"="red", "Correlation_CV5"="yellow", "Correlation_CV10"="green")) +
  labs(title = "diabet pearson",
       x = "#train/#data ratio",
       y = "pearson corr") +
  scale_x_continuous(breaks = r) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title


print(plot2)


#TIME COMPLEXITY
calculate_time = function(nested_list) {
  estimations = sapply(nested_list, function(x) x$time)
  time = mean(estimations)
  return(time)
}

Jeff_time = unlist(lapply(jeff_results2, calculate_time))
CV2_time = unlist(lapply(CV2_results2, calculate_time))
CV5_time = unlist(lapply(CV5_results2, calculate_time))
CV10_time = unlist(lapply(CV10_results2, calculate_time))





# Create a matrix for the proportions and times
time_matrix = rbind(Jeff_time, CV2_time, CV5_time, CV10_time)
colnames(time_matrix) = c("0.1", "0.3", "0.5", "0.7", "0.9")



# Create a bar plot
barplot(time_matrix, beside = TRUE, 
        legend.text = rownames(time_matrix),
        names.arg = colnames(time_matrix),
        args.legend = list(cex = 0.75, x="topright", inset = c(0, -0.35)))







######################################################################
#####################################################################
#NEW DF LETTER







df_letters = read.table("letter-recognition.data", sep = ",", header = TRUE)
dataframe = df_letters


num_splits = 50
#THESE USED FOR THE DEPTHS PLOTS
jeff_results = list()
CV2_results=list()
CV5_results=list()
CV10_results=list()
test_results=list()


depths= seq.int(1, 12, by = 1)


r = c(0.1, 0.3, 0.5, 0.9)
plots_list = list()
for (prop in r) {
  for (depth in depths) {
    splits_jeff = list()
    splits_test=list()
    splits_CV2=list()
    splits_CV5=list()
    splits_CV10=list()
    for (i in 1:num_splits) {
      # Generate random indices for the current split proportion
      num_samples = floor(prop * nrow(dataframe))
      indices = sample(1:nrow(dataframe), num_samples)
      
      # Create the training dataset for the current split
      train_data = dataframe[indices, ]
      test_data = dataframe[-indices, ]
      
      
      cv2_estimation = CV(train_data, 2, depth)
      cv5_estimation = CV(train_data, 5, depth)
      cv10_estimation = CV(train_data, 10, depth)
      
      splits_CV2[[i]]=list(estimation=cv2_estimation)
      splits_CV5[[i]]=list(estimation=cv5_estimation)
      splits_CV10[[i]]=list(estimation=cv10_estimation)
      
      # Perform classification using decision tree
      model = rpart(Class ~ ., data = train_data, method = "class", parms=list(split="gini"),
                    control = rpart.control(minsplit=2, minbucket=1, cp=0, maxdepth = depth))
      
      
      
      
      predictions = vector("character", length(nrow(test_data)))    
      
      depth_model = max(rpart:::tree.depth(as.numeric(rownames(model$frame))))
      
      if (depth_model == depth){
        predictions=predict(model, newdata=test_data, type = "class")
        estimation=Gen_error(model, train_data)
      } else {
        predictions = NA
        estimation = NA
      }
      test_estimation = mean(predictions != test_data$Class)
      
      splits_jeff[[i]] = list(estimation=estimation)    
      splits_test[[i]] = list(estimation=test_estimation)
    }
    jeff_results[[as.character(depth)]] = splits_jeff
    test_results[[as.character(depth)]] = splits_test
    CV2_results[[as.character(depth)]]= splits_CV2
    CV5_results[[as.character(depth)]]= splits_CV5
    CV10_results[[as.character(depth)]]= splits_CV10
    
  }
  
  
  mean_performance_jeff = sapply(jeff_results, function(metrics) mean(sapply(metrics, function(split) split$estimation), na.rm = TRUE))
  mean_performance_test = sapply(test_results, function(metrics) mean(sapply(metrics, function(split) split$estimation), na.rm = TRUE))
  mean_performance_CV2 = sapply(CV2_results, function(metrics) mean(sapply(metrics, function(split) split$estimation), na.rm = TRUE))
  mean_performance_CV5 = sapply(CV5_results, function(metrics) mean(sapply(metrics, function(split) split$estimation), na.rm = TRUE))
  mean_performance_CV10 = sapply(CV10_results, function(metrics) mean(sapply(metrics, function(split) split$estimation), na.rm = TRUE))
  
  
  
  data = data.frame(depths = depths,
                     Jeff = mean_performance_jeff,
                     Test = mean_performance_test,
                     CV2 = mean_performance_CV2,
                     CV5 = mean_performance_CV5,
                     CV10 = mean_performance_CV10)
  #remove NaN
  
  max_depth = max(max(data$depths[complete.cases(data$Jeff)]),
                   max(data$depths[complete.cases(data$Test)]),
                   max(data$depths[complete.cases(data$CV2)]),
                   max(data$depths[complete.cases(data$CV5)]),
                   max(data$depths[complete.cases(data$CV10)]))
  
  
  
  # Melt the data frame for easier plotting
  melted_data = melt(data, id.vars = "depths", variable.name = "Method")
  
  # Create the plot
  current_plot = ggplot(melted_data, aes(x = depths, y = value, color = Method)) +
    geom_line() +
    geom_point() +
    labs(x = "Tree depth", y = "Generalization Error", color = "Method") +
    ggtitle(paste("letter (r:", prop, ")")) +
    scale_x_continuous(limits = c(1, max_depth)) +  # Set x-axis limits
    theme_minimal()
  
  plots_list[[as.character(prop)]] = current_plot
}

# Arrange the plots in a grid
grid.arrange(grobs = plots_list, ncol = 2)  # Adjust ncol as needed




#THESE USED FOR THE MSE/PEARSON PLOTS
r = c(0.1, 0.3, 0.5, 0.7, 0.9)
jeff_results2 = list()
CV2_results2=list()
CV5_results2=list()
CV10_results2=list()
test_results2=list()

pearson_jeff=c()
pearson_CV2=c()
pearson_CV5=c()
pearson_CV10=c()

for (prop in r) {
  splits_jeff_ = list()
  splits_CV2_=list()
  splits_CV5_=list()
  splits_CV10_=list()
  
  test_for_pearson=c()
  p_pearson_jeff=c()
  p_pearson_CV2=c()
  p_pearson_CV5=c()
  p_pearson_CV10=c()
  
  
  for (i in 1:num_splits) {
    # Generate random indices for the current split proportion
    num_samples = floor(prop * nrow(dataframe))
    indices = sample(1:nrow(dataframe), num_samples)
    
    # Create the training dataset for the current split
    train_data = dataframe[indices, ]
    test_data = dataframe[-indices, ]
    
    #model for ground truth
    model = rpart(Class ~ ., data = train_data, 
                  method = "class", parms=list(split="gini"))
    
    predictions=predict(model, newdata=test_data, type = "class")
    test_estimation = mean(predictions != test_data$Class)
    test_for_pearson=append(test_for_pearson, test_estimation)
    
    # Estimate cross-validation error using trainControl and train from caret
    for (n in c(2,5,10)){
      
      train_ctrl = trainControl(method = "cv", number = n)  
      start_time_cv=Sys.time()
      cv_model = train(Class ~ ., data = train_data, method = "rpart",
                       trControl = train_ctrl, parms = list(split = "gini"))
      cv_estimation = 1-mean(cv_model$results$Accuracy)
      end_time_cv=Sys.time()
      cv_mse_obs=(test_estimation-cv_estimation)^2
      
      time_cv=end_time_cv-start_time_cv
      if (n==2){
        splits_CV2_[[i]]=list(estimation=cv_mse_obs, time=time_cv)
        p_pearson_CV2=append(p_pearson_CV2, cv_estimation)
      }else if(n==5){
        splits_CV5_[[i]]=list(estimation=cv_mse_obs, time=time_cv)
        p_pearson_CV5=append(p_pearson_CV5, cv_estimation)
      }else {
        splits_CV10_[[i]]=list(estimation=cv_mse_obs, time=time_cv)
        p_pearson_CV10=append(p_pearson_CV10, cv_estimation)}
    }
    
    
    #estimate generalization error
    start_time_jeff=Sys.time()
    estimation_jeff=Gen_error(model, train_data)
    end_time_jeff=Sys.time()
    
    p_pearson_jeff=append(p_pearson_jeff, estimation_jeff)
    mse_jeff_obs=(test_estimation-estimation_jeff)^2
    time_jeff= end_time_jeff - start_time_jeff
    
    splits_jeff_[[i]] = list(estimation=mse_jeff_obs, time=time_jeff)
    
  }
  jeff_results2[[as.character(prop)]] = splits_jeff_
  CV2_results2[[as.character(prop)]]= splits_CV2_
  CV5_results2[[as.character(prop)]]= splits_CV5_
  CV10_results2[[as.character(prop)]]= splits_CV10_
  
  
  cor_jeff=cor(test_for_pearson, p_pearson_jeff, method = "pearson")
  pearson_jeff=append(pearson_jeff, cor_jeff)
  cor_CV2=cor(test_for_pearson, p_pearson_CV2, method="pearson")
  pearson_CV2=append(pearson_CV2, cor_CV2)
  cor_CV5=cor(test_for_pearson, p_pearson_CV5, method="pearson")
  pearson_CV5=append(pearson_CV5, cor_CV5)
  cor_CV10=cor(test_for_pearson, p_pearson_CV10, method="pearson")
  pearson_CV10=append(pearson_CV10, cor_CV10)
  
}


# Function to extract MSE estimations for each proportion
calculate_average = function(nested_list) {
  estimations = sapply(nested_list, function(x) x$estimation)
  average = mean(estimations)
  return(average)
}

calculate_std = function(nested_list){
  estimations = sapply(nested_list, function(x) x$estimation)
  std=sd(estimations)
  return (std)
}

# Apply the function to each nested list within the main list
jeff = lapply(jeff_results2, calculate_average)
mse_jeff= unlist(jeff)

CV2 = lapply(CV2_results2, calculate_average)
mse_CV2=unlist(CV2)

CV5 = lapply(CV5_results2, calculate_average)
mse_CV5=unlist(CV5)

CV10 = lapply(CV10_results2, calculate_average)
mse_CV10 = unlist(CV10)

#MSE standard deviation
std_jeff=lapply(jeff_results2, calculate_std)
std_dev_mse_jeff=unlist(std_jeff)

std_CV2 = lapply(CV2_results2, calculate_std)
std_dev_mse_CV2=unlist(std_CV2)

std_CV5 = lapply(CV5_results2, calculate_std)
std_dev_mse_CV5=unlist(std_CV5)

std_CV10 = lapply(CV10_results2, calculate_std)
std_dev_mse_CV10=unlist(std_CV10)

# Create a data frame for the mse plot
plot_data = data.frame(
  Proportion = r,
  MSE_Jeff = mse_jeff,
  SD_Jeff = std_dev_mse_jeff,
  MSE_CV2 = mse_CV2,
  SD_CV2 = std_dev_mse_CV2,
  MSE_CV5 = mse_CV5,
  SD_CV5 = std_dev_mse_CV5,
  MSE_CV10 = mse_CV10,
  SD_CV10 = std_dev_mse_CV10
)

# Create the plot 
plot0 = ggplot(plot_data, aes(x = Proportion)) +
  geom_line(aes(y = MSE_Jeff, color="MSE_Jeff"), linetype = "solid", size = 1) +
  geom_errorbar(aes(ymin = MSE_Jeff - SD_Jeff, ymax = MSE_Jeff + SD_Jeff, color="MSE_Jeff"),
                width = 0.1, size = 0.5) +
  geom_line(aes(y = MSE_CV2, color="MSE_CV2"), linetype = "solid", size = 1) +
  geom_errorbar(aes(ymin = MSE_CV2 - SD_CV2, ymax = MSE_CV2 + SD_CV2, color="MSE_CV2"),
                width = 0.1, size = 0.5) +
  geom_line(aes(y = MSE_CV5, color="MSE_CV5"), linetype = "solid", size = 1) +
  geom_errorbar(aes(ymin = MSE_CV5 - SD_CV5, ymax = MSE_CV5 + SD_CV5, color="MSE_CV5"),
                width = 0.1, size = 0.5) +
  geom_line(aes(y = MSE_CV10, color="MSE_CV10"), linetype = "solid", size = 1) +
  geom_errorbar(aes(ymin = MSE_CV10 - SD_CV10, ymax = MSE_CV10 + SD_CV10, color="MSE_CV10"),
                width = 0.1, size = 0.5) +
  scale_color_manual(values=c("MSE_Jeff"="blue", "MSE_CV2"="red", "MSE_CV5"="yellow", "MSE_CV10"="green")) +
  labs(title = "letter mse",
       x = "#train/#data ratio",
       y = "MSE vs measured error") +
  scale_x_continuous(breaks = r) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title


print(plot0)

# Create a data frame for the pearson plot
plot_pearson = data.frame(
  Proportion = r,
  Correlation_Jeff = pearson_jeff,
  Correlation_CV2 = pearson_CV2,
  Correlation_CV5 = pearson_CV5,
  Correlation_CV10 = pearson_CV10
)

# Create the plot using ggplot2
plot2 = ggplot(plot_pearson, aes(x = Proportion)) +
  geom_line(aes(y = Correlation_Jeff, color="Correlation_Jeff"), linetype = "solid", size = 1) +
  geom_line(aes(y = Correlation_CV2, color="Correlation_CV2"), linetype = "solid", size = 1) +
  geom_line(aes(y = Correlation_CV5, color="Correlation_CV5"), linetype = "solid", size = 1) +
  geom_line(aes(y = Correlation_CV10, color="Correlation_CV10"), linetype = "solid", size = 1) +
  scale_color_manual(values=c("Correlation_Jeff"="blue", "Correlation_CV2"="red", "Correlation_CV5"="yellow", "Correlation_CV10"="green")) +
  labs(title = "letter pearson correlation",
       x = "#train/#data ratio",
       y = "pearson corr") +
  scale_x_continuous(breaks = r) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title


print(plot2)


#TIME COMPLEXITY
calculate_time = function(nested_list) {
  estimations = sapply(nested_list, function(x) x$time)
  time = mean(estimations)
  return(time)
}

Jeff_time = unlist(lapply(jeff_results2, calculate_time))
CV2_time = unlist(lapply(CV2_results2, calculate_time))
CV5_time = unlist(lapply(CV5_results2, calculate_time))
CV10_time = unlist(lapply(CV10_results2, calculate_time))





# Create a matrix for the proportions and times
time_matrix <- rbind(Jeff_time, CV2_time, CV5_time, CV10_time)
colnames(time_matrix) <- c("0.1", "0.3", "0.5", "0.7", "0.9")



# Create a bar plot
barplot(time_matrix, beside = TRUE, 
        legend.text = rownames(time_matrix),
        names.arg = colnames(time_matrix),
        args.legend = list(cex = 0.75, x="topright", inset = c(0, -0.35)))
