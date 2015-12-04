# == PACKAGES ==============================================
#
# neural net package
#if (! require(neuralnet, quietly=TRUE)) {
#  install.packages("neuralnet",
#                   repos="http://cran.us.r-project.org")
#  library(neuralnet)
#}
# deep neural net with H2O
if (! require(h2o, quietly=TRUE)) {
  install.packages("h2o",
                   repos="http://cran.us.r-project.org")
  library(h2o)
}
# deep neural net with deepnet
if (! require(deepnet, quietly=TRUE)) {
  install.packages("deepnet",
                   repos="http://cran.us.r-project.org")
  library(deepnet)
}

# ==========================================================

inputs <- as.matrix(read.csv("./playing_around/combiFeaturesAv.csv", header=TRUE, stringsAsFactors = FALSE)[,-1])
targets <- as.matrix(read.csv("./playing_around/drugSynergies.csv", header=TRUE, stringsAsFactors = FALSE)[,4])
targetRange <- max(targets) - min(targets)
ntargets <- targets + abs(min(targets))
ntargets <- targets/targetRange

# Splitting into training and x-validating set

idx <- sample(1:nrow(inputs), round(nrow(inputs)*3/4))

inputsTrain <- inputs[idx,]
targetsTrain <- ntargets[idx,]

inputsVal <- inputs[-idx,]
targetsVal <- ntargets[-idx,]

# ==========================================================
# H2O
# if you get a RCurl missing dependency problem, install it:
# sudo apt-get install libcurl4-openssl-dev
# set up a local cluster with 1GB RAM
localH2o = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE)
#dataTrain <- cbind(targetsTrain, inputsTrain)
#dataTrainH2o <- as.h2o(localH2o, destination_frame = "dataTrainH2o")
dataH2o <- h2o.importFile("./playing_around/dataH2O_temp.csv")
dataTrainH2o <- dataH2o[idx,]
dataXValH2o <- dataH2o[-idx,1:ncol(dataH2o)]
h2oModel <- h2o.deeplearning(x = 2: ncol(dataH2o),
                             y = 1,
                             training_frame = dataTrainH2o,
                             activation = "TanhWithDropout",
                             input_dropout_ratio = 0.2,
                             hidden_dropout_ratios = c(0.5),
                             hidden = c(50), 
                             epochs = 100)
h2oPredictions <- h2o.predict(h2oModel, dataXValH2o)

randForPredictions <- h2o.randomForest(x = 2:ncol(dataH2o),
                      y = 1,
                      training_frame = dataTrainH2o,
                      validation_frame = dataXValH2o,
                      sample_rate = 0.2,
                      max_depth = 100)



# ==========================================================
# deepnet

deepnetModel <- dbn.dnn.train(inputsTrain, targetsTrain,
                              hidden = c(10),
                              activationfun = "sigm",
                              learningrate = 0.01,
                              output = "sigm",
                              numepochs = 5,
                              momentum = 0,
                              batchsize = 30, # train networks with batches of this many samples, this param affects the training time and not so much the performance
                              )
predictions <- nn.predict(deepnetModel, inputsVal)
predictions <- predictions * targetRange - abs(min(targets))