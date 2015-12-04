## Start a local cluster with 1GB RAM (default)
library(h2o)
localH2O <- h2o.init(ip = "localhost", port = 54321, startH2O = TRUE)

## Start a local cluster with 2GB RAM
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, 
                    Xmx = '2g')

## Convert Breast Cancer into H2O
#dat <- BreastCancer[, -1]  # remove the ID column
#dat_h2o <- as.h2o(localH2O, dat, key = 'dat')

## Import MNIST CSV as H2O
#dat_h2o <- h2o.importFile(localH2O, path = ".../mnist_train.csv")