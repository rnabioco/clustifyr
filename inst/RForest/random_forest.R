
# install.packages('devtools')
devtools::install_github('NCBI-Hackathons/RClusterCT')
library(RClusterCT)
# Load data manually
# load(file = "/home/YueYvetteHao/data/pbmc4k_vargenes.rda")
# load(file = "/home/YueYvetteHao/data/pbmc4k_matrix.rda")
# load(file = "/home/YueYvetteHao/data/pbmc4k_meta.rda")

# install and load packages
# install.packages("randomForest")
# install.packages("caret")
# install.packages("mlbench")
# install.packages("e1071")
# install.packages("rpart")
library(randomForest)
library(caret)
library(mlbench)
library(e1071)
# library(rpart)

source('reduce_matrix.R')

# Input for random forest
pbmc4 <- reduce_matrix(pbmc4k_matrix, pbmc4k_meta, pbmc4k_vargenes)
pbmc5 <- reduce_matrix(pbmc5_matrix, pbmc5_meta, pbmc4k_vargenes)

names(pbmc5)[names(pbmc5) == "cluster"] <- "classified"

pbmc4.sub <- pbmc4[, colnames(pbmc4) %in% colnames(pbmc5)]
pbmc5.sub <- pbmc5[, colnames(pbmc5) %in% colnames(pbmc4.sub)]


seed <- 100
# Algorithm Tune (tuneRF)
set.seed(seed)
# mtry <- sqrt(ncol(pbmc4))
x <- pbmc4[,2:ncol(pbmc4)]
y <- pbmc4[,1]
bestmtry <- tuneRF(x, y, stepFactor=1.5, improve=1e-5, ntree=500)
print(bestmtry)
# bestmtry is 256 when ntree is 500

# Tune model using caret
# Custom tuning
# https://machinelearningmastery.com/tune-machine-learning-algorithms-in-r/
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes


metric <- "Accuracy"
control <- trainControl(method="repeatedcv", number=10, repeats=3)
tunegrid <- expand.grid(.mtry=c(250:260), .ntree=c(500, 1000, 1500, 2000))
set.seed(seed)
rf_tune <- train(classified~., data=pbmc4, method=customRF, metric=metric, tuneGrid=tunegrid, trControl=control)
summary(rf_tune)
plot(rf_tune)


# Optimal parameters
bmtry<- 256
bntree<- 1000
#bntree<-1500
#bntree<-2000

# Build the model
rFmodel <- randomForest(classified ~ ., data = pbmc4, ntree = bntree, mtry = bmtry, importance = TRUE)
rFmodel.sub <- randomForest(classified ~ ., data = pbmc4.sub, ntree = bntree, mtry = bmtry, importance = TRUE)

save(rFmodel.sub,file = "rFmodelsub.RData")

# To check important variables
importance(rFmodel)        
varImpPlot(rFmodel) 
save(rFmodel,file = "rFmodel.RData")
load("rFmodel.RData")

# Fit the model with new data
pred <- predict(rFmodel.sub, pbmc5.sub, type = "class")
# Checking classification accuracy
table(pred, pbmc5.sub$classified) 
