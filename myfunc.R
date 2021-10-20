cv_penLogistic <- function(y,X,method=c("vanilla",
                                        "fowardAIC",
                                        "forwardBIC",
                                        "stepwiseAIC",
                                        "stepwiseBIC",
                                        "ridge",
                                        "lasso",
                                        "firth"),
                           folds=10,
                           repeats=20,
                           seed=1)
{
  set.seed(seed)
  n <- nrow(X)
  
  error_mat <- matrix(NA,n,repeats)
  dat <- data.frame(y=y,X=X)
  # repeats是最大的循环, 每个repeat都要执行k-fold CV
  for (r in 1:repeats) 
  {
    sets <- sample(rep(1:folds,n)[1:n],n)
    # 划分 训练集和 测试集
    for(i in 1:folds){
      testSet  <- which(sets==i)
      trainSet <- (1:n)[-testSet] 
      testData    <- dat[testSet, ]
      trainData   <- dat[trainSet, ]
      # 针对各种方法设置相应的训练模型:
      if (method=="vanilla") {
        model       <- glm(y~.,family=binomial,data=trainData)
      }
      
      if (method%in%c("fowardAIC","forwardBIC","stepwiseAIC","stepwiseBIC")) 
      {
        null = glm(y~1,data=trainData,family=binomial)
        full = glm(y~.,data=trainData,family=binomial)
      }
      
      if (method=="fowardAIC") {
        model <- step(null,scope=list(lower=null,upper=full),k=2,trace=0)
      }
      
      if (method=="forwardBIC") {
        model <- step(null,scope=list(lower=null,upper=full),k=log(n),trace=0)
      }     
      
      if (method=="stepwiseAIC") {
        model <- step(full,k=2,trace=0)
      }           
      
      if (method=="stepwiseBIC") {
        model <- step(full,k=log(n),trace=0)
      }
      
      if (method=="ridge") {
        
        model <- cv.glmnet(data.matrix(trainData[,-1]), trainData$y, 
                           alpha = 0, family = "binomial")
      }
      
      if (method=="lasso") {
        
        model <- cv.glmnet(data.matrix(trainData[,-1]), trainData$y, 
                           alpha = 1, family = "binomial")
      }
      
      if (method=="firth") {
        model <- brglm(y~., data=trainData)
      }      
      # 针对各个方法进行test
      if (!(method%in%c("ridge","lasso"))) {
        testProb <- predict(model, testData, type="response")
      }
      
      if (method%in%c("ridge","lasso") ) {
        res <- predict(model, newx=data.matrix(testData[,-1]), s="lambda.1se")
        testProb <- 1/(1 + exp(-res))
      }
      # 记录test的误差
      testPred <- round(testProb)
      errs <- as.numeric(testData$y!=testPred)
      error_mat[testSet,r] <- errs
    }
  }
  result<-list(error_mat, model)
  return(result)
}