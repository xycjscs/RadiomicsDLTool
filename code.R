library(tidyverse)
library(caret)
library(pROC)
library(glmnet)
# library(DMwR)
library(rmda)
library(ggpubr)
library(ModelGood)
library(rms)
library(mRMRe)
library(DescTools)
library(Publish)
library(pheatmap)


rm(list=ls())
setwd("C:/Users/212768837/Downloads/demo code (1)/Radiomics")

# fpath <- choose.files()
# data <- read.csv(fpath)

source('assist_file.R')

set.seed(123456)

clinics_column <- 2

fpath <- choose.files()

dt_all <- read.csv(fpath)

dt_all$Label <- factor(dt_all$Label, ordered = T)

dt <- dt_all[, -c(2:clinics_column)]

if(!is_empty(nearZeroVar(dt)))
{
  dt <- dt[, -nearZeroVar(dt)]
}


# clinical data

dt_cli <- dt_all[, c(1:clinics_column)]

acc_val_train <- numeric(0)
acc_val_test <- numeric(0)

all_idx <- list()
train_list <- list()
test_list <- list()
cvfit_list <- list()
fit_list <- list()
out_list <- list()

for(numIt  in c(1:10))
{
  s_pre <- preProcess(dt, method = 'medianImpute')
  dt <- predict(s_pre, dt)
  idx_train <- createDataPartition(dt$Label, p = 0.7, list = F)
  dt_train <- dt[idx_train, ]
  dt_test <- dt[-idx_train, ]
  
  # browser()
  
  all_idx[[numIt]] <- idx_train
  step_pre <- preProcess(dt_train, method = c('center', 'scale'))
  dt_train_pre <- predict(step_pre, dt_train) %>% as.data.frame
  
  dt_train_pre0 <- dt_train_pre
  # dt_train_pre <- SMOTE(Label~., data = dt_train_pre)
  dt_test_pre <- predict(step_pre, dt_test)
  
  dt_mrmr <- mRMR.data(dt_train_pre)
  f_sel <- mRMR.classic(data = dt_mrmr, target_indices = c(1), feature_count = 20)
  
  # browser()
  
  dt_train_pre <- select(dt_train_pre, c('Label', featureNames(f_sel)[unlist(solutions(f_sel))]))
  dt_train_pre0 <- select(dt_train_pre0, c('Label', featureNames(f_sel)[unlist(solutions(f_sel))]))
  dt_test_pre <- select(dt_test_pre, c('Label', featureNames(f_sel)[unlist(solutions(f_sel))]))
  
  x <- as.matrix(dt_train_pre[, -1])
  y <- dt_train_pre$Label
  
  cv.fit <- cv.glmnet(x, y, family = 'binomial')
  fit <- glmnet(x, y, family = 'binomial')
  
  
  train_list[[numIt]] <- dt_train_pre0
  test_list[[numIt]] <- dt_test_pre
  cvfit_list[[numIt]] <- cv.fit
  fit_list[[numIt]] <- fit
  
  # browser()
  pre_res_test <- as.vector(predict(fit, newx = as.matrix(dt_test_pre[, -1]), s = cv.fit$lambda.min))
  roc_res_test <- pROC::roc(dt_test_pre$Label, pre_res_test)
  
  
  
  pre_res_train <- as.vector(predict(fit, newx = x, s = cv.fit$lambda.min))
  roc_res_train <- pROC::roc(dt_train_pre$Label, pre_res_train)
  
  
  dir_sign_test <- roc_res_test$direction
  dir_sign_train <- roc_res_train$direction
  
  if(dir_sign_test == dir_sign_train)
  {
    acc_val_test <- c(acc_val_test, pROC::auc(roc_res_test))
    acc_val_train <- c(acc_val_train, pROC::auc(roc_res_train))
  }
  else
  {
    acc_val_test <- c(acc_val_test, 0)
    acc_val_train <- c(acc_val_train, 0)
  }
}

idx_vec <- c(1:length(acc_val_test))
idx_vec <- idx_vec[acc_val_train > acc_val_test]
acc_val <- acc_val_test[acc_val_train > acc_val_test]
init_idx <- which.max(acc_val)
sel_idx <- idx_vec[init_idx]

idx_train <- all_idx[[sel_idx]]
grp_info <- tibble(Label = dt$Label, Group = 'Test')
grp_info$Radscore <- 0
grp_info$Group[idx_train] <- 'Training'


dt_train_final <- train_list[[sel_idx]]
dt_test_final <- test_list[[sel_idx]]
cvfit <- cvfit_list[[sel_idx]]
fit <- fit_list[[sel_idx]]

s = cvfit$lambda.min


pre_res_test <- as.vector(predict(fit, newx = as.matrix(dt_test_final[, -1]), s = s))
pre_res_test_prob <- as.vector(predict(fit, newx = as.matrix(dt_test_final[, -1]), s = s, 
                                       type = 'link'))
roc_res_test <- pROC::roc(dt_test_final$Label, pre_res_test, ci = T)

out_res_test <- ifelse(pre_res_test > coords(roc_res_test, x = 'best')[[1]], 1, 0)
conf_mat_test <- confusionMatrix(as.factor(out_res_test),as.factor(dt_test_final$Label))
rec_test <- c(conf_mat_test$overall[c(1, 3, 4)], conf_mat_test$byClass[c(1:4)])

pre_res_train <- as.vector(predict(fit, newx = as.matrix(dt_train_final[, -1]), s = s))
pre_res_train_prob <- as.vector(predict(fit, newx = as.matrix(dt_train_final[, -1]), s = s, 
                                        type = 'link'))
roc_res_train <- pROC::roc(dt_train_final$Label, pre_res_train, ci = T)

out_res_train <- ifelse(pre_res_train > coords(roc_res_train, x = 'best')[[1]], 1, 0)
conf_mat_train <- confusionMatrix(as.factor(out_res_train), as.factor(dt_train_final$Label))
rec_train <- c(conf_mat_train$overall[c(1, 3, 4)], conf_mat_train$byClass[c(1:4)])

rec_rad <- data.frame(rbind(rec_train, rec_test), row.names = c('Train', 'Test'))

write.csv(rec_rad, file = 'res_radiomics.csv')

grp_info$Radscore[idx_train] <- pre_res_train
grp_info$Radscore[-idx_train] <- pre_res_test

write_csv(grp_info, 'group_info.csv')

cutoff_radiomics <- coords(roc_res_train, x = 'best')

## rad score
dt_final_test <- tibble(Label = dt_test_final$Label, rad_score = pre_res_test)
dt_final_arr <- arrange(dt_final_test, rad_score)
dt_final_arr$x <- 1:nrow(dt_final_arr)

dt_final_train <- tibble(Label = dt_train_final$Label, rad_score = pre_res_train)

p_train <- ggboxplot(x = 'Label', y = 'rad_score', data = dt_final_train,
                     add = 'jitter', color = 'Label', palette = 'jco') + 
  ylim(-3, 3) + 
  stat_compare_means(method = 'wilcox.test') + 
  geom_hline(yintercept = coords(roc_res_train, x = 'best')[[1]]) + theme_bw()
p_test <- ggboxplot(x = 'Label', y = 'rad_score', data = dt_final_test,
                    add = 'jitter', color = 'Label', palette = 'jco') +
  ylim(-3, 3) +
  stat_compare_means(method = 'wilcox.test') + 
  geom_hline(yintercept = coords(roc_res_train, x = 'best')[[1]]) + theme_bw()




coefs <- coefficients(fit, s = s)
useful_feature <- unlist(coefs@Dimnames)[coefs@i + 1]
useful_feature <- useful_feature[-1]

dt_coef <- data.frame(Feature = useful_feature, Coef = coefs@x[-1])
dt_coef <- arrange(dt_coef, desc(Coef))
dt_coef$Feature <- factor(dt_coef$Feature, 
                          levels = as.character(dt_coef$Feature))

p_coef <- ggplot(aes(x = Feature, y = Coef), data = dt_coef)
p_coef <- p_coef + geom_col(fill = 'blue', width = 0.7) + coord_flip() + 
  theme_bw() + ylab('Coefficients')

final_data_test <- add_column(dt_test_final, radscore = pre_res_test)
final_data_test <- select(final_data_test, c('Label', useful_feature, 'radscore'))
write_csv(final_data_test, path = 'dataset_test.csv')
final_data_train <- add_column(dt_train_final, radscore = pre_res_train)
final_data_train <- select(final_data_train, c('Label', useful_feature, 'radscore'))
write_csv(final_data_train, path = 'dataset_train.csv')



fit_train <- glm(Label~rad_score, data = dt_final_train, family = 'binomial')

dt_dca_train <- dt_final_train
dt_dca_train$Label <- as.numeric(dt_dca_train$Label) - 1
dca_curve <- decision_curve(Label~rad_score, data = dt_dca_train)

radscore <- paste('Radscore = ', paste(round(coefs@x[-1], 3), 
                                       useful_feature, sep = '*', collapse = '+'), '+', round(coefs@x[1], 3))
print(radscore)
write_file(radscore, path = 'radscore.txt')

# figure1
oldpar <- par(mfrow = c(2, 1))
plot(cvfit)
# figure2
plot(fit, s = s, xvar = 'lambda')
abline(v = log(cvfit$lambda.min), lty = 2)
par(oldpar)
# figure3

oldpar <- par(mfrow = c(1, 2))
plot(roc_res_train, print.auc = T, 
     print.auc.pattern = 'AUC: %.2f(%.2f-%.2f)', legacy.axes = T)

# figure4
plot(roc_res_test, print.auc = T, legacy.axes = T,
     print.auc.pattern = 'AUC: %.2f(%.2f-%.2f)')

par(oldpar)
# figure5
ggarrange(p_train, p_test, ncol = 2)



#figure 6
p_coef + theme(axis.text.y = element_text(size = 12))

