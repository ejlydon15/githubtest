#Poisson LASSO predictions for 2021 and Poisson GBM for 2021




    

   



library(car)
library(MASS)
library(leaps)
library(ggplot2)
library(olsrr)
library(openxlsx)
library(caret)
library(coefplot)
library(randomForest)
library(vip)
library(cowplot)



#Date of data in yyyy-mm-dd format
datadate<- '2020-12-17'
#moddate<- '2020.11.18RiskProfile'
moddate<- '2020.12.20'
i<- 'KnownPAS'

#writing l2 loss and abs loss functions
l2<- function(x,y){
  se<- (x-y)^2
  # sse<- sum(se)
  l2<- mean(se)
  return(l2)
}

l1<- function(x,y){
  se<- abs(x-y)
  # sse<- sum(se)
  l1<- mean(se)
  return(l1)
}

pl1<- function(x,y){
  se<- (abs(x-y))
  # sse<- sum(se)
  l1<- mean(se)
  pl1<- l1/mean(x)
  return(pl1)
}


#Read in these two for Known PAS Data
raw.2018<- read.csv(paste("/Users/elydon/Desktop/IRS-Filing Season/Data/Analysis Sets/", datadate, "KnownPASAnalysisSet2018.csv", sep=''),
                    header=T)
raw.2019<- read.csv(paste("/Users/elydon/Desktop/IRS-Filing Season/Data/Analysis Sets/", datadate, "KnownPASAnalysisSet2019.csv", sep=''),
                    header=T)

#raw.2020<- read.csv(paste("/Users/elydon/Desktop/IRS-Filing Season/Data/Analysis Sets/", datadate, "KnownPASAnalysisSet2020.csv", sep=''),
#                    header=T)
raw.2020<- read.csv(paste("/Users/elydon/Desktop/IRS-Filing Season/Data/Analysis Sets/", '2020-12-20', "KnownPASAnalysisSet2020.csv", sep=''),
                    header=T)

#Cleaning variable names
names(raw.2018)[!(names(raw.2018) %in% c('Mapped.System'))]<- substr(names(raw.2018)[!(names(raw.2018) %in% c('Mapped.System'))], 2, 50)
names(raw.2019)[!(names(raw.2019) %in% c('Mapped.System'))]<- substr(names(raw.2019)[!(names(raw.2019) %in% c('Mapped.System'))], 2, 50)
names(raw.2020)[!(names(raw.2020) %in% c('Mapped.System'))]<- substr(names(raw.2020)[!(names(raw.2020) %in% c('Mapped.System'))], 2, 50)  

#Creating working sets then stripping their variable names
work2018<- raw.2018[order(raw.2018$`2018_fs_total_production_issues`, decreasing=T), ]
work2019<- raw.2019[order(raw.2019$`2019_fs_total_production_issues`, decreasing= T), ]
work2020<- raw.2020[order(raw.2020$`2020_fs_total_production_issues`, decreasing= T), ]

names(work2018)[!(names(work2018) %in% c('Mapped.System'))]<- substr(names(work2018)[!(names(work2018) %in% c('Mapped.System'))], 6, 50)
names(work2019)[!(names(work2019) %in% c('Mapped.System'))]<- substr(names(work2019)[!(names(work2019) %in% c('Mapped.System'))], 6, 50)
names(work2020)[!(names(work2020) %in% c('Mapped.System'))]<- substr(names(work2020)[!(names(work2020) %in% c('Mapped.System'))], 6, 50) 

#Grabbing only systems on the PSL list 
top2020<- work2020$Mapped.System[1:35]

set4<- c("septdec_test_incidents",
         "cy_ad_open_by_bsm", "cy_ad_open_by_seid",
         "cy_eops_open_by_bsm", "cy_eops_open_by_seid",
         "cy_uns_open_by_bsm", "cy_uns_open_by_seid",
         "cy_ad_tech_srts", "cy_eops_tech_srts",
         "cy_med_time_to_resolve", "fs_app_changes_icads",
         "fs_cyber_icads",
         "fs_data_extract_icads","fs_dbservices_icads",
         "fs_engineernet_icads",  "fs_infra_icads",
         "fs_new_ad_icads", "fs_penalty_int_icads",
         "fs_server_enviro_icads", "fs_software_lic_icads",
         'octdec_prod_normal_changes', 'octdec_prod_emergency_changes',
         'octdec_nonprod_normal_changes', 'octdec_nonprod_emergency_changes')

set<- set4

mods<- list()

#using 2020 regressed on 2019 for training 
train<- work2019[work2019$Mapped.System %in% top2020, ]
train<- train[order(train$Mapped.System), ]
out1<- work2020[(work2020$Mapped.System %in% train$Mapped.System), c('Mapped.System', 'fs_total_production_issues')] #2020 incidents
out1<- out1[order(out1$Mapped.System), ]

actualvsfit2020.top30<- data.frame(Mapped.System= out1$Mapped.System, 
                                   Actual2020Incidents= out1$fs_total_production_issues)

out1<- out1$fs_total_production_issues




scales<- build_scales(train[ ,set])
covars<- data.frame(fastScale(dataSet = train[ ,set], scales=scales))
rid<- which(colSums(is.na(covars)) == nrow(covars))
if(!is.na(rid[1])){covars<- covars[, -rid]}



set.seed(13)
cvlasso<- cv.glmnet(as.matrix(covars), out1,
                    type.measure = 'mse', 
                    family='poisson', 
                    lambda= c(.001, .005, .01, .05, .1, .5, 1, 5, 10), 
                    maxit=50000,
                    #intercept = F, 
                    standardize = T)
lasso<- glmnet(as.matrix(covars), out1,
               family='poisson', 
               lambda= cvlasso$lambda.min,
               standardize= T)

coeffs<- data.frame(Vars= as.vector(row.names(coef(lasso))[c(1,which(lasso$beta!=0)+1)]), Values=as.vector(c(round(lasso$a0, 5), round(as.vector(lasso$beta[lasso$beta!=0]), 5))))
coeffs$Vars[grep('icads', coeffs$Vars, invert = T)]<- paste('2019', coeffs$Vars[grep('icads', coeffs$Vars, invert = T)])
coeffs$Vars[grep('icads', coeffs$Vars, invert = F)]<- paste('2020', coeffs$Vars[grep('icads', coeffs$Vars, invert = F)])



mods[[paste('Set 4 Poisson Lasso Results (Outcome All FS 2020 Incidents)- Top 30 Data', i, 'Data', sep='')]]<- coeffs




#Running 2020 through to get 2021 Predictions
#now predicting 2021 with this model
#test<- work2020[work2020$Mapped.System %in% top2020, ]
testA<-work2020
test<- data.frame(fastScale(dataSet = testA[ ,set], scales=scales))
if(!is.na(rid[1])){test<- test[ ,!(names(test) %in% names(rid))]}

fit2021from2019<- data.frame(Mapped.System=testA$Mapped.System, Fitted=round(predict(lasso,  type='response', s='lambda.min', newx=as.matrix(test))))
names(fit2021from2019)[2]<- 'Fitted'


 

fit2021from2019$lb<- round(fit2021from2019$Fitted - 1.96*sd(fit2021from2019$Fitted)/sqrt(162))
fit2021from2019$ub<- round(fit2021from2019$Fitted + 1.96*sd(fit2021from2019$Fitted)/sqrt(162))
fit2021from2019<- fit2021from2019[order(fit2021from2019$Fitted, decreasing = T), ]


#PREDICTIONS
actual<- fit2021from2019[order(fit2021from2019$Fitted, decreasing = T), ]
actual<- actual[1:23, ]
actual$Mapped.System<- factor(actual$Mapped.System, levels=actual$Mapped.System)

g<- ggplot(data=actual, aes(x=Mapped.System, y=Fitted)) +
  geom_bar(stat='identity',  fill='gray',  alpha=.01) +
  geom_line(aes(x=as.numeric(Mapped.System), y=Fitted), col='blue', size=1.5) + 
  #geom_ribbon(data=actual, aes(x=as.numeric(Mapped.System), ymin=lb, ymax=ub), fill='gray', alpha=.6) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 40, hjust=1, size=24)) + ylab('Number of Incidents') + xlab('Mapped System') + 
  ggtitle('Predicted FS 2021 Incidents from Poisson LASSO Models Trained on 2019 Data') + ylim(c(0, 25)) + theme(text = element_text(size=28), axis.text.x = element_text(size=28))

jpeg(paste("/Users/elydon/Desktop/IRS-Filing Season/Output/Model Results/", moddate, "/PoissLassPred2021from2019.jpg", sep=""), width = 2200, height=1500)
print(g)
dev.off()

rm(actual, g )












rm(lasso, coeffs, scales, cvlasso)








