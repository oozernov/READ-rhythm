library(reshape)
library(dplyr)
library(plotrix)
library(gdata)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(psych)
library(DescTools)
library(corrplot)
library(ggm)
setwd ("~/Dropbox/Ola/Data Files")

####Import all the data

l_ILF=read.csv("l_ilf_fa2.txt.csv", header = T, na.strings = c("", "NA"))
r_ILF=read.csv("r_ilf_fa.csv", header = T, na.strings = c("", "NA"))
l_AF=read.csv("l_af_fa.csv", header = T, na.strings = c("", "NA"))
r_AF=read.csv("r_af_fa.csv", header = T, na.strings = c("", "NA"))
l_SLF=read.csv("l_slf_fa.csv", header = T, na.strings = c("", "NA"))
r_SLF=read.csv("r_slf_fa.csv", header = T, na.strings = c("", "NA"))
maj_cc=read.csv("f_maj_fa.csv",header = T, na.strings = c("", "NA"))
min_cc=read.csv("f_min_fa.csv",header = T, na.strings = c("", "NA"))
QC=read.csv("READ_DTIQC.csv", header = T, na.strings = c("", "NA"))
d=read.csv("Read_Final_data_Dec20_2016.csv")
d_rhythm=read.csv("~/Dropbox/Ola/READ_rhythm_data.csv", header = T, na.strings = c("", "NA")) #for behavioral


#Merge and filter
d_rhythm=merge(d_rhythm, QC, "READ", all=FALSE)
d_rhythm = d_rhythm %>% filter(d_rhythm$num_rem<10)  #for brain
#d_rhythm = d_rhythm %>% filter(d_rhythm$Wave<4) #for beh
d_rhythm$mean <-rowMeans(d_rhythm[3:4], na.rm=TRUE) #for beh
d_rhythm = d_rhythm %>% filter(d_rhythm$mean>0) # for beh


#####Behavioral analysis#######
d_rhythm$PA3<-(rowMeans(d_rhythm[c('CTELss_T4','CTBWss_T4')], na.rm=TRUE))
d_rhythm$RAN3<-(rowMeans(d_rhythm[c('RANOss_T4','RAN2ss_T4','RANLss_T4')], na.rm=TRUE))
d_rhythm$WM3<-(rowMeans(d_rhythm[c('CTNRss_T4','CTMDss_T4')], na.rm=TRUE))


#### correlation matrix ####
d_rhythm_beh = d_rhythm %>%dplyr::select(mean,PA3,RAN3,
                                  WM3,KBITss_T4,W3WIss_T4)
x = cor(d_rhythm_beh, use = "pairwise.complete.obs")
corrplot(x, method="circle", tl.cex = 0.5, type = "lower")
corrplot(cor(x[,1:6])[1:6,1, drop=FALSE], cl.pos='n')
pairs(d_rhythm_beh)

cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

res1 <- cor.mtest(d_rhythm_beh, 0.95)
## specialized the insignificant value according to the significant level
corrplot(x, p.mat = res1[[1]], sig.level=0.05, insig = "blank", method="circle", tl.cex = 0.5, type = "lower")

####Correlation Plot####
cor.test(d_rhythm$W3WIss_T4, d_rhythm$mean)

ggplot(d_rhythm, aes(x=W3WIss_T4, y=mean))+
  geom_point()+
  labs(x="Word ID (ss)", y="Rhythm Discrimination (d')")+
  stat_smooth(method="lm", col="red")

####Regression Models####
m1<-lm(W3WIss_T4~T2_AGEmos+GenderCoded+KBITss+PA3+mean,data=d_rhythm)
anova(m1)
eta_sq(m1)
calc.relimp(m1, type = c("lmg"),
            rela = FALSE)
m2<-boot.relimp(
  m1,
  b = 1000,
  typesel = c("lmg"),
  rank = TRUE,
  diff = TRUE,
  rela = FALSE
)
booteval.relimp(m2) # print result
results <- (booteval.relimp(m2, sort = TRUE)) # plot result
plot(results, level = 0.8, names.abbrev = 10,cex=1.5,
     main = "Relative Importance")

#####AFQ#####
#READ_Rhythm_Data <- rename(READ_Rhythm_Data, c(ID="READ"))
d_raf=merge(d_rhythm, r_AF, "READ",all=FALSE)
d_laf=merge(d_rhythm, l_AF, "READ",all=FALSE)
d_lilf=merge(d_rhythm, l_ILF, "READ",all=FALSE)
d_rilf=merge(d_rhythm, r_ILF, "READ",all=FALSE)
d_lslf=merge(d_rhythm, l_SLF, "READ",all=FALSE)
d_rslf=merge(d_rhythm, r_SLF, "READ",all=FALSE)
d_maj=merge(d_rhythm, maj_cc, "READ",all=FALSE)
d_min=merge(d_rhythm, min_cc, "READ",all=FALSE)
#####################Functions################################

AFQRhythmanalysis <- function(df, dfArr,Range){
  for (i in Range){ 
    fit <- lm(df[,c(i)]~GenderCoded,data=df,na.action=na.exclude)
    df$tmp <- resid(fit)
    norm=shapiro.test(df$tmp)
    if (norm$p.value <=0.05) {
      x <- cor.test(df$tmp, dfArr, use="pairwise", 
                    adjust="none",method="spearman", alpha=.05)
    }  else {
      x <- cor.test(df$tmp, dfArr, use="pairwise", 
                    adjust="none",method="pearson", alpha=.05) 
    }
    
    if (x[3] < 0.05) {  
      print(names(df[i]))
      print(x$p.value)
      print(x$estimate)
    }	
  }
}

#######################################

AFQRhythmanalysis(d_raf, d_raf$mean,310:359)
AFQRhythmanalysis(d_laf, d_laf$mean,310:359)
AFQRhythmanalysis(d_lslf, d_lslf$mean,310:359)
AFQRhythmanalysis(d_rslf, d_rslf$mean,310:359)
AFQRhythmanalysis(d_lilf, d_lilf$mean,310:359)
AFQRhythmanalysis(d_rilf, d_rilf$mean,310:359)
AFQRhythmanalysis(d_maj, d_maj$mean,310:359)
AFQRhythmanalysis(d_min, d_min$mean,310:359)


##################Linear Models##################################

## mean of significant nodes

#d_laf$meanFA <- rowMeans(d_laf[329:333], na.rm=TRUE) #mean of all significant columns
d_rslf$meanFA <- rowMeans(d_rslf[343:357], na.rm=TRUE) #mean of all significant columns


r_slf2<-d_rslf%>%dplyr::select('mean','Gender','T2_AGEmos','KBITss','CTELss','CTmeanSS0s','meanFA')
r_slf2<-na.omit(r_slf2)
lm=lm(mean~Gender+T2_AGEmos+KBITss+CTELss+meanFA, data=r_slf2) #SLF predicting rhythm
eta_sq(lm)
r_slf2$Gender<-as.factor(r_slf2$Gender)


calc.relimp(lm, type = c("lmg"),
            rela = TRUE)
boot<-boot.relimp(
  lm,
  b = 1000,
  typesel = c("lmg"),
  rank = TRUE,
  diff = TRUE,
  rela = TRUE
)
booteval.relimp(boot) # print result
results <- (booteval.relimp(boot, sort = TRUE)) # plot result
plot(results, level = 0.8, names.abbrev = 10,
     main="Relative Importance for Rhythm",xlim=c(0,6),)



#find what correlates with righ SLF

large_rslf<-merge(d,r_SLF)
large_rslf$meanFA <- rowMeans(large_rslf[141:155], na.rm=TRUE) #mean of all significant columns


for (i in 4:24){
  x<-cor.test(large_rslf[,c(i)], large_rslf$meanFA, method="pearson",use="pairwise")
    if (x$p.value < 0.05) {  
  print(names(large_rslf[i]))
  print(x$p.value)
  print(x$statistic)
}
}

cor.test(large_rslf$meanFA,large_rslf$CTELss)
cor.test(large_rslf$meanFA,large_rslf$CTMDss)

ggplot(d_rslf, aes(x=W3WIraw, y=meanFA))+
  geom_point()+
  labs(x="Phonological Awareness (ss)", y="Mean FA in Right SLF (nodes 34-50)")+
  stat_smooth(method="lm", col="red")

########Plots#####################
#save residuals

fit <- lm(meanFA~T2_AGEmos+GenderCoded,data=d_rslf,na.action=na.exclude)
d_rslf$resid <- resid(fit)


ggplot(d_rslf, aes(x=mean, y=resid))+
  geom_point()+
  labs(x="Rhythm Performance (d')", y="Mean FA in Right SLF (nodes 34-50)")+
  stat_smooth(method="lm", col="red")


######
