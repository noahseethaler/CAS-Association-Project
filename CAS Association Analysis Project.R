#BS858 Project
library(tidyverse)

data <- read.csv("C:/Users/noahs/OneDrive/Desktop/Fall 2024/BS858 Statistical Genetics/Project/project_sib_pheno_and_RV_data.csv")
view(data)
sum(is.na(data))

#1a Heritability
ldl = c(data$STLDL1,data$STLDL2)
meanldl = mean(ldl,na.rm=T)
sdldl = sd(ldl,na.rm=T)
N = nrow(data)

adjkid1 = data$STLDL1 - meanldl
adjkid2 = data$STLDL2 - meanldl
icc = sum(adjkid1*adjkid2)/sdldl^2/(N-1)

h2 = 2*icc

print(h2)

#1b Familial
#Recurrence Risk Ratio
summary(data$age2)

proband.affected <- data %>% filter(CAS1 == 2)
sib.affected <- proband.affected %>% filter (CAS2 ==2 )
b <- nrow(proband.affected)
a <- nrow(sib.affected)


#2 Allele Frequencies for the SNPs
project_freqx <- read.delim("C:/Users/noahs/OneDrive/Desktop/Fall 2024/BS858 Statistical Genetics/Plink/project_freqx.frqx")
#view(project_freqx)
data2 <- project_freqx %>% select (1:7) %>% mutate(A1.FREQ = ((C.HOM.A1.*2) + C.HET.) / (6000*2),
                                                                       A2.FREQ = ((C.HOM.A2.*2) + C.HET.) / (6000*2)) %>% view()


#3 HWE
project_hwe <- read.csv("C:/Users/noahs/OneDrive/Desktop/Fall 2024/BS858 Statistical Genetics/Plink/project_hwe.hwe", sep="")
data3 <- project_hwe
view(data3)
hwe <- data3 %>% filter(P < 0.05) %>% view()

#4 FWER
project_ld <- read.csv("C:/Users/noahs/OneDrive/Desktop/Fall 2024/BS858 Statistical Genetics/Plink/project_ld.ld", sep="")
data4 <- project_ld
view(data4)


#5 Power


#6
#SNP rs116418597 and rs61849569 


#7 

#8 CAS Association
project.logisitc.adj.assoc <- read.csv("C:/Users/noahs/OneDrive/Desktop/Fall 2024/BS858 Statistical Genetics/Plink/project.logisitc.adj.assoc.logistic", sep="")
data8 <- project.logisitc.adj.assoc
view(data8)
thresh8 <- data8b %>% filter(TEST == "ADD") %>% filter(P < 0.002) %>% view()

project.logisitc.adj5.assoc <- read.csv("C:/Users/noahs/OneDrive/Desktop/Fall 2024/BS858 Statistical Genetics/Plink/project.logisitc.adj5.assoc.logistic", sep="")
data8b <- project.logisitc.adj5.assoc
view(data8b)
data8b$beta <- log(data8b$OR)

data.compare <- data8b %>% filter(TEST == "ADD") %>%select(SNP, A1, TEST, OR, SE, L95, U95, beta, P)
view(data.compare)

write.csv(data.compare, "CAS Compare.csv", row.names = FALSE)

#9 LDL-C Association 
project.linear.adj.assoc <- read.csv("C:/Users/noahs/OneDrive/Desktop/Fall 2024/BS858 Statistical Genetics/Plink/project.linear.adj.assoc.linear", sep="")
data9 <- project.linear.adj.assoc
view(data9)

thresh9 <- data9 %>% filter(TEST == "ADD") %>% filter(P < 0.002) %>% view()

#10 Rare Variant Association
view(data10)

data10 <- data %>% select(-7, -8, -9, -10)
geno<- data10[7:18]
geno$CMC <- apply(geno,1,sum)
CMC <- geno$CMC
geno$CAST <- ifelse(geno$CMC>0,1,0)
CAST <- geno$CAST

data10$CAS1[data10$CAS1 == 1] <- 0
data10$CAS1[data10$CAS1 == 2] <- 1
data10 <- as.data.frame(cbind(data10,CAST))
data10 <- as.data.frame(cbind(data10,CMC))
view(data10)

data10$rv1maf <- sum(data10$RV.1)/ 12000
data10$rv2maf <- sum(data10$RV.2)/ 12000
data10$rv3maf <- sum(data10$RV.3)/ 12000
data10$rv4maf <- sum(data10$RV.4)/ 12000
data10$rv5maf <- sum(data10$RV.5)/ 12000
data10$rv6maf <- sum(data10$RV.6)/ 12000
data10$rv7maf <- sum(data10$RV.7)/ 12000
data10$rv8maf <- sum(data10$RV.8)/ 12000
data10$rv9maf <- sum(data10$RV.9)/ 12000
data10$rv10maf <- sum(data10$RV.10)/ 12000
data10$rv11maf <- sum(data10$RV.11)/ 12000
data10$rv12maf <- sum(data10$RV.12)/ 12000

maf <- as.numeric(data10[1,21:32] )
weights <- as.vector(1 / sqrt(6000*maf*(1-maf)))
weights
geno2<-geno[1:12]
geno2<-as.matrix(geno2)
xxx <- sapply(1:12, function(i)geno2[,i]*weights[i])
xxx
view(xxx)
data10$MB <- apply(xxx,1,sum)
view(data10)
model10.CAST <- glm(CAS1 ~  sex1 + age1 + hospital + CAST, binomial("logit"), data = data10)
summary(model10.CAST)
#hospital not significant, CAST not significant

model10.CMC <- glm(CAS1 ~  sex1 + age1 + hospital + CMC, binomial("logit"), data = data10)
summary(model10.CMC)
# hospital not significant, CMC not significant

model10.MB <- glm(CAS1 ~  sex1 + age1 + hospital + MB, binomial("logit"), data = data10)
summary(model10.MB)
# hospital not significant, MB not significant

