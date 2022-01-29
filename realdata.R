############# Real data analysis
rm(list=ls())
library(rTensor)
source('FTSVD2.R')
library(MASS)
library(fdapace)
library(latex2exp)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)

setwd('D:/Github/functional_tensor_svd')




##############
# 1. ECAM data
##############
counts0 <- read.csv("genus_count_cleaned.csv", row.names=1)
metadata <- read.csv("genus_metadata_cleaned.csv", row.names=1)
table(rownames(counts0)==rownames(metadata))
metauni <- unique(metadata[,c('studyid', 'delivery', 'diet')])
study.ID <- metauni$studyid
nsub <- nrow(metauni)

microbiome.data <- format_tfpca(counts0, metadata$day_of_life, metadata$studyid, threshold=0.9, pseudo_count=0.5)
# summarize the available month data...
count = rep(0,25)
p1 = nrow(metauni)
p2 = dim(microbiome.data[[1]])[1] - 1
for (i in 1:p1) {
  tn = round(microbiome.data[[i]][1,]/30)
  tn[tn==25] = 24
  tn = unique(tn)
  count[tn+1] = count[tn+1] + 1
}

# transform the data to a 42-50-25 tensor with missing values.

n = 19
microbiome.tensor = array(NA, dim= c(42,50,19))
for (i in 1:p1) {
  tn = round(microbiome.data[[i]][1,]/30)
  tn[tn==25] = 24
  for (k in 1:length(tn)) {
    if(tn[k]<=12)
      microbiome.tensor[i,,tn[k]+1] = microbiome.data[[i]][2:(p2+1),k]
    else
      microbiome.tensor[i,,ceiling(tn[k]/2)+7] = microbiome.data[[i]][2:(p2+1),k]
  }
}

data.interpolate = microbiome.tensor
for (i in 1:42) {
  for (j in 1:50) {
    if(all(is.na(microbiome.tensor[i,j,]))) data.interpolate[i,j,] = 0
    else data.interpolate[i,j,] = interpolate(microbiome.tensor[i,j,])
  }
}
data.interpolate = as.tensor(data.interpolate)

tn = c(0:12,14,16,18,20,22,24)

#res = FTSVD(data.interpolate, f.grid = c(list(NULL),list(NULL),list(tn/24)), 
#                  rank=3, penalty = 1e-5)
data.interpolate = data.interpolate - mean(data.interpolate@data)

lambda_max = 0
for (i in 1:3){
  lambda_max = max(lambda_max, svd(k_unfold(data.interpolate, i)@data)$d[1])
}
data.interpolate = data.interpolate / (lambda_max / sqrt(dim(data.interpolate)[3]))

# compare different models
r.sq = rep(0,10)
bic = r.sq
for (r in 1:10){
  res = FTSVD(data.interpolate, f.grid = c(list(NULL),list(NULL),list(tn/24)), 
              rank=r, Hbound = 5)
  res.assess = FTSVD.model.assess(data.interpolate, res, r)
  r.sq[r] = res.assess[1]
  bic[r] = res.assess[2]
}
plot(r.sq)
plot(bic)



# loading analysis.
res = FTSVD(data.interpolate, f.grid = c(list(NULL),list(NULL),list(tn/24)), 
            rank=4, Hbound = 5)
A.PC = res[[2]][,2:4]
A.PC[,1] = -A.PC[,1]
colnames(A.PC) = c("Component 1","Component 2","Component 3")
A.data <- metauni
rownames(A.data) <- A.data$studyid
A.data <- cbind(A.PC, A.data)
npc <- ncol(A.PC)
p_deliv_list <- vector("list", npc*(npc-1)/2)
ij <- 1
colnames(A.data)[5] = "Delivery"
for (i in 1:(npc-1)){
  for (j in (i+1):npc){
    p_deliv_list[[ij]] <- local({
      ij <- ij
      i <- i
      j <- j
      ptmp <- ggplot(data=A.data, aes(x=A.data[,i], y=A.data[,j], color=Delivery)) + 
        geom_point(size=3) + theme(text = element_text(size = 15)) +
        labs(x=colnames(A.data)[i], y=colnames(A.data)[j]) + scale_fill_discrete(name = "Delivery Mode")
    })
    ij <- ij+1
  }
}
print(ggarrange(grobs=p_deliv_list,nrow=1,common.legend=TRUE))

ggarrange(p_deliv_list[[1]],p_deliv_list[[2]],p_deliv_list[[3]],nrow=1,common.legend=TRUE,legend = 'bottom')
# save 3.5-by-by-8 inch


# functional analysis
Phi.data <- res[[4]]*10
Phi.data = - Phi.data
t(Phi.data)%*%Phi.data

npc <- ncol(Phi.data)
Phi.data <- data.frame( Month=seq(0,1,0.01)*24,value = as.vector(Phi.data),
                        Component = paste0('Comp ', as.vector(t(matrix(rep(1:npc,length(seq(0,1,0.01))),npc,)))))

ptime <- ggplot(data=Phi.data, aes(x=Month, y=value, color=Component)) + geom_line(size=1) +
  theme(text = element_text(size = 15), legend.title = element_blank(), legend.position="bottom", axis.title.y = element_blank()) +
  ggtitle("Functional loading estimations")
ptime


# p_time_list = vector("list", npc)
# for (i in 1:npc){
#     p_time_list[[i]] <- local({
#       i <- i
#       ptmp <- ggplot(data=Phi.data, aes(x=Month, y=Phi.data[,i+1])) + ylim(-2,2)+
#         geom_line(size=1.5) + theme(text = element_text(size = 15)) + xlim(0,24)+
#         labs(x="Month", y='') + ggtitle(colnames(Phi.data)[i+1])
#     })
# }
# ggarrange(p_time_list[[1]],p_time_list[[2]],p_time_list[[3]],nrow=1,common.legend=TRUE,legend = 'bottom')




# observed curves aggregated by PC.
# transform the result...

agg_feature_2 = ttm(data.interpolate, t(res[[3]][,2]),2)@data[,1,]
agg_feature_3 = ttm(data.interpolate, t(res[[3]][,3]),2)@data[,1,]
agg_feature_4 = ttm(data.interpolate, -t(res[[3]][,4]),2)@data[,1,]

mean = c(apply(agg_feature_2[which(metauni$delivery=="Vaginal"),], 2, mean),
           apply(agg_feature_2[which(metauni$delivery=="Cesarean"),], 2, mean),
         apply(agg_feature_3[which(metauni$delivery=="Vaginal"),], 2, mean),
         apply(agg_feature_3[which(metauni$delivery=="Cesarean"),], 2, mean),
         apply(agg_feature_4[which(metauni$delivery=="Vaginal"),], 2, mean),
         apply(agg_feature_4[which(metauni$delivery=="Cesarean"),], 2, mean))
sd = c(apply(agg_feature_2[which(metauni$delivery=="Vaginal"),], 2, FUN=function(x){sd(x)/sqrt(length(x))}),
      apply(agg_feature_2[which(metauni$delivery=="Cesarean"),], 2, FUN=function(x){sd(x)/sqrt(length(x))}),
        apply(agg_feature_3[which(metauni$delivery=="Vaginal"),], 2, FUN=function(x){sd(x)/sqrt(length(x))}),
apply(agg_feature_3[which(metauni$delivery=="Cesarean"),], 2, FUN=function(x){sd(x)/sqrt(length(x))}),
apply(agg_feature_4[which(metauni$delivery=="Vaginal"),], 2, FUN=function(x){sd(x)/sqrt(length(x))}),
apply(agg_feature_4[which(metauni$delivery=="Cesarean"),], 2, FUN=function(x){sd(x)/sqrt(length(x))})
)

agg_obs = data.frame(mean = mean, se = sd, Month = rep(tn, 6), 
                     Component = rep(c("Component 2", "Component 3", "Component 4"), each = length(tn)*2),
                       Delivery = as.factor(rep(c("Vaginal","Cesarean","Vaginal","Cesarean","Vaginal","Cesarean"),each=length(tn))))

p_feat_obs <- ggplot(data = agg_obs,
                     aes(x=Month, y=mean, group=Delivery, color=Delivery)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=mean-qnorm(0.95)*se, ymax=mean+qnorm(0.95)*se, 
                  color=Delivery, fill=Delivery), linetype=2, alpha=0.3) +
  theme(text = element_text(size = 15), axis.title.y = element_blank(), legend.title = element_blank(), legend.position="bottom") +
  ggtitle('Observed aggregated trajectory') + 
  facet_grid(Component~.) 
p_feat_obs
## save as 6-by-4 inch






# investigate the top quantile
# B.data <- as.data.frame(res[[3]])
# colnames(B.data) = c("PC1","PC2","PC3")
# B.data[,1] = B.data[,1]
# rownames(B.data) = rownames(microbiome.data[[1]])[2:123]
# taxon_sel <- rownames(B.data)[order(-abs(B.data[,2]))[1:2]]
# taxon_sel
# tab_B_obs <- NULL
# for (i in 1:length(taxon_sel)){
#   value <- unlist(sapply(microbiome.data, function(x){x[taxon_sel[i],]}))
#   time_point <- unlist(sapply(microbiome.data, function(x){x['time_point',]}))
#   nobs <- sapply(microbiome.data, function(x){ncol(x)})
#   subID <- unlist(mapply(function(i){rep(names(microbiome.data)[i], nobs[i])}, 
#                          1:length(nobs)))
#   tmp <- data.frame(value=value, time=time_point, bacname=taxon_sel[i],
#                     studyid=subID)
#   tab_B_obs <- rbind(tab_B_obs, merge(tmp, metauni))
# }
# p_ftsel_obs2 <- ggplot(data=tab_B_obs, aes(x=time, y=value, group=studyid, color=delivery)) + 
#   geom_point() + ggtitle('PC2 leading bacteria') + 
#   facet_wrap(as.factor(tab_B_obs$bacname), nrow=3) + xlab('days')
# p_ftsel_obs2











##################
# 2. food harvest data
##################
rm(list=ls())
setwd('D:/Github/functional_tensor_svd')
source("FTSVD2.R")
raw.data = read.csv("crops/Production_Crops_E_All_Data_NOFLAG.csv")
raw.data = raw.data[which(raw.data$Element=="Production"),]
Item.set = unique(raw.data$Item)
Area.set = unique(raw.data$Area)


# only consider areas instead of country
Area.set = Area.set[213:236]
Area.set = Area.set[c(-1,-7,-12,-18,-24)]


# first convert it to a tensor with NA.
crops = array(NA, dim = c(length(Area.set),length(Item.set),59))
for (i in 1:nrow(raw.data)) {
  item.ind = which(Item.set == raw.data[i,"Item"])
  area.ind = which(Area.set == raw.data[i,"Area"])
  crops[area.ind, item.ind, ] = as.numeric(raw.data[i,8:66])
}

# 1st-round thresholding 
flag = !is.na(crops)
Item.set = Item.set[apply(flag, 2, sum) > 900]
# remove some items that leads bias.
Item.set = Item.set[c(1:2,4:7,10:11,13,15,20:21,36:41,43:48,50:51)]
crops = array(NA, dim = c(length(Area.set),length(Item.set),59))
for (i in 1:nrow(raw.data)) {
  item.ind = which(Item.set == raw.data[i,"Item"])
  area.ind = which(Area.set == raw.data[i,"Area"])
  crops[area.ind, item.ind, ] = as.numeric(raw.data[i,8:66])
}

# Interpolate the missing NAs.


crops.interpolate = crops
for (i in 1:length(Area.set)) {
  for (j in 1:length(Item.set)) {
    if(all(is.na(crops[i,j,]))) crops.interpolate[i,j,] = 0
    else crops.interpolate[i,j,] = interpolate(crops[i,j,])
  }
}

crops.log = log(crops.interpolate + 1/2)

# Calculate differential...
crops.diff = array(0, dim = c(length(Area.set), length(Item.set), 58))
for (i in 1:length(Area.set)) {
  for (j in 1:length(Item.set)) {
    crops.diff[i,j,] = crops.log[i,j,2:59] - crops.log[i,j,1:58]
  }
}


save(crops.log, Area.set, Item.set, file='Crop_by_area.Rdata')

# Apply tensor methods

crops.tensor = as.tensor(crops.interpolate)
# rescale the data
crops.tensor = crops.tensor - mean(crops.tensor@data)
lambda_max = 0
for (i in 1:3){
  lambda_max = max(lambda_max, svd(k_unfold(crops.tensor, i)@data)$d[1])
}
crops.tensor = crops.tensor / (lambda_max / sqrt(dim(crops.tensor)[3]))
crops.tensor = crops.tensor - mean(crops.tensor@data)
r.sq = rep(0,10)
bic = r.sq
for (r in 1:10){
  res = FTSVD(crops.tensor, f.grid = c(list(NULL),list(NULL),list((1:59)/59)), 
              rank=r, Hbound = 5)
  res.assess = FTSVD.model.assess(crops.tensor, res, r, (1:59)/59)
  r.sq[r] = res.assess[1]
  bic[r] = res.assess[2]
}
plot(r.sq)
plot(bic)




crops.est = FTSVD(crops.tensor, f.grid = c(list(NULL),list(NULL),list((1:59)/59)), 
                  rank=2, Hbound = 2)

A.est = data.frame(PC1 = crops.est[[2]][,1], PC2 = -crops.est[[2]][,2])
row.names(A.est) = Area.set
A.est$name = row.names(A.est)
sp1 <- ggplot(A.est, aes(PC1, PC2, label = rownames(A.est)))+
  geom_point(size = 4) + geom_text(data=subset(A.est,A.est[,1]>0.18),aes(PC1,PC2,label=name),size=4,vjust=1.5) +
  xlim(0,0.7) + ylim(-0.8,0.5) + theme(text = element_text(size = 15)) + ggtitle("Area") +
xlab('Component 1') + ylab('Component 2')

B.est = data.frame(PC1 = -crops.est[[3]][,1], PC2 = crops.est[[3]][,2])
row.names(B.est) = Item.set
B.est$name = row.names(B.est)
sp2 <- ggplot(B.est, aes(PC1, PC2, label = rownames(B.est)))+
  geom_point(size = 4) + geom_text(data=subset(B.est,B.est[,1]>0.1),aes(PC1,PC2,label=name),size=4,vjust=1.5) +
  xlim(0,0.7) + ylim(-0.8,0.8) + theme(text = element_text(size = 15)) + ggtitle("Item") + xlab('Component 1') + ylab('Component 2')
plot(PC2~PC1, col="lightblue", pch=19, cex=2,data=B.est)
text(PC2~PC1, labels=rownames(B.est[which(B.est[,1]>0.1),]),data=subset(B.est,B.est[,1]>0.1), cex=0.9, font=2)

t = seq(1961,2019,length.out = 101)

Phi.est = data.frame(Year = t, PC = rep(c("Component 1","Component 2"),each = 101),
                     value = c(-crops.est[[4]][,1], -crops.est[[4]][,2]))
sp3 <- ggplot(Phi.est, aes(x=Year,y=value)) + geom_line(aes(color = PC,linetype=PC),size=1.5) +
  scale_color_manual(values = c("darkred", "steelblue")) + theme(text = element_text(size = 15), legend.position = c(0.7, 0.2),legend.title = element_blank()) +
  labs(y='') + xlim(1961,2019) + ggtitle("Time")

ggarrange(sp1,sp2,sp3,nrow=1)


plot(crops.est[[4]][,1])
plot(-crops.est[[4]][,2])
plot(crops.est[[4]][,3])
plot(crops.est[[4]][,4])

C.est = data.frame(PC1 = crops.est[[3]][,1], PC2 = crops.est[[3]][,2])
row.names(B.est) = Item.set
plot(PC2~PC1, col="lightblue", pch=19, cex=2,data=B.est)
text(PC2~PC1, labels=rownames(B.est),data=B.est, cex=0.9, font=2)