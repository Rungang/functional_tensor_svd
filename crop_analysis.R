rm(list=ls())
library(rTensor)
source('FTSVD.R')
library(MASS)
library(fdapace)
library(latex2exp)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
raw.data = read.csv("Production_Crops_E_All_Data_NOFLAG.csv")
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

# Apply tensor methods

crops.tensor = as.tensor(crops.interpolate)
# rescale the data
crops.tensor = crops.tensor - mean(crops.tensor@data)
Y = crops.tensor
a.est = svd(k_unfold(Y,1)@data)$u[,1]
Y2 = ttm(Y, t(as.matrix(a.est)), 1)
lambda_lower = svd(k_unfold(Y2,2)@data)$d[1] /sqrt(59)

crops.est = FTSVD(crops.tensor, f.grid = c(list(NULL),list(NULL),list((1:59)/59)), 
                  rank=2, alpha = 200*lambda_lower)

A.est = data.frame(PC1 = crops.est[[2]][,1], PC2 = crops.est[[2]][,2])
row.names(A.est) = Area.set
A.est$name = row.names(A.est)
sp1 <- ggplot(A.est, aes(PC1, PC2, label = rownames(A.est)))+
  geom_point(size = 4) + geom_text(data=subset(A.est,A.est[,1]>0.15),aes(PC1,PC2,label=name),size=4,vjust=1.5) +
  xlim(-0.1,0.7) + ylim(-0.8,0.5) + theme(text = element_text(size = 15)) + ggtitle("Area") +
  xlab('Component 1') + ylab('Component 2')

B.est = data.frame(PC1 = -crops.est[[3]][,1], PC2 = -crops.est[[3]][,2])
row.names(B.est) = Item.set
B.est$name = row.names(B.est)
sp2 <- ggplot(B.est, aes(PC1, PC2, label = rownames(B.est)))+
  geom_point(size = 4) + geom_text(data=subset(B.est,B.est[,1]>0.1),aes(PC1,PC2,label=name),size=4,vjust=1.5) +
  xlim(-0.1,0.7) + ylim(-0.8,0.8) + theme(text = element_text(size = 15)) + ggtitle("Item") + xlab('Component 1') + ylab('Component 2')

t = seq(1961,2019,length.out = 101)

Phi.est = data.frame(Year = t, PC = rep(c("Component 1","Component 2"),each = 101),
                     value = c(-crops.est[[4]][,1], -crops.est[[4]][,2]))
sp3 <- ggplot(Phi.est, aes(x=Year,y=value)) + geom_line(aes(color = PC,linetype=PC),size=1.5) +
  scale_color_manual(values = c("darkred", "steelblue")) + theme(text = element_text(size = 15), legend.position = c(0.7, 0.2),legend.title = element_blank()) +
  labs(y='') + xlim(1961,2019) + ggtitle("Time")

ggarrange(sp1,sp2,sp3,nrow=1)