############# ECAM data analysis
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


counts0 <- read.csv("genus_count_cleaned.csv", row.names=1)
metadata <- read.csv("genus_metadata_cleaned.csv", row.names=1)
table(rownames(counts0)==rownames(metadata))
metauni <- unique(metadata[,c('studyid', 'delivery', 'diet')])
study.ID <- metauni$studyid
nsub <- nrow(metauni)

microbiome.data <- format_tfpca(counts0, metadata$day_of_life, metadata$studyid, threshold=0.9, pseudo_count=0.5)
# summarize the available month data
count = rep(0,25)
p1 = nrow(metauni)
p2 = dim(microbiome.data[[1]])[1] - 1
for (i in 1:p1) {
  tn = round(microbiome.data[[i]][1,]/30)
  tn[tn==25] = 24
  tn = unique(tn)
  count[tn+1] = count[tn+1] + 1
}

# transform the data to a 42-50-19 tensor with missing values.
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
data.interpolate = data.interpolate - mean(data.interpolate@data)
Y = data.interpolate
#data.interpolate = data.interpolate / (norm(k_unfold(data.interpolate,1)@data,'F') / sqrt(42*50))
a.est = svd(k_unfold(Y,1)@data)$u[,1]
Y2 = ttm(Y, t(as.matrix(a.est)), 1)
lambda_lower = svd(k_unfold(Y2,2)@data)$d[1] /sqrt(length(tn))

# apply FTSVD
res = FTSVD(data.interpolate, f.grid = c(list(NULL),list(NULL),list(tn/24)), 
            rank=3, alpha = 200 * lambda_lower)

# loading analysis.
A.PC = res[[2]][,1:3]
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


agg_feature_1 = ttm(data.interpolate, t(res[[3]][,1]),2)@data[,1,]
agg_feature_2 = ttm(data.interpolate, t(res[[3]][,2]),2)@data[,1,]
agg_feature_3 = ttm(data.interpolate, -t(res[[3]][,3]),2)@data[,1,]

mean = c(apply(agg_feature_1[which(metauni$delivery=="Vaginal"),], 2, mean),
         apply(agg_feature_1[which(metauni$delivery=="Cesarean"),], 2, mean),
         apply(agg_feature_2[which(metauni$delivery=="Vaginal"),], 2, mean),
         apply(agg_feature_2[which(metauni$delivery=="Cesarean"),], 2, mean),
         apply(agg_feature_3[which(metauni$delivery=="Vaginal"),], 2, mean),
         apply(agg_feature_3[which(metauni$delivery=="Cesarean"),], 2, mean))
sd = c(apply(agg_feature_1[which(metauni$delivery=="Vaginal"),], 2, FUN=function(x){sd(x)/sqrt(length(x))}),
       apply(agg_feature_1[which(metauni$delivery=="Cesarean"),], 2, FUN=function(x){sd(x)/sqrt(length(x))}),
       apply(agg_feature_2[which(metauni$delivery=="Vaginal"),], 2, FUN=function(x){sd(x)/sqrt(length(x))}),
       apply(agg_feature_2[which(metauni$delivery=="Cesarean"),], 2, FUN=function(x){sd(x)/sqrt(length(x))}),
       apply(agg_feature_3[which(metauni$delivery=="Vaginal"),], 2, FUN=function(x){sd(x)/sqrt(length(x))}),
       apply(agg_feature_3[which(metauni$delivery=="Cesarean"),], 2, FUN=function(x){sd(x)/sqrt(length(x))})
)

agg_obs = data.frame(mean = mean, se = sd, Month = rep(tn, 6), 
                     Component = rep(c("Component 1", "Component 2", "Component 3"), each = length(tn)*2),
                     Delivery = as.factor(rep(c("Vaginal","Cesarean","Vaginal","Cesarean","Vaginal","Cesarean"),each=length(tn))))

p_feat_obs <- ggplot(data = subset(agg_obs, Component!="Component 1"),
                     aes(x=Month, y=mean, group=Delivery, color=Delivery)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=mean-qnorm(0.95)*se, ymax=mean+qnorm(0.95)*se, 
                  color=Delivery, fill=Delivery), linetype=2, alpha=0.3) +
  theme(text = element_text(size = 15), axis.title.y = element_blank(), legend.title = element_blank(), legend.position="bottom") +
  ggtitle('Observed aggregated trajectory') + 
  facet_grid(Component~.) 
p_feat_obs