


## produce vertical distribution plot for gobies--------------

library("ggplot2")
library("plyr")
library("reshape")

goby <- read.csv("Goby_Vert.csv")

goby$depth <- (goby$MinZ+goby$MaxZ)/2
goby$depth <- round_any(goby$depth,10)

goby$number <- goby$Pre + goby$Post



# 
# ggplot(goby, aes(x=as.factor(depth),y=number)) + 
#   geom_boxplot() +
#   ylim(c(0,10)) 
# # scale_y_reverse( lim=c(90,0)) +
# xlab("Mean Abundance") + ylab("Depth (m)") + theme_bw() + 
#   theme(axis.title.x = element_text(face="bold", size=14)) + 
#   theme(axis.title.y = element_text( face="bold", size=14)) +
#   theme(axis.text.x = element_text(size=8, vjust=.05)) + 
#   theme(axis.text.y = element_text(size=8, vjust=.05))
# 
# 
# 
# 
# gobysum <- ddply(goby,~depth, summarise,
#                  mean=mean(number),
#                  sum=sum(number))
# 
# 
# 
# 
# ggplot(gobysum, aes(y=mean,group=depth)) + 
#   geom_boxplot() 
#   scale_y_reverse( lim=c(90,0)) +
#   xlab("Mean Abundance") + ylab("Depth (m)") + theme_bw() + 
#   theme(axis.title.x = element_text(face="bold", size=14)) + 
#   theme(axis.title.y = element_text( face="bold", size=14)) +
#   theme(axis.text.x = element_text(size=8, vjust=.05)) + 
#   theme(axis.text.y = element_text(size=8, vjust=.05))
#   
#   
#   
  

    
    
goby2 <- goby[goby$number !=0,]
goby2 <- subset(goby2, depth!=20 & depth!=40 & depth!=60)
goby2 <- subset(goby2, ROT <= 4)

goby_expand1<- as.data.frame(goby2[rep(row.names(goby2), goby2$Pre), c("depth")])
colnames(goby_expand1) <- "depth"
goby_expand1$flex <- "pre"

goby_expand2<- as.data.frame(goby2[rep(row.names(goby2), goby2$Post), c("depth")])
colnames(goby_expand2) <- "depth"
goby_expand2$flex <- "post"

goby_expand <- rbind(goby_expand1, goby_expand2)


ggplot(goby_expand,aes(x=depth)) + 
  # geom_density(size=1.5,color="black") + 
  geom_histogram(binwidth=10) +
  coord_flip() + scale_x_reverse(lim=c(90,0),breaks=c(90,70,50,30,10)) +
  # ylim(c(0.01,0.06)) +
  xlab("Depth (m)") + ylab("Proportion of Larvae") + theme_bw() + 
  theme(axis.title.x = element_text(face="bold", size=14)) + 
  theme(axis.title.y = element_text( face="bold", size=14)) +
  theme(axis.text.y = element_text(size=8, vjust=.05)) +
  theme(axis.text.x = element_blank())


# ggplot(goby,aes(y=number,x=depth)) + 
#   geom_smooth(method="loess") +
#   coord_flip() + scale_x_reverse(lim=c(90,0),breaks=c(90,70,50,30,10)) +
#   # ylim(c(0.01,0.06)) +
#   xlab("Depth (m)") + ylab("Proportion of Larvae") + theme_bw() + 
#   theme(axis.title.x = element_text(face="bold", size=14)) + 
#   theme(axis.title.y = element_text( face="bold", size=14)) +
#   theme(axis.text.y = element_text(size=8, vjust=.05)) +
#   theme(axis.text.x = element_blank())


## produce a vertical distribution for mahi mahi -------------



library("ggplot2")
library("plyr")
library("reshape")

setwd("C:/Users/mattf/Desktop")
mahi <- read.csv("coryphaena_bdos.csv")

mahi$depth <- (mahi$mindepth+mahi$maxdepth)/2
mahi$depth <- round_any(mahi$depth,10)

mahi$number <- mahi$preflexion + mahi$postflexion



# 
# ggplot(mahi, aes(x=as.factor(depth),y=number)) + 
#   geom_boxplot() +
#   ylim(c(0,10)) 
# # scale_y_reverse( lim=c(90,0)) +
# xlab("Mean Abundance") + ylab("Depth (m)") + theme_bw() + 
#   theme(axis.title.x = element_text(face="bold", size=14)) + 
#   theme(axis.title.y = element_text( face="bold", size=14)) +
#   theme(axis.text.x = element_text(size=8, vjust=.05)) + 
#   theme(axis.text.y = element_text(size=8, vjust=.05))
# 
# 
# 
# 
# mahisum <- ddply(mahi,~depth, summarise,
#                  mean=mean(number),
#                  sum=sum(number))
# 
# 
# 
# 
# ggplot(mahisum, aes(y=mean,group=depth)) + 
#   geom_boxplot() 
#   scale_y_reverse( lim=c(90,0)) +
#   xlab("Mean Abundance") + ylab("Depth (m)") + theme_bw() + 
#   theme(axis.title.x = element_text(face="bold", size=14)) + 
#   theme(axis.title.y = element_text( face="bold", size=14)) +
#   theme(axis.text.x = element_text(size=8, vjust=.05)) + 
#   theme(axis.text.y = element_text(size=8, vjust=.05))
#   
#   
#   


mahi2 <- mahi[mahi$number !=0 & !is.na(mahi$depth),]
mahi2 <- subset(mahi2, depth!=20 & depth!=40 & depth!=60)
mahi2 <- subset(mahi2, ROT <= 4)

mahi_expand1<- as.data.frame(mahi2[rep(row.names(mahi2), mahi2$preflexion), c("depth")])
colnames(mahi_expand1) <- "depth"
mahi_expand12<- as.data.frame(mahi2[rep(row.names(mahi2), mahi2$preflexion), c("rotation")])
colnames(mahi_expand12) <- "rotation"
mahi_expand_pre <- cbind(mahi_expand1,mahi_expand12)
mahi_expand_pre$flex <- "preflexion"

mahi_expand2<- as.data.frame(mahi2[rep(row.names(mahi2), mahi2$postflexion), c("depth")])
colnames(mahi_expand2) <- "depth"
mahi_expand22<- as.data.frame(mahi2[rep(row.names(mahi2), mahi2$postflexion), c("rotation")])
colnames(mahi_expand22) <- "rotation"
mahi_expand_post <- cbind(mahi_expand2,mahi_expand22)
mahi_expand_post$flex <- "postflexion"

mahi_expand <- rbind(mahi_expand_pre, mahi_expand_post)


ggplot(mahi_expand,aes(y=..density..,x=depth)) + 
  # geom_density(size=1.5,color="black") + 
  geom_histogram(binwidth=10) +
  coord_flip() + scale_x_reverse(lim=c(90,0),breaks=c(90,70,50,30,10)) +
  # ylim(c(0.01,0.06)) +
  facet_wrap(~flex) +
  xlab("Depth (m)") + ylab("Proportion of Larvae") + theme_bw() + 
  theme(axis.title.x = element_text(face="bold", size=14)) + 
  theme(axis.title.y = element_text( face="bold", size=14)) +
  theme(axis.text.y = element_text(size=8, vjust=.05)) 


# ggplot(mahi,aes(y=number,x=depth)) + 
#   geom_smooth(method="loess") +
#   coord_flip() + scale_x_reverse(lim=c(90,0),breaks=c(90,70,50,30,10)) +
#   # ylim(c(0.01,0.06)) +
#   xlab("Depth (m)") + ylab("Proportion of Larvae") + theme_bw() + 
#   theme(axis.title.x = element_text(face="bold", size=14)) + 
#   theme(axis.title.y = element_text( face="bold", size=14)) +
#   theme(axis.text.y = element_text(size=8, vjust=.05)) +
#   theme(axis.text.x = element_blank())

mahi_expand$group <- "nope"

for (i in 1:nrow(mahi_expand)){
  if (mahi_expand[i,c("rotation")] <= 4){
    mahi_expand[i,c("group")] <- "ocean"
  } else {
    mahi_expand[i,c("group")] <- "intrusion"
  }
}

ggplot(mahi_expand,aes(y=..density..,x=depth)) + 
  # geom_density(size=1.5,color="black") + 
  geom_histogram(binwidth=10) +
  coord_flip() + scale_x_reverse(lim=c(90,0),breaks=c(90,70,50,30,10)) +
  # ylim(c(0.01,0.06)) +
  facet_grid(group~flex) +
  xlab("Depth (m)") + ylab("Proportion of Larvae") + theme_bw() + 
  theme(axis.title.x = element_text(face="bold", size=14)) + 
  theme(axis.title.y = element_text( face="bold", size=14)) +
  theme(axis.text.y = element_text(size=8, vjust=.05)) 


