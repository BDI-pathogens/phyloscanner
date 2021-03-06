plot.age.pyramid <- function(){
  indir = "~/Documents/PhD/Xiaoyue"
  
  load(file.path(indir,"data",'190327_sampling_by_gender_age.rda'))	
  des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGEYRS, SEX, INMIGRANT))	
  
  # restrict age
  des = subset(des, AGEYRS < 50)
  
  make.plots("fishing")
  make.plots("island")
}

make.plots = function(commun){
  
  library(ggplot2)
  library(data.table)
  library("stringr")
  
  if(commun == "fishing")
  { 
    des = subset(des, substr(COMM_NUM_A, 1, 1) == "f")
  } else
  {
    des =subset(des, !(substr(COMM_NUM_A, 1, 1) == "f"))
  }
  
  #	aggregate participation
  dp		<- des[, list(ELIG=length(PERM_ID), PART=length(which(PARTICIPATED==1))), by=c('AGEYRS','SEX')]
  
  opt.exclude.onART.from.denominator	<- 1
  
  dh		<- des[, list(HIVNEG=length(which(HIV_1517==0 & PARTICIPATED==1))), by=c('AGEYRS','SEX')]
  
  des		<- subset(des, HIV_1517==1& PARTICIPATED==1)
  if(opt.exclude.onART.from.denominator)
  {
    stopifnot(des[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
    des	<- subset(des, SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))		
  }
  
  ds<-des[,list(ARTNAIVE=length(PERM_ID), SEQ=length(which(MIN_PNG_OUTPUT==1))),
          by=c('AGEYRS','SEX')]
  
  dsamp<-merge(merge(dp,ds,by=c('AGEYRS','SEX')), dh, by=c('AGEYRS','SEX'))
  dsamp[,NOT_SEQ:=ARTNAIVE-SEQ]
  dsamp[,ART:=PART-ARTNAIVE-HIVNEG]# PART = HIVPOS (:= ARTNAIVE + ART) + HIVNEG 
  dsamp[,HIVNEG:=PART-ARTNAIVE-ART]
  dsamp[,NOT_PART:= ELIG-PART]
  
  dsamp<-subset(dsamp,select = c('AGEYRS','SEX','NOT_SEQ','ART', 'HIVNEG', 'NOT_PART','SEQ'))
  setnames(dsamp,c('SEQ','NOT_SEQ','ART', 'HIVNEG','NOT_PART'),
           c('sequenced','ART-naive and not sequenced','on ART', 'HIV negative','not participated'))
  
  dsamp<-melt(dsamp,id.vars = c('AGEYRS','SEX'))
  
  cbPalette <- c("#999999", "#8E97DA", "#468B00", "#FF4500", "#FFD700")
  
  dsamp$variable = factor(dsamp$variable,levels = c('not participated','HIV negative','on ART','ART-naive and not sequenced','sequenced'))
  labels_facet = c(`not participated` = "not participated",
                   `HIV negative` = 'HIV negative',
                   `on ART or HIV negative` = 'on ART',
                   `ART-naive and not sequenced` = "ART-naive and \n not sequenced",
                   `sequenced` = "sequenced")
  
  ggplot(dsamp, aes(x=AGEYRS,  fill=variable)) +
    geom_bar(data=dsamp[SEX=='M'], stat="identity", aes(y=value)) +
    geom_bar(data=dsamp[SEX=='F'], stat="identity", aes(y=-value)) +
    geom_hline(yintercept=0, colour="white", lwd=1) +
    coord_flip() +
    theme_bw()+
    scale_fill_manual(values=cbPalette)+
    theme(text = element_text(size = 20)) +
    scale_y_continuous(limits = c(-1000,1000), breaks = c(-500,0,500), labels = abs(c(-500,0,500))) +
    labs(y="\n Count", x="\n Age\n",fill='Number of individuals in each category', title =paste("Age distribution to participation in", commun, "community"), 
         subtitle =  paste("Female", str_pad("", width = 50, side = "both"), "Male")) +
    theme(legend.position="bottom",legend.direction="vertical", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  ggsave(device = "pdf", w = 10, h = 11, file = file.path(indir, "figures", paste0("pyramid", commun, ".pdf")))
  
  ggplot(dsamp, aes(x=AGEYRS, fill=variable)) +
    geom_bar(data=dsamp[SEX=='M'], stat="identity", aes(y=value)) +
    geom_bar(data=dsamp[SEX=='F'], stat="identity", aes(y=-value)) +
    geom_hline(yintercept=0, colour="white", lwd=1) +
    coord_flip() +
    theme_bw()+
    scale_fill_manual(values=cbPalette)+
    facet_grid(~variable, labeller = as_labeller(labels_facet))+
    theme(text = element_text(size = 20))+
    scale_y_continuous(limits = c(-600,600), breaks = c(-300,0,300), labels = abs(c(-300,0,300))) +
    labs(y="\n Count", x="\n Age\n",fill='Number of individuals in each category', title ="Age distribution to participation \n") +
    theme(legend.position="bottom",legend.direction="vertical", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    annotate("text", x = 52, y = 300, label = "Male")+
    annotate("text", x = 52, y = -300, label = "Female")
  
  
  dsamp.rate<-merge(merge(dp,ds,by=c('AGEYRS','SEX')), dh, by=c('AGEYRS','SEX'))
  dsamp.rate[,NOT_SEQ:=ARTNAIVE-SEQ]
  dsamp.rate[,ART:=PART-ARTNAIVE-HIVNEG]# PART = HIVPOS (:= ARTNAIVE + ART) + HIVNEG 
  dsamp.rate[,HIVNEG:=PART-ARTNAIVE-ART]
  dsamp.rate[,NOT_PART:= ELIG-PART]
  dsamp.rate[,TOTAL:=SEQ+NOT_SEQ+ART+ HIVNEG+NOT_PART]
  dsamp.rate[,NOT_SEQ:=NOT_SEQ/TOTAL]
  dsamp.rate[,ART:=ART/TOTAL]
  dsamp.rate[,HIVNEG:=HIVNEG/TOTAL]
  dsamp.rate[,NOT_PART:=NOT_PART/TOTAL]
  dsamp.rate[,SEQ:=SEQ/TOTAL]
  
  dsamp.rate<-subset(dsamp.rate,select = c('AGEYRS','SEX','NOT_SEQ','ART','HIVNEG','NOT_PART','SEQ'))
  setnames(dsamp.rate,c('SEQ','NOT_SEQ','HIVNEG', 'ART','NOT_PART'),
           c('sequenced','ART-naive and not sequenced','HIV negative', 'on ART','not participated'))
  
  dsamp.rate<-melt(dsamp.rate,id.vars = c('AGEYRS','SEX'))
  
  dsamp.rate$variable = factor(dsamp.rate$variable,levels = c('not participated','HIV negative','on ART', 'ART-naive and not sequenced','sequenced'))
  
  ggplot(dsamp.rate, aes(x=AGEYRS,  fill=variable)) +
    geom_bar(data=dsamp.rate[SEX=='M'], stat="identity", aes(y=value)) +
    geom_bar(data=dsamp.rate[SEX=='F'], stat="identity", aes(y=-value)) +
    geom_hline(yintercept=0, colour="white", lwd=1) +
    coord_flip() +
    theme_bw()+
    scale_fill_manual(values=cbPalette)+
    theme(text = element_text(size = 20))+
    scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.5), labels = abs(seq(-1,1,0.5))) +
    labs(y="\n Count", x="\n Age\n",fill='Number of individuals in each category', title =paste("Age distribution to participation in", commun, "community \n (proportion)"), 
         subtitle =  paste("Female", str_pad("", width = 50, side = "both"), "Male")) +
    theme(legend.position="bottom",legend.direction="vertical", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  ggsave(device = "pdf", w = 10, h = 11 , file = file.path(indir, "figures", paste0("pyramid_prop", commun, ".pdf")))
  
  ggplot(dsamp.rate, aes(x=AGEYRS, fill=variable)) +
    geom_bar(data=dsamp.rate[SEX=='M'], stat="identity", aes(y=value)) +
    geom_bar(data=dsamp.rate[SEX=='F'], stat="identity", aes(y=-value)) +
    geom_hline(yintercept=0, colour="white", lwd=1) +
    coord_flip() +
    theme_bw()+
    scale_fill_manual(values=cbPalette)+
    facet_grid(~variable, labeller = as_labeller(labels_facet))+
    theme(text = element_text(size = 20))+
    scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.5), labels = abs(seq(-1,1,0.5))) +
    labs(y="\n Count", x="\n Age\n",fill='Number of individuals in each category', title ="Age distribution to participation (proportion) \n") +
    theme(legend.position="bottom",legend.direction="vertical", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    annotate("text", x = 52, y = 0.5, label = "Male")+
    annotate("text", x = 52, y = -0.5, label = "Female")
}
