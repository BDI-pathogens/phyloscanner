plot.age.pyramid <- function(){
  library(ggplot2)
  library(data.table)
  load('180322_sampling_by_gender_age.rda')	
  des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))	
  
  # remove individuals in the communities where no sequences were obtained successfully
  tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
  des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')
  
  # age group
  des[,AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("15-19","20-24","25-29","30-34","35-39","40-44","45+"))]

  #	aggregate participation
  dp		<- des[, list(ELIG=length(PERM_ID), PART=length(which(PARTICIPATED==1))), by=c('AGE_AT_MID_C','SEX')]
  
  opt.exclude.onART.from.denominator	<- 1
  load("180322_sampling_by_gender_age.rda")
  des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))		
  
  # remove individuals in the communities where no sequences were obtained successfully
  tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
  des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')
  
  # age group
  des[,AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("15-19","20-24","25-29","30-34","35-39","40-44","45+"))]
 
  #	aggregate sequencing
  des		<- subset(des, HIV_1517==1)
  if(opt.exclude.onART.from.denominator)
  {
    stopifnot(des[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
    des	<- subset(des, SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))		
  }
  
  ds<-des[,list(ARTNAIVE=length(PERM_ID), SEQ=length(which(MIN_PNG_OUTPUT==1))),
          by=c('AGE_AT_MID_C','SEX')]
  
  dsamp<-merge(dp,ds,by=c('AGE_AT_MID_C','SEX'))
  dsamp[,NOT_SEQ:=ARTNAIVE-SEQ]
  dsamp[,ART_OR_HIVNEG:=PART-ARTNAIVE]
  dsamp[,NOT_PART:= ELIG-PART]
  dsamp<-subset(dsamp,select = c('AGE_AT_MID_C','SEX','NOT_SEQ','ART_OR_HIVNEG','NOT_PART','SEQ'))
  setnames(dsamp,c('SEQ','NOT_SEQ','ART_OR_HIVNEG','NOT_PART'),
           c('sequenced','ART-naive and not sequenced','on ART or HIV negative','not participated'))
  
  dsamp<-melt(dsamp,id.vars = c('AGE_AT_MID_C','SEX'))
  
  cbPalette <- c("#999999", "#E69F00", "#0072B2", "#D55E00")
  
  ggplot(dsamp, aes(x=AGE_AT_MID_C)) +
    geom_bar(data=dsamp[SEX=='M'], aes(y=value, fill=factor(variable,levels = c('not participated','on ART or HIV negative','ART-naive and not sequenced','sequenced'))), stat="identity") +
    geom_bar(data=dsamp[SEX=='F'], aes(y=-value, fill=factor(variable,levels = c('not participated','on ART or HIV negative','ART-naive and not sequenced','sequenced'))), stat="identity") +
    geom_hline(yintercept=0, colour="white", lwd=1) +
    coord_flip() + 
    theme_bw()+
    scale_fill_manual(values=cbPalette)+
    theme(text = element_text(size = 20))+
    scale_x_discrete(expand = c(0,0))+ 
    scale_y_continuous(limits = c(-5100,5100))+ 
    labs(y="\n Count", x="\n Age groups\n",fill='Number of individuals in each category') +
    theme(legend.position="bottom",legend.direction="vertical")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    ggtitle("Female                                                     Male \n")
  
  dsamp.rate<-merge(dp,ds,by=c('AGE_AT_MID_C','SEX'))
  dsamp.rate[,NOT_SEQ:=ARTNAIVE-SEQ]
  dsamp.rate[,ART_OR_HIVNEG:=PART-ARTNAIVE]
  dsamp.rate[,NOT_PART:= ELIG-PART]
  dsamp.rate[,TOTAL:=SEQ+NOT_SEQ+ART_OR_HIVNEG+NOT_PART]
  dsamp.rate[,NOT_SEQ:=NOT_SEQ/TOTAL]
  dsamp.rate[,ART_OR_HIVNEG:=ART_OR_HIVNEG/TOTAL]
  dsamp.rate[,NOT_PART:=NOT_PART/TOTAL]
  dsamp.rate[,SEQ:=SEQ/TOTAL]
  
  dsamp.rate<-subset(dsamp.rate,select = c('AGE_AT_MID_C','SEX','NOT_SEQ','ART_OR_HIVNEG','NOT_PART','SEQ'))
  setnames(dsamp.rate,c('SEQ','NOT_SEQ','ART_OR_HIVNEG','NOT_PART'),
           c('sequenced','ART-naive and not sequenced','on ART or HIV negative','not participated'))
  
  dsamp.rate<-melt(dsamp.rate,id.vars = c('AGE_AT_MID_C','SEX'))
  
  ggplot(dsamp.rate, aes(x=AGE_AT_MID_C)) +
    geom_bar(data=dsamp.rate[SEX=='M'], aes(y=value, fill=factor(variable,levels = c('not participated','on ART or HIV negative','ART-naive and not sequenced','sequenced'))), stat="identity") +
    geom_bar(data=dsamp.rate[SEX=='F'], aes(y=-value, fill=factor(variable,levels = c('not participated','on ART or HIV negative','ART-naive and not sequenced','sequenced'))), stat="identity") +
    geom_hline(yintercept=0, colour="white", lwd=1) +
    coord_flip() + 
    theme_bw()+
    scale_fill_manual(values=cbPalette)+
    theme(text = element_text(size = 20))+
    scale_x_discrete(expand = c(0,0))+ 
    scale_y_continuous(limits = c(-1,1))+ 
    labs(y="\n Count", x="\n Age groups\n",fill='Proportion of individuals in each categories') +
    theme(legend.position="bottom",legend.direction="vertical")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    ggtitle("Female                                                     Male \n")
}