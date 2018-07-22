
Downsizing female IBS
------------------------

**Upload file from github respository**

```all<-read.delim("LLD_1135_basics_IBS_menstruating", row.names = 1, header = T, sep = "\t")

females<-all[all$Sex=="Female",]
femaleibs<-females[females$IBS==1,]

randomsel = function(df,n){
  return(df[sample(nrow(df),n),])
}


random1<-as.data.frame(randomsel(femaleibs,24))
random2<-as.data.frame(randomsel(femaleibs,24))
random3<-as.data.frame(randomsel(femaleibs,24))
random4<-as.data.frame(randomsel(femaleibs,24))
random5<-as.data.frame(randomsel(femaleibs,24))
```

