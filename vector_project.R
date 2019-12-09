library(ggplot2)
t<-read.table('~/Desktop/projected.csv', sep='\t',header=F)
scale<-10000
names<-as.vector(t[,1])
names<-sapply(strsplit(names, "_"), head, 1)
names<-gsub("GM", "NA", names)
v1<-as.vector(t[,2])/scale
v2<-as.vector(t[,3])/scale
v3<-as.vector(t[,4])/scale
v4<-as.vector(t[,5])/scale
allvecs<-list(v1,v2,v3,v4)
vars<-c(var(v1),var(v2),var(v3),var(v4))
ethnicity<-read.table('~/Desktop/1000pop.csv', sep=',',header=F) 

groups<-as.vector(unlist(ethnicity[match(names,ethnicity[,1]),2]))

m<-data.frame(v1,v2,v3,v4,names,groups)
plot<-ggplot(data=m,aes(v1,v3), size = 1)+geom_point(aes(colour = factor(groups)))


allgroups<-unique(groups)
allgroups<-allgroups[which(allgroups!='NA')]

centers<-list()
sds<-list()
for (group in allgroups) {
	
	center<-c()
	sd<-c()
	index<-which(groups==group)
	
	for (vec in allvecs) {
		
		center<-c(center,mean(vec[index]))
		sd<-c(sd,sd(vec[index]))
	}

	centers[[group]]<-center
	
	sds[[group]]<-sd
}


alldistances<-list()
for (i in 1:length(rownames(t))){
	
	sample<-as.vector(unlist(t[i,2:(1+length(allvecs))]))/scale
	
	distances<-c()
	relatives<-c()
	
	for (j in 1:length(centers)) {
		
		sd<-as.vector(unlist(sds[j]))
		center<-as.vector(unlist(centers[j]))
		
		distance<-sum(((sample-center))^2)
		
		distances<-c(distances,distance)
		
		
	}
	
	distances_rela<-c()
	for (j in 1:length(centers)) {
		
		sd<-as.vector(unlist(sds[j]))
		center<-as.vector(unlist(centers[j]))
		
		distance_rela<-sum(((sample-center)/sd)^2)
		
		distances_rela<-c(distances_rela,distance_rela)
		
	}
	
	
	within<-allgroups[which(distances_rela<10)]
	
	if (length(within)==0){
		
		within<-'Unknown'
		
	}
	
	else if (length(within)==1){
		
		within<-within[1]
		
	}
	else {
		
		within<-paste(within, collapse = '&')
	}
	
	
	name<-names[i]
	
	alldistances[[name]]<-c(distances,distances_rela,within)
		
}

df<-as.data.frame(do.call(rbind, alldistances))
df<-cbind(df, groups)
colnames(df)<-c(allgroups,paste(allgroups, '_relavtive'),'reported','determined')
df























