afterlife=matrix(c(375,134,435,147),nrow=2,byrow=TRUE)
dimnames(afterlife)=list(c("Laki-laki","Perempuan"),c("Percaya","Tidak"))
afterlife
names(dimnames(afterlife))=c("JK","Kepercayaan")
afterlife

JK=c("Perempuan","Perempuan","Laki-laki","Laki-laki")
Kepercayaan=c("percaya","Tidak","Percaya","Tidak")
Jumlah=c(435,147,375,134)
afterlife1=data.frame(JK,Kepercayaan,Jumlah)
afterlife1

attach(afterlife1)
data=tapply(Jumlah,list(JK,Kepercayaan),c)
data
names(dimnames(data))=c("JK","Kepercayaan")
data

data1=data[,c(2,1)]
data1

#Cara1
total=sum(afterlife)
afterlife/total

###Ynta,fungsi : apply(x,MARGIN,FUN,...)
total_baris=apply(afterlife,1,sum)
total_baris
total_kolom=apply(afterlife,2,sum)
total_kolom

#Proporsi 
#
prop.baris=sweep(afterlife,1,total_baris,"/"))
round(prop.baris,3)
prop.kolom=sweep(afterlife,2,total_kolom,"/"))
round(prop.kolom,3)

#Cara2
afterlife.test=prop.test(afterlife)
afterlife.test
afterlife.test=prop.test(afterlife,correct=F)
afterlife.test1

afterlife.test$estimate
afterlife.test$conf.int
round(afterlife.test$conf.int,3)

#Asosiasi 2x2
##beda peluang
bp=afterlife.test$estimate[1]-afterlife.test$estimate[2]
bp

#RR
rr=afterlife.test$estimate[1]/afterlife.test$estimate[2]
rr

#Odds Rasio
afterlife.test$estimate
odds=afterlife.test$estimate/(1-afterlife.test$estimate)
odds
or=odds[1]/odds[2]
or

atau
(afterlife[1,1]*afterlife[2,2])/afterlife[2,1]*afterlife[1,2])

#Confident Interval
theta=odds[1]/odds[2]
ASE=sqrt(sum(1/afterlife))
logtheta.CI=log(theta)+c(-1,1)*1.96*ASE
logtheta.CI
exp(logtheta.CI)


odds.ratio <- function(x, pad.zeros=FALSE, conf.level=0.95){ 
if (pad.zeros){ 
if (any(x==0)) x <- x + 0.5} 
theta<-x[1,1]*x[2,2]/(x[2,1]*x[1,2]) 
ASE<-sqrt(sum(1/x)) 
CI<-exp(log(theta) 
+ c(-1,1) * qnorm(0.5*(1+conf.level)) *ASE ) 
list(estimator=theta, 
ASE=ASE, 
conf.interval=CI, 
conf.level=conf.level)}
odds.ratio(afterlife)


uk.as=function(x,alpha=0,05)
#Penaksir peluang bersyarat
n1.<-x[1,1]+x[1,2]
n2.<-x[2,1]+x[2,2]
p1.1<-x[1,1]/n1.
p1.2<-x[2,1]/n2.
#Ukuran Asosiasi
bp<-p1.1-p1.2
rr<-p1.1/p1.2
or<-(x[1,1]*x[2,2])/(x[1,2]*x[2,1])
#Standard error
se.bp<-sqrt(p1.1*(1-p1.1)/n1.+p1.2*(1-p1.2)/n2.)
se.rr<-sqrt((1-p1.1)/(p1.1*n1.)+(1-p1.2)/(p1.2*n2.))
se.or<-sqrt(1/x[1,1]+1/x[1,2]+1/x[2,1]+1/x[2,2])
#Interval Confident
lb.bp<-bp-qnorm(1-alpha/2)*se.bp
ub.bp<-bp+qnorm(1-alpha/2)*se.bp
lb.rr<-log(rr)-qnorm(1-alpha/2)*se.rr
ub.rr<-log(rr)+qnorm(1-alpha/2)*se.rr
lb.or<-log(or)-qnorm(1-alpha/2)*se.or
ub.or<-log(or)+qnorm(1-alpha/2)*se.or
#output
output<-data.frame(Penaksir=bp,Std.Err=se.bp,B.Bawah=lb.bp,B.Atas=ub.bp)
output<-rbind(output,c(log(rr),se.rr,lb.rr,ub.rr),c(rr,NA,exp(lb.rr),exp(ub.rr)))
output<-rbind(output, c(log(or),se.or,lb.or,ub.or),c(or,NA,exp(lb.or),exp(ub.or)))
rownames(output)<-c("Beda.Pel","log.RR","Res.Relatif","log.OR","Odds.Ratio")
output=t(output)
print(output)
}
uk.as(afterlife)

#Uji Independensi
##Chi-square
chisq.test(afterlife)
chisq.test(afterlife,simulate.p.value=TRUE,B=10000)

#Fungsi yang menampilkan uji independensi untuk tabe; kontingensi
ind.test <-function(x, digits=4) { 
total <- sum(x) 
rowsum <- apply(x,1,sum) 
colsum <- apply(x,2,sum) 
prop <- x/total
rowprop <- sweep(x,1,rowsum,"/") 
colprop <- sweep(x,2,colsum,"/") 
expected <- (matrix(rowsum) %*% t(matrix(colsum))) / total 
dimnames(expected) <- dimnames(x) 
resid <- (x-expected)/sqrt(expected) 
adj.resid <- resid / 
sqrt((1-matrix(rowsum)/total) %*% t(1-matrix(colsum)/total)) 
df <- prod(dim(x)-1) 
X2 <- sum(residˆ2) 
attr(X2,"P-value") <- 1-pchisq(X2,df)
# Perhatikan Nilai frekuensi yang nol. 
tmp <- x*log(x/expected) 
tmp[x==0] <- 0 
G2 <- 2 * sum(tmp) 
attr(G2,"P-value") <- 1-pchisq(G2,df) 
list(sample.size=total,row.totals=rowsum, 
col.totals=colsum, overall.proportions=prop, 
row.proportions=rowprop, col.proportions=colprop, 
expected.freqs=expected, residuals=resid, 
adjusted.residuals=adj.resid, chi.square=X2, 
likelihood.ratio.stat=G2, df=df)
}
ind.test(afterlife)

#Data Ordinal#
cacat=matrix(c(17066,48,14464,38,788,5,126,1,37,1),5,2,byrow=T)
dimnames(cacat)=

tren.test <- function(x, u, v){ 
prop.tot<-x/sum(x) 
prop.row<-apply(prop.tot,1,sum) 
prop.col<-apply(prop.tot,2,sum) 
u.rat<-sum(u*prop.col) 
v.rat<-sum(v*prop.row) 
r<-sum(matrix(v-v.rat)%*%t(matrix(u-u.rat))*prop.tot)/sqrt(sum((u-u.rat)^2*prop.col)*sum((v-v.rat)^2*prop.row))
M2<-(sum(x)-1)*r^2 
attr(M2,"P-value")<-1-pchisq(M2,1)
##Pendekatan Normal untuk Sampel Besar 
M<-sqrt(M2) 
attr(M,"P-value")<-1-pnorm(M) 
list(korelasi=r, 
stat.uji.M2=M2, 
stat.uji.M=M) 
}
tren.test(cacat,u,v)
tren.test(cacat,c(0,1),c(0,0.5,1.5,4,7))

#EksakFisher
install.packages("ctest")
library(ctest)
teh=matrix(x(3,1,1,3),ncol=2)
dimnames(teh)=list(cangkir=c)"susu","teh"),Tamu-c("susu","teh"))
teh
fisher.test(teh)
fisher.test(teh,alternative="greater")


#TIGA arah
dp<-c(53,414,11,37,0,16,4,139) 
dp<-array(dp, dim=c(2,2,2))
dimnames(dp)<-list(DeathPen=c("yes","no"),Defendant=c("white","black"),Victim=c("white","black"))
dp

Menggunakan fungsi “xtabs” 
datalabel<-list(defendant=c("white","black"),death=c("yes","no") , victim=c("white","black")) 
data.dp<-expand.grid(defendant=c("white","black"),death=c("yes","no"),victim=c("white","black"))
data.dp<-cbind(data.dp,count=c(53,414,11,37,0,16,4,139) )
xtabs(count~death+defendant+victim,data=data.dp)
dp
ftable(dp,row.vars=c("Victim","Defendant"),col.vars="DeathPen")

odds_ratio_parsial = function(x){
return ((x[1,1]+0.5)*(x[2,2]+0.5)/((x[2,1]+0.5)*(x[1,2]+0.5)))}
parsial=apply(dp,3,odds_ratio_parsial)
parsial

library(vcd) 
summary(oddsratio(dp,log=F,stratum=3))

marginal_data=c(53+0,414+16,11+4,37+139)
marginal_data=array(marginal_data,dim=c(2,2))
dimnames(marginal_data)=list(DeathPen=c("yes","no"),Defendant=c("white","black"))
marginal_data

odds_ratio_marginal = function(x){
 return ((x[1,1]+0.5)*(x[2,2]+0.5)/((x[2,1]+0.5)*(x[1,2]+0.5)))}
marginal = odds_ratio_marginal(marginal_data)
marginal