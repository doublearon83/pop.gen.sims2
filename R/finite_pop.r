#' Simulating single locus, two allele Mendelian populations
#' with finite population sizes over many generations
#'
#' @description \code{finite_pop} simulates Mendelian populations over many generations. Specify initial allele
#' frequencies, population size, and the number of generations.  The simulation produces each generation of
#' the population by randomly selecting 2N (ND) alleles for N (ND) individuals from a gamete pool based on 
#' allele frequencies.
#' @param p1g2 Allele 1 frequency
#' @param p2g2 Allele 2 frequency
#' @param alleles allele types (use 1 and 2)
#' @param ND Population size
#' @param generations number of generations
#' @details \code{finite_pop} returns three plots: 
#' 1) Allele 1 frequency versus number of generations, 2) Observed and expected heterozygosity versus number
#' of generations, and 3) Wright's F with a \link{loess} regression best fit line.
#' It also returns Fisher's exact tests (\link{fisher.test}) comparing the observed genotype counts to the expected
#' and initial values.
#' @export
#' @examples
#' finite_pop(0.3,0.7,c(1,2),100,100)

finite_pop<-function (p1g2,p2g2,alleles,ND,generations) {
  
  #genotype simulation
  Ho_f<-numeric(generations)
  He_f<-numeric(generations)
  Fstatisticf<-numeric(generations)
  ffA1<-numeric(generations)
  ffA2<-numeric(generations)
  p1g21<-p1g2;p2g21<-p2g2
  for (i in 1:generations) {
    g2A1<-sample(alleles,ND,replace=TRUE,prob=c(p1g2,p2g2))
    g2A2<-sample(alleles,ND,replace=TRUE,prob=c(p1g2,p2g2))
    g2<-cbind(g2A1,g2A2)
    g2<-data.frame(g2)
    p1g2<-length(g2[g2==1])/(2*ND)
    p2g2<-length(g2[g2==2])/(2*ND)
    
    ##calculate F-statistic and test for HW 
    #expected Hardy-Weinberg genotype frequencies given allele frequency
    HWf11<-(((sum(length(g2[g2[,1]==1,][,1]),length(g2[g2[,2]==1,][,2])))/(2*ND))^2)
    HWf22<-(((sum(length(g2[g2[,1]==2,][,1]),length(g2[g2[,2]==2,][,2])))/(2*ND))^2)
    HWf12<-(((sum(length(g2[g2[,1]==1,][,1]),length(g2[g2[,2]==1,][,2])))/(2*ND))*((sum(length(g2[g2[,1]==2,][,1]),length(g2[g2[,2]==2,][,2])))/(2*ND))*2)
    
    #observed genotype frequencies
    ff11<-length(g2[g2[,1]==1 & g2[,2]==1,][,1])/ND
    ff22<-length(g2[g2[,1]==2 & g2[,2]==2,][,1])/ND
    ff12<-sum(length(g2[g2[,1]==1 & g2[,2]==2,][,1]),length(g2[g2[,1]==2 & g2[,2]==1,][,2]))/ND
    
    #Wright's F
    Ho_f[i]<-ff12
    He_f[i]<-HWf12
    Fstatisticf[i]<-1-(ff12/HWf12)
    #Ho_f;He_f;Fstatisticf
    
    #allele frequencies
    ffA1[i]<-length(g2[g2==1])/(2*ND);ffA1
    ffA2[i]<-length(g2[g2==2])/(2*ND);ffA2
  }
  
  #Chi-squared test to determine if observed genotype frequencies differ from expected and initial
  print("Fisher's exact test for difference in genotype counts from expected")
  print(suppressWarnings(fisher.test(matrix(c(ff11*ND,ff22*ND,ff12*ND,HWf11*ND,HWf22*ND,HWf12*ND),nrow=3))))
  
  print("Fisher's exact test for difference in genotype counts from initial")
  print(suppressWarnings(fisher.test(matrix(c(ff11*ND,ff22*ND,ff12*ND,p1g21^2*ND,p2g21^2*ND,2*p1g21*p2g21*ND),nrow=3))))
  
  par( mfrow = c( 2, 2 ) )
  plot(c(1:generations),ffA1,type="l",lwd=2,ylab="A1 freq.", xlab="Generations")#plot of allele 1 over time
  plot(c(1:generations),Ho_f,type="l",lwd=2, ylab="Ho and He", xlab="Generations")#plot of He over time
  points(c(1:generations),He_f,type="l",lwd=2,lty=1, col="lightgray")#plot of Ho over time
  legend("topright",c("Ho","He"),lty=c(1,1),col=c("black","lightgray"),cex=1)
  plot(c(1:generations),Fstatisticf,type="l",lty=1,col="lightgray",ylab="Wright's F", xlab="Generations")#plot of F over time
  points(predict(loess(Fstatisticf~c(1:generations))),type="l",col="red",lwd=2)
  legend("topright",c("F","Loess line"),lty=c(1,1),col=c("lightgray","red"),cex=1)
  
}