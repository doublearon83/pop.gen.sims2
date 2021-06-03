#' Simulating (large) Mendelian populations
#'
#' @description \code{lp} simulates Mendelian populations. Specify allele
#' frequencies and population size.  The simulation generates a population with N individuals 
#' and 2N alleles by randomly selecting 2N alleles for N individuals from a gamete pool based on 
#' allele frequencies. 
#' @param p1 Allele 1 frequency
#' @param p2 Allele 2 frequency
#' @param alleles allele types (use 1 and 2)
#' @param N Population size
#' @details \code{lp} returns the allele frequencies, observed and expected heterozygosities,
#' Wright's F, and Fisher's exact tests (\link{fisher.test}) comparing the observed genotype counts to the expected
#' and initial values.
#' @export
#' @examples
#' lp(0.3,0.7,c(1,2),100)

lp<-function(p1,p2,alleles,N){
  
  #Genotypes 
  probs<-c(p1,p2) #vector of allele frequencies, the probability that an allele is represented in a genotype depends on the specified allele frequencies.
  allele1<-sample(alleles,N,replace=TRUE,prob=probs) #selecting the first allele
  allele2<-sample(alleles,N,replace=TRUE,prob=probs) #and the second allele
  genotypes<-data.frame(allele1,allele2) #creating object of genotypes from alleles
  genotypes #printing genotypes
  p11<-length(genotypes[genotypes[,1]==1 & genotypes[,2]==1,][,1])/N #calculating genotype 11 freq.
  p12<-sum(length(genotypes[genotypes[,1]==1 & genotypes[,2]==2,][,1]),length(genotypes[genotypes[,1]==2 & genotypes[,2]==1,][,1]))/N #calculating genotype 12 (heterozygote) freq.)
  p22<-length(genotypes[genotypes[,1]==2 & genotypes[,2]==2,][,1])/N #calculating genotype 22 freq.
  sum(p11,p12,p22) #seeing that genotype frequencies sum to one.
  
  #####*****RUN ONLY TO THIS POINT FIRST*****#####
  ##calculate F-statistic and test for HW
  #expected Hardy-Weinberg genotype frequencies given allele frequency
  HWlp11<-(((sum(length(genotypes[genotypes[,1]==1,][,1]),length(genotypes[genotypes[,2]==1,][,2])))/(2*N))^2)
  HWlp22<-(((sum(length(genotypes[genotypes[,1]==2,][,1]),length(genotypes[genotypes[,2]==2,][,2])))/(2*N))^2)
  HWlp12<-(((sum(length(genotypes[genotypes[,1]==1,][,1]),length(genotypes[genotypes[,2]==1,][,2])))/(2*N))*((sum(length(genotypes[genotypes[,1]==2,][,1]),length(genotypes[genotypes[,2]==2,][,2])))/(2*N))*2)
  
  #observed genotype frequencies
  flp11<-length(genotypes[genotypes[,1]==1 & genotypes[,2]==1,][,1])/N
  flp22<-length(genotypes[genotypes[,1]==2 & genotypes[,2]==2,][,1])/N
  flp12<-sum(length(genotypes[genotypes[,1]==1 & genotypes[,2]==2,][,1]),length(genotypes[genotypes[,1]==2 & genotypes[,2]==1,][,2]))/N
  
  #Genetic diversity
  print("Heterozygosity")
  Fstatisticnr<-1-(flp12/HWlp12)# Wright's F
  diversity<-c(flp12,HWlp12,Fstatisticnr);names(diversity)<-c("Ho","He","Wright's F");print(diversity)
  
  #Allele frequency
  freqs<-c(length(genotypes[genotypes==1])/(2*N),length(genotypes[genotypes==2])/(2*N));names(freqs)<-c("p1","p2")
  print("Allele frequencies")
  print(freqs)
  
  #Chi-squared test to determine if observed genotype frequencies differ from expected and initial 
  print("Fisher's exact test for difference in genotype counts from expected")
  print(suppressWarnings(fisher.test(matrix(c(flp11*N,flp22*N,flp12*N,HWlp11*N,HWlp22*N,HWlp12*N),nrow=3))))
  
  print("Fisher's exact test for difference in genotype counts from initial")
  print(suppressWarnings(fisher.test(matrix(c(flp11*N,flp22*N,flp12*N,p1^2*N,p2^2*N,2*p1*p2*N),nrow=3))))
  
}