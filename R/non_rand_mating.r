#' Simulating non-random mating Mendelian populations 
#'
#' @description \code{non_rand_mating} simulates non-random mating Mendelian populations by specifying allele
#' frequencies, population size, and outcrossing rate.  The simulation generates a parent population with N individuals 
#' with 2N alleles by randomly sampling from a gamete pool based on the specified allele frequencies.  Then parents are randomly selected to mate where the
#' outcrossing rate determines whether they mate with themselves or another randomly selected individual. 
#' @param outcrossingrate probability an individual mates with a randomly selected individual from the parent population
#' @param p1 Allele 1 frequency
#' @param p2 Allele 2 frequency
#' @param alleles allele types (use 1 and 2)
#' @param N Population size
#' @details \code{non_rand_mating} returns the allele frequencies, observed and expected heterozygosities,
#' Wright's F, and Fisher's exact (\link{fisher.test}) tests comparing the observed genotype counts to the expected
#' and initial values.
#' @export
#' @examples
#' non_rand_mating(0.5,0.3,0.7,c(1,2),1000)

non_rand_mating<-function (outcrossingrate,p1,p2,alleles,N){
  
  #Genotypes 
  probs<-c(p1,p2) #vector of allele frequencies, the probability that an allele is represented in a genotype depends on the specified allele frequencies.
  allele1<-sample(alleles,N,replace=TRUE,prob=probs) #selecting the first allele
  allele2<-sample(alleles,N,replace=TRUE,prob=probs) #and the second allele
  genotypes<-data.frame(allele1,allele2) #creating object of genotypes from alleles
  p11<-length(genotypes[genotypes[,1]==1 & genotypes[,2]==1,][,1])/N #calculating genotype 11 freq.
  p12<-sum(length(genotypes[genotypes[,1]==1 & genotypes[,2]==2,][,1]),length(genotypes[genotypes[,1]==2 & genotypes[,2]==1,][,1]))/N #calculating genotype 12 (heterozygote) freq.)
  p22<-length(genotypes[genotypes[,1]==2 & genotypes[,2]==2,][,1])/N #calculating genotype 22 freq.
  
  offspringallele1<-numeric(N)
  offspringallele2<-numeric(N)
  matingevent<-character(N)
  for (i in 1:N) {
    parent<-genotypes[sample(1:nrow(genotypes),1),]
    matingprob<-runif(1)
    matingevent<-if (matingprob>outcrossingrate) {"self"} else {"out"}
    offspringallele1[i]<-sample(parent,1)
    offspringallele2[i]<-if (matingprob>outcrossingrate) {
      sample(parent,1)
    } else {
      genotypes[sample(1:nrow(genotypes),1,replace=TRUE),1]
    }
  }
  OA1<-matrix(offspringallele1)
  OA2<-matrix(offspringallele2)
  output<-data.frame(OA1,OA2,matingevent)
  output$OA1<-as.numeric(output$OA1)
  output$OA2<-as.numeric(output$OA2)
  
  ##calculate F-statistic and test for HW
  #expected Hardy-Weinberg genotype frequencies given allele frequency
  HWnr11<-(((sum(length(output[output[,1]==1,][,1]),length(output[output[,2]==1,][,2])))/(2*N))^2)
  HWnr22<-(((sum(length(output[output[,1]==2,][,1]),length(output[output[,2]==2,][,2])))/(2*N))^2)
  HWnr12<-(((sum(length(output[output[,1]==1,][,1]),length(output[output[,2]==1,][,2])))/(2*N))*((sum(length(output[output[,1]==2,][,1]),length(output[output[,2]==2,][,2])))/(2*N))*2)
  
  #observed genotype frequencies
  fnr11<-length(output[output[,1]==1 & output[,2]==1,][,1])/N
  fnr22<-length(output[output[,1]==2 & output[,2]==2,][,1])/N
  fnr12<-sum(length(output[output[,1]==1 & output[,2]==2,][,1]),length(output[output[,1]==2 & output[,2]==1,][,2]))/N
  
  #####*****RUN ONLY TO THIS POINT FIRST*****#####
  #Genetic diversity
  print("Heterozygosity")
  Fstatisticnr<-1-(fnr12/HWnr12)# Wright's F
  diversity<-c(fnr12,HWnr12,Fstatisticnr);names(diversity)<-c("Ho","He","Wright's F");print(diversity)
  
  #Allele frequency
  freqs<-c(length(output[output==1])/(2*N),length(output[output==2])/(2*N));names(freqs)<-c("p1","p2")
  print("Allele frequencies")
  print(freqs)
  
  #Chi-squared test to determine if observed genotype frequencies differ from expected and initial
  
  print("Fisher's exact test for difference in genotype counts from expected")
  print(suppressWarnings(fisher.test(matrix(c(fnr11*N,fnr22*N,fnr12*N,HWnr11*N,HWnr22*N,HWnr12*N),nrow=3))))
  print("Fisher's exact test for difference in genotype counts from initial")
  print(suppressWarnings(fisher.test(matrix(c(fnr11*N,fnr22*N,fnr12*N,p1^2*N,p2^2*N,2*p1*p2*N),nrow=3))))
  
}