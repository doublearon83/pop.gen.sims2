#' Simulating single locus, two allele Mendelian populations with
#' finite population sizes and natural selection over many generations
#'
#' @description \code{finite_selection} simulates finite Mendelian populations with selection over many generations.
#' Specify initial allele frequencies, population size, number of generations, selection coefficient, dominance coefficient, and number of replicates.
#' The simulation produces each generation of the population by randomly selecting 2N (ND) alleles for N (ND) individuals from a gamete pool based on
#' the allele frequencies as determined by the relative fitnesses of the genotypes.
#' @param s Selection coefficient
#' @param h Dominance coefficient
#' @param p1g2_nr Allele 1 frequency
#' @param p2g2_nr Allele 2 frequency
#' @param alleles Allele types (use 1 and 2)
#' @param ND_nr Population size
#' @param generations_nr Number of generations
#' @param n_reps Number simulation replicates
#' @details \code{finite_selection} returns three plots:
#' 1) Allele 1 frequency versus number of generations, 2) Observed and expected heterozygosity (expected based on
#' allele frequencies and drift) versus number of generations, and 3) Wright's F with a \link{loess} regression best fit line.
#' It also returns Fisher's exact tests (\link{fisher.test}) comparing the observed genotype counts to the expected
#' and initial values.
#' @export
#' @examples
#' finite_selection(0.5,0.5,0.3,0.7,c(1,2),25,100,100)

finite_selection<-function(s,h,p1g2_nr,p2g2_nr,alleles,ND_nr,generations_nr,n_rep){

  #genotype simulation averages over replicates
  Ho_f_nr_rep<-matrix(0,nrow=generations_nr,ncol=n_rep)
  He_f_nr_rep<-matrix(0,nrow=generations_nr,ncol=n_rep)
  Fstatisticf_nr_rep<-matrix(0,nrow=generations_nr,ncol=n_rep)
  ffA1_nr_rep<-matrix(0,nrow=generations_nr,ncol=n_rep)
  ffA2_nr_rep<-matrix(0,nrow=generations_nr,ncol=n_rep)
  p1g2_nr2<-p1g2_nr;p2g2_nr2<-p2g2_nr

  for (jj in 1:n_rep) {

    #genotype simulation
    p1g2_nr1<-p1g2_nr
    p2g2_nr1<-p2g2_nr
    Ho_f_nr<-numeric(generations_nr)
    He_f_nr<-numeric(generations_nr)
    Fstatisticf_nr<-numeric(generations_nr)
    ffA1_nr<-numeric(generations_nr)
    ffA2_nr<-numeric(generations_nr)
    for (i in 1:generations_nr) {

      g2A1_nr<-numeric(ND_nr)
      g2A2_nr<-numeric(ND_nr)
      p1g2_nr1<-p1g2_nr1+((s*p1g2_nr1*p2g2_nr1)*(h*(1-2*p2g2_nr1)+p2g2_nr1))/(1-(2*p1g2_nr1*p2g2_nr1*h*s)-(s*p2g2_nr1^2))
      p2g2_nr1<-1-p1g2_nr1
      for (j in 1:ND_nr) {
        g2A1_nr[j]<-sample(alleles,1,replace=TRUE,prob=c(p1g2_nr1,p2g2_nr1))
        g2A2_nr[j]<-sample(alleles,1,replace=TRUE,prob=c(p1g2_nr1,p2g2_nr1))
      }
      g2AO1_nr<-matrix(g2A1_nr)
      g2AO2_nr<-matrix(g2A2_nr)
      g2_nr<-data.frame(g2AO1_nr,g2AO2_nr)
      p1g2_nr1<-length(g2_nr[g2_nr==1])/(2*ND_nr)
      p2g2_nr1<-length(g2_nr[g2_nr==2])/(2*ND_nr)

      ##calculate F-statistic and test for HW
      #expected Hardy-Weinberg genotype frequencies given allele frequency
      HWf11_nr<-(((sum(length(g2_nr[g2_nr[,1]==1,][,1]),length(g2_nr[g2_nr[,2]==1,][,2])))/(2*ND_nr))^2)
      HWf22_nr<-(((sum(length(g2_nr[g2_nr[,1]==2,][,1]),length(g2_nr[g2_nr[,2]==2,][,2])))/(2*ND_nr))^2)
      HWf12_nr<-(((sum(length(g2_nr[g2_nr[,1]==1,][,1]),length(g2_nr[g2_nr[,2]==1,][,2])))/(2*ND_nr))*((sum(length(g2_nr[g2_nr[,1]==2,][,1]),length(g2_nr[g2_nr[,2]==2,][,2])))/(2*ND_nr))*2)

      #observed genotype frequencies
      ff11_nr<-length(g2_nr[g2_nr[,1]==1 & g2_nr[,2]==1,][,1])/ND_nr
      ff22_nr<-length(g2_nr[g2_nr[,1]==2 & g2_nr[,2]==2,][,1])/ND_nr
      ff12_nr<-sum(length(g2_nr[g2_nr[,1]==1 & g2_nr[,2]==2,][,1]),length(g2_nr[g2_nr[,1]==2 & g2_nr[,2]==1,][,2]))/ND_nr

      #Wright's F
      Ho_f_nr[i]<-ff12_nr
      He_f_nr[i]<-HWf12_nr
      Fstatisticf_nr[i]<-1-(ff12_nr/HWf12_nr)
      #Ho_f_nr;He_f_nr;Fstatisticf_nr

      #allele frequencies
      ffA1_nr[i]<-length(g2_nr[g2_nr==1])/(2*ND_nr)
      ffA2_nr[i]<-length(g2_nr[g2_nr==2])/(2*ND_nr)
    }

    Ho_f_nr_rep[,jj]<-Ho_f_nr
    He_f_nr_rep[,jj]<-He_f_nr
    Fstatisticf_nr_rep[,jj]<-Fstatisticf_nr
    ffA1_nr_rep[,jj]<-ffA1_nr
    ffA2_nr_rep[,jj]<-ffA2_nr
    print(jj)
    if (jj==n_rep) {print("Simulation complete")}
  }

  #####*****RUN ONLY TO THIS POINT FIRST*****#####
  #Chi-squared test to determine if observed genotype frequencies differ from expected and initial
  print("Fisher's exact test for difference in genotype counts from expected")
  print(suppressWarnings(fisher.test(matrix(c(ff11_nr*ND_nr,ff22_nr*ND_nr,ff12_nr*ND_nr,HWf11_nr*ND_nr,HWf22_nr*ND_nr,HWf12_nr*ND_nr)
                                            ,nrow=3))))

  print("Fisher's exact test for difference in genotype counts from initial")
  print(suppressWarnings(fisher.test(matrix(c(ff11_nr*ND_nr,ff22_nr*ND_nr,ff12_nr*ND_nr,p1g2_nr2^2*ND_nr,p2g2_nr2^2*ND_nr,2*p1g2_nr2*p2g2_nr2*ND_nr)
                                            ,nrow=3))))

  #Plots
  par( mfrow = c( 2, 2 ) )
  plot(c(1:generations_nr),rowMeans(ffA1_nr_rep,na.rm=T),type="l",ylab="A1 freq.",xlab="Generations",lwd=2)#plot of allele 1 over time
  plot(c(1:generations_nr),rowMeans(Ho_f_nr_rep,na.rm=T),type="l",ylab="He and Ho",xlab="Generations",lwd=2, ylim=c(0,max(rowMeans(He_f_nr_rep,na.rm=T))*2))#plot of He over time
  points(c(1:generations_nr),rowMeans(He_f_nr_rep,na.rm=T),type="l",lwd=2,lty=1, col="lightgray")#plot of Ho over time
  points(c(1:generations_nr),(p1g2_nr2*p2g2_nr2*2)*(1-(1/(2*ND_nr)))^(1:generations_nr),type="l",lwd=2,lty=2,col="red")
  legend("topright",c("Ho","He","He - drift"),lty=c(1,1,2),col=c("black","lightgray", "red"),cex=1)
  plot(c(1:generations_nr),rowMeans(Fstatisticf_nr_rep,na.rm=T),type="l",ylab="Wright's F",xlab="Generations",lwd=2)#plot of F over time
  points(predict(loess(rowMeans(Fstatisticf_nr_rep,na.rm=T)~c(1:generations_nr))),type="l",col="red",lwd=2)
  legend("bottomright",c("F","Loess line"),lty=c(1,1),col=c("black","red"),cex=1)

}
