#Use file.remove if you would like to delete the previously made .csv doc (so you can use the same name for new file)
# file.remove("C:/Users/tdtibbet/Desktop/RMediationErrorsStudy1.csv")
# file.remove("C:/Users/tdtibbet/Desktop/MacKinnon2004SummaryInformationStudy1.csv")
# file.remove("C:/Users/tdtibbet/Desktop/MacKinnon2004AccuracySummaryStudy1.csv")

#Load required packages
library(RMediation) #RMediation is required for medci() function

start_time <- Sys.time()
# set.seed(134)
set.seed(119)

nsims=10000
alpha=.05
ns=c(50, 100, 200, 500, 1000)
as=c(0, .14, .39, .59, 0, 0, 0, .14, .14, .14, .39, .39, .39, .59, .59, .59)
bs=c(0, 0, 0, 0, .14, .39, .59, .14, .39, .59, .14, .39, .59, .14, .39, .59)

#Make some empty matrices that will fill up later
CIlimits=matrix(nrow=2, ncol=2)
truebelow=matrix(nrow=nsims, ncol=2)
trueabove=matrix(nrow=nsims, ncol=2)
suminfo=matrix(nrow=nsims, ncol=4)
iteration=matrix(nrow=length(ns)*length(as)*nsims, ncol=1)
summaryinformation=matrix(nrow=length(ns)*length(as), ncol=9)
accursum=matrix(nrow=length(ns)*length(as), ncol=8)

#Loop through all combinations of a and b paths for indirect effect
for (h in 1: length(ns)) {
  for (k in 1:length(as)) {
    amat=t(cbind(0, matrix(as[k], nrow=1, ncol=1)))
    cprimebmat=t(cbind(0, 0, matrix(bs[k], nrow=1, ncol=1)))
    #Loop through number of simulations for each combination of a and b paths
    for (i in 1:nsims) {
      #Generate vector of random numbers from standard normal distribution for values of X variable and construct design matrix
      Xs=rnorm(ns[h], 0, 1)
      Xsdm=cbind(matrix(1, nrow=ns[h], ncol=1), Xs)
      #Calculate M for each X variable using corresponding regression equation adding random error and create combined design matrix for M and X
      Ms=(Xsdm%*%amat+rnorm(ns[h], 0, 1))
      XsMsdm=cbind(Xsdm, Ms)
      #Calculate Y for each X and M variable using corresponding regression equation adding random error
      Ys=XsMsdm%*%cprimebmat+rnorm(ns[h], 0, 1)
      #Combine X, M, and Y vectors to form data frame of observations (IV, mediator, and DV for each)
      sampledat=data.frame(X=Xs, M=Ms, Y=Ys)
      #Calculate estimate of a and its std error (by calculating var-cov matrix, Sigma, and taking sqrt() of bottom right entry)
      Mcoefs=solve(t(Xsdm)%*%Xsdm)%*%t(Xsdm)%*%Ms
      a=Mcoefs[2,]
      Mhats=Xsdm%*%Mcoefs
      MSSresid=t(Ms-Mhats)%*%(Ms-Mhats)
      MSigma=(MSSresid/(ns[h]-nrow(Mcoefs)))[1,1]*(solve(t(Xsdm)%*%Xsdm))
      astderr=sqrt(MSigma[2,2])
      #Calculate estimate of b and its std error (by calculating var-cov matrix, Sigma, and taking sqrt() of bottom right entry)
      Ycoefs=solve(t(XsMsdm)%*%XsMsdm)%*%t(XsMsdm)%*%Ys
      b=Ycoefs[3,]
      Yhats=XsMsdm%*%Ycoefs
      YSSresid=t(Ys-Yhats)%*%(Ys-Yhats)
      YSigma=(YSSresid/(ns[h]-nrow(Ycoefs)))[1,1]*(solve(t(XsMsdm)%*%XsMsdm))
      bstderr=sqrt(YSigma[3,3])
      
      #DISTRIB OF THE PRODUCT M CONFIDENCE INTERVAL FOR INDIRECT EFFECT
      #Generate CI using RMediation package
      rmedCI=medci(mu.x=a, se.x=astderr, mu.y=b, se.y=bstderr, rho=0, alpha=alpha, type="dop")
      #If RMediation cannot calculate CI, output values that generated the error to .csv file and redo the iteration (otherwise proceed with sim)
      if (NA %in% rmedCI$`97.5% CI`) {
        write.table(data.frame(n=ns[h], mu.x=a, se.x=astderr, mu.y=b, se.y=bstderr, rho=0, alpha=alpha), 
                    file = "C:/Users/tdtibbet/Desktop/RMediationErrorsStudy1.csv",
                    append=TRUE, sep = ",",
                    quote = FALSE,
                    row.names=FALSE,
                    col.names=FALSE)
        i=i-1
      } else {
        #Calculate estimate of ab and its std error (based on multivariate delta method)
        ab=a*b
        abstderr=sqrt(a^2*bstderr^2+b^2*astderr^2)
        #Reverse engineer critical values using RMediation CIs and plug them into CI formulas using correct std errors
        #Place lower limit in first row of CI matrix and upper limit in second row
        CIlimits[1,2]=ab+((rmedCI$`97.5% CI`[1]-ab)/rmedCI$SE)*abstderr
        CIlimits[2,2]=ab+((rmedCI$`97.5% CI`[2]-ab)/rmedCI$SE)*abstderr
        
        #Z TEST CONFIDENCE INTERVAL FOR INDIRECT EFFECT (use 1.96 for z critical value as original sim did)
        #Place lower limit in first row of CI matrix and upper limit in second row
        CIlimits[1,1]=ab-1.96*abstderr
        CIlimits[2,1]=ab+1.96*abstderr

        #Assess accuracy of CIs by seeing if true indirect effect is below lower limit or above upper limit
        for (j in 1:ncol(CIlimits)) {
          if (CIlimits[1,j]>as[k]*bs[k]) {
            truebelow[i,j]=1
            trueabove[i,j]=0
          } else if (CIlimits[2,j]<as[k]*bs[k]) {
            truebelow[i,j]=0
            trueabove[i,j]=1
          } else {
            truebelow[i,j]=0
            trueabove[i,j]=0
          }
        }
        suminfo[i,1:4]=c(a, b, astderr, bstderr)
        element2=(h-1)*length(as)*nsims+(k-1)*nsims+i
        iteration[element2]=element2
      }
      
    }
    #Calculate average a and b estimates, their empirical std errors, and the average std error estimates of a and b to check for sim errors, then output
    element=(h-1)*length(as)+k
    summaryinformation[element,1]=ns[h]
    summaryinformation[element,2]=as[k]
    summaryinformation[element,3]=bs[k]
    summaryinformation[element,4:5]=colSums(suminfo[,1:2])/nsims
    summaryinformation[element,6:7]=apply(suminfo[,1:2], 2, sd)
    summaryinformation[element,8:9]=colSums(suminfo[,3:4])/nsims
    #Calculate proportion of times true indirect effect falls below/above CIs
    propbelow=colSums(truebelow)/nsims
    propabove=colSums(trueabove)/nsims
    #Output proportions to results matrix
    accursum[element,1]=ns[h]
    accursum[element,2]=as[k]
    accursum[element,3]=bs[k]
    accursum[element,4]=as[k]*bs[k]
    accursum[element,5:6]=propbelow
    accursum[element,7:8]=propabove
  }
}

#Convert summary information matrix to dataframe and write to .csv file
summaryinformation=as.data.frame(summaryinformation)
colnames(summaryinformation)=c("n", "a", "b", "average sample a", "average sample b", 
                               "empirical standard error of sample a", "empirical standard error of sample b",
                               "average sample std error of a", "average sample std error of b")
write.table(summaryinformation, file = "C:/Users/tdtibbet/Desktop/MacKinnon2004SummaryInformationStudy1.csv",
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)

#Convert results matrix to dataframe and write to .csv file
accursum=as.data.frame(accursum)
colnames(accursum)=c("n", "a", "b", "ab", "z proportion below", "M proportion below", 
                          "z proportion above", "M proportion above")
write.table(accursum, file = "C:/Users/tdtibbet/Desktop/MacKinnon2004AccuracySummaryStudy1.csv",
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)

end_time <- Sys.time()

end_time - start_time