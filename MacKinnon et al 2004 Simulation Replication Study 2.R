#Create file path to use throughout code where you would like to output .csv files to
##Get path directory where this script is saved using rstudioapi package
filepath <- dirname(rstudioapi::getActiveDocumentContext()$path)
#Use file.remove if you would like to delete the previously made .csv doc (so you can use the same name for new file)
# file.remove(paste(filepath,"/RMediationErrorsStudy2.csv", sep =""))
# file.remove(paste(filepath,"/MacKinnon2004SummaryInformationStudy2.csv", sep =""))
# file.remove(paste(filepath,"/MacKinnon200495PercentAccuracySummaryStudy2.csv", sep =""))
# file.remove(paste(filepath,"/MacKinnon200490PercentAccuracySummaryStudy2.csv", sep =""))
# file.remove(paste(filepath,"/MacKinnon200480PercentAccuracySummaryStudy2.csv", sep =""))
# file.remove(paste(filepath,"/MacKinnon2004RejectRateSummaryStudy2.csv", sep =""))

#Load required packages
library(RMediation) #RMediation is required for medci() function

start_time <- Sys.time()
#set.seed(134)
set.seed(1200)

nsims=1000
nboot=1000
nMC=1000
#Set alpha level (used in calculation of M test and jackknife confidence limits)
alpha=c(.05, .10, .20)
#Set confidence level (used in calculation of bootstrap and Monte Carlo confidence limits; see Med, Mod, Condit Process Analysis by Hayes)
confidence=c(95, 90, 80) 
ns=c(25, 50, 100, 200)
as=c(0, 0, 0, 0, .14, .39, .59, .14, .14, .39)
bs=c(0, .14, .39, .59, .14, .39, .59, .39, .59, .59)

#Make some empty matrices that will fill up later
CIlimits95=matrix(nrow=2, ncol=8)
CIlimits90=matrix(nrow=2, ncol=8)
CIlimits80=matrix(nrow=2, ncol=8)
truebelow95=matrix(nrow=nsims, ncol=8)
truebelow90=matrix(nrow=nsims, ncol=8)
truebelow80=matrix(nrow=nsims, ncol=8)
trueabove95=matrix(nrow=nsims, ncol=8)
trueabove90=matrix(nrow=nsims, ncol=8)
trueabove80=matrix(nrow=nsims, ncol=8)

rejectrate=matrix(nrow=nsims, ncol=8)

bootdistrib=matrix(nrow=nboot, ncol=1)
bootTdistrib=matrix(nrow=nboot, ncol=1)

suminfo=matrix(nrow=nsims, ncol=4)
iteration=matrix(nrow=length(ns)*length(as)*nsims, ncol=1)
summaryinformation=matrix(nrow=length(ns)*length(as), ncol=9)
accursum95=matrix(nrow=length(ns)*length(as), ncol=20)
accursum90=matrix(nrow=length(ns)*length(as), ncol=20)
accursum80=matrix(nrow=length(ns)*length(as), ncol=20)
rejectratesum=matrix(nrow=length(ns)*length(as), ncol=12)

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
      #Generate 95% CI using RMediation package
      rmedCI95=medci(mu.x=a, se.x=astderr, mu.y=b, se.y=bstderr, rho=0, alpha=alpha[1], type="dop")
      #If RMediation cannot calculate CI, output values that generated the error to .csv file and redo the iteration
      if (NA %in% rmedCI95$`97.5% CI`) {
        write.table(data.frame(mu.x=a, se.x=astderr, mu.y=b, se.y=bstderr, rho=0, alpha=alpha[1]), 
                    file = paste(filepath,"/RMediationErrorsStudy2.csv", sep =""),
                    append=TRUE, sep = ",",
                    quote = FALSE,
                    row.names=FALSE,
                    col.names=FALSE)
        i=i-1
        #Else, continue the simulation as normal
      } else {
        #Generate 90% and 80% CIs using RMediation package
        rmedCI90=medci(mu.x=a, se.x=astderr, mu.y=b, se.y=bstderr, rho=0, alpha=alpha[2], type="dop")
        rmedCI80=medci(mu.x=a, se.x=astderr, mu.y=b, se.y=bstderr, rho=0, alpha=alpha[3], type="dop")
        #Calculate estimate of ab and its std error (based on multivariate delta method)
        ab=a*b
        abstderr=sqrt(a^2*bstderr^2+b^2*astderr^2)
        #Reverse engineer critical values using RMediation CIs and plug them into CI formulas using correct std errors
        #Place lower limit in first row of CI matrix and upper limit in second row
        CIlimits95[1,2]=ab+((rmedCI95$`97.5% CI`[1]-ab)/rmedCI95$SE)*abstderr
        CIlimits95[2,2]=ab+((rmedCI95$`97.5% CI`[2]-ab)/rmedCI95$SE)*abstderr
        CIlimits90[1,2]=ab+((rmedCI90$`95% CI`[1]-ab)/rmedCI90$SE)*abstderr
        CIlimits90[2,2]=ab+((rmedCI90$`95% CI`[2]-ab)/rmedCI90$SE)*abstderr
        CIlimits80[1,2]=ab+((rmedCI80$`90% CI`[1]-ab)/rmedCI80$SE)*abstderr
        CIlimits80[2,2]=ab+((rmedCI80$`90% CI`[2]-ab)/rmedCI80$SE)*abstderr
        
        #Z TEST CONFIDENCE INTERVAL FOR INDIRECT EFFECT (use 1.96 for z critical value as original sim did)
        #Place lower limit in first row of CI matrix and upper limit in second row
        CIlimits95[1,1]=ab-round(qnorm(p=.975, mean=0, sd=1), digits=3)*abstderr
        CIlimits95[2,1]=ab+round(qnorm(p=.975, mean=0, sd=1), digits=3)*abstderr
        CIlimits90[1,1]=ab-round(qnorm(p=.95, mean=0, sd=1), digits=3)*abstderr
        CIlimits90[2,1]=ab+round(qnorm(p=.95, mean=0, sd=1), digits=3)*abstderr
        CIlimits80[1,1]=ab-round(qnorm(p=.90, mean=0, sd=1), digits=3)*abstderr
        CIlimits80[2,1]=ab+round(qnorm(p=.90, mean=0, sd=1), digits=3)*abstderr
        
        #Begin JACKKNIFE procedure by removing one observation from the sample at a time & calculating ab for each
        jackknifedistrib=matrix(nrow=nrow(sampledat), ncol=1)
        for (j in 1:nrow(sampledat)) {
          jacksampledat=sampledat[-j,]
          #Calculate jackknife sample estimate of a
          jackmedmodel=summary(lm(M~X, data=jacksampledat))
          jacka=jackmedmodel$coefficients[2,1]
          #Calculate jackknife sample estimate of b
          jackdirmodel=summary(lm(Y~X+M, data=jacksampledat))
          jackb=jackdirmodel$coefficients[3,1]
          #Calculate jackknife sample estimate of ab and add to matrix
          jackab=jacka*jackb
          jackknifedistrib[j,1]=jackab
        }
        #Calculate the average jackknife estimate of ab across jackknife samples
        jackavgab=sum(jackknifedistrib)/nrow(jackknifedistrib)
        #Calculate jackknife std error
        jackabstderr=sqrt(((ns[h]-1)/ns[h])*sum((jackknifedistrib-jackavgab)^2))
        #Place lower limit in first row of CI matrix and upper limit in second row
        CIlimits95[1,3]=jackavgab-qnorm(1-alpha[1]/2)*jackabstderr
        CIlimits95[2,3]=jackavgab+qnorm(1-alpha[1]/2)*jackabstderr
        CIlimits90[1,3]=jackavgab-qnorm(1-alpha[2]/2)*jackabstderr
        CIlimits90[2,3]=jackavgab+qnorm(1-alpha[2]/2)*jackabstderr
        CIlimits80[1,3]=jackavgab-qnorm(1-alpha[3]/2)*jackabstderr
        CIlimits80[2,3]=jackavgab+qnorm(1-alpha[3]/2)*jackabstderr
        
        #Begin BOOTSTRAP procedure by resampling with replacement
        for (j in 1:nboot) {
          bootdatamatrix=as.matrix(sampledat[sample(nrow(sampledat), ns[h], replace=TRUE), ])
          #Create design matrices and vectors of Ms and Ys for calculation of bootstrap a and b
          bootXsdm=cbind(matrix(1, nrow=ns[h], ncol=1), bootdatamatrix[,1])
          bootMs=bootdatamatrix[,2]
          bootXsMsdm=cbind(matrix(1, nrow=ns[h], ncol=1), bootdatamatrix[,1:2])
          bootYs=bootdatamatrix[,3]
          #Calculate bootstrap a and its std error (by calculating var-cov matrix, Sigma, and taking sqrt() of bottom right entry)
          bootMcoefs=solve(t(bootXsdm)%*%bootXsdm)%*%t(bootXsdm)%*%bootMs
          boota=bootMcoefs[2,]
          bootMhats=bootXsdm%*%bootMcoefs
          bootMSSresid=t(bootMs-bootMhats)%*%(bootMs-bootMhats)
          bootMSigma=(bootMSSresid/(ns[h]-nrow(bootMcoefs)))[1,1]*(solve(t(bootXsdm)%*%bootXsdm))
          bootastderr=sqrt(bootMSigma[2,2])
          #Calculate bootstrap b and its std error (by calculating var-cov matrix, Sigma, and taking sqrt() of bottom right entry)
          bootYcoefs=solve(t(bootXsMsdm)%*%bootXsMsdm)%*%t(bootXsMsdm)%*%bootYs
          bootb=bootYcoefs[3,]
          bootYhats=bootXsMsdm%*%bootYcoefs
          bootYSSresid=t(bootYs-bootYhats)%*%(bootYs-bootYhats)
          bootYSigma=(bootYSSresid/(ns[h]-nrow(bootYcoefs)))[1,1]*(solve(t(bootXsMsdm)%*%bootXsMsdm))
          bootbstderr=sqrt(bootYSigma[3,3])
          #Calculate bootstrap ab and std error and add to matrix
          bootab=boota*bootb
          bootabstderr=sqrt(boota^2*bootbstderr^2+bootb^2*bootastderr^2)
          bootdistrib[j,1]=bootab
          #Calculate bootstrap t and add to matrix for bootstrap t method
          bootTdistrib[j,1]=(bootab-ab)/bootabstderr
        }
        #Sort the bootstrap sampling distrib from smallest to largest indirect effects
        sortbootdistrib=sort(bootdistrib)
        
        #PERCENTILE BOOTSTRAP
        #Place lower limit in first row of CI matrix and upper limit in second row
        CIlimits95[1,4]=sortbootdistrib[floor(0.005*nboot*(100-confidence[1]))]
        CIlimits95[2,4]=sortbootdistrib[ceiling(nboot*(1-0.005*(100-confidence[1]))+1)]
        CIlimits90[1,4]=sortbootdistrib[floor(0.005*nboot*(100-confidence[2]))]
        CIlimits90[2,4]=sortbootdistrib[ceiling(nboot*(1-0.005*(100-confidence[2]))+1)]
        CIlimits80[1,4]=sortbootdistrib[floor(0.005*nboot*(100-confidence[3]))]
        CIlimits80[2,4]=sortbootdistrib[ceiling(nboot*(1-0.005*(100-confidence[3]))+1)]
        
        #Begin BIAS-CORRECTION procedure
        #Find proportion of bootstrap estimates that fall below original estimate and convert proportion to z-score
        zadj=qnorm(sum(sortbootdistrib<ab)/nboot)
        #Add 2 times proportion's z-score to upper and lower limits of bootstrap CI z-scores, then convert back to percentiles
        adjlo95=pnorm(qnorm((100-confidence[1])/200)+2*zadj)
        adjhi95=pnorm(qnorm(1-((100-confidence[1])/200))+2*zadj)
        adjlo90=pnorm(qnorm((100-confidence[2])/200)+2*zadj)
        adjhi90=pnorm(qnorm(1-((100-confidence[2])/200))+2*zadj)
        adjlo80=pnorm(qnorm((100-confidence[3])/200)+2*zadj)
        adjhi80=pnorm(qnorm(1-((100-confidence[3])/200))+2*zadj)
        #Place lower limit in first row of CI matrix and upper limit in second row
        CIlimits95[1,5]=sortbootdistrib[max(floor(nboot*adjlo95), 1)]
        CIlimits95[2,5]=sortbootdistrib[min(ceiling(nboot*adjhi95), nboot)]
        CIlimits90[1,5]=sortbootdistrib[max(floor(nboot*adjlo90), 1)]
        CIlimits90[2,5]=sortbootdistrib[min(ceiling(nboot*adjhi90), nboot)]
        CIlimits80[1,5]=sortbootdistrib[max(floor(nboot*adjlo80), 1)]
        CIlimits80[2,5]=sortbootdistrib[min(ceiling(nboot*adjhi80), nboot)]
        
        #Begin BOOTSTRAP-T procedure
        #Sort the bootstrap t distrib from smallest to largest
        sortbootTdistrib=sort(bootTdistrib)
        #Place lower limit in first row of CI matrix and upper limit in second row
        CIlimits95[1,6]=ab-sortbootTdistrib[ceiling(nboot*(1-0.005*(100-confidence[1]))+1)]*abstderr
        CIlimits95[2,6]=ab-sortbootTdistrib[floor(0.005*nboot*(100-confidence[1]))]*abstderr
        CIlimits90[1,6]=ab-sortbootTdistrib[ceiling(nboot*(1-0.005*(100-confidence[2]))+1)]*abstderr
        CIlimits90[2,6]=ab-sortbootTdistrib[floor(0.005*nboot*(100-confidence[2]))]*abstderr
        CIlimits80[1,6]=ab-sortbootTdistrib[ceiling(nboot*(1-0.005*(100-confidence[3]))+1)]*abstderr
        CIlimits80[2,6]=ab-sortbootTdistrib[floor(0.005*nboot*(100-confidence[3]))]*abstderr
        
        #Begin BOOTSTRAP-Q procedure
        #Calculate skewness of bootstrap t distribution using formulae in Manly (1997) on pg. 59
        bootTvariance=sum((bootTdistrib-sum(bootTdistrib)/nboot)^2)/nboot
        bootTskewness=sum((bootTdistrib-sum(bootTdistrib)/nboot)^3)/(bootTvariance*sqrt(bootTvariance))
        #Calculate the Q statistic using formula in Manly (1997) on pg. 59
        bootQdistrib=bootTdistrib+((bootTskewness*bootTdistrib^2)/3)+((bootTskewness^2*bootTdistrib^3)/27)+bootTskewness/(6*nboot)
        #Sort the bootstrap Q distrib from smallest to largest
        sortbootQdistrib=sort(bootQdistrib)
        #Find and then transform critical values of Q distribution using formula in Manly (1997) on pg. 60 
        transformQCVlower95=(3*(sign(1+bootTskewness*(sortbootQdistrib[floor(0.005*nboot*(100-confidence[1]))]-bootTskewness/(6*nboot)))*abs(1+bootTskewness*(sortbootQdistrib[floor(0.005*nboot*(100-confidence[1]))]-bootTskewness/(6*nboot)))^(1/3)-1))/bootTskewness
        transformQCVupper95=(3*(sign(1+bootTskewness*(sortbootQdistrib[ceiling(nboot*(1-0.005*(100-confidence[1]))+1)]-bootTskewness/(6*nboot)))*abs(1+bootTskewness*(sortbootQdistrib[ceiling(nboot*(1-0.005*(100-confidence[1]))+1)]-bootTskewness/(6*nboot)))^(1/3)-1))/bootTskewness
        transformQCVlower90=(3*(sign(1+bootTskewness*(sortbootQdistrib[floor(0.005*nboot*(100-confidence[2]))]-bootTskewness/(6*nboot)))*abs(1+bootTskewness*(sortbootQdistrib[floor(0.005*nboot*(100-confidence[2]))]-bootTskewness/(6*nboot)))^(1/3)-1))/bootTskewness
        transformQCVupper90=(3*(sign(1+bootTskewness*(sortbootQdistrib[ceiling(nboot*(1-0.005*(100-confidence[2]))+1)]-bootTskewness/(6*nboot)))*abs(1+bootTskewness*(sortbootQdistrib[ceiling(nboot*(1-0.005*(100-confidence[2]))+1)]-bootTskewness/(6*nboot)))^(1/3)-1))/bootTskewness
        transformQCVlower80=(3*(sign(1+bootTskewness*(sortbootQdistrib[floor(0.005*nboot*(100-confidence[3]))]-bootTskewness/(6*nboot)))*abs(1+bootTskewness*(sortbootQdistrib[floor(0.005*nboot*(100-confidence[3]))]-bootTskewness/(6*nboot)))^(1/3)-1))/bootTskewness
        transformQCVupper80=(3*(sign(1+bootTskewness*(sortbootQdistrib[ceiling(nboot*(1-0.005*(100-confidence[3]))+1)]-bootTskewness/(6*nboot)))*abs(1+bootTskewness*(sortbootQdistrib[ceiling(nboot*(1-0.005*(100-confidence[3]))+1)]-bootTskewness/(6*nboot)))^(1/3)-1))/bootTskewness
        #Place lower limit in first row of CI matrix and upper limit in second row
        CIlimits95[1,7]=ab-transformQCVupper95*abstderr
        CIlimits95[2,7]=ab-transformQCVlower95*abstderr
        CIlimits90[1,7]=ab-transformQCVupper90*abstderr
        CIlimits90[2,7]=ab-transformQCVlower90*abstderr
        CIlimits80[1,7]=ab-transformQCVupper80*abstderr
        CIlimits80[2,7]=ab-transformQCVlower80*abstderr
        
        #Begin MONTE CARLO procedure
        #Generate random samples of a and b from normal distributions with mean and std dev equal to sample values
        MCa=rnorm(n=nMC, mean=a, sd=astderr)
        MCb=rnorm(n=nMC, mean=b, sd=bstderr)
        #Calculate MC distribution of ab 
        MCdistrib=MCa*MCb
        #Sort the MC distrib from smallest to largest
        sortMCdistrib=sort(MCdistrib)
        #Place lower limit in first row of CI matrix and upper limit in second row
        CIlimits95[1,8]=sortMCdistrib[floor(0.005*nMC*(100-confidence[1]))]
        CIlimits95[2,8]=sortMCdistrib[ceiling(nMC*(1-0.005*(100-confidence[1]))+1)]
        CIlimits90[1,8]=sortMCdistrib[floor(0.005*nMC*(100-confidence[2]))]
        CIlimits90[2,8]=sortMCdistrib[ceiling(nMC*(1-0.005*(100-confidence[2]))+1)]
        CIlimits80[1,8]=sortMCdistrib[floor(0.005*nMC*(100-confidence[3]))]
        CIlimits80[2,8]=sortMCdistrib[ceiling(nMC*(1-0.005*(100-confidence[3]))+1)]
        
        #Assess accuracy of CIs (by seeing if true indirect effect is below lower limit or above upper limit)
        ##and assess rejection rate of CIs (by seeing if lower limit is greater than 0 or upper limit is less than 0)
        for (j in 1:ncol(CIlimits95)) {
          if (CIlimits95[1,j]>as[k]*bs[k]) {
            truebelow95[i,j]=1
            trueabove95[i,j]=0
            truebelow90[i,j]=1
            trueabove90[i,j]=0
            truebelow80[i,j]=1
            trueabove80[i,j]=0
          } else if (CIlimits95[2,j]<as[k]*bs[k]) {
            truebelow95[i,j]=0
            trueabove95[i,j]=1
            truebelow90[i,j]=0
            trueabove90[i,j]=1
            truebelow80[i,j]=0
            trueabove80[i,j]=1
          } else if (CIlimits90[1,j]>as[k]*bs[k]) {
            truebelow95[i,j]=0
            trueabove95[i,j]=0
            truebelow90[i,j]=1
            trueabove90[i,j]=0
            truebelow80[i,j]=1
            trueabove80[i,j]=0
          } else if (CIlimits90[2,j]<as[k]*bs[k]) {
            truebelow95[i,j]=0
            trueabove95[i,j]=0
            truebelow90[i,j]=0
            trueabove90[i,j]=1
            truebelow80[i,j]=0
            trueabove80[i,j]=1
          } else if (CIlimits80[1,j]>as[k]*bs[k]) {
            truebelow95[i,j]=0
            trueabove95[i,j]=0
            truebelow90[i,j]=0
            trueabove90[i,j]=0
            truebelow80[i,j]=1
            trueabove80[i,j]=0
          } else if (CIlimits80[2,j]<as[k]*bs[k]) {
            truebelow95[i,j]=0
            trueabove95[i,j]=0
            truebelow90[i,j]=0
            trueabove90[i,j]=0
            truebelow80[i,j]=0
            trueabove80[i,j]=1
          } else {
            truebelow95[i,j]=0
            trueabove95[i,j]=0
            truebelow90[i,j]=0
            trueabove90[i,j]=0
            truebelow80[i,j]=0
            trueabove80[i,j]=0
          }
          if (CIlimits95[1,j]>0 | CIlimits95[2,j]<0) {
            rejectrate[i,j]=1
          } else {
            rejectrate[i,j]=0
          }
        }
        #Store summary information to make sure sim is running properly
        suminfo[i,1:4]=c(a, b, astderr, bstderr)
        #Store each iteration number to make sure no iterations are skipped/repeated due to RMediation loop 
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
    propbelow95=colSums(truebelow95)/nsims
    propabove95=colSums(trueabove95)/nsims
    propbelow90=colSums(truebelow90)/nsims
    propabove90=colSums(trueabove90)/nsims
    propbelow80=colSums(truebelow80)/nsims
    propabove80=colSums(trueabove80)/nsims
    #Output proportions to results matrices
    element=(h-1)*length(as)+k
    accursum95[element,1]=ns[h]
    accursum95[element,2]=as[k]
    accursum95[element,3]=bs[k]
    accursum95[element,4]=as[k]*bs[k]
    accursum95[element,5:12]=propbelow95
    accursum95[element,13:20]=propabove95
    accursum90[element,1]=ns[h]
    accursum90[element,2]=as[k]
    accursum90[element,3]=bs[k]
    accursum90[element,4]=as[k]*bs[k]
    accursum90[element,5:12]=propbelow90
    accursum90[element,13:20]=propabove90
    accursum80[element,1]=ns[h]
    accursum80[element,2]=as[k]
    accursum80[element,3]=bs[k]
    accursum80[element,4]=as[k]*bs[k]
    accursum80[element,5:12]=propbelow80
    accursum80[element,13:20]=propabove80
    
    #Calculate rejection rates
    rejectionrate=colSums(rejectrate)/nsims
    #Output a effect, b effect, average ab effect, and rejection rates to results matrices
    rejectratesum[element,1]=ns[h]
    rejectratesum[element,2]=as[k]
    rejectratesum[element,3]=bs[k]
    rejectratesum[element,4]=as[k]*bs[k]
    rejectratesum[element,5:12]=rejectionrate
  }
}

#Convert summary information matrix to dataframe and write to .csv file
summaryinformation=as.data.frame(summaryinformation)
colnames(summaryinformation)=c("n", "a", "b", "average sample a", "average sample b", 
                               "empirical standard error of sample a", "empirical standard error of sample b",
                               "average sample std error of a", "average sample std error of b")
write.table(summaryinformation, file = paste(filepath,"/MacKinnon2004SummaryInformationStudy2.csv", sep =""),
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)

#Convert results matrices to dataframes and write to .csv files
accursum95=as.data.frame(accursum95)
colnames(accursum95)=c("n", "a", "b", "ab", 
                       "z proportion below", "M proportion below", "jackknife proportion below", "percentile boot proportion below",
                       "bias-corrected boot proportion below", "boot-t proportion below", "boot-Q proportion below", "Monte Carlo proportion below",
                       "z proportion above", "M proportion above", "jackknife proportion above", "percentile boot proportion above",
                       "bias-corrected boot proportion above", "boot-t proportion above", "boot-Q proportion above", "Monte Carlo proportion above")
write.table(accursum95, file = paste(filepath,"/MacKinnon200495PercentAccuracySummaryStudy2.csv", sep =""),
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)

accursum90=as.data.frame(accursum90)
colnames(accursum90)=c("n", "a", "b", "ab", 
                       "z proportion below", "M proportion below", "jackknife proportion below", "percentile boot proportion below",
                       "bias-corrected boot proportion below", "boot-t proportion below", "boot-Q proportion below", "Monte Carlo proportion below",
                       "z proportion above", "M proportion above", "jackknife proportion above", "percentile boot proportion above",
                       "bias-corrected boot proportion above", "boot-t proportion above", "boot-Q proportion above", "Monte Carlo proportion above")
write.table(accursum90, file = paste(filepath,"/MacKinnon200490PercentAccuracySummaryStudy2.csv", sep =""),
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)

accursum95=as.data.frame(accursum80)
colnames(accursum80)=c("n", "a", "b", "ab", 
                       "z proportion below", "M proportion below", "jackknife proportion below", "percentile boot proportion below",
                       "bias-corrected boot proportion below", "boot-t proportion below", "boot-Q proportion below", "Monte Carlo proportion below",
                       "z proportion above", "M proportion above", "jackknife proportion above", "percentile boot proportion above",
                       "bias-corrected boot proportion above", "boot-t proportion above", "boot-Q proportion above", "Monte Carlo proportion above")
write.table(accursum80, file = paste(filepath,"/MacKinnon200480PercentAccuracySummaryStudy2.csv", sep =""),
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)

rejectratesum=as.data.frame(rejectratesum)
colnames(rejectratesum)=c("n", "a", "b", "ab", 
                          "z reject rate", "M reject rate", "jackknife reject rate", "percentile boot reject rate",
                          "bias-corrected boot reject rate", "boot-t reject rate", "boot-Q reject rate", "Monte Carlo reject rate")
write.table(rejectratesum, file = paste(filepath,"/MacKinnon2004RejectRateSummaryStudy2.csv", sep =""),
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)


end_time <- Sys.time()

end_time - start_time
