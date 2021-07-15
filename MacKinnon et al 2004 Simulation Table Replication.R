##############################Table 1##############################
#Read in .csv file
study1=read.csv("C:/Users/tdtibbet/Desktop/MacKinnon2004AccuracySummaryStudy1.csv")
#Keep only cases where a path is zero
study1out=study1[study1$a==0,]
#Remove ab column
study1out=study1out[,-4]

library(tidyverse)
#Combine proportions above and below z CI into a single var
wide = unite(study1out,  col = z, 4, 6, sep=" ", remove =TRUE) 
#Combine proportions above and below M CI into a single var
wide = unite(wide,  col = M, 5, 6, sep=" ", remove =TRUE) 
#Convert to wide table format to get z and M vars for each sample size
wide = pivot_wider(wide, names_from = n, values_from = c(z,M))
#Convert back to long table format and split vars to get z-M indicator var
long = pivot_longer(wide, cols = z_50:M_1000,
                          names_to = c("type","n"), 
                          names_sep = "_", 
                          values_to = "bounds")
#Convert one last time to wide table format for final format
wide2=pivot_wider(long, names_from=n, values_from=bounds)
#Split lower and upper proportions into separate columns
wide2=separate(wide2,
               col='50',
               sep=' ',
               into=c("50 lower", "50 upper"))
wide2=separate(wide2,
               col='100',
               sep=' ',
               into=c("100 lower", "100 upper"))
wide2=separate(wide2,
               col='200',
               sep=' ',
               into=c("200 lower", "200 upper"))
wide2=separate(wide2,
               col='500',
               sep=' ',
               into=c("500 lower", "500 upper"))
wide2=separate(wide2,
               col='1000',
               sep=' ',
               into=c("1000 lower", "1000 upper"))
#Calculate average rows for z and M at bottom of table
wide3=rbind(wide2, setNames(data.frame('ave', 'ave', 'z',
                                       sum(study1out[study1out$n==50,4])/nrow(study1out[study1out$n==50,]),
                                       sum(study1out[study1out$n==50,6])/nrow(study1out[study1out$n==50,]),
                                       sum(study1out[study1out$n==100,4])/nrow(study1out[study1out$n==100,]),
                                       sum(study1out[study1out$n==100,6])/nrow(study1out[study1out$n==100,]),
                                       sum(study1out[study1out$n==200,4])/nrow(study1out[study1out$n==200,]),
                                       sum(study1out[study1out$n==200,6])/nrow(study1out[study1out$n==200,]),
                                       sum(study1out[study1out$n==500,4])/nrow(study1out[study1out$n==500,]),
                                       sum(study1out[study1out$n==500,6])/nrow(study1out[study1out$n==500,]),
                                       sum(study1out[study1out$n==1000,4])/nrow(study1out[study1out$n==1000,]),
                                       sum(study1out[study1out$n==1000,6])/nrow(study1out[study1out$n==1000,])),
                                       names(wide2))
                                        )
wide3=rbind(wide3, setNames(data.frame('ave', 'ave', 'M',
                                       sum(study1out[study1out$n==50,5])/nrow(study1out[study1out$n==50,]),
                                       sum(study1out[study1out$n==50,7])/nrow(study1out[study1out$n==50,]),
                                       sum(study1out[study1out$n==100,5])/nrow(study1out[study1out$n==100,]),
                                       sum(study1out[study1out$n==100,7])/nrow(study1out[study1out$n==100,]),
                                       sum(study1out[study1out$n==200,5])/nrow(study1out[study1out$n==200,]),
                                       sum(study1out[study1out$n==200,7])/nrow(study1out[study1out$n==200,]),
                                       sum(study1out[study1out$n==500,5])/nrow(study1out[study1out$n==500,]),
                                       sum(study1out[study1out$n==500,7])/nrow(study1out[study1out$n==500,]),
                                       sum(study1out[study1out$n==1000,5])/nrow(study1out[study1out$n==1000,]),
                                       sum(study1out[study1out$n==1000,7])/nrow(study1out[study1out$n==1000,])),
                            names(wide2))
)
#Write table to .csv file
write.table(wide3, file = "C:/Users/tdtibbet/Desktop/MacKinnon2004Table1.csv",
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)








##############################Table 2##############################
#Read in .csv file
study1=read.csv("C:/Users/tdtibbet/Desktop/MacKinnon2004AccuracySummaryStudy1.csv")
#Keep only nonzero cases from table in paper
study1out=study1[study1$a!=0 & study1$b!=0 & study1$a <= study1$b,]
#Remove ab column
study1out=study1out[,-4]

library(tidyverse)
#Combine proportions above and below z CI into a single var
wide = unite(study1out,  col = z, 4, 6, sep=" ", remove =TRUE) 
#Combine proportions above and below M CI into a single var
wide = unite(wide,  col = M, 5, 6, sep=" ", remove =TRUE) 
#Convert to wide table format to get z and M vars for each sample size
wide = pivot_wider(wide, names_from = n, values_from = c(z,M))
#Convert back to long table format and split vars to get z-M indicator var
long = pivot_longer(wide, cols = z_50:M_1000,
                    names_to = c("type","n"), 
                    names_sep = "_", 
                    values_to = "bounds")
#Convert one last time to wide table format for final format
wide2=pivot_wider(long, names_from=n, values_from=bounds)
#Split lower and upper proportions into separate columns
wide2=separate(wide2,
               col='50',
               sep=' ',
               into=c("50 lower", "50 upper"))
wide2=separate(wide2,
               col='100',
               sep=' ',
               into=c("100 lower", "100 upper"))
wide2=separate(wide2,
               col='200',
               sep=' ',
               into=c("200 lower", "200 upper"))
wide2=separate(wide2,
               col='500',
               sep=' ',
               into=c("500 lower", "500 upper"))
wide2=separate(wide2,
               col='1000',
               sep=' ',
               into=c("1000 lower", "1000 upper"))
#Calculate average rows for z and M at bottom of table
wide3=rbind(wide2, setNames(data.frame('ave', 'ave', 'z',
                                       sum(study1out[study1out$n==50,4])/nrow(study1out[study1out$n==50,]),
                                       sum(study1out[study1out$n==50,6])/nrow(study1out[study1out$n==50,]),
                                       sum(study1out[study1out$n==100,4])/nrow(study1out[study1out$n==100,]),
                                       sum(study1out[study1out$n==100,6])/nrow(study1out[study1out$n==100,]),
                                       sum(study1out[study1out$n==200,4])/nrow(study1out[study1out$n==200,]),
                                       sum(study1out[study1out$n==200,6])/nrow(study1out[study1out$n==200,]),
                                       sum(study1out[study1out$n==500,4])/nrow(study1out[study1out$n==500,]),
                                       sum(study1out[study1out$n==500,6])/nrow(study1out[study1out$n==500,]),
                                       sum(study1out[study1out$n==1000,4])/nrow(study1out[study1out$n==1000,]),
                                       sum(study1out[study1out$n==1000,6])/nrow(study1out[study1out$n==1000,])),
                            names(wide2))
)
wide3=rbind(wide3, setNames(data.frame('ave', 'ave', 'M',
                                       sum(study1out[study1out$n==50,5])/nrow(study1out[study1out$n==50,]),
                                       sum(study1out[study1out$n==50,7])/nrow(study1out[study1out$n==50,]),
                                       sum(study1out[study1out$n==100,5])/nrow(study1out[study1out$n==100,]),
                                       sum(study1out[study1out$n==100,7])/nrow(study1out[study1out$n==100,]),
                                       sum(study1out[study1out$n==200,5])/nrow(study1out[study1out$n==200,]),
                                       sum(study1out[study1out$n==200,7])/nrow(study1out[study1out$n==200,]),
                                       sum(study1out[study1out$n==500,5])/nrow(study1out[study1out$n==500,]),
                                       sum(study1out[study1out$n==500,7])/nrow(study1out[study1out$n==500,]),
                                       sum(study1out[study1out$n==1000,5])/nrow(study1out[study1out$n==1000,]),
                                       sum(study1out[study1out$n==1000,7])/nrow(study1out[study1out$n==1000,])),
                            names(wide2))
)
#Write table to .csv file
write.table(wide3, file = "C:/Users/tdtibbet/Desktop/MacKinnon2004Table2.csv",
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)










##############################Table 3##############################
#Read in .csv file
study2=read.csv("C:/Users/tdtibbet/Desktop/MacKinnon200495PercentCIAccuracySummaryStudy2.csv")
#Create dataframe that contains indicator var (=0 if ab is 0 and =1 if ab is nonzero)
study2out=study2
study2out$abind=ifelse(study2out$ab==0,0,1)
#Remove a, b, and ab columns
study2out=study2out[,-2:-4]

library(tidyverse)
#Average lower proportions and upper proportions across combinations of sample size and null/nonzero indirect effects
study2out1=study2out
study2out1=unite(study2out,  col = abind_n, 18, 1, sep=" ", remove =TRUE) 
study2out1 =  study2out1 %>% 
              group_by(abind_n) %>% 
              summarise(z.below.ave = mean(z.proportion.below),
                        M.below.ave = mean(M.proportion.below),
                        jackknife.below.ave = mean(jackknife.proportion.below),
                        percentile.below.ave = mean(percentile.boot.proportion.below),
                        bias.corrected.below.ave = mean(bias.corrected.boot.proportion.below),
                        boot.t.below.ave = mean(boot.t.proportion.below),
                        boot.Q.below.ave = mean(boot.Q.proportion.below),
                        Monte.Carlo.below.ave = mean(Monte.Carlo.proportion.below),
                        z.above.ave = mean(z.proportion.above),
                        M.above.ave = mean(M.proportion.above),
                        jackknife.above.ave = mean(jackknife.proportion.above),
                        percentile.above.ave = mean(percentile.boot.proportion.above),
                        bias.corrected.above.ave = mean(bias.corrected.boot.proportion.above),
                        boot.t.above.ave = mean(boot.t.proportion.above),
                        boot.Q.above.ave = mean(boot.Q.proportion.above),
                        Monte.Carlo.above.ave = mean(Monte.Carlo.proportion.above)
                        )

#Combine proportions above and below z CI into a single var
wide = unite(study2out1,  col = z, 2, 10, sep=" ", remove =TRUE) 
#Combine proportions above and below M CI into a single var
wide = unite(wide,  col = M, 3, 10, sep=" ", remove =TRUE) 
#Combine proportions above and below Jackknife CI into a single var
wide = unite(wide,  col = jackknife, 4, 10, sep=" ", remove =TRUE) 
#Combine proportions above and below Percentile Bootstrap CI into a single var
wide = unite(wide,  col = percentile, 5, 10, sep=" ", remove =TRUE) 
#Combine proportions above and below Bias-Corrected Boostrap CI into a single var
wide = unite(wide,  col = bias.corrected, 6, 10, sep=" ", remove =TRUE) 
#Combine proportions above and below Bootstrap T CI into a single var
wide = unite(wide,  col = boot.t, 7, 10, sep=" ", remove =TRUE) 
#Combine proportions above and below Bootstrap Q CI into a single var
wide = unite(wide,  col = boot.Q, 8, 10, sep=" ", remove =TRUE)
#Combine proportions above and below Monte Carlo CI into a single var
wide = unite(wide,  col = Monte.Carlo, 9, 10, sep=" ", remove =TRUE) 

#Separate the ab indicator/sample size variable into two variables
wide=separate(wide,
               col='abind_n',
               sep=' ',
               into=c("abind", "n"))

#Convert to wide table format to get z, M, etc. vars for each sample size
wide = pivot_wider(wide, names_from = n, values_from = c(z, M, jackknife, percentile, bias.corrected, boot.t, boot.Q, Monte.Carlo))
#Convert back to long table format and split vars to get z, M, etc. indicator var
long = pivot_longer(wide, cols = 2:33,
                    names_to = c("type","n"), 
                    names_sep = "_", 
                    values_to = "bounds")
#Convert one last time to wide table format for final format
wide2=pivot_wider(long, names_from=n, values_from=bounds)
#Split lower and upper proportions into separate columns
wide2=separate(wide2,
               col='25',
               sep=' ',
               into=c("25 lower", "25 upper"))
wide2=separate(wide2,
               col='50',
               sep=' ',
               into=c("50 lower", "50 upper"))
wide2=separate(wide2,
               col='100',
               sep=' ',
               into=c("100 lower", "100 upper"))
wide2=separate(wide2,
               col='200',
               sep=' ',
               into=c("200 lower", "200 upper"))

#Reorder columns so smaller sample sizes come first
wide2=wide2 %>% 
  select(abind, type, `25 lower`, `25 upper`, `50 lower`, `50 upper`, `100 lower`, `100 upper`, `200 lower`, `200 upper`)

#Write table to .csv file
write.table(wide2, file = "C:/Users/tdtibbet/Desktop/MacKinnon2004Table3.csv",
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)









##############################Table 4##############################
#Read in .csv files for 80%, 90%, and 95% CIs
study2.80=read.csv("C:/Users/tdtibbet/Desktop/MacKinnon200480PercentCIAccuracySummaryStudy2.csv")
study2.90=read.csv("C:/Users/tdtibbet/Desktop/MacKinnon200490PercentCIAccuracySummaryStudy2.csv")
study2.95=read.csv("C:/Users/tdtibbet/Desktop/MacKinnon200495PercentCIAccuracySummaryStudy2.csv")


#Remove n, a, b, and ab columns
study2.80=study2.80[,-1:-4]
study2.90=study2.90[,-1:-4]
study2.95=study2.95[,-1:-4]

#count the cases in which proportion outside confidence limits exceeds Bradley's liberal robustness criterion
numout.80=matrix(nrow=nrow(study2.80), ncol=ncol(study2.80))
for (j in 1:ncol(study2.80)) {
  for (i in 1:nrow(study2.80)) {
    if (study2.80[i,j]<.05 | study2.80[i,j]>.15) {
      numout.80[i,j]=1
    } else {
      numout.80[i,j]=0
    }
  }
}

numoutsum.80=t(colSums(numout.80))
colnames(numoutsum.80)=colnames(study2.80)
numoutsum.80=as.data.frame(numoutsum.80)

numout.90=matrix(nrow=nrow(study2.90), ncol=ncol(study2.90))
for (j in 1:ncol(study2.90)) {
  for (i in 1:nrow(study2.90)) {
    if (study2.90[i,j]<.025 | study2.90[i,j]>.075) {
      numout.90[i,j]=1
    } else {
      numout.90[i,j]=0
    }
  }
}

numoutsum.90=t(colSums(numout.90))
colnames(numoutsum.90)=colnames(study2.90)
numoutsum.90=as.data.frame(numoutsum.90)

numout.95=matrix(nrow=nrow(study2.95), ncol=ncol(study2.95))
for (j in 1:ncol(study2.95)) {
  for (i in 1:nrow(study2.95)) {
    if (study2.95[i,j]<.0125 | study2.95[i,j]>.0375) {
      numout.95[i,j]=1
    } else {
      numout.95[i,j]=0
    }
  }
}

numoutsum.95=t(colSums(numout.95))
colnames(numoutsum.95)=colnames(study2.95)
numoutsum.95=as.data.frame(numoutsum.95)



#Combine proportions above and below z CI into a single var
wide.80 = unite(numoutsum.80,  col = z, 1, 9, sep=" ", remove =TRUE)
wide.90 = unite(numoutsum.90,  col = z, 1, 9, sep=" ", remove =TRUE)
wide.95 = unite(numoutsum.95,  col = z, 1, 9, sep=" ", remove =TRUE)
#Combine proportions above and below M CI into a single var
wide.80 = unite(wide.80,  col = M, 2, 9, sep=" ", remove =TRUE) 
wide.90 = unite(wide.90,  col = M, 2, 9, sep=" ", remove =TRUE) 
wide.95 = unite(wide.95,  col = M, 2, 9, sep=" ", remove =TRUE) 
#Combine proportions above and below Jackknife CI into a single var
wide.80 = unite(wide.80,  col = jackknife, 3, 9, sep=" ", remove =TRUE) 
wide.90 = unite(wide.90,  col = jackknife, 3, 9, sep=" ", remove =TRUE) 
wide.95 = unite(wide.95,  col = jackknife, 3, 9, sep=" ", remove =TRUE) 
#Combine proportions above and below Percentile Bootstrap CI into a single var
wide.80 = unite(wide.80,  col = percentile, 4, 9, sep=" ", remove =TRUE)
wide.90 = unite(wide.90,  col = percentile, 4, 9, sep=" ", remove =TRUE)
wide.95 = unite(wide.95,  col = percentile, 4, 9, sep=" ", remove =TRUE)
#Combine proportions above and below Bias-Corrected Boostrap CI into a single var
wide.80 = unite(wide.80,  col = bias.corrected, 5, 9, sep=" ", remove =TRUE) 
wide.90 = unite(wide.90,  col = bias.corrected, 5, 9, sep=" ", remove =TRUE) 
wide.95 = unite(wide.95,  col = bias.corrected, 5, 9, sep=" ", remove =TRUE) 
#Combine proportions above and below Bootstrap T CI into a single var
wide.80 = unite(wide.80,  col = boot.t, 6, 9, sep=" ", remove =TRUE) 
wide.90 = unite(wide.90,  col = boot.t, 6, 9, sep=" ", remove =TRUE) 
wide.95 = unite(wide.95,  col = boot.t, 6, 9, sep=" ", remove =TRUE) 
#Combine proportions above and below Bootstrap Q CI into a single var
wide.80 = unite(wide.80,  col = boot.Q, 7, 9, sep=" ", remove =TRUE)
wide.90 = unite(wide.90,  col = boot.Q, 7, 9, sep=" ", remove =TRUE)
wide.95 = unite(wide.95,  col = boot.Q, 7, 9, sep=" ", remove =TRUE)
#Combine proportions above and below Monte Carlo CI into a single var
wide.80 = unite(wide.80,  col = Monte.Carlo, 8, 9, sep=" ", remove =TRUE) 
wide.90 = unite(wide.90,  col = Monte.Carlo, 8, 9, sep=" ", remove =TRUE) 
wide.95 = unite(wide.95,  col = Monte.Carlo, 8, 9, sep=" ", remove =TRUE) 


#Convert to long table format and split vars to get z, M, etc. indicator var
long.80 = pivot_longer(wide.80, cols = 1:8,
                    names_to = c("type"), 
                    values_to = "bounds")
long.90 = pivot_longer(wide.90, cols = 1:8,
                       names_to = c("type"), 
                       values_to = "bounds")
long.95 = pivot_longer(wide.95, cols = 1:8,
                       names_to = c("type"), 
                       values_to = "bounds")

#Split lower and upper proportions into separate columns
long2.80=separate(long.80,
               col='bounds',
               sep=' ',
               into=c("Left80", "Right80"))
long2.90=separate(long.90,
                  col='bounds',
                  sep=' ',
                  into=c("Left90", "Right90"))
long2.95=separate(long.95,
                  col='bounds',
                  sep=' ',
                  into=c("Left95", "Right95"))

#Convert Left and Right columns back to numeric type (instead of character) and then sum to get Total var for each CI type
long2.80=transform(long2.80, Left80=as.numeric(Left80), Right80=as.numeric(Right80))
long2.80$Total80=apply(long2.80[,2:3], 1, sum)
long2.90=transform(long2.90, Left90=as.numeric(Left90), Right90=as.numeric(Right90))
long2.90$Total90=apply(long2.90[,2:3], 1, sum)
long2.95=transform(long2.95, Left95=as.numeric(Left95), Right95=as.numeric(Right95))
long2.95$Total95=apply(long2.95[,2:3], 1, sum)

#attach a Total row at bottom of table
long3.80=rbind(long2.80, setNames(data.frame('Total',
                                             sum(long2.80[,2]),
                                             sum(long2.80[,3]),
                                             sum(long2.80[,4])),
                                  names(long2.80))
                            
)
long3.90=rbind(long2.90, setNames(data.frame('Total',
                                             sum(long2.90[,2]),
                                             sum(long2.90[,3]),
                                             sum(long2.90[,4])),
                                  names(long2.90))
               
)
long3.95=rbind(long2.95, setNames(data.frame('Total',
                                             sum(long2.95[,2]),
                                             sum(long2.95[,3]),
                                             sum(long2.95[,4])),
                                  names(long2.95))
               
)

#Combine the 80% CI, 90% CI, and 95% CI matrices to form final table
FinalTable4=cbind(long3.80, long3.90[,-1], long3.95[,-1])
#Add Overall Left, Right, and Total columns to final table
FinalTable4$LeftOverall=apply(FinalTable4[,c(2, 5, 8)], 1, sum)
FinalTable4$RightOverall=apply(FinalTable4[,c(3, 6, 9)], 1, sum)
FinalTable4$TotalOverall=apply(FinalTable4[,c(4, 7, 10)], 1, sum)

#Write table to .csv file
write.table(FinalTable4, file = "C:/Users/tdtibbet/Desktop/MacKinnon2004Table4.csv",
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)









##############################Table 5##############################
#Read in .csv file
study2reject=read.csv("C:/Users/tdtibbet/Desktop/MacKinnon2004RejectRateSummaryStudy2.csv")
#Create dataframe that contains indicator var (=0 if ab is 0 and =1 if ab is nonzero)
study2rejectout=study2reject
study2rejectout$abind=ifelse(study2rejectout$ab==0,0,1)
#Remove a, b, and ab columns
study2rejectout=study2rejectout[,-2:-4]

library(tidyverse)
#Average rejection rates across combinations of sample size and null/nonzero indirect effects
study2rejectout1=study2rejectout
study2rejectout1=unite(study2rejectout,  col = abind_n, 10, 1, sep=" ", remove =TRUE) 
study2rejectout1 =  study2rejectout1 %>% 
  group_by(abind_n) %>% 
  summarise(z = mean(z.reject.rate),
            M = mean(M.reject.rate),
            jackknife = mean(jackknife.reject.rate),
            percentile = mean(percentile.boot.reject.rate),
            bias.corrected = mean(bias.corrected.boot.reject.rate),
            boot.t = mean(boot.t.reject.rate),
            boot.Q = mean(boot.Q.reject.rate),
            Monte.Carlo = mean(Monte.Carlo.reject.rate)
  )

#Separate the ab indicator/sample size variable into two variables
study2rejectout1=separate(study2rejectout1,
              col='abind_n',
              sep=' ',
              into=c("abind", "n"))

#Convert to wide table format to get z, M, etc. vars for each sample size
wide = pivot_wider(study2rejectout1, names_from = n, values_from = c(z, M, jackknife, percentile, bias.corrected, boot.t, boot.Q, Monte.Carlo))
#Convert back to long table format and split vars to get z, M, etc. indicator var
long = pivot_longer(wide, cols = 2:33,
                    names_to = c("type","n"), 
                    names_sep = "_", 
                    values_to = "Reject Rate")
#Convert one last time to wide table format for final format
wide2=pivot_wider(long, names_from=n, values_from="Reject Rate")


#Reorder columns so smaller sample sizes come first
wide2=wide2 %>% 
  select(abind, type, `25`, `50`, `100`, `200`)

#Write table to .csv file
write.table(wide2, file = "C:/Users/tdtibbet/Desktop/MacKinnon2004Table5.csv",
            append=TRUE, sep = ",",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE)

##########ADDITIONAL FORMATTING FOR ALL TABLES WAS DONE IN EXCEL##########