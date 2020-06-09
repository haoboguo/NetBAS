rm(list=ls())
start = Sys.time()
debug = 0
#library(igraph)
library(foreach)
library(doMC)
registerDoMC(8) #too many cores may takes more RAM, 10K = 5 core X 11 min

# Yeast PIN 
pairs= read.csv("Data/yeast.pin.csv", stringsAsFactors = FALSE)
if (debug > 9) {
  pairs = pairs[1:1E3,]
}

#First, define a function to calculate the difference
 diff.value = function( inpairs, invalues ) {
   names(inpairs) = c("id1", "id2")
   names(invalues) = c("id", "value")
   inpairs$value1= invalues$value[match(inpairs$id1, invalues$id)];
   inpairs$value2= invalues$value[match(inpairs$id2, invalues$id)];
   ret = mean( abs( inpairs$value1 - inpairs$value2 ), na.rm=T );
 } 

# Janssens data

#list.files(path='Data/Dang_Existing_Data')
#dang.tb = read.csv("Data/Dang_Existing_Data/2012-12-17-UVW.tab", sep='\t', stringsAsFactors = FALSE)
#list.files(path='Data/')
jtb.ori = read.csv("Data/e15.prot.csv", stringsAsFactors=FALSE)
summary(jtb.ori)

columns_ = c("X7.8","X10.7" , "X14","X17.8","X22","X26.8","X32.2","X38.4","X45.4","X53.3","X62.2","X72.3");

jtb.ori[jtb.ori< 0] = 0.1 # two entries have negative values (noises)

jtb.ratios = jtb.ori[, columns_ ] / jtb.ori$X7.8

# 
jtb.ratios = log2(jtb.ratios)
jtb.ratios$ORF = jtb.ori$ORF
summary(jtb.ratios)

# permutation

#for( i in 2:2) {
for( i in 2:12) {
col_ = i
exp_ = names(jtb.ratios)[col_]
sub_ = jtb.ratios[, c(13,col_)]

hist.name <- paste("eLife15.all/", "eLife15.", "prot.hist.", i, ".pdf", sep="")
pdf(hist.name, width=5, height=4, paper='special')
hist(sub_[,2], main=paste('Jassens ratios ', exp_))
dev.off()

# Simiilar to Fraser02 Fig 3A. Use permutation to calculate p-value of difference

 # calculate the observed difference in measurement of interest
 diff.value.obs = diff.value ( pairs, sub_ );
 paste( "Observed delta value = ", diff.value.obs); 

####loop over MS02 null networks

ms02files = list.files(path='/scr/hguo/R/permutations/rls.ratio/bg163/ms02/', pattern = ".csv")
if (debug > 0 ) {ms02files = ms02files[1: 10000] }
permutated.diff.values = c();

# file = "ms02.1.csv" #debug
permutated.diff.values= foreach( fi =1:length(ms02files) ) %dopar% {  #this may takes large RAM
#for (fi in 1 : length(ms02files)){
  if ( debug > 0 ) { print(ms02files[fi]) }
  ms02_pairs= read.csv(paste("/scr/hguo/R/permutations/rls.ratio/bg163/ms02/", ms02files[fi], sep=''),
                       colClasses = c("character", "character"))
  ms02_pairs = ms02_pairs[,1:2]
  #permutated.diff.values = c( permutated.diff.values, diff.value ( ms02_pairs, dang.sub ))
  diff.value ( ms02_pairs, sub_ )
}

 permutated.diff.values = unlist(permutated.diff.values)

 # calulate p-value
 sub = permutated.diff.values[ permutated.diff.values <= diff.value.obs ]
 p.value = length(sub) / length(permutated.diff.values);

 #Calculate Z-score here 
 Z = ( diff.value.obs - mean(permutated.diff.values, na.rm = T) ) / sd(permutated.diff.values)
 
 print(paste("p-value of observed delta value =", p.value, " Z-score=", Z));

 xmin = min(c( permutated.diff.values,diff.value.obs) )*0.98
 xmax = max(c( permutated.diff.values,diff.value.obs) )*1.04
 
 # generate a figure
 file.name <- paste("eLife15.all/", "eLife15.", "prot.", i, ".png", sep="")
 png(file.name, width=5, height=4, res=1200, unit="in")
 hist( permutated.diff.values,  br=20, main = paste(exp_, ", p=", round(p.value,3), ", z=", round(Z, 3), sep=""), xlim = c(xmin, xmax), xlab=expression(paste(Delta, "Ratio",sep="")), ylab="Frequency, 10k Null Models" );
 arrows( diff.value.obs, 1000, diff.value.obs, 10, col="red", lwd=1 );
 text( diff.value.obs, 1200, paste("obs=", "\n", round(diff.value.obs,3)), cex=0.6);
 dev.off()
 
} #end of i-th column

system("convert -delay 80 eLife15.all/eLife15.prot.*.png eLife15.all.eLife15.pro.gif")

end = Sys.time()
print( end - start)

