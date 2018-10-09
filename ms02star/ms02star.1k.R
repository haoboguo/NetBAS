# 20180614 error corrected
# 20180522 revise to deal with highly connected networks

rm(list=ls())

#R -f file --args input_network_csv_filename output_csv_filename in_ncycles debug

options(echo=TRUE) # if you want see commands in output file 
args <- commandArgs(trailingOnly = TRUE)
#print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandsArgs(trailingOnly=FALSE))

for (p in 1:1000) {

infile = args[1]
outdir = args[2]
outfile = paste(outdir, "/", "ms02.", p, ".csv", sep="")
in_ncycles = 100
debug = 9

#set.seed(100)
print(paste("infile=[",infile, "]\n"));

f_match_degree = function (inTb) {
  inTb$DegreePercentile1 = degree$DegreePercentile[match(inTb$id1, degree$ORF)];
  inTb$DegreePercentile2 = degree$DegreePercentile[match(inTb$id2, degree$ORF)];
  inTb$DegreePercentileProduct = inTb$DegreePercentile1 * inTb$DegreePercentile2; 
  inTb = inTb[order(inTb$DegreePercentileProduct, decreasing = FALSE), ];
  return (inTb); 
}

f_redo_rest.tb = function ( inTb ) {
  redo.tb = inTb[ inTb$selfpairs==1 | inTb$tag_counts>1, ]
  rest.tb = inTb[ inTb$selfpairs==0 & inTb$tag_counts==1, ]
  return(list("redo.tb"=redo.tb, "rest.tb" = rest.tb))
}

# recursive_permutation_no_selfpairing_v0.01 rpns
#This function only perform permutation, but does NOT check for correctness of permutation. 
#Correctness of permutation should be check outside of this function!
recursive_permutation_no_selfpairing_v0.01 = function( inpairs,  ncycles=5, indebug=0, preserve_rate = 0.1 ) { 
    if(indebug>0) {
      cat(paste( "\n**", '(rpns) Start ncycles=', ncycles, ", preserve_rate=", preserve_rate, 
                 ", length(inpairs[,1])=",length(inpairs[,1]), "**\n\n" )); 
    }
    longids = c(as.character(inpairs[,1]), as.character(inpairs[,2]) )
    longids = sample(longids)
    len = length(inpairs[,1])
    
    newpairs0 = data.frame( cbind( longids[1:len], longids[(len+1): (2*len)]) )
    newpairs = t(apply(newpairs0, 1, sort)); #oder id1 and id2
    newpairs = data.frame(newpairs); 
    if(indebug > 2) {#check ids ordering results
      print(paste("===Before:Aftersort  ===="),NULL); 
      print(cbind( newpairs0, newpairs));
    }   
    names(newpairs) = c('id1', 'id2')
    newpairs$id1 = as.character( newpairs$id1)
    newpairs$id2 = as.character( newpairs$id2)    
    
    newpairs$tag =  paste(newpairs[,1], newpairs[,2], sep="_")
    counts = table( newpairs$tag )
    newpairs$tag_counts = counts[newpairs$tag]
    # if(indebug>8) {    counts;    }
    newpairs$selfpairs = ifelse( newpairs$id1 == newpairs$id2, 1, 0 )
    
    newpairs = f_match_degree(newpairs); 
    
    redo.tb = newpairs[ newpairs$selfpairs==1 | newpairs$tag_counts>1, ]
    rest.tb = newpairs[ newpairs$selfpairs==0 & newpairs$tag_counts==1, ]
    
    if(indebug>8) {
      print(paste('(rpns)ncycles=', ncycles, "===redopairs===="),NULL); 
      print (redo.tb);
      print(paste('(rpns)ncycles=', ncycles, "===restpairs===="),NULL); 
      print(rest.tb);
      print(paste("================="),NULL)
    }
    
    if( length(redo.tb[,1])>=1 ) {
       if ( ( ncycles == 0) | (length(rest.tb[,1]) < 1) ) { 
        print(paste("ncycles reached zero, ncycles, OR, not enough data in rest.tb for randomization", ncycles) );
        if (indebug > 8 ) {
          write.csv(redo.tb, "tmp/redo_tb.csv");
          write.csv(rest.tb, "tmp/rest_tb.csv");
        }
        return( newpairs ); # no more randomization inside of this function
      } else  {
        degreeProduct.cutoff = quantile(rest.tb$DegreePercentileProduct, prob = (1 - preserve_rate)); 
        unchangedpairs = rest.tb[ rest.tb$DegreePercentileProduct >degreeProduct.cutoff, ]
        selectedpairs = rbind(redo.tb,rest.tb[ rest.tb$DegreePercentileProduct <=degreeProduct.cutoff , ] )   
        
        if (indebug > 0) {
          print(paste('(rpns) ncycles=', ncycles, ", degreeProduct.cutoff=",degreeProduct.cutoff, 
                      ", length(redo.tb[,1])=",length(redo.tb[,1]),
                      ", length(rest.tb[,1])=",length(rest.tb[,1]) ),NULL); 
        }
        ncycles = ncycles - 1; 
        return( rbind(unchangedpairs, 
                      recursive_permutation_no_selfpairing_v0.01(selectedpairs, ncycles,indebug = indebug, preserve_rate = preserve_rate )))#20180522, recursive trap bug
      }
    } else if (length(redo.tb[,1])==0) { #20180522
      if(indebug>0) {
        print(paste('(rpns) SUCCESS, ncycles=', ncycles, 
                    ", length(redo.tb[,1])=",length(redo.tb[,1]),
                    ", length(rest.tb[,1])=",length(rest.tb[,1]) ),NULL); 
      }      
      return (newpairs )
    } 

}#end of function


success_flag = 0; 
global_cycles = 5;

net = read.csv( infile, header = T, colClasses = c("character", "character") )
head(net)

longids = c( net[,1], net[,2] );
degree = sort( table( longids ), decreasing = TRUE ); 
degree = data.frame(degree)
names(degree) = c("ORF", "degree"); 
degree$ORF = as.character( degree$ORF); 
# Pick highly connected hub nodes, and prioritize their permutations
# Notes of caution: This means degree-degree profile are not random in permuted networks
degree$DegreePercentile = 0;
for ( i in 1:length(degree[,1])) {
  degree$DegreePercentile[i] = 1 - length(which(degree$degree> degree$degree[i]))/ length(degree[,1]);
}

  perNet0 = recursive_permutation_no_selfpairing_v0.01(net, ncycles = in_ncycles, indebug = debug, preserve_rate = runif(1,min=0.5,max=0.8)); 
  perNet0 = f_match_degree( perNet0 );
  
  x = f_redo_rest.tb( perNet0 );
  redo.PerNet0 = x$redo.tb; 
  rest.PerNet0 = x$rest.tb;
  
  if ( length(redo.PerNet0[,1])==0 ) {
    write.csv(perNet0, outfile, quote=F, row.names=F)
  } else {
    system(paste("touch ", outfile, "_err")); 
  }
 } 
