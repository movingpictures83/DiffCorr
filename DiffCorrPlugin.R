# Title     : TODO
# Objective : TODO
# Created by: stebliankin
# Created on: 10/14/19

##########################################################
# Set up remote session
#library("ssh")
#session <- ssh_connect("vsteb002@castalia.cis.fiu.edu")
##########################################################

library("DiffCorr")

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
        pfix = prefix()
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
##########################################################
# Constants setup
##########################################################
CONTROL_CORR_FILE <<- paste(pfix, parameters["control", 2], sep="/")
DISEASE_CORR_FILE <<- paste(pfix, parameters["disease", 2], sep="/")
N_SAMPLES_CONTROL <<- as.integer(parameters["numcontrol", 2])
N_SAMPLES_DISEASE <<- as.integer(parameters["numdisease", 2])
}

run <- function() {}

output <- function(outputfile) {

##########################################################
# Load Matrixes
##########################################################

#g1 - Non Users (CONTROL)
g1_corr <- read.table(CONTROL_CORR_FILE, sep=",", header = TRUE, row.names = 1)
#g1_corr <- matrix(g1_corr)

#g2 - users (DISEASE)
g2_corr <- read.table(DISEASE_CORR_FILE, sep=",", header = TRUE, row.names = 1)


#g2_corr <- matrix(g2_corr)
#comp.2.cc.fdr(output.file = "/Users/stebliankin/Desktop/SabrinaProject/scripts/DifferentialNetwork/test_out.txt", g1_corr, g2_corr, method="spearman", threshold = 0.05, p.adjust.methods = "none")

#g1 <- read.table("/Users/stebliankin/Desktop/SabrinaProject/PlumaPipline/MultiOmics/multiomics.users.csv", sep=",", header=TRUE, row.names = 1)
#g1 <- t(g1)

#g2 <- read.table("/Users/stebliankin/Desktop/SabrinaProject/PlumaPipline/MultiOmics/multiomics.nonUsers.csv", sep=",", header=TRUE, row.names = 1)
#g2 <- t(g2)

#comp.2.cc.fdr(output.file = "/Users/stebliankin/Desktop/SabrinaProject/scripts/DifferentialNetwork/test_out.txt", g1, g2, method="spearman", threshold = 0.05)


# Initialize output differential correlation
diff_corr_df <- data.frame(matrix(ncol=ncol(g1_corr), nrow=nrow(g1_corr)))
colnames(diff_corr_df) <- colnames(g1_corr)
rownames(diff_corr_df) <- rownames(g1_corr)

# Initialize a matrix with sign from wich the correlation is changing 
sign_corr_df <- data.frame(matrix(ncol=ncol(g1_corr), nrow=nrow(g1_corr)))
colnames(sign_corr_df) <- colnames(g1_corr)
rownames(sign_corr_df) <- rownames(g1_corr)



for(row in 1:nrow(g1_corr)) {
  for(col in 1:ncol(g1_corr)) {
    p_val <- compcorr(N_SAMPLES_DISEASE,g2_corr[row, col], N_SAMPLES_CONTROL, g1_corr[row, col])$pval
    if (p_val < 0.05) {
      diff_corr_df[row, col] <- g2_corr[row, col] -  g1_corr[row, col]
      
      if (g1_corr[row, col]>0){# change from +
        from_sign <- "+"
      }
      else if (g1_corr[row, col]<0){# change from -
        from_sign <- "-"
      }
      else if (g1_corr[row, col]==0){# change from 0
        from_sign <- "0"
      }
      
      if (g2_corr[row, col]>0){# change from +
        to_sign <- "+"
      }
      else if (g2_corr[row, col]<0){# change from -
        to_sign <- "-"
      }
      else if (g2_corr[row, col]==0){# change from 0
        to_sign <- "0"
      }
      sign_corr_df[row, col] <- paste(from_sign, to_sign, sep="/")
      
    }
    else{
      diff_corr_df[row, col]<-0
      sign_corr_df[row, col] <- paste("0", "0", sep="/")
    }
  }
}

row <- 30
col <- 4

p_val <- compcorr(25,g2_corr[row, col], 25, g1_corr[row, col])$pval
if (p_val < 0.05) {
  diff <- g2_corr[row, col] -  g1_corr[row, col]
}


row <- 30
col <- 4

p_val <- compcorr(25,g2_corr[row, col], 25, g1_corr[row, col])$pval
if (p_val < 0.05) {
  diff_corr_df[row, col] <- g2_corr[row, col] -  g1_corr[row, col]
  
  if (g1_corr[row, col]>0){# change from +
    from_sign <- "+"
  }
  else if (g1_corr[row, col]<0){# change from -
    from_sign <- "-"
  }
  else if (g1_corr[row, col]==0){# change from 0
    from_sign <- "0"
  }
  
  if (g2_corr[row, col]>0){# change from +
    to_sign <- "+"
  }
  else if (g2_corr[row, col]<0){# change from -
    to_sign <- "-"
  }
  else if (g2_corr[row, col]==0){# change from 0
    to_sign <- "0"
  }
  sign_corr_df[row, col] <- paste(from_sign, to_sign, sep="/")
  
}


write.csv(diff_corr_df, paste(outputfile, "diff.csv", sep="."))
write.csv(sign_corr_df, paste(outputfile, "sign.csv", sep="."))
}
