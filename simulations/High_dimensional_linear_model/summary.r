library(xtable)
library(mvtnorm)
source("statistics.R")

main <- function(number)
{
	file_names = list.files()
	tran = "lm"
	n_args = c(3,2)
	jobs = assig(n_args)
	
	result1 = vector("list", 6)
	result2 = vector("list", 6)
	result3 = vector("list", 6)
	result4 = vector("list", 6)
	
	for (nb in 1:6) {
	    id = jobs[, nb]
	    n=200
	    p=500
	    tran = "lm"
	    error = c("norm", "contaminate", "t2")[id[1]]
	    rho = c(0.3, 0.7)[id[2]]
	    case_name = paste(tran, n, p, rho, error, sep="-")
	    case = file_names[grep(case_name, file_names)]
	    result = fits(case, 3, 2, n,  p,  rho, error, tran)
	    result1[[nb]] = result[[1]]
	    result2[[nb]] = result[[2]]
	    result3[[nb]] = result[[3]]
	    result4[[nb]] = result[[4]]
	}

  result1[[1]] = rbind(result1[[1]], result1[[2]], result1[[3]])
  result1[[4]] = rbind(result1[[4]], result1[[5]], result1[[6]])
  result2[[1]] = rbind(result2[[1]], result2[[2]], result2[[3]])
  result2[[4]] = rbind(result2[[4]], result2[[5]], result2[[6]])
  result3[[1]] = rbind(result3[[1]], result3[[2]], result3[[3]])
  result3[[4]] = rbind(result3[[4]], result3[[5]], result3[[6]])
    
  result1 = cbind(result1[[1]], result1[[4]])
  result2 = cbind(result2[[1]], result2[[4]])
  result3 = cbind(result3[[1]], result3[[4]])
    
  for (j in 2:6) result4[[1]] = rbind(result4[[1]], result4[[j]])
    
  FITS  = result1
  NNAVC = result2
  NNANZ = result3
  IDTFALL = result4[[1]]
  
  tex = xtable(NNAVC,digits=c(0))
  print(tex, file = "NNAVC.tex")
  tex = xtable(NNANZ,digits=c(0))
  print(tex, file = "NNANZ.tex")
  tex = xtable(FITS)
  print(tex, file = "FITS.tex")
  tex = xtable(IDTFALL)
  print(tex, file = "IDTF.tex")
}
  
args <- commandArgs(TRUE)
main(args[1])
