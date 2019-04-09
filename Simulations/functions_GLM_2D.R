###################################################
##   Functions to compute Beta, SE(Betas) and T  ##
###################################################

tbeta <- function(volume,design){

 
  time <- length(design) 
  dimens <- dim(volume) 

  varact <- c(var(design))
  meanact <- c(mean(design))
  b1 <- apply(volume,c(1,2),cov,y=design)/varact 
}

tsebeta <- function(b1, volume, design){

  time <- length(design) 
  dimens <- dim(volume) 

  varact <- c(var(design)) 
  meanact <- c(mean(design))

  b0 <- apply(volume,c(1,2),mean)-b1*meanact 
  predact <- array(array(b1,dim=c(dimens[1],dimens[2],1))%*%array(design,dim=c(1,time)),dim=dimens) 
  pred <- array(rep(b0,time),dim=c(dimens[1],dimens[2],time)) + predact 
  help <- (pred-volume)^2 
  se2 <- apply(help,c(1,2),sum)/(time-2) 
  help2 <- sum((design-mean(design))^2) 
  sb1 <- sqrt(se2/help2)

}

tvolume <- function(b1, sb1) {  

  tmap <- b1/sb1

}

glmresid <- function(b1, volume, design){

  time <- length(design) 
  dimens <- dim(volume)

  varact <- c(var(design)) 
  meanact <- c(mean(design)) 

  b0 <- apply(volume,c(1,2,3),mean)-b1*meanact 
  predact <- array(array(b1,dim=c(dimens[1],dimens[2],dimens[3],1))%*%array(design,dim=c(1,time)),dim=dimens) 
  pred <- array(rep(b0,time),dim=c(dimens[1],dimens[2],dimens[3],time)) + predact 
  help <- (volume-pred) 

}

