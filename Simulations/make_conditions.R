##########################################
##      SIMULATION CONDITIONS           ##
##########################################


deltaval <- c(0.01,0.02,0.03)
scanval <- c(64,128,320,680,900)
percval <- c(0.75,0.80,0.85,0.90,0.95,0.99)


delta <- rep(deltaval, each=length(scanval))
scans <- rep(scanval, times=length(deltaval))
perc <- rep(percval, times=length(deltaval)*length(scanval))


conditions <- as.data.frame(cbind(delta,scans))
colnames(conditions) <- c("delta","scans")
head(conditions)
conditions

write.csv(conditions,file="pathwhereyourscriptsarestored/conditions_hpc.txt",row.names=FALSE)

