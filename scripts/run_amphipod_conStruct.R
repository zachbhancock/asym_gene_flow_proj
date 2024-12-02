################################################################
################################################################
#
# running conStruct analyses on amphipod data
#
################################################################
################################################################

library(conStruct)
library(parallel)
library(foreach)

################################
# run conStruct
################################

set.seed(123)

load("amphipod_conStruct_infile.Robj")

cl <- parallel::makeCluster(12)
doParallel::registerDoParallel(cl)

outs <- foreach::foreach(i = 1:4,.packages="conStruct") %dopar% {
	conStruct(spatial=TRUE,
			  K=i,
			  freqs=genData_new,
			  geoDist=geoDist_neighbor_new,
			  coords=coords_new,
			  prefix=sprintf("shovelbugs_K%s_wypt",i),
			  n.chains=5,
			  n.iter=2e3,
			  make.figs=TRUE,
			  save.files=TRUE)
}

doParallel::stopImplicitCluster()

if(FALSE){
	
layer.contributions <- matrix(NA,nrow=4,ncol=4)
load("shovelbugs_K1_conStruct.results.Robj")
load("shovelbugs_K1_data.block.Robj")

# calculate layer contributions
bestChain <- which.max(unlist(lapply(conStruct.results,function(x){x$MAP$lpd})))
layer.contributions[,1] <- c(calculate.layer.contribution(conStruct.results[[bestChain]],data.block),rep(0,3))
tmp <- conStruct.results[[1]]$MAP$admix.proportions

for(i in 2:4){
    # load the conStruct.results.Robj and data.block.Robj
    #   files saved at the end of a conStruct run
    load(sprintf("shovelbugs_K%s_conStruct.results.Robj",i))
    load(sprintf("shovelbugs_K%s_data.block.Robj",i))
    
    # match layers up across runs to keep plotting colors consistent
    #   for the same layers in different runs
	bestChain <- which.max(unlist(lapply(conStruct.results,function(x){x$MAP$lpd})))
    tmp.order <- match.layers.x.runs(tmp,conStruct.results[[bestChain]]$MAP$admix.proportions)  

    # calculate layer contributions
    layer.contributions[,i] <- c(calculate.layer.contribution(conStruct.results=conStruct.results[[bestChain]],
                                                             data.block=data.block,
                                                             layer.order=tmp.order),
                                    rep(0,4-i))
    tmp <- conStruct.results[[bestChain]]$MAP$admix.proportions[,tmp.order]
}

pdf(file="amphipod_layer_contributions.pdf")
	barplot(layer.contributions,
	        col=c("blue", "red", "goldenrod1","forestgreen"),
	        xlab="",
	        ylab="layer contributions",
	        names.arg=paste0("K=",1:4))
dev.off()
}