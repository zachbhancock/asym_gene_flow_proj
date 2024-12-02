################################################################
################################################################
#
# prepping conStruct analyses on amphipod data
#
################################################################
################################################################

library(conStruct)

################################
# load in data
################################

# geographic metadata
metadata <- read.table("new_indv_coords.txt",header=TRUE)
coords <- cbind(metadata$long,metadata$lat)
row.names(coords) <- metadata[,1]
geoDist <- fields::rdist.earth(coords,miles=FALSE)

# genetic data
genData <- conStruct::structure2conStruct(infile="feems_filt_struct.txt",
										  onerowperind=FALSE,
										  start.loci=3,
										  start.samples=2,
										  missing.datum=-9,
										  outfile="amphipod_conStruct_data")

################
# deal with missing data 
################

# assess missing data
#	missing data by sample
md.s <- apply(genData,1,function(x){length(which(is.na(x)))})/ncol(genData)
#	missing data by locus (capped at 25%)
md.l <- apply(genData,2,function(x){length(which(is.na(x)))})/nrow(genData)

# drop all samples genotyped at <50% of loci
#	and regenerate data objects for running conStruct
toDrop <- which(md.s > 0.3)
coords_new <- coords[-toDrop,]
geoDist_new <- fields::rdist.earth(cbind(coords_new))
genData_new <- genData[-toDrop,]

# create a geographic distance matrix 
#	that's the sum of nearest neighbor distances
unique_coords <- cbind(unique(coords[,1]),unique(coords[,2]))
row.names(unique_coords) <- unique(unlist(lapply(strsplit(metadata$id,"_"),function(x){paste0(x[-1],collapse="_")})))
#plot(coords, type='n') ; text(coords,labels=row.names(coords))
#plot(unique_coords,type='n') ; text(unique_coords,labels=row.names(unique_coords))
#plot(unique_coords,type='n') ; text(unique_coords,labels=1:18)

geoDist_neighbor <- matrix(NA,nrow=nrow(coords),ncol=nrow(coords))
neighbor_dists <- unlist(
					lapply(1:17,
						function(i){
							fields::rdist.earth(unique_coords[i+1,,drop=FALSE],
												unique_coords[i,,drop=FALSE])}))
for(i in 1:nrow(coords)){
	for(j in 1:nrow(coords)){
		p1 <- metadata$popNumber[i]
		p2 <- metadata$popNumber[j]
		if(p1==p2){
			geoDist_neighbor[i,j] <- 0	
		} else {
			if(p2 > p1){
				geoDist_neighbor[i,j] <- sum(neighbor_dists[p1:(p2-1)])
			} else if(p1 > p2){
				geoDist_neighbor[i,j] <- sum(neighbor_dists[p2:(p1-1)])	
			}
		}
	}
}

geoDist_neighbor_new <- geoDist_neighbor[-toDrop,-toDrop]

################
# save data output objects
################

save(coords_new,geoDist_new,genData_new,geoDist_neighbor_new,file="amphipod_conStruct_infile.Robj")

