#Rscript to convert .ibd file to .sims in preparation for running MAPS
#Edit lines 37-38 according to the length IBD segment you want to filter for
#This script was written by Hussein Al-Asadi and is available in its original format here: https://github.com/halasadi/MAPS/tree/master/convert_beagle

meta_info <- data.frame(fid = meta_info_Ethiopia$ID, iid = meta_info_Ethiopia$ID, long = meta_info_Ethiopia$longitude,
			lat = meta_info_Ethiopia$latitude, stringsAsFactors = FALSE)

hapibd <- read.table("/share/hennlab/users/espless/PhasingAndIBD/hapibd_Pagani_Scheinfeldt/IBDOut/Pagani_Scheinfeldt_ibd_all_merged_cM.ibd", header = TRUE, stringsAsFactors = FALSE)

id1_relabel <- paste0(hapibd$id1, "_", hapibd$id1, "_", "1")

id2_relabel <- paste0(hapibd$id2, "_", hapibd$id2, "_", "2")

ibd_data <- data.frame(id1_relabel, id2_relabel, hapibd$start, hapibd$end, stringsAsFactors = FALSE)
#added stringsAsFactors

names(ibd_data) <- c('id1', 'id2', 'start', 'end')

ids <- paste0(meta_info$iid, "_", meta_info$fid)

locs <- paste0(meta_info$lat, " ", meta_info$long)

ids <- c(paste0(ids, "_1"), paste0(ids, "_2"))
locs <- c(locs, locs)

n <- length(ids)

# set up the matrix
ibd_summary <- matrix(nrow = n, ncol = n, 0)
rownames(ibd_summary) <- ids
colnames(ibd_summary) <- ids

# compute the lengths of the PSC segments
lengths <- as.numeric(ibd_data$end)-as.numeric(ibd_data$start)

# highlight PSC segments greater than 2cM
lowerBnd <- 6
upperBnd <- Inf
selected_inds <- which(lengths > lowerBnd & lengths < upperBnd)

for (i in 1:length(selected_inds)){
	id1 <- ibd_data$id1[selected_inds[i]]
	id2 <- ibd_data$id2[selected_inds[i]]
	if (id1 %in% ids & id2 %in% ids){
	       ibd_summary[id1, id2] = ibd_summary[id1, id2] + 1
	       ibd_summary[id2, id1] = ibd_summary[id1, id2]
	   }
    }

# write similarity matrix and coordinates to a .sims and .coord file respectively
write.table(ibd_summary, file = paste0("/share/hennlab/users/espless/MAPS/Pagani_Scheinfeldt/6_Inf/maps_", lowerBnd, "_", upperBnd, ".sims"), quote=FALSE, sep = " ", row.names = FALSE, col.names=FALSE)
