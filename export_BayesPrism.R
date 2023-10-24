#export visium

#work under /data/peer/tinyi/sawyers/visium
#conda activate r_env

library(Seurat)
library(ggplot2)
library(BayesPrism)



load("/data/peer/tinyi/sawyers/visium/bp.res.week10.singulatorRef_cellBender_4CTMarker.rdata")


ct_obj <- CreateAssayObject(counts = t(get.fraction(bp.res, which.theta="final", state.or.type="type")))
dat.merged[['BP.frac']] <- ct_obj


Znk_obj <- CreateAssayObject(counts = t(apply(bp.res@posterior.initial.cellType@Z,c(1,3),sum)))
dat.merged[['BP.Znk']] <- Znk_obj


cord <- do.call(rbind.data.frame,lapply(dat.merged@images,function(i) i@coordinates))
dat.merged@meta.data <- cbind.data.frame(dat.merged@meta.data, cord)


dat.merged@meta.data$Barcode <- substr(rownames(dat.merged@meta.data), 27, 44)


#load expression


#exclued #17, a normal epithelial cell

ct.names <- colnames(get.fraction(bp.res, which.theta="final", state.or.type="type"))

tumor_cells <- ct.names[!ct.names %in% c("stromal","myeloid","normal","endo")]


dat.merged[['tum.exp']] <- CreateAssayObject(counts = t(rowSums(bp.res@posterior.initial.cellType@Z[,, tumor_cells], dims=2)))


library(readr)


df <- as.data.frame(dat.merged[["tum.exp"]]@counts)
df <- cbind.data.frame(gene_ids =rownames(df), df)
write_csv( df, file="week10Tissue.bp.tum_snrnaSigulatorCellbender_4CTMarker.csv")




