#compare Hotspot analysis (deconvolved vs undeconvolved)

#work under /data/peer/tinyi/sawyers/scanpy

library(Seurat)
library(BayesPrism)

#load joe's data
load("/data/peer/tinyi/sawyers/visium/ref.dat.sig.rdata")

cell.type.labels.new <- cell.type.labels

cell.type.labels.new[cell.state.labels=="endothelial"] <- "endothelial"
cell.type.labels.new[cell.state.labels=="lymphatic"] <- "lymphatic"
cell.type.labels.new[cell.state.labels=="pericyte"] <- "pericyte"
cell.type.labels.new[cell.state.labels=="glial"] <- "glial"
cell.type.labels.new[cell.state.labels %in% c("mse-1", "mse-2")] <- "stromal"
cell.type.labels.new[cell.state.labels %in% c("TFF3", "Mutant_L1","Mutant_L2","Mutant_B1","NEPC-P","NEPC") ] <- "tumor"
cell.type.labels.new[cell.type.labels %in% c("Macrophages", "Neutrophils","DC") ] <- "myeloid"


selected.ct <- c("endothelial", "stromal", "tumor", "myeloid")
selected.idx <- cell.type.labels.new %in% selected.ct

diff.exp.stat <- get.exp.stat(sc.dat= ref.dat[selected.idx,],# filter genes to reduce memory use
                              cell.type.labels= cell.type.labels.new[selected.idx],
                              cell.state.labels= cell.type.labels.new[selected.idx],
                              psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=20 #number of threads
                                          )

save(diff.exp.stat, file="diff.exp.stat_JoeData_4CellTypes_forBenchmark.rdata")


#map back to gene symbol 
#also subset on genes selected by hotspot analysis
autoCor_all <- read.csv("all_hotspot_autoCorr_snrnaCB_4CTMarker.csv")
gene_annot <- read.csv("/lila/data/peer/tinyi/sawyers/scanpy/snrna_singulator_adata_var_cellBender.csv")

gene_annot <- gene_annot[gene_annot$X %in% autoCor_all$Gene, ]

diff.exp.stat.symbol <- lapply(diff.exp.stat, function(i) {
	i.mappable <- i[rownames(i) %in% gene_annot$gene_ids,]
	rownames(i.mappable) <- gene_annot[match(rownames(i.mappable), gene_annot$gene_ids), "X"]
	i.mappable
	})



#get markers for t.test 

stat <- diff.exp.stat.symbol
pval.max=0.01
lfc.min=0.1
topN <-100
stat.filtered <- lapply(1:length(stat), function(i) {
        stat.i <- stat[[i]]
        stat.i.filtered <- stat.i[stat.i$pval.up.min < pval.max &
            stat.i$min.lfc > lfc.min, , drop = F]
        
        stat.i.filtered <- stat.i.filtered[order(stat.i.filtered$min.lfc,decreasing=TRUE)[1: min(topN, nrow(stat.i.filtered))],]   
        cat(names(stat)[i], ": ", nrow(stat.i.filtered), "\n")
        return(stat.i.filtered)
    })
    

markers.all.top100 <- lapply(stat.filtered, rownames)
names(markers.all.top100) <- names(stat)




#now lets check the correlation between hotspot runs

autoCor_all <- read.csv("all_hotspot_autoCorr_snrnaCB_4CTMarker.csv")

autoCor_tum <- read.csv("tum_hotspot_autoCorr_snrnaCB_4CTMarker.csv")

autoCor_tum <- autoCor_tum[match(autoCor_all$Gene, autoCor_tum$Gene),]



check.autoCor.over.genes <- function(autoCor_all, autoCor_tum, genes){
	
	stopifnot(all(autoCor_all$Gene== autoCor_tum $Gene))
	
	gene.idx <- autoCor_all$Gene %in% genes
	cat("correlation=", cor(autoCor_all[gene.idx, "Z"], autoCor_tum[gene.idx, "Z"], use="complete"),"\n")
	print(t.test (autoCor_all[gene.idx, "Z"]-autoCor_tum[gene.idx, "Z"],alternative="greater"))
	print(t.test (autoCor_all[gene.idx, "Z"]-autoCor_tum[gene.idx, "Z"]))
}


check.autoCor.over.genes(autoCor_all, autoCor_tum, markers.all.top100 $tumor)
# correlation= 0.9915758

        # One Sample t-test

# data:  autoCor_all[gene.idx, "Z"] - autoCor_tum[gene.idx, "Z"]
# t = -2.963, df = 99, p-value = 0.9981
# alternative hypothesis: true mean is greater than 0
# 95 percent confidence interval:
 # -2.615014       Inf
# sample estimates:
# mean of x
# -1.675889


        # One Sample t-test

# data:  autoCor_all[gene.idx, "Z"] - autoCor_tum[gene.idx, "Z"]
# t = -2.963, df = 99, p-value = 0.003815
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
 # -2.7981710 -0.5536069
# sample estimates:
# mean of x
# -1.675889

  
check.autoCor.over.genes(autoCor_all, autoCor_tum, markers.all.top100 $myeloid)
# correlation= 0.9792372

        # One Sample t-test

# data:  autoCor_all[gene.idx, "Z"] - autoCor_tum[gene.idx, "Z"]
# t = 5.0441, df = 99, p-value = 1.033e-06
# alternative hypothesis: true mean is greater than 0
# 95 percent confidence interval:
 # 1.824774      Inf
# sample estimates:
# mean of x
 # 2.720197


        # One Sample t-test

# data:  autoCor_all[gene.idx, "Z"] - autoCor_tum[gene.idx, "Z"]
# t = 5.0441, df = 99, p-value = 2.067e-06
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
 # 1.650140 3.790254
# sample estimates:
# mean of x
 # 2.720197 

check.autoCor.over.genes(autoCor_all, autoCor_tum, markers.all.top100 $stromal)
# correlation= 0.8959966

        # One Sample t-test

# data:  autoCor_all[gene.idx, "Z"] - autoCor_tum[gene.idx, "Z"]
# t = 6.3686, df = 99, p-value = 3.011e-09
# alternative hypothesis: true mean is greater than 0
# 95 percent confidence interval:
 # 7.621387      Inf
# sample estimates:
# mean of x
 # 10.30911


        # One Sample t-test

# data:  autoCor_all[gene.idx, "Z"] - autoCor_tum[gene.idx, "Z"]
# t = 6.3686, df = 99, p-value = 6.022e-09
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
  # 7.097201 13.521017
# sample estimates:
# mean of x
 # 10.30911

check.autoCor.over.genes(autoCor_all, autoCor_tum, markers.all.top100 $endothelial)
# correlation= 0.9760018

        # One Sample t-test

# data:  autoCor_all[gene.idx, "Z"] - autoCor_tum[gene.idx, "Z"]
# t = 2.9488, df = 99, p-value = 0.00199
# alternative hypothesis: true mean is greater than 0
# 95 percent confidence interval:
 # 0.5037446       Inf
# sample estimates:
# mean of x
 # 1.152914


        # One Sample t-test

# data:  autoCor_all[gene.idx, "Z"] - autoCor_tum[gene.idx, "Z"]
# t = 2.9488, df = 99, p-value = 0.00398
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
 # 0.3771372 1.9286905
# sample estimates:
# mean of x
 # 1.152914
 

#make violin plot

library(ggplot2)
library(ggallin)


 data_summary <- function(x) {
    m <- median(x)
    ymin <- as.numeric(quantile(x)[2])
    ymax <-as.numeric(quantile(x)[4])
    return(c(y=m,ymin=ymin,ymax=ymax))
 }


#make violin plot (two groups)
autoCor.df <- do.call(rbind.data.frame, lapply(1:length(markers.all.top100),function(i){
	marker.gene.idx <- autoCor_all$Gene %in% markers.all.top100 [[i]]
	all.df <- data.frame(marker_cell_type= names(markers.all.top100)[i], Z= autoCor_all[marker.gene.idx, "Z"], exp.type="all")
	tum.df <- data.frame(marker_cell_type= names(markers.all.top100)[i], Z= autoCor_tum[marker.gene.idx, "Z"], exp.type="tum")
	rbind.data.frame(all.df, tum.df )
}))


p <- ggplot(autoCor.df, aes(x= marker_cell_type, y= Z, fill= exp.type)) +
			geom_violin(trim=FALSE, scale="width")+ 
			stat_summary(fun.data= data_summary,geom="pointrange", color="white", position=position_dodge(0.9)) +
			theme_classic() + ggtitle("Hotspot auto-correlation") 

pdf("autoCor.violin.scrnaCB_4CTMarker.pdf",pointsize=8,useDingbats=FALSE)
print(p)
dev.off()





#make violin plot (diff between all and tumor)
autoCor.df <- do.call(rbind.data.frame, lapply(1:length(markers.all.top100),function(i){
	marker.gene.idx <- autoCor_all$Gene %in% markers.all.top100 [[i]]
	delta_Z <- autoCor_tum[marker.gene.idx, "Z"] - autoCor_all[marker.gene.idx, "Z"]
	dat.df <- data.frame(marker_cell_type= names(markers.all.top100)[i], delta_Z = delta_Z)
	dat.df
}))


# pdf("autoCor.delta_violin.scrnaCB_4CTMarker.pdf",pointsize=8,useDingbats=FALSE)

# p <- ggplot(autoCor.df, aes(x= marker_cell_type, y= delta_Z, fill= marker_cell_type)) +
			# geom_violin(trim=FALSE, scale="width")+ 
			# stat_summary(fun.data= data_summary,geom="pointrange", color="white", position=position_dodge(0.9)) +
			# theme_classic() + ggtitle("Hotspot auto-correlation deconovolved-undeconvolved") + scale_y_continuous(trans = 'pseudo_log',breaks = seq(-100, 50, by = 25))
# print(p)

# p <- ggplot(autoCor.df, aes(x= marker_cell_type, y= delta_Z, fill= marker_cell_type)) +
			# geom_violin(trim=FALSE, scale="width")+ 
			# stat_summary(fun.data= data_summary,geom="pointrange", color="white", position=position_dodge(0.9)) +
			# theme_classic() + ggtitle("Hotspot auto-correlation deconovolved-undeconvolved") 
# print(p)

# dev.off()

p <- ggplot(autoCor.df, aes(x= marker_cell_type, y= delta_Z, fill= marker_cell_type)) +
			geom_violin(trim=FALSE, scale="width")+ 
			stat_summary(fun.data= data_summary,geom="pointrange", color="white", position=position_dodge(0.9)) +
			theme_classic() + ggtitle("Hotspot auto-correlation deconovolved-undeconvolved") + scale_y_continuous(trans = 'pseudo_log',breaks = c(-100,-50, -25, -10, 0,10,25,50) ) +
  coord_flip() + geom_hline(yintercept = 0, linetype="dashed", color = "red")
  
  
ggsave(filename = "autoCor.delta_violin.scrnaCB_4CTMarker.pdf", plot = p, width = 10, height = 6)









#how about genes not specific to tumor?

stat <- diff.exp.stat.symbol
pval.max=0.05

nonmarkers <- rownames(stat[[1]])[apply(do.call(cbind,lapply(stat, "[[", "pval.up.min")),1,min)  > pval.max]

 
check.autoCor.over.genes(autoCor_all, autoCor_tum, nonmarkers)
# correlation= 0.990215

        # One Sample t-test

# data:  autoCor_all[gene.idx, "Z"] - autoCor_tum[gene.idx, "Z"]
# t = -0.43095, df = 403, p-value = 0.6666
# alternative hypothesis: true mean is greater than 0
# 95 percent confidence interval:
 # -0.7970751        Inf
# sample estimates:
 # mean of x
# -0.1651763


        # One Sample t-test

# data:  autoCor_all[gene.idx, "Z"] - autoCor_tum[gene.idx, "Z"]
# t = -0.43095, df = 403, p-value = 0.6667
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
 # -0.9186622  0.5883095
# sample estimates:
 # mean of x
# -0.1651763



#draw density plot
densScatterplot <- function(x1, x2, uselog=FALSE, n=256, range.values=NULL, ...) {
  df <- data.frame(x1, x2)

  if(uselog) {
    x <- densCols(log(x1,10),log(x2,10), colramp=colorRampPalette(c("black", "white")))
  } else {
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  }
  df$dens <- col2rgb(x)[1,] + 1L

  ## Map densities to colors
#  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  cols <- colorRampPalette(c("light gray", "#000099", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(n)
  df$col <- cols[df$dens]
	
	if(is.null(range.values)) range.values <- c(x1,x2)
	range.values <-range.values[is.finite(range.values)]
	.limits<-range(range.values)
	
  ## Plot it, reordering rows so that densest points are plotted on top
  if(uselog) {
    plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, log="xy", xlim=.limits, ylim=.limits,asp=1,...)
  } else {
    plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2,asp=1, xlim=.limits, ylim=.limits, ...)
  }

  abline(a=0,b=1,col="red",lty=2)
  
  text(x=quantile(range.values,prob=0.99) , 
  	   y=quantile(range.values,prob=0.75),
  	   labels=paste("cor=",round(cor(x1,x2,use="complete"),3),sep=""))
  
}


#get markers for density plot
stat <- diff.exp.stat.symbol
pval.max=0.01
lfc.min=0.1
stat.filtered <- lapply(1:length(stat), function(i) {
        stat.i <- stat[[i]]
        stat.i.filtered <- stat.i[stat.i$pval.up.min < pval.max &
            stat.i$min.lfc > lfc.min, , drop = F]
        
        cat(names(stat)[i], ": ", nrow(stat.i.filtered), "\n")
        return(stat.i.filtered)
    })
    

markers.all <- lapply(stat.filtered, rownames)
names(markers.all) <- names(stat)

# endothelial :  158
# stromal :  310
# tumor :  2583
# myeloid :  199


pdf("autoCor_dens.scrnaCB_4CTMarker.pdf", useDingbats=FALSE)

#range.values<-c(autoCor_all[,"Z"],  autoCor_tum[,"Z"])
range.values <- c(0,200)

selected.idx <- autoCor_all$Gene %in% markers.all $tumor
densScatterplot(x1= autoCor_all[selected.idx,"Z"], 
				x2= autoCor_tum[selected.idx,"Z"], 
				main="tumor markers",
				xlab="autoCor, all", ylab="autoCor, deconvolved",
				range.values= range.values)

selected.idx <- autoCor_all$Gene %in% unlist(markers.all[c('myeloid','stromal','endo')])
densScatterplot(x1= autoCor_all[selected.idx,"Z"], 
				x2= autoCor_tum[selected.idx,"Z"], 
				main="non-tumor markers",
				xlab="autoCor, all", ylab="autoCor, deconvolved",
				range.values= range.values)

selected.idx <- autoCor_all$Gene %in% nonmarkers
densScatterplot(x1= autoCor_all[selected.idx,"Z"], 
				x2= autoCor_tum[selected.idx,"Z"], 
				main="non-markers",
				xlab="autoCor, all", ylab="autoCor, deconvolved",
				range.values= range.values)


dev.off()





#test cross-correlation


pwCor_all <- read.csv("all_hotspot_pairWiseCorr_snrnaCB_4CTMarker.csv",row.names=1,check.names=F)

pwCor_tum <- read.csv("tum_hotspot_pairWiseCorr_snrnaCB_4CTMarker.csv",row.names=1,check.names=F)

pwCor_tum <- pwCor_tum[match(rownames(pwCor_all), rownames(pwCor_tum)), match(colnames(pwCor_all), colnames(pwCor_tum))]

pwCor_all <- as.matrix(pwCor_all)
pwCor_tum <- as.matrix(pwCor_tum)

#mask out diagonal
diag(pwCor_all) <- NA
diag(pwCor_tum) <- NA



#check pairwise correlation


check.pwCor.over.genes <- function(pwCor_all, pwCor_tum, genes_A, genes_B=NULL){
	
	stopifnot(all(rownames(pwCor_all)== rownames(pwCor_tum)))
	
	if(is.null(genes_B)) genes_B <- genes_A

	gene.idx.A <- rownames(pwCor_all) %in% genes_A
	gene.idx.B <- rownames(pwCor_all) %in% genes_B
	
	pwCor_all_sub <- pwCor_all[gene.idx.A, gene.idx.B]
	pwCor_tum_sub <- pwCor_tum[gene.idx.A, gene.idx.B]
	
	if(nrow(pwCor_all_sub) == ncol(pwCor_all_sub)){
		#if it is a square matrix 
		pwCor_all_vec <- pwCor_all_sub[upper.tri(pwCor_all_sub)]
		pwCor_tum_vec <- pwCor_tum_sub[upper.tri(pwCor_tum_sub)]
	}
	else{
		pwCor_all_vec <- as.vector(pwCor_all_sub)
		pwCor_tum_vec <- as.vector(pwCor_tum_sub)
	}
	
	cat("correlation=", 
		cor(pwCor_all_vec, pwCor_tum_vec, use="complete"),
		"\n")
		
	print(t.test (abs(pwCor_all_vec)-abs(pwCor_tum_vec),alternative="greater"))
	print(t.test (abs(pwCor_all_vec)-abs(pwCor_tum_vec)))

	
	return(list(pwCor_all_vec= pwCor_all_vec, pwCor_tum_vec= pwCor_tum_vec))
}


dat <- check.pwCor.over.genes (pwCor_all, pwCor_tum, markers.all.top100$tumor)

# correlation= 0.9746704

        # One Sample t-test

# data:  abs(pwCor_all_vec) - abs(pwCor_tum_vec)
# t = -4.1562, df = 4949, p-value = 1
# alternative hypothesis: true mean is greater than 0
# 95 percent confidence interval:
 # -0.1938161        Inf
# sample estimates:
 # mean of x
# -0.1388531


        # One Sample t-test

# data:  abs(pwCor_all_vec) - abs(pwCor_tum_vec)
# t = -4.1562, df = 4949, p-value = 3.291e-05
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
 # -0.20434927 -0.07335698
# sample estimates:
 # mean of x
# -0.1388531


dat <- check.pwCor.over.genes (pwCor_all, pwCor_tum, markers.all.top100$tumor,
					   unlist(markers.all.top100[c('myeloid','stromal','endo')]))
# correlation= 0.8715495

        # One Sample t-test

# data:  abs(pwCor_all_vec) - abs(pwCor_tum_vec)
# t = 15.146, df = 19999, p-value < 2.2e-16
# alternative hypothesis: true mean is greater than 0
# 95 percent confidence interval:
 # 0.3093338       Inf
# sample estimates:
# mean of x
# 0.3470234


        # One Sample t-test

# data:  abs(pwCor_all_vec) - abs(pwCor_tum_vec)
# t = 15.146, df = 19999, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
 # 0.3021129 0.3919339
# sample estimates:
# mean of x
# 0.3470234

dat <- check.pwCor.over.genes (pwCor_all, pwCor_tum, 
					   unlist(markers.all.top100[c('myeloid','stromal','endo')]),
					   unlist(markers.all.top100[c('myeloid','stromal','endo')]))
					   
# correlation= 0.7195284

        # One Sample t-test

# data:  abs(pwCor_all_vec) - abs(pwCor_tum_vec)
# t = 94.495, df = 19899, p-value < 2.2e-16
# alternative hypothesis: true mean is greater than 0
# 95 percent confidence interval:
 # 3.266412      Inf
# sample estimates:
# mean of x
  # 3.32428


        # One Sample t-test

# data:  abs(pwCor_all_vec) - abs(pwCor_tum_vec)
# t = 94.495, df = 19899, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
 # 3.255325 3.393235
# sample estimates:
# mean of x
  # 3.32428









pdf("pwCor_dens.snrnaCB_4CTMarker.pdf", useDingbats=FALSE)

my.marker <- markers.all.top100

range.values <- c(as.vector(pwCor_all), as.vector(pwCor_tum))							
range.values <- range.values[!is.na(range.values)]
range.values <- c(-50,50)

dat <- check.pwCor.over.genes (pwCor_all, pwCor_tum, 
								my.marker $tumor)
								
plot.df <- rbind.data.frame( data.frame(Z= dat[[1]], decon="all", type="tumor.vs.tumor"),
								data.frame(Z= dat[[2]], decon="tum", type="tumor.vs.tumor"))
								
densScatterplot(x1= dat[[1]], 
				x2= dat[[2]], 
				main="tumor vs tumor",
				xlab="pairwiseCor, all", ylab="pairwiseCor, deconvolved",
				range.values= range.values)

dat <- check.pwCor.over.genes (pwCor_all, pwCor_tum, 
								my.marker $tumor, 
								unlist(my.marker[c('myeloid','stromal','endo')]))
plot.df <- rbind.data.frame(plot.df,
			rbind.data.frame( data.frame(Z= dat[[1]], decon="all", type="tumor.vs.non-tumor"),
								data.frame(Z= dat[[2]], decon="tum", type="tumor.vs.non-tumor")))


densScatterplot(x1= dat[[1]], 
				x2= dat[[2]], 
				main="tumor vs non-tumor",
				xlab="pairwiseCor, all", ylab="pairwiseCor, deconvolved",
				range.values= range.values)

dat <- check.pwCor.over.genes (pwCor_all, pwCor_tum, 
								unlist(my.marker[c('myeloid','stromal','endo')]))
plot.df <- rbind.data.frame(plot.df,
			rbind.data.frame( data.frame(Z= dat[[1]], decon="all", type="non-tumor.vs.non-tumor"),
								data.frame(Z= dat[[2]], decon="tum", type="non-tumor.vs.non-tumor")))

densScatterplot(x1= dat[[1]], 
				x2= dat[[2]], 
				main="non-tumor vs non-tumor",
				xlab="pairwiseCor, all", ylab="pairwiseCor, deconvolved",
				range.values= range.values)

dev.off()


plot.df$Z <- abs(plot.df$Z)

#violin plots

p <- ggplot(plot.df, aes(x= type, y= Z, fill= decon)) +
			geom_violin(trim=FALSE, scale="width")+ 
			stat_summary(fun.data= data_summary,geom="pointrange", color="white", position=position_dodge(0.9)) +
			theme_classic() + ggtitle("Hotspot pairwise correlation") 

pdf("pwCor.violin.snranCB_4CTMarker.pdf",pointsize=8,useDingbats=FALSE)
print(p)
dev.off()


#plot delta 

dat <- check.pwCor.over.genes (pwCor_all, pwCor_tum, 
								my.marker $tumor)

plot.df <-  data.frame(delta_Z= abs(dat[[2]])-abs(dat[[1]]), type="tumor.vs.tumor")
								
			

dat <- check.pwCor.over.genes (pwCor_all, pwCor_tum, 
								my.marker $tumor, 
								unlist(my.marker[c('myeloid','stromal','endo')]))
plot.df <- rbind.data.frame(plot.df, 
			data.frame(delta_Z= abs(dat[[2]])-abs(dat[[1]]), type="tumor.vs.non-tumor")) 

dat <- check.pwCor.over.genes (pwCor_all, pwCor_tum, 
								unlist(my.marker[c('myeloid','stromal','endo')]))

plot.df <- rbind.data.frame(plot.df, 
			data.frame(delta_Z= abs(dat[[2]])-abs(dat[[1]]), type="non-tumor.vs.non-tumor")) 





p <- ggplot(plot.df, aes(x= type, y= delta_Z, fill= type)) +
			geom_violin(trim=FALSE, scale="width")+ 
			stat_summary(fun.data= data_summary,geom="pointrange", color="white", position=position_dodge(0.9)) +
			theme_classic() + ggtitle("Abs, Hotspot auto-correlation deconovolved-undeconvolved") + scale_y_continuous(trans = 'pseudo_log',breaks = c(-50, -25, -10, 0,10,25) ) +
  coord_flip() + geom_hline(yintercept = 0, linetype="dashed", color = "red")
  
  
ggsave(filename = "pwCor.delta_violin.snranCB_4CTMarker.pdf", plot = p, width = 10, height = 6)


# pdf("pwCor.delta_violin.snranCB_4CTMarker.pdf",pointsize=8,useDingbats=FALSE,paper="A4r")

# p <- ggplot(plot.df, aes(x= type, y= delta_Z, fill= type)) +
			# geom_violin(trim=FALSE, scale="width")+ 
			# stat_summary(fun.data= data_summary,geom="pointrange", color="white", position=position_dodge(0.9)) +
			# theme_classic() + ggtitle("Abs, Hotspot auto-correlation deconovolved-undeconvolved") 
# print(p) + geom_hline(yintercept = 0, linetype="dashed", color = "red")

# dev.off()







# # 
# nontumor.idx <- rownames(pwCor_all) %in% c(markers.all.symbol$myeloid, markers.all.symbol$stromal, markers.all.symbol$endo)

# selected.idx <- tumor.idx | nontumor.idx

# pdf("all.vs.bp_pw.pdf", useDingbats=FALSE)
# plot(x=as.vector(pwCor_all[selected.idx, selected.idx]), y=as.vector(pwCor_tum[selected.idx, selected.idx]), cex=0.1)
# plot(x=as.vector(pwCor_all[tumor.idx, tumor.idx]), y=as.vector(pwCor_tum[tumor.idx, tumor.idx]), cex=0.1, pch= 1,col="red")
# plot(x=as.vector(pwCor_all[tumor.idx, nontumor.idx]), y=as.vector(pwCor_tum[tumor.idx, nontumor.idx]), cex=0.1, pch= 1,col="purple")
# plot(x=as.vector(pwCor_all[nontumor.idx, nontumor.idx]), y=as.vector(pwCor_tum[nontumor.idx, nontumor.idx]), cex=0.1, pch= 1,col="blue")

# dev.off()





#compare sc vs st









