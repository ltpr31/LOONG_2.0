

library(Seurat)
library(monocle)
library(ggpubr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

setwd(/path/to/work/path)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
options(future.globals.maxSize = 500 * 1024^3)
height <- 7
width <- 6


for (tissue in c("d0", "d14")) {
stRNA <- readRDS(paste("./",tissue,".rds", sep = ""))
####################### add ALBST values:
##Ignore this step, if you have added ALBST to your RDS file.
meta <- NULL
meta <- stRNA@meta.data
head(meta)
meta$ID <- row.names(meta)
albst <- read.table(paste("/path/to","/",tissue,"/","out/albst.txt", sep = ""), header = T, sep = " ")  ## input ALBST values
meta <- merge(meta, albst, by = "ID", all.x = T, sort = F)
row.names(meta) <- meta$ID
stRNA@meta.data <- meta[colnames(stRNA),]
head(stRNA@meta.data)



####################### subset spots with Width (optional)
##Ignore this step, if it is not necessary for your tissue.
meta <- stRNA@meta.data
subset_spots <- meta[meta$Width > 0,]
stRNA <- subset(stRNA, cells = rownames(subset_spots)  )



####################### create monocle dataset
DefaultAssay(stRNA) <- "RNA"
message(DefaultAssay(stRNA))
assay="RNA"
head(stRNA@meta.data)
expr_matrix <- stRNA@assays[[assay]]@counts
p_data <- stRNA@meta.data
head(p_data)
f_data <- data.frame(gene_short_name = row.names(stRNA),row.names = row.names(stRNA))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
		phenoData = pd,
		featureData = fd,
		lowerDetectionLimit = 0.5,
		expressionFamily = negbinomial.size())
cds<-estimateSizeFactors(cds)
cds <- estimateDispersions(cds)



####################### infer pseudotime
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
high_var_genes_info <- disp_table[disp_table$gene_id %in% disp.genes, ]
write.table(high_var_genes_info, file=paste("./",tissue,".monocle_highvarible.Gene.info.tsv",sep=""), row.names = F, col.names = T, quote = F)
p1<-plot_ordering_genes(cds)
pdf(file=paste("./",tissue,".selected_order_gene.pdf",sep=""), height=height, width=width)
print(p1)
dev.off()
png(filename=paste("./",tissue,".selected_order_gene.png",sep=""), height=height*300, width=width*300, res=300, units="px")
print(p1)
dev.off()
rm(stRNA)
cds <- reduceDimension(cds, max_components = 2,
		method = 'DDRTree')
cds <- orderCells(cds)



################### plot cell with ALBST and Width
reduction <- p_data[,c("ALBST","Width")]
reduction <- t(reduction)
cds@reducedDimS <- reduction
cds@reducedDimK <- reduction
saveRDS(cds, file = paste("./",tissue,".spatial_trajectory.rds",sep=""))
df<-data.frame(barcode=row.names(pData(cds)),pData(cds))
write.table(df,file=paste("./",tissue,".Pseudotime.metadata.tsv",sep=""),sep="\t",quote=F,row.names=F,col.names = T)
p=plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.3)
pdf(file=paste("./",tissue,".Pseudotime.pdf",sep=""), height=height*2/7, width=width)
p
dev.off()
png(filename=paste("./",tissue,".Pseudotime.png",sep=""), height=height*300*2/7, width=width*300, res=300, units="px")
p
dev.off()
p=plot_cell_trajectory(cds, color_by = "State", cell_size = 0.3)+scale_colour_manual(values =  col_vector)
pdf(file=paste("./",tissue,".State.pdf",sep=""), height=height*2/7, width=width)
p
dev.off()
png(filename=paste("./",tissue,".State.png",sep=""), height=height*300*2/7, width=width*300, res=300, units="px")
p
dev.off()



######################### replace pseudotime with ALBST：
pData(cds)$Pseudotime <- pData(cds)$ALBST



######################### identify ALBST related genes:
Time_diff <- differentialGeneTest(cds, cores = 1,fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff=arrange(Time_diff,qval)
write.table(data.frame(ID=rownames(Time_diff),Time_diff),file=paste("./",tissue,".albst_related_all.tsv",sep=""),sep="\t",quote=F,row.names=F,col.names = T)
Time_genes <- top_n(Time_diff, n = 2000, desc(qval)) %>% filter(qval < 0.01) %>% pull(gene_short_name) %>% as.character()
pdf(file=paste("./",tissue,".topALBST_related_genes_heatmap.pdf",sep=""), height=height*2, width=width)
p=plot_pseudotime_heatmap(cds[Time_genes,],
		num_clusters = 5,
		cores = 1,
		show_rownames = T, return_heatmap=T)
dev.off()
png(filename=paste("./",tissue,".topALBST_related_genes_heatmap.png",sep=""), height=height*300*4, width=width*300, res=300, units="px")
plot_pseudotime_heatmap(cds[Time_genes,],
		num_clusters = 5,
		cores = 1,
		show_rownames = T)
dev.off()
clusters <- cutree(p$tree_row, k = 5)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
write.table(data.frame(ID=rownames(clustering),clustering),file=paste("./",tissue,".topALBST_related_genes_heatmap_cluster.tsv",sep=""),sep="\t",quote=F,row.names=F,col.names = T)


}


