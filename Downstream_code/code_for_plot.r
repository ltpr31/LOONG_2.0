
library(dplyr)
library(ggplot2)
library(Seurat)
library(Matrix)
library(dplyr)
library(gridExtra)
library(magick)
library(jsonlite)
library(patchwork)
library(tidyverse)
library(grid)

setwd("path/to/work/path")
obj <- readRDS("./object.rds")


############################## add ALBST values：
meta <- NULL
meta <- obj@meta.data
head(meta)
meta$ID <- row.names(meta)
albst <- read.table("/path/to/out/albst.txt", header = T, sep = " ")  ## input ALBST result
meta <- merge(meta, albst, by = "ID", all.x = T, sort = F)
row.names(meta) <- meta$ID

start <- 0
end <- 10
step <- 0.1
thresholds <- seq(start, end, by = step)
for (i in 1:(length(thresholds) - 1)) {
  if (i == 1){
    meta$pos[meta$ALBST >= thresholds[i] & meta$ALBST <= thresholds[i + 1]] <- thresholds[i]
  }else{
    meta$pos[meta$ALBST > thresholds[i] & meta$ALBST <= thresholds[i + 1]] <- thresholds[i]
  }
}
unique(meta$pos)




############################# mark cell (optional)
mark <- read.table("/path/to/mark/result/select.txt",header = F, sep = " ")
mark$mark <- "target"
meta <- merge(meta, mark, by.x = "ID", by.y = "V1", all.x = T, sort = F)
rownames(meta) <- meta$ID
meta$mark <- ifelse(is.na(meta$mark), "Other", meta$mark)



#########################################
obj@meta.data <- meta[colnames(obj),]
head(obj@meta.data)
saveRDS(obj, file = "object_rearranged.rds")




##########################################################  plot part：
##########################################################  plot part：
obj <- readRDS("object_rearranged.rds")
DefaultAssay(obj) <- "RNA"  #or DefaultAssay(obj) <- "SCT"
meta <- obj@meta.data
meta$ID <- row.names(meta)




####################### ALBST feature plot:
p <- ggplot(meta, aes(x = col, y = row, colour = ALBST)) +
  geom_point(size = 0.25, stroke = 0) +
  theme_light()+
  scale_colour_gradientn( colours = c("#000080", "#4169E1", "#FFD700", "#FF8C00", "#B22222"),
                          values = scales::rescale(c(0, 0.3, 0.5, 0.7, 1))  )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())+
  coord_fixed(ratio = 1) +
  scale_y_reverse() +
  labs(x = "x", y = "y" )

w <- 7
h <- 6
pdf(file= "./albst_featureplot.pdf", width = w, height = h)
print(p)
dev.off()
png(filename="./albst_featureplot.png", width = w*300, height = h*300, res=300, units="px")
print(p)
dev.off()




################## ALBST feature plot with tissue image:
img <- as.raster(image_read("/path/to/tissue.png")) %>% 
  inset_element(left = 0, bottom = 0, right = 1,top = 1, 
                on_top = FALSE)
inf <- image_read("/path/to/tissue.png") %>% 
  image_info()
width <- inf$width
height <- inf$height
x_scale_factor <- 1  # adjust it to fit your practical condition
y_scale_factor <- 1  # adjust it to fit your practical condition
## Plot the coordinates on tissue image
p <- ggplot() +
  geom_point(coord, mapping = aes(x = col * x_scale_factor, y = row * y_scale_factor, color = ALBST),  size = 0.7) + 
  scale_x_continuous(limits = c(0, inf$width), expand = c(0, 0)) + 
  scale_y_reverse(limits = c(inf$height, 0), expand = c(0, 0)) +
  scale_colour_gradientn( colours = c("#000080", "#4169E1", "#FFD700", "#FF8C00", "#B22222"),
                          values = scales::rescale(c(0, 0.3, 0.5, 0.7, 1)))+
  coord_fixed(ratio = 1) + theme_void() + img
w <- 8
h <- 7
pdf(file="./albst_featureplot.pdf", width = w, height = h)
print(p)
dev.off()
png(filename="./albst_featureplot.png", width = w*300, height = h*300, res=300, units="px")
print(p)
dev.off()




#################### barplot for gene expression:
obj <- readRDS("object_rearranged.rds")
meta <- obj@meta.data
####scale data###
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 1e6)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = rownames(obj), verbose = TRUE, vars.to.regress = NULL)
## SCTransform result was used in macaque cortex, DefaultAssay(obj) <- "SCT"
idlist <- c("gene1", "gene2", "gene3")
##idlist <- c("Ang4", "B4galt1", "Hmgcs2")
##idlist <- c("LOC-11g13890")
##idlist <- c("TF")
for (gene in idlist) {
## get normalized expression level:
  data <- obj@assays[["RNA"]]@data  ## or: data <- obj@assays[["SCT"]]@data 
  data <- data[gene,]
  data <- as.data.frame(data)
  data$ID <- row.names(data)
  data <- merge(data, meta[,c("ID","ALBST","Width")], by = "ID", all.x = T, sort = F)
## segment data according to ALBST:
  start <- 0
  end <- 10
  step <- 0.2
  thresholds <- seq(start, end, by = step)
  for (i in 1:(length(thresholds) - 1)) {
    test <- NULL
    test <- meta[meta$ALBST > thresholds[i]  & meta$ALBST <= thresholds[i + 1],]
    if(dim(test)[1] > 0){
## subset data
      subset_obj <- NULL
      if (i == 1){
        subset_obj <- subset(obj, ALBST >= thresholds[i] & ALBST <= thresholds[i + 1])
      }else{
        subset_obj <- subset(obj, ALBST > thresholds[i] & ALBST <= thresholds[i + 1])
      }
      
      data$group[data$ID %in% colnames(subset_obj)] <- thresholds[i]
      avg_expression <- NULL
      avg_expression = AverageExpression(subset_obj,
                                         features = gene,
                                         slot = 'data')
      avg_expression = avg_expression[[1]]
      data$ave_exp[data$ID %in% colnames(subset_obj)] <- avg_expression[1,1]
      
    }else{
      next
    }
  }
  df_bar <- NULL
  df_bar <- data[,c("group","ave_exp")]
  df_bar <- df_bar %>% distinct(group, .keep_all = TRUE)
## fit values
  fit <- mgcv::gam(ave_exp ~ s(group, bs = "cs"), data = df_bar)
  df_bar$fit <- predict(fit, newdata = df_bar)
  df_bar <- arrange(df_bar, group)
  write.table(df_bar, paste("./", gene, "_data.txt", sep=""), col.names = T, row.names = F, sep = "\t", quote = F)
## barplot
p1 <- ggplot(df_bar, aes(x = as.numeric(as.character(group)), y = ave_exp)) +
    geom_col(aes(fill = ave_exp), width = 0.1) +
    theme_light()+
    scale_fill_gradient(low = "white", high = "#336699") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_line(aes(y = fit, colour = fit), size = 1)+
    scale_colour_gradient(high = "#FF9900", low = "#0000FF")+
    scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10))+
    theme(aspect.ratio = 0.2)+
    ggtitle(gene)+
    labs(x = "ALBST", y = "Average expression" )
## scatter plot
p2 <- ggplot(data, aes(x = as.numeric(as.character(ALBST)), y = Width, colour = data)) +
    geom_point(size = 0.8, stroke = 0) +
    theme_light()+
    scale_color_gradientn( colours = c('black','#330066',"#0000FF","#FF9900") )+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())+
    scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10))+
    coord_fixed(ratio = 1)+
    labs(x = "ALBST", y = "Width" )
  w <- 7
  h <- 6
  pdf(file=paste("./", gene, "_albst_barplot.pdf", sep=""), width = w, height = h)
  grid.arrange(p1, p2, ncol = 1)
  dev.off()
  png(filename=paste("./", gene, "_albst_barplot.png", sep=""), width = w*300, height = h*300, res=300, units="px")
  grid.arrange(p1, p2, ncol = 1)
  dev.off()
}




########################### Stacked barplot for cell proportion:
obj <- readRDS("object_rearranged.rds")
meta <- obj@meta.data
start <- 0
end <- 10
step <- 0.2
thresholds <- seq(start, end, by = step)
df_bar <- data.frame(
  group = numeric(),
  celltype = character(),
  ratio = numeric(),
  stringsAsFactors = FALSE
)
for (i in 1:(length(thresholds) - 1)) {
## subset cells
  if (i == 1) {
    subset_cells <- meta[meta$ALBST >= thresholds[i] & meta$ALBST <= thresholds[i + 1], ]
  } else {
    subset_cells <- meta[meta$ALBST > thresholds[i] & meta$ALBST <= thresholds[i + 1], ]
  }
## cell proportion
  if (nrow(subset_cells) > 0) {
    celltype_counts <- subset_cells %>%
      group_by(SubClass) %>%  ## select a class according to your requirement
      summarise(count = n(), .groups = 'drop')
    total_count <- sum(celltype_counts$count)
    celltype_counts$ratio <- celltype_counts$count / total_count
    df_bar <- rbind(df_bar, data.frame(
      group = thresholds[i],
      celltype = celltype_counts$SubClass,
      ratio = celltype_counts$ratio
    ))
  }
}
write.table(df_bar, "./Cell_proportion.txt", col.names = T, row.names = F, sep = "\t", quote = F)
df_bar$celltype <- factor(df_bar$celltype, levels = unique(df_bar$celltype) )
require(RColorBrewer)
mycol=c(brewer.pal(8,"Set1"),brewer.pal(7,"Set2"),brewer.pal(8,"Set3"),brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(7,"Accent"))
p <- ggplot(df_bar, aes(x = group, y = ratio, fill = celltype)) +
  geom_col(width = 0.2, position = "stack") + 
  scale_fill_manual(values = mycol) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + 
  scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10)) +
  theme_light() +
  labs(x = "ALBST", y = "Cell Proportion", title =  "SubClass Proportion") +
  ggtitle("cell proportion") +
  theme(panel.grid = element_blank(), aspect.ratio = 0.5)
w <- 7
h <- 8
pdf(file = "Subclass_ratio.pdf", width = w, height = h)
print(p)
dev.off()
png(filename= "Subclass_ratio.png", width = w*600, height = h*600, res=600, units="px")
print(p)
dev.off()




############################# Barplot for specific cell type:
df_specific <- df_bar[df_bar$celltype == "OLG",]  ## select a cell type according to your requirement
fit <- mgcv::gam(ratio ~ s(group, bs = "cs"), data = df_specific)
df_specific$fit <- predict(fit, newdata = df_specific)
df_specific <- arrange(df_specific, group)
p <- ggplot(df_specific, aes(x = as.numeric(as.character(group)), y = ratio)) +
  geom_col(aes(fill = ratio), width = 0.2) +
  scale_fill_gradient(low = "white", high = "#336699") +
  geom_line(data = df_specific, aes(x = group, y = fit, colour = fit), size = 1, inherit.aes = FALSE) +
  scale_colour_gradient(high = "#FF9900", low = "#0000FF") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10)) +
  theme_light() +
  labs(x = "ALBST", y = "Cell Proportion") +
  ggtitle("Specific proportion") +
  theme(panel.grid = element_blank(), aspect.ratio = 0.2)
w <- 7
h <- 8
pdf(file = "Specific_ratio.pdf", width = w, height = h)
print(p)
dev.off()
png(filename= "Specific_ratio.png", width = w*600, height = h*600, res=600, units="px")
print(p)
dev.off()




########################### Stacked barplot for subtypes:
start <- 0
end <- 10
step <- 0.2
thresholds <- seq(start, end, by = step)
df_bar <- data.frame(
  group = numeric(),
  celltype = character(),
  ratio = numeric(),
  stringsAsFactors = FALSE
)
meta_OLG <- meta[meta$SubClass == "OLG",]  ## select a cell group according to your requirement
for (i in 1:(length(thresholds) - 1)) {
## subset cells
  if (i == 1) {
    subset_cells <- meta_OLG[meta_OLG$ALBST >= thresholds[i] & meta_OLG$ALBST <= thresholds[i + 1], ]
  } else {
    subset_cells <- meta_OLG[meta_OLG$ALBST > thresholds[i] & meta_OLG$ALBST <= thresholds[i + 1], ]
  }
## cell proportion
  if (nrow(subset_cells) > 0) {
    celltype_counts <- subset_cells %>%
      group_by(celltype) %>%  ## select a class according to your requirement
      summarise(count = n(), .groups = 'drop')
    total_count <- sum(celltype_counts$count)
    celltype_counts$ratio <- celltype_counts$count / total_count
    df_bar <- rbind(df_bar, data.frame(
      group = thresholds[i],
      celltype = celltype_counts$celltype,
      ratio = celltype_counts$ratio
    ))
  }
}
write.table(df_bar, "./Subtypes_proportion.txt", col.names = T, row.names = F, sep = "\t", quote = F)
df_bar$celltype <- factor(df_bar$celltype, levels = c("OLG.7", "OLG.6", "OLG.4", "OLG.9","OLG.10", "OLG.5", "OLG.3","OLG.8" ))  ## adjust it to fit your practical condition
mycol <- c("#FF7F00","#FFFF33","#FDC086",  "#A6D854",  "#F781BF", "#377EB8", "#D21B1B","#A65628",   "#984EA3", "#4DAF4A", "#FC8D62","#8DA0CB")
p <- ggplot(df_bar, aes(x = group, y = ratio, fill = celltype)) +
  geom_col(width = 0.2, position = "stack") + 
  scale_fill_manual(values = mycol) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + 
  scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10)) +
  theme_light() +
  labs(x = "ALBST", y = "Cell Proportion" ) +
  ggtitle("Subtypes proportion") +
  theme(panel.grid = element_blank(), aspect.ratio = 0.5)
w <- 7
h <- 8
pdf(file = "Subtypes_ratio.pdf", width = w, height = h)
print(p)
dev.off()
png(filename= "Subtypes_ratio.png", width = w*600, height = h*600, res=600, units="px")
print(p)
dev.off()

