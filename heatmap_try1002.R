library(readr) 
library(readxl) 
library(tidyverse) 
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ComplexHeatmap) #Heatmap
library(circlize) #colours
library(RColorBrewer) # colours

# DATEN EINLESEN

rnaseq_data <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 1) 
genelist <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 4)
EGFR_patients <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 3)
TCGA_clinicaldata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-HNSC_Master_clinical data_v1.xlsx", sheet = 1)
progenydata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA_HNSCC_progeny.xlsx", sheet = 1)

# OBJEKTE VORBEREITEN

# NR.1 Heatmapdaten

heatmap <- rnaseq_data %>% 
  filter(rnaseq_data$Entrez_Gene_ID %in% genelist$GeneID_all)
colnames(heatmap) <- gsub("-01","",colnames(heatmap), fixed=TRUE)
rownames(heatmap) <- heatmap$Entrez_Gene_ID 
heatmap <- as.data.frame(heatmap)
heatmap <- heatmap[-1]
heatmap <- as.matrix.data.frame(heatmap)

heatmap <- log(heatmap+1)
heatmap <- scale(t(heatmap), center = T, scale = T)#subtracting mean and dividing with SD
heatmap <- t(heatmap)

clustering <- hclust(dist(t(heatmap), method = "euclidean"), method = "ward.D2") #cluster distance und cluster methode
t <- cutree(clustering, 4) #wie viele cluster
cluster <- as.data.frame(t) 
cluster$Patient_ID <- as.factor(rownames(cluster))
cluster <- cluster %>%
  dplyr::rename(`Cluster` = `t`)

# NR.2 TOP ANNOTATIONS
progdata <- progenydata
progdata <- as.data.frame(progdata)
rownames(progdata) <- progdata$...1
progdata$...1 <- NULL
progdata <- t(progdata)
progdata <- as.data.frame(progdata)
progdata <- progdata %>%
  dplyr::rename(`PROGENy` = `EGFR`)
progdata$Patient_ID <- rownames(progdata)
rownames(progdata) <- NULL
progdata <- progdata %>%
  mutate(progdata, PROGENy = case_when(PROGENy >= 0.5 ~ "high",
                                           PROGENy <= -0.5 ~ "low",
                                       (PROGENy < 0.5 & PROGENy > -0.5) ~ "moderate"))


EGFR_patients <- as.data.frame(EGFR_patients)
rownames(EGFR_patients) <- EGFR_patients$Patient_ID
EGFR_patients$Patient_ID <- NULL

TCGA_clinicaldata <- as.data.frame(TCGA_clinicaldata)
EGFR_patients$Patient_ID <- c(rownames(EGFR_patients))
TCGA_clinicaldata <- left_join(EGFR_patients, TCGA_clinicaldata, by = "Patient_ID")
TCGA_clinicaldata <- left_join(TCGA_clinicaldata, cluster, by = "Patient_ID")
TCGA_clinicaldata <- left_join(TCGA_clinicaldata, progdata, by = "Patient_ID")
rownames(TCGA_clinicaldata) <- TCGA_clinicaldata$Patient_ID
TCGA_clinicaldata$Patient_ID <- NULL

clinical <- TCGA_clinicaldata %>%
  dplyr::select(Cluster, Smoking, Alcohol, HPV16, OS_5y_event, Gender , PROGENy, EGFR) %>%
  mutate(OS_5y_event = as.factor(ifelse(OS_5y_event == "0", "Yes", "No"))) %>%
  as.data.frame(.)
clinical <- clinical %>%
  mutate(clinical, Cluster = case_when(Cluster == 1 ~ "B1",
                                        Cluster == 2 ~ "B2b",
                                       Cluster == 3 ~ "A1",
                                       Cluster == 4 ~ "B2a"))
clinical <- clinical %>%  
  mutate(clinical, OS_5y_event = case_when(OS_5y_event == "No" ~ "Yes",
                                           OS_5y_event == "Yes" ~ "No"))



# NR. 3 LEFT ANNOTATIONS

genelist <- as.data.frame(genelist)
genelist <- genelist[match(rownames(heatmap),genelist$GeneID_all),]
genelist$Group[is.na(genelist$Group)] <- "up"
rownames(genelist) <- genelist$GeneID_all
genelist$GeneID_all <- NULL
genelist <- genelist %>%
  dplyr::rename(`Regulation` = `Group`)

# FARBEN BESTIMMEN
mypal <- brewer.pal(3, "Set1")
mypal2 <- brewer.pal(9, "Set1")

annocolours <- list(
  EGFR = c(high=mypal[2],low="purple", moderate=mypal[3]),
  PROGENy = c(high=mypal[2],low="purple", moderate=mypal[3]),
  `OS_5y_event` = c("Yes" =mypal[2], "No"=mypal[1]),
  Smoking = c(Yes=mypal[2],No =mypal[1]),
  HPV16 = c(Positive =mypal2[2], Negative= mypal[1]),
  Alcohol= c(Yes =mypal[2],No=mypal[1]),
  Gender = c(Male =mypal[2], Female =mypal[1]),
  Cluster = c(A1 = "green", B1 = "red1", B2a = "purple", B2b = "blue"))

col_rna = colorRamp2(c(-4, 0, 4), c("blue", "#EEEEEE", "red"), space = "LAB", transparency = 0)

# ANNOTATIONEN DEFINIEREN

gene_annotation = HeatmapAnnotation(df = genelist, which = "row", col = list(Regulation = c(up=mypal[2], down=mypal[1])),annotation_name_side = "top",annotation_name_gp = gpar(fontsize = 12))
clinical_annotations = HeatmapAnnotation(df = clinical, which = "column", col = annocolours, annotation_name_side = "right",annotation_name_gp = gpar(fontsize = 12))

# HEATMAP ERSTELLEN

HM <- ComplexHeatmap::Heatmap(heatmap,                              
                              col = col_rna,
                              name = "Expression",
                              clustering_distance_rows = "euclidean",
                              clustering_method_rows = "ward.D2",
                              #cluster_rows = TRUE,
                              #cluster_row_slices = TRUE,
                              #column_split =  EGFR_patients$EGFR,
                              clustering_distance_columns = "euclidean",
                              clustering_method_columns = "ward.D2",
                              column_title =  c("", "", "", ""),
                              column_title_side = "bottom",
                              column_title_gp = gpar(fontsize = 8),
                              row_split = 2,
                              column_split= 4,
                              #   row_split = genelist$Group,
                              #column_split = EGFR_patients$EGFR,
                              row_title = c("", ""),
                              row_title_rot = 0,
                              #row_title_gp = gpar(fontsize = 8),
                              show_column_names = FALSE,
                              show_row_names = FALSE, 
                              #use_raster = TRUE, 
                              #raster_device = "png",
                              top_annotation = clinical_annotations,
                              left_annotation = gene_annotation
)
draw(HM)
