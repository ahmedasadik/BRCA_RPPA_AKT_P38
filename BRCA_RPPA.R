library(purrr)
library(ggplot2)
library(gridExtra)
library(reshape)
library(miceadds)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(ggbeeswarm)

## group functions
group_no_z_FUN <- function(x,y){
  x_med <- median(x)
  x_mad <- mad(x)
  group <- vector("character", length = length(x))
  group[which(x >= (x_med + (y*x_mad)))] <- "high"
  group[which(x < (x_med - (y*x_mad)))] <- "low"
  group
}

## load TCGA BRCA RNA_seq dataset
load.Rdata("./TCGA-BRCATranscriptome_ProfilingWed_Oct__4_171633_2017.RData", "brca")

## RPPA results
breast_rpp <- read.csv("./TCGA-BRCA-L4.csv", sep = ",", stringsAsFactors = F)
breast_rpp$new_sample_ID <- unlist(map(breast_rpp$Sample_ID, substring, first=1, last=16))

## Extract dataframes from the Summarized experiment
BRCA_coldata <- SummarizedExperiment::colData(brca)
cnt_data <- SummarizedExperiment::assay(brca)

## Filtering count data based on total counts per feature
samplesTP <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(cnt_data),typesample = c("TP")) 
counts_TP <- cnt_data[,match(samplesTP, colnames(cnt_data))]
counts_cpm <- edgeR::cpm(counts_TP, log = T)

## new BRCA colData matching only the primary tumor patients
new_BRCA_ColData <- BRCA_coldata[match(colnames(counts_TP),BRCA_coldata$barcode),]

## Targets
tor_targets <- c("AKT_pT308","P38_pT180Y182","P38MAPK","AKT_pS473","P70S6K_pT389", "MTOR_pS2448",
                     "AKT","P70S6K1","MTOR")

## RPPA_data_frame
ALL_df <- data.frame(pt_sample_ID=new_BRCA_ColData$sample, stringsAsFactors = F)
ALL_df[,tor_targets] <- matrix(nrow = dim(ALL_df)[1], ncol = length(tor_targets))
ALL_df[,tor_targets] <- breast_rpp[match(ALL_df$pt_sample_ID, breast_rpp$new_sample_ID),c(tor_targets)]
rppa_df <- ALL_df[!is.na(ALL_df$AKT_pT308),]

metadata_rppa <- new_BRCA_ColData[match(rppa_df$pt_sample_ID, new_BRCA_ColData$sample),]

rppa_df <- data.frame(rppa_df, sex=metadata_rppa$gender,
                      birth_Y=metadata_rppa$year_of_birth,
                      age=metadata_rppa$subtype_Age.at.Initial.Pathologic.Diagnosis,
                      days_to_death = metadata_rppa$days_to_death,
                      days_to_last = metadata_rppa$days_to_last_follow_up,
                      surv_days = metadata_rppa$subtype_OS.Time,
                      vital_status = metadata_rppa$vital_status,
                      vital_status_sb = metadata_rppa$subtype_Vital.Status,
                      surv_doa = metadata_rppa$subtype_OS.event,
                      DoA=metadata_rppa$subtype_OS.event,
                      stage=metadata_rppa$subtype_AJCC.Stage,
                      pam50=metadata_rppa$subtype_PAM50.mRNA,
                      ER=metadata_rppa$subtype_ER.Status,
                      PR=metadata_rppa$subtype_PR.Status,
                      Tnm=metadata_rppa$subtype_Tumor,
                      tNm=metadata_rppa$subtype_Node,
                      tnM=metadata_rppa$subtype_Metastasis,stringsAsFactors = F)

rppa_df$days_to_death2 <- rppa_df$days_to_death
rppa_df$days_to_death2 <- ifelse(!is.na(rppa_df$days_to_death), rppa_df$days_to_death, 
                                 ifelse(is.na(rppa_df$days_to_death) & !is.na(rppa_df$days_to_last), rppa_df$days_to_last,
                                        ifelse(is.na(rppa_df$days_to_last), rppa_df$surv_days, NA)))

## Combined rppa dataset with targets and groups
groups_AKT_T308 <- map2(list(rppa_df$AKT_pT308),c(0,0.5,1,1.5), group_no_z_FUN) %>% do.call(cbind,.)
colnames(groups_AKT_T308) <- paste("AKT_T308_", c(0,05,1,15), sep = "")
rppa_all <- data.frame(rppa_df, groups_AKT_T308, stringsAsFactors = F)
write.table(rppa_all,"./rppa_all.txt", sep = "\t")


## P38 group assignments in the AKT groups
rppa_all$AKT_pT308_hi_lo_p38_grps <- rppa_all$AKT_T308_0

rppa_all$AKT_pT308_hi_lo_p38_grps[which(rppa_all$AKT_pT308_hi_lo_p38_grps=="high")] <- paste(rppa_all$AKT_pT308_hi_lo_p38_grps[which(rppa_all$AKT_pT308_hi_lo_p38_grps=="high")],
                                                                                             group_no_z_FUN(rppa_all$P38_pT180Y182[which(rppa_all$AKT_pT308_hi_lo_p38_grps=="high")],0), sep = "_")

rppa_all$AKT_pT308_hi_lo_p38_grps[which(rppa_all$AKT_pT308_hi_lo_p38_grps=="low")] <- paste(rppa_all$AKT_pT308_hi_lo_p38_grps[which(rppa_all$AKT_pT308_hi_lo_p38_grps=="low")],
                                                                                            group_no_z_FUN(rppa_all$P38_pT180Y182[which(rppa_all$AKT_pT308_hi_lo_p38_grps=="low")],0), sep = "_")
## Relevel groups
new_levels <- c(high_high="high_high", high_low="high_low", low_high="low_high", low_low="low_low")
rppa_all$AKT_pT308_hi_lo_p38_grps <- factor(rppa_all$AKT_pT308_hi_lo_p38_grps, levels = new_levels, labels = names(new_levels))

## Custom labels and colors
new_labs <- c('Akt-pT308', 'p38-pT180/Y182', 'p38 (MAPK)', 'Akt-pS473', 'p70-S6K-pT389', 'mTOR-pS2448',
              'Akt', 'p70-S6K', 'mTOR')

new_cols <- c("5 126 151", "188 188 188", "242 97 37", "109 109 109")
new_cols_hex <- sapply(strsplit(new_cols, " "), function(x)
  rgb(x[1], x[2], x[3], maxColorValue=255))

## Generate plots
rppa_to_plot <- rppa_all[,c("pt_sample_ID","AKT_pT308_hi_lo_p38_grps",tor_targets)] %>% melt()

facet_plot <- ggboxplot(rppa_to_plot, x = "AKT_pT308_hi_lo_p38_grps", y = "value", legend="right",
                        font.xtick=c(10,"plain","black"), color = "AKT_pT308_hi_lo_p38_grps",
                        palette = new_cols_hex, x.text.angle = 45, xlab = FALSE, ylab="z - score",
                        outlier.size=0, facet.by = "variable", scales="free_y",
                        panel.labs = list(variable=new_labs))+
  ggbeeswarm::geom_quasirandom(aes_string(color="AKT_pT308_hi_lo_p38_grps"), size=0.2)+
  stat_compare_means(comparisons = list(c("high_high" ,"high_low"), c("low_high", "low_low")),
                     label = "p.signif", method = "wilcox")

pdf("./AKT_groups_targets_leg_top.pdf", width = 9, height = 12, compress = F, pointsize = 28)
ggpar(facet_plot, legend.title = "AKT_P38_groups", legend = "top")
dev.off()

