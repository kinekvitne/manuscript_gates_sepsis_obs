# set wd
setwd("~/Desktop/manuscript_gates_sepsis_obs")

library(tidyverse)
library(mixOmics)
library(ggpubr)
library(vegan)
library(caret)
library(limma) 
library(patchwork)
library(rstatix)
library(Spectra) 
library(MsBackendMgf) 
library(homologueDiscoverer)
library(lme4)
library(lmerTest)
library(ComplexUpset)
library(viridis)

palette_del <- c("#6A66A3", "#84A9C0")
palette_hmo <- c("#80B6A3", "#E9946E")
palette_hmo_rev <- rev(palette_hmo)
palette_wat <- c("#4D6A6D", "#563440")
palette_wat_rev <- rev(palette_wat)
palette_stage <- c("Colostrum" = "#441b57", "Transitional" = "#78c25d", "Mature" = "#fce621")


###########################################################################################
# Function for PEG removal
plotAnnotatedStatic <- function(annotated, legend_setting = "bottom"){
  annotated <- mutate(annotated,
                      homologue_id = if_else(is.na(homologue_id), 0L, homologue_id),
                      homologue_id = as.factor(homologue_id))
  colourCount = length(unique(annotated$homologue_id))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  ncolor = length(unique(annotated$homologue_id))
  annotated <- arrange(annotated, desc(homologue_id))
  g <- ggplot(annotated, aes(group = homologue_id)) +
    geom_point(aes(x = rt, y = mz, color = homologue_id, size = homologue_id,
                   shape = homologue_id, alpha = homologue_id)) +
    geom_line(data = filter(annotated, homologue_id != 0),
              aes(x = rt, y = mz, group = homologue_id), alpha = 0.4, size = 0.1) +
    ggtitle("Annotated Peak Table") +
    scale_colour_manual(values = c("black", getPalette(colourCount-1))) +
    scale_alpha_manual(values =  c(0.2, rep(0.85, times = ncolor))) +
    scale_size_manual(values = c(0.1, rep(1.5, times = ncolor))) +
    scale_shape_manual(values = c(1, rep(19, times = ncolor))) +
    xlab("Retention Time (s)") +
    ylab("Mass to Charge Ratio") +
    theme(legend.position=legend_setting, text = element_text(family="mono"))
  return(g)
}

# Read in data
feature_table <- read_csv("data/obs_others_iimn_fbmn_quant.csv")
colnames(feature_table)[3] <- "RT"
metadata <- read_csv("data/qiita_metadata.csv") 
sequence_fecal <- read_csv("data/sequence_stool_plasma.csv")
sequence_milk <- read_csv("data/sequence_milk.csv")
sequence <- rbind(sequence_fecal, sequence_milk)

annotations <- read.delim("data/merged_results_with_gnps_all_libs_obs.tsv") %>%
  dplyr::filter(!str_detect(pattern = "REFRAME", LibraryName)) # remove drug library
annotations$X.Scan. <- as.character(annotations$X.Scan.)

info_feature <- feature_table %>% dplyr::select(1:3,7)
colnames(info_feature) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature$Feature <- as.character(info_feature$Feature)

info_feature_complete <- info_feature %>% 
  left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  dplyr::select(1:4,18,24) %>% 
  dplyr::filter(RT >= 0.20) 

# Data table
data <- feature_table %>% 
  dplyr::filter(RT >= 0.20) %>%
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE)

data$SampleID <- gsub(".mzML Peak area", "", data$SampleID)

# Metadata
metadata_metabolomics <- tibble(SampleID = data$SampleID) %>% 
  left_join(metadata %>% mutate(SampleID = as.character(SampleID)), by = "SampleID") %>% 
  left_join(sequence, by = c("SampleID" = "File Name")) %>%
  mutate(Stage = case_when(
    infant_age_days < 5 ~ "Colostrum",
    infant_age_days > 15 ~ "Mature",
    TRUE ~ "Transitional"
  )) %>% 
  dplyr::filter(!str_detect(SampleID, "Blank|QC")) %>% 
  mutate(water_treatment = case_when(
    drink_water_safe == "nothing" ~ "untreated",
    TRUE ~ "treated"
  )) %>% 
  mutate(AssetIndex2 = case_when(
    AssetIndexQuintile == "Quintile 1" ~ "Q12",
    AssetIndexQuintile == "Quintile 2" ~ "Q12",
    AssetIndexQuintile == "Quintile 3" ~ "Q35",
    AssetIndexQuintile == "Quintile 4" ~ "Q35",
    AssetIndexQuintile == "Quintile 5"~ "Q35",
  )) %>%  
  mutate(water_boil = case_when(
    str_detect(pattern = "boil", drink_water_safe) ~ "Yes_boil", TRUE ~ "No")) %>%
  mutate(water_filter = ifelse(water_filter == 1, "Yes_filter", "No")) 

### BASELINE CHARACTERISTICS ###
baseline <- metadata_metabolomics %>% 
  group_by(part_id) %>% slice(1) %>%
  ungroup() %>% 
  dplyr::filter(!SampleID == "41875") # Remove sample without metadata

table(baseline$sex)
mean(baseline$birthweight)
sd(baseline$birthweight)
mean(baseline$gestational_age_birth)/7
sd(baseline$gestational_age_birth)/7
mean(baseline$mother_age_years)
sd(baseline$mother_age_years)
table(baseline$mode_delivery)
table(baseline$hmo_Secretor)
table(baseline$maternal_antibiotics)
table(baseline$AssetIndexQuintile)
table(baseline$water_treatment)


# Investigate total peak area
data_TIC <- data.frame(TIC = rowSums(data %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

data_TIC %>% 
  ggscatter("Order", "TIC", add = "reg.line") +
  stat_cor() 

# Check sample type
sample_stool_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "feces", sample_type))
sample_plasma_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "blood", sample_type))
sample_milk_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "milk", sample_type))

six_stool_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Mix", SampleID))
blank_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Blank", SampleID))

# Check TIC overtime
sample_stool_tic %>%
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 2e9) +
  stat_cor() # remove one sample (113485)

sample_plasma_tic %>%
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 1e9) +
  stat_cor()

sample_milk_tic %>%
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 1e9) +
  stat_cor()

# Check internal standard and remove sample where there was a shift
fbmn_IS <- annotations %>% dplyr::filter(str_detect(Compound_Name, regex("Sulfam|Sulfad", ignore_case = TRUE))) %>% 
  distinct(X.Scan., .keep_all = TRUE) %>% dplyr::filter(Organism != "BILELIB19")

# Extract IS features for the table
table_IS <- data %>% column_to_rownames("SampleID") %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% 
  dplyr::filter(ID %in% fbmn_IS$X.Scan.) %>% column_to_rownames("ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("SampleID") %>% dplyr::filter(!(str_detect(SampleID, "Blank|Mix")))
# IS is missing in most of the samples. cannot check CV

# Check features per sample type
data_blank <- data %>% dplyr::filter(str_detect(pattern = "Blank", SampleID))
data_sixmix <- data %>% dplyr::filter(str_detect(pattern = "Mix", SampleID))
data_sample <- data %>% dplyr::filter(!(str_detect(pattern = "Mix|Blank", SampleID)))

# Blank
blanks_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                  Mean_blank = data_blank %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_blank =  data_blank %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_blank, SD_blank, CV_blank) %>% 
  dplyr::filter(Mean_blank > 0) %>% arrange(desc(Mean_blank))

# Six mix
sixmix_feature_info <- data.frame(Feature = colnames(data_sixmix)[-1],
                                  Mean_sixmix = data_sixmix %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_sixmix = data_sixmix %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sixmix = SD_sixmix/Mean_sixmix) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_sixmix, SD_sixmix, CV_sixmix) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% arrange(desc(Mean_sixmix))
# there is a shift in RT of standards. 

# Sample
sample_feature_info <- data.frame(Feature = colnames(data_sample)[-1],
                                  Mean_sample = data_sample %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_sample =  data_sample %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_sample, SD_sample, CV_sample) %>% 
  dplyr::filter(Mean_sample > 0) %>% arrange(desc(Mean_sample))


# Features to be removed Sample/Blank < 5
feature_to_remove <- blanks_feature_info %>% left_join(sample_feature_info) %>%
  dplyr::filter(Mean_blank > 0) %>% 
  dplyr::mutate(sample_Blank = Mean_sample/Mean_blank) %>% 
  dplyr::filter(sample_Blank < 5 | is.na(sample_Blank))

# Data with blank removal
data_clean <- data %>% 
  dplyr::filter(SampleID != "113485") %>% # remove sample with high TIC
  dplyr::select(-c(feature_to_remove$Feature)) %>%
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")

# Features to be removed Sample/QCmix < 5
feature_to_remove_mix <- sixmix_feature_info %>% left_join(sample_feature_info) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% 
  dplyr::mutate(sample_Mix = Mean_sample/Mean_sixmix) %>% 
  dplyr::filter(sample_Mix < 5 | is.na(sample_Mix)) %>% 
  dplyr::filter(!(Feature %in% feature_to_remove$Feature))

# Data with QCmix removal
data_clean2 <- data_clean %>% dplyr::select(-c(feature_to_remove_mix$Feature))

# Remove feature before 0.2 minutes and after 8 minutes
feature_to_remove_rt <- info_feature_complete %>% dplyr::filter(RT < 0.2 | RT > 8) %>%
  dplyr::filter(!(Feature %in% feature_to_remove$Feature))  %>%
  dplyr::filter(!(Feature %in% feature_to_remove_mix$Feature))

# Final cleaned table
data_clean3 <- data_clean2 %>% dplyr::select(-c(feature_to_remove_rt$Feature))


# PCA raw data
PCA_raw <- mixOmics::pca(data_clean3 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

i <- "sample_type"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Whole", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()


#######################################################################################
# Remove plastic contamination
feature_check <- info_feature %>% dplyr::select(1:3) %>% 
  dplyr::filter(Feature %in% colnames(data_clean3)) %>% as_tibble()
colnames(feature_check) <- c("peak_id", "mz", "rt")
feature_check$mz <- as.double(feature_check$mz)
feature_check$rt <- as.double(feature_check$rt)

feature_blank <- data_clean3 %>% 
  dplyr::filter(SampleID == "Blank_Solvent20") %>%
  column_to_rownames("SampleID") %>%
  t() %>% as.data.frame() %>% rownames_to_column("peak_id")

feature_check_blank <- feature_check %>% left_join(feature_blank)

feature_check_blank$peak_id <- as.integer(feature_check_blank$peak_id)
colnames(feature_check_blank)[4] <- "intensity"

# Remove PEGs
peak_table_PEG <- detectHomologues(feature_check_blank, mz_steps = c(44.0262),
                                   rt_min = 0.01, rt_max = 100, ppm_tolerance = 5, 
                                   min_series_length = 5, search_mode = "targeted", 
                                   step_mode = "increment", verbose = TRUE)

plotAnnotatedStatic(peak_table_PEG)

sdb <- sdbCreate(peak_table_PEG, sample_origin = "data")
sdb1_info <- info_feature_complete %>% dplyr::filter(Feature %in% sdb$sample_peak_id)
sdb$sample_peak_id <- as.character(sdb$sample_peak_id)

features_contaminat <- info_feature_complete %>% 
  dplyr::filter(Feature %in% sdb$sample_peak_id) %>%
  dplyr::filter(!Feature == "13518")

# Check one PEG
peg <- data_clean3 %>% dplyr::select(SampleID, `1075`) %>% 
  left_join(metadata_metabolomics, by = "SampleID")
colnames(peg)[2] <- "PEG"

peg %>% dplyr::mutate(LogPEG = log2(PEG+1)) %>% 
  ggscatter(x = "SampleID", y = "PEG") + scale_color_viridis_c()

# Clean data from contaminants
data_clean4 <- data_clean3 %>% 
  dplyr::select(-c(features_contaminat$Feature)) 

other_PEGs <- info_feature_complete %>% 
  dplyr::filter(str_detect(pattern = "glycol", Compound_Name))

data_clean5 <- data_clean4 %>%
  dplyr::select(-any_of(as.character(other_PEGs$Feature)))
#######################################################################################


# Keep only samples
data_sample_final <- data_clean5 %>% 
  dplyr::filter(!(str_detect(pattern = "Blank|Mix", SampleID)))

# RCLR transformation
data_sample_clr <- decostand(data_sample_final %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

i <- "sample_type"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = "All Samples",
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))


# Remove sample with NA description and the fecal clustering with blood
data_sample_filter <- data_sample_final %>%
  dplyr::filter(!(SampleID %in% c(41875, 33194)))


# RCLR transformation
data_sample_clr <- decostand(data_sample_filter %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

i <- "sample_type"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = "All Samples",
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic(), legend = "None") +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 6),
        axis.text = element_text(size = 8))

#ggsave(plot = PCA_plot, filename = "whole_metabolome.svg", device = "svg", dpi = "retina", height = 2, width = 2)

# PERMANOVA
dist_metabolites <- vegdist(data_sample_clr, method = "euclidean")
disper_type <- betadisper(dist_metabolites, PCA_whole_scores$sample_type)
anova(disper_type)
permanova <- adonis2(dist_metabolites ~ sample_type, PCA_whole_scores, na.action = na.omit, by = "terms")


# Separate the different sample types
sample_stool <- metadata_metabolomics %>% dplyr::filter(str_detect(pattern = "feces", sample_type))
sample_plasma <- metadata_metabolomics %>% dplyr::filter(str_detect(pattern = "blood", sample_type))
sample_milk <- metadata_metabolomics %>% dplyr::filter(str_detect(pattern = "milk", sample_type))


##### BILE ACIDS VALIDATION #####

BA_query <- read.delim("data/Gates_postFBMNvalidation_bileacids_2percent.tsv")
names(BA_query)[names(BA_query) == "X.Scan."] <- "Feature"
BA_query_filter <- BA_query %>% 
  dplyr::filter(str_detect(pattern = "Did not pass", query_validation)) %>% 
  dplyr::filter(str_detect(pattern = "GNPS-BILE-ACID-MODIFICATIONS", LibraryName))

info_feature_complete_filter <- info_feature_complete %>%
  mutate(Compound_Name = if_else(Feature %in% BA_query_filter$Feature, NA_character_, Compound_Name))

#write_csv(x = info_feature_complete_filter, file = "info_feature_complete_filter.csv")






####################
# MILK #
###################
data_milk <- data_sample_filter %>% 
  dplyr::filter(SampleID %in% sample_milk$SampleID) %>% 
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")

# RCLR transformation
data_milk_clr <- decostand(data_milk %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_milk <- mixOmics::pca(data_milk_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                          ncomp = 2, center = TRUE, scale = TRUE)
PCA_milk_scores <- data.frame(PCA_milk$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_milk_plots <- list()

for (i in c("mother_age_years", "gestational_age_birth", "mode_delivery", "water_filter",
            "sex", "hmo_Secretor")) {
  
  PCA_plot <- PCA_milk_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - milk", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_milk$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_milk$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_milk_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_milk_plots[[i]] <- PCA_plot
  
}

PCA_milk_plots_final <- wrap_plots(PCA_milk_plots, nrow = 3)

i <- "infant_age_days"

PCA_plot <- PCA_milk_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - milk", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_milk$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_milk$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""), legend = "none",
            ggtheme = theme_classic()) + scale_color_viridis_c() +
  theme(plot.title = element_text(size = 9), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(plot = PCA_plot, filename = "PCA_Milk_Age_leg.svg", device = "svg", dpi = "retina", width = 2, height = 2)

# PERMANOVA
dist_metabolites <- vegdist(data_milk_clr, method = "euclidean")
disper_delivery <- betadisper(dist_metabolites, PCA_milk_scores$infant_age_days)
anova(disper_delivery)
permanova <- adonis2(dist_metabolites ~ infant_age_days + hmo_Secretor + water_treatment + mode_delivery + host_subject_id,
                     PCA_milk_scores, na.action = na.omit, by = "terms")
# Shows significant effect of time; supports looking at lactation stages. Mode of delivery is significant, water_treatement and secretor status are borderline (p<0.10)


# PLSDA milk - maternal secretor status
PCA_milk_scores_filter <- PCA_milk_scores  %>% 
  dplyr::filter(hmo_Secretor != "unknown")

PLSDA_milk_secretor <- mixOmics::plsda(data_milk_clr %>% rownames_to_column("SampleID") %>%
                                         dplyr::filter(SampleID %in% PCA_milk_scores_filter$SampleID)%>%
                                         column_to_rownames("SampleID") %>%
                                         select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                                       PCA_milk_scores_filter$hmo_Secretor, ncomp = 2, scale = TRUE)
PLSDA_milk_secretor_scores <- data.frame(PLSDA_milk_secretor$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics) %>%
  mutate(hmo_Secretor = as.factor(hmo_Secretor))

PLSDA_milk_secretor_plot <- PLSDA_milk_secretor_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "hmo_Secretor", alpha = 0.6, title = "PLSDA milk secretor",
            xlab = paste("Component 1 (", round(PLSDA_milk_secretor$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_milk_secretor$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = palette_hmo_rev) +
  geom_point(data = PLSDA_milk_secretor_scores %>% group_by(hmo_Secretor) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = hmo_Secretor), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(plot = PLSDA_milk_secretor_plot, filename = "PLSDA_milk_secretor_plot.svg", device = "svg", dpi = "retina", width = 5, height = 2.8)

Loadings_milk_secretor <- plotLoadings(PLSDA_milk_secretor, plot = FALSE, contrib = "max")$X %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_milk_secretor <- perf(PLSDA_milk_secretor, validation = "Mfold", folds = 5, nrepeat = 99, progressBar = TRUE) 
#plot(perf_plsda_milk_secretor, legend = FALSE)

VIPs_milk_secretor <- as.data.frame(mixOmics::vip(PLSDA_milk_secretor))
VIPs_milk_secretor_filter <- dplyr::filter(VIPs_milk_secretor, VIPs_milk_secretor$comp1 > 1)
VIPs_milk_secretor_filter$ID <- rownames(VIPs_milk_secretor_filter)
VIPs_milk_secretor_select <- VIPs_milk_secretor_filter %>% dplyr::select(ID, comp1)
VIPs_milk_secretor_Load <- VIPs_milk_secretor_select %>% 
  left_join(Loadings_milk_secretor, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete_filter, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#VIPs_milk_secretor_Load2 <- VIPs_milk_secretor_Load
#VIPs_milk_secretor_Load2$GroupContrib <- ifelse(VIPs_milk_secretor_Load2$GroupContrib == 1, "secretor", "non-secretor")
#write_csv(x = VIPs_milk_secretor_Load2, file = "Human_milk_Secretor_status_VIP.csv")

# Check top most relevant features
VIPs_milk_seretoryes <- VIPs_milk_secretor_Load %>% dplyr::filter(GroupContrib == "1") %>% head(70)
VIPs_milk_seretorno <- VIPs_milk_secretor_Load %>% dplyr::filter(GroupContrib == "0") %>% head(70)

data_milk_ratio_hmo <- data_milk %>%
  dplyr::select("SampleID", VIPs_milk_seretoryes$ID, VIPs_milk_seretorno$ID) %>%
  dplyr::mutate(yes = rowSums(select(., VIPs_milk_seretoryes$ID))) %>%
  dplyr::mutate(no = rowSums(select(., VIPs_milk_seretorno$ID))) %>%
  dplyr::mutate(Ratio = log(yes/no + 1)) %>%
  left_join(metadata_metabolomics)

# Plot ratio
plot_ratio_milk_secretor <- data_milk_ratio_hmo %>% 
  dplyr::filter(hmo_Secretor %in% c(0,1)) %>%
  dplyr::mutate(hmo_Secretor = factor(hmo_Secretor, levels = c(0, 1))) %>%
  ggboxplot(x = "hmo_Secretor", y = "Ratio", add = "jitter",
            add.params = list(color = "hmo_Secretor", alpha = 0.6), color = "hmo_Secretor", 
            palette = palette_hmo_rev,
            xlab = "Secretor status",
            title = "Differential features from PLS-DA model") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

plot_ratio_milk_secretor_time <- data_milk_ratio_hmo %>% 
  dplyr::filter(hmo_Secretor %in% c(0,1)) %>%
  dplyr::mutate(hmo_Secretor = factor(hmo_Secretor, levels = c(1, 0))) %>%
  ggscatter(x = "infant_age_days", y = "Ratio", add = "loess", ylab = "Ln(yes/no)", alpha = 0.5,
            add.params = list(color = "hmo_Secretor", alpha = 0.6), color = "hmo_Secretor", 
            palette = palette_hmo, xlab = "Infant age (days)", legend = "none",
            title = "Milk Ratio - Secretor status") +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Linear mixed effect model
model <- data_milk_ratio_hmo %>% 
  dplyr::filter(hmo_Secretor %in% c(0,1)) %>% 
  lmer(formula = Ratio ~ hmo_Secretor + infant_age_days + (1|host_subject_id))
summary(model)

##### CREATE BOXPLOT FOR FEAT 583 (2FL M/Z 471.1709) AND 769 (LDFT, M/Z 657.2211)

# 2FL
data_milk_feat_2FL <- data_milk %>% 
  dplyr::select(SampleID, `583`) %>%
  dplyr::mutate(Log_feat = log2(`583` + 1)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

data_milk_feat_2FL_colostrum <- data_milk_feat_2FL %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_feat, hmo_Secretor) %>%
  dplyr::filter(infant_age_days < 5) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_feat == max(Log_feat, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_2FL_colostrum_plot <- data_milk_feat_2FL_colostrum %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_feat", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev, 
            xlab = "Colostrum", ylab = "Log2(peak area)") + stat_compare_means() + ylim (-0.1, 20) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) 

data_milk_feat_2FL_trans <- data_milk_feat_2FL %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_feat, hmo_Secretor) %>%
  dplyr::filter(infant_age_days >= 5  & infant_age_days < 16) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_feat == max(Log_feat, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_2FL_trans_plot <- data_milk_feat_2FL_trans %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_feat", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Transitional", ylab = FALSE) + stat_compare_means() + ylim (-0.1, 20) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_milk_feat_2FL_mature <- data_milk_feat_2FL %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_feat, hmo_Secretor) %>%
  dplyr::filter(infant_age_days >= 16 & infant_age_days < 90) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_feat == max(Log_feat, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_2FL_mature_plot <- data_milk_feat_2FL_mature %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_feat", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Mature", ylab = FALSE) + stat_compare_means() +  ylim (-0.1, 20) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))


combined_2FL_lactation_stages <- ggarrange(data_milk_feat_2FL_colostrum_plot, data_milk_feat_2FL_trans_plot,
                                           data_milk_feat_2FL_mature_plot, nrow = 1)

p1 <- compare_means(Log_feat ~ hmo_Secretor, data = data_milk_feat_2FL_colostrum)$p
p2 <- compare_means(Log_feat ~ hmo_Secretor, data = data_milk_feat_2FL_trans)$p
p3 <- compare_means(Log_feat ~ hmo_Secretor, data = data_milk_feat_2FL_mature)$p
p_raw <- c(p1, p2, p3)
p_adj <- p.adjust(p_raw, method = "BH") 
#ggsave(plot = combined_2FL_lactation_stages, filename = "combined_2FL_lactation_stages.svg", device = "svg", dpi = "retina", width = 2.6, height = 2.3)

# LDFT
data_milk_feat_LDFT <- data_milk %>% 
  dplyr::select(SampleID, `769`, `583`) %>%
  dplyr::mutate(Log_feat = log2(`769` + 1)) %>%
  dplyr::mutate(Log_2FL = log2(`583` + 1)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

data_milk_feat_LDFT_colostrum <- data_milk_feat_LDFT %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_feat, Log_2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days < 5) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_2FL == max(Log_2FL, na.rm = TRUE)) %>% 
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_LDFT_colostrum_plot <- data_milk_feat_LDFT_colostrum %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_feat", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev, 
            xlab = "Colostrum", ylab = "Log2(peak area)") + stat_compare_means() + ylim (-0.1, 22) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_milk_feat_LDFT_trans <- data_milk_feat_LDFT %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_feat, Log_2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days >= 5 & infant_age_days < 16) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_2FL == max(Log_2FL, na.rm = TRUE)) %>% 
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_LDFT_trans_plot <- data_milk_feat_LDFT_trans %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_feat", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Transitional", ylab = FALSE) + stat_compare_means() + ylim (-0.1, 22) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_milk_feat_LDFT_mature <- data_milk_feat_LDFT %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_feat, Log_2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days >= 16 & infant_age_days < 90) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_2FL == max(Log_2FL, na.rm = TRUE)) %>% 
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_LDFT_mature_plot <- data_milk_feat_LDFT_mature %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_feat", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Mature", ylab = FALSE) + stat_compare_means() +  ylim (-0.1, 22) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))


combined_LDFT_lactation_stages <- ggarrange(data_milk_feat_LDFT_colostrum_plot, data_milk_feat_LDFT_trans_plot,
                                           data_milk_feat_LDFT_mature_plot, nrow = 1)

p1 <- compare_means(Log_feat ~ hmo_Secretor, data = data_milk_feat_LDFT_colostrum)$p
p2 <- compare_means(Log_feat ~ hmo_Secretor, data = data_milk_feat_LDFT_trans)$p
p3 <- compare_means(Log_feat ~ hmo_Secretor, data = data_milk_feat_LDFT_mature)$p
p_raw <- c(p1, p2, p3)
p_adj <- p.adjust(p_raw, method = "BH") 
#ggsave(plot = combined_LDFT_lactation_stages, filename = "combined_LDFT_lactation_stages.svg", device = "svg", dpi = "retina", width = 2.6, height = 2.3)


### Targeted HMO ###
hmo_data <- read_csv("data/SEPSIS_HMO_Report_Obs_20221201.csv") %>% 
  dplyr::mutate(SampleID = as.character(SampleID)) %>%
  dplyr::left_join(metadata_metabolomics) %>% 
  dplyr::mutate(hmo_Secretor = as.factor(hmo_Secretor))

# Targeted 2FL
data_milk_feat_2FL_colostrum <- hmo_data %>% 
  dplyr::select(host_subject_id, infant_age_days, x2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days < 5) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(x2FL == max(x2FL, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_2FL_colostrum_plot <- data_milk_feat_2FL_colostrum %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "x2FL", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev, 
            xlab = "Colostrum", ylab = "nmol/mL") + stat_compare_means() + ylim (-0.1, 17000) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) 

data_milk_feat_2FL_trans <- hmo_data %>% 
  dplyr::select(host_subject_id, infant_age_days, x2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days >= 5 & infant_age_days < 16) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(x2FL == max(x2FL, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_2FL_trans_plot <- data_milk_feat_2FL_trans %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "x2FL", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Transitional", ylab = FALSE) + stat_compare_means() + ylim (-0.1, 17000) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_milk_feat_2FL_mature <- hmo_data %>% 
  dplyr::select(host_subject_id, infant_age_days, x2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days >= 16 & infant_age_days < 90) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(x2FL == max(x2FL, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_2FL_mature_plot <- data_milk_feat_2FL_mature %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "x2FL", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Mature", ylab = FALSE) + stat_compare_means() +  ylim (-0.1, 17000) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))


combined_2FL_lactation_stages <- ggarrange(data_milk_feat_2FL_colostrum_plot, data_milk_feat_2FL_trans_plot,
                                           data_milk_feat_2FL_mature_plot, nrow = 1)

p1 <- compare_means(x2FL ~ hmo_Secretor, data = data_milk_feat_2FL_colostrum)$p
p2 <- compare_means(x2FL ~ hmo_Secretor, data = data_milk_feat_2FL_trans)$p
p3 <- compare_means(x2FL ~ hmo_Secretor, data = data_milk_feat_2FL_mature)$p
p_raw <- c(p1, p2, p3)
p_adj <- p.adjust(p_raw, method = "BH") 
#ggsave(plot = combined_2FL_lactation_stages, filename = "combined_2FL_lactation_stages.svg", device = "svg", dpi = "retina", width = 2.8, height = 2.3)

# Targeted DFLac
data_milk_feat_DFLac_colostrum <- hmo_data %>% 
  dplyr::select(host_subject_id, infant_age_days, DFLac, x2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days < 5) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(x2FL == max(x2FL, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_DFLac_colostrum_plot <- data_milk_feat_DFLac_colostrum %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "DFLac", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev, 
            xlab = "Colostrum", ylab = "nmol/mL") + stat_compare_means() + ylim (-0.1, 6000) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) 

data_milk_feat_DFLac_trans <- hmo_data %>% 
  dplyr::select(host_subject_id, infant_age_days, DFLac, x2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days >= 5 & infant_age_days < 16) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(x2FL == max(x2FL, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_DFLac_trans_plot <- data_milk_feat_DFLac_trans %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "DFLac", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Transitional", ylab = FALSE) + stat_compare_means() + ylim (-0.1, 6000) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_milk_feat_DFLac_mature <- hmo_data %>% 
  dplyr::select(host_subject_id, infant_age_days, DFLac, x2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days >= 16 & infant_age_days < 90) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(x2FL == max(x2FL, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_milk_feat_DFLac_mature_plot <- data_milk_feat_DFLac_mature %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "DFLac", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Mature", ylab = FALSE) + stat_compare_means() +  ylim (-0.1, 6000) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))


combined_DFLac_lactation_stages <- ggarrange(data_milk_feat_DFLac_colostrum_plot, data_milk_feat_DFLac_trans_plot,
                                           data_milk_feat_DFLac_mature_plot, nrow = 1)

p1 <- compare_means(DFLac ~ hmo_Secretor, data = data_milk_feat_DFLac_colostrum)$p
p2 <- compare_means(DFLac ~ hmo_Secretor, data = data_milk_feat_DFLac_trans)$p
p3 <- compare_means(DFLac ~ hmo_Secretor, data = data_milk_feat_DFLac_mature)$p
p_raw <- c(p1, p2, p3)
p_adj <- p.adjust(p_raw, method = "BH") 
#ggsave(plot = combined_DFLac_lactation_stages, filename = "combined_DFLac_lactation_stages.svg", device = "svg", dpi = "retina", width = 2.8, height = 2.3)


# PLSDA milk - water treatment
PCA_milk_scores_water <- PCA_milk_scores  %>% 
  dplyr::filter(water_treatment != "NA")

PLSDA_milk_water <- mixOmics::plsda(data_milk_clr %>% rownames_to_column("SampleID") %>%
                                      dplyr::filter(SampleID %in% PCA_milk_scores_water$SampleID)%>%
                                      column_to_rownames("SampleID") %>%
                                      select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                    PCA_milk_scores_water$water_treatment, ncomp = 2, scale = TRUE)

PLSDA_milk_water_scores <- data.frame(PLSDA_milk_water$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics) %>% 
  mutate(water_treatment = as.factor(water_treatment))

PLSDA_milk_water_plot <- PLSDA_milk_water_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "water_treatment", alpha = 0.6, title = "PLSDA milk water",
            xlab = paste("Component 1 (", round(PLSDA_milk_water$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_milk_water$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = palette_wat) +
  geom_point(data = PLSDA_milk_water_scores %>% group_by(water_treatment) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = water_treatment), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(plot = PLSDA_milk_water_plot, filename = "PLSDA_milk_water_plot.svg", device = "svg", dpi = "retina", width = 5, height = 2.4)

Loadings_milk_water <- plotLoadings(PLSDA_milk_water, plot = FALSE, contrib = "max")$X %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_milk_water <- perf(PLSDA_milk_water, validation = "Mfold", folds = 5, nrepeat = 99, progressBar = TRUE) 
#plot(perf_plsda_milk_water, legend = FALSE)

VIPs_milk_water <- as.data.frame(mixOmics::vip(PLSDA_milk_water))
VIPs_milk_water_treatment <- dplyr::filter(VIPs_milk_water, VIPs_milk_water$comp1 > 1)
VIPs_milk_water_treatment$ID <- rownames(VIPs_milk_water_treatment)
VIPs_milk_water_select <- VIPs_milk_water_treatment %>% dplyr::select(ID, comp1)
VIPs_milk_water_Load <- VIPs_milk_water_select %>% 
  left_join(Loadings_milk_water, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete_filter, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#write_csv(x = VIPs_milk_water_Load, file = "Human_milk_Water_treatment_VIP.csv")

# Check top most relevant features
VIPs_milk_treated <- VIPs_milk_water_Load %>% dplyr::filter(GroupContrib == "treated") %>% head(70)
VIPs_milk_untreated <- VIPs_milk_water_Load %>% dplyr::filter(GroupContrib == "untreated") %>% head(70)

data_milk_ratio_wat <- data_milk %>%
  dplyr::select("SampleID", VIPs_milk_treated$ID, VIPs_milk_untreated$ID) %>%
  dplyr::mutate(treated = rowSums(select(., VIPs_milk_treated$ID))) %>%
  dplyr::mutate(untreated = rowSums(select(., VIPs_milk_untreated$ID))) %>%
  dplyr::mutate(Ratio = log(treated/untreated + 1)) %>%
  left_join(metadata_metabolomics)

# Plot ratio
plot_ratio_milk_water <- data_milk_ratio_wat %>% 
  dplyr::filter(infant_age_days <30) %>% 
  ggboxplot(x = "water_treatment", y = "Ratio", add = "jitter", ylab = "Ln(treated/untreated)",
            add.params = list(color = "water_treatment", alpha = 0.6), legend = "none",
            palette = palette_wat_rev, xlab = "Water treatment",
            title = "Differential features from PLS-DA model") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

plot_ratio_milk_water_time <- data_milk_ratio_wat %>% 
  dplyr::filter(water_treatment != "NA") %>% 
  ggscatter(x = "infant_age_days", y = "Ratio", add = "loess", ylab = "Ln(treated/untreated)", alpha = 0.5,
            add.params = list(color = "water_treatment", alpha = 0.6), color = "water_treatment",
            palette = palette_wat, xlab = "Infant age (days)", legend = "none",
            title = "Milk Ratio - Water treatment") +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Linear mixed effect model
model <- data_milk_ratio_wat %>% 
  lmer(formula = Ratio ~ water_treatment + infant_age_days + (1|host_subject_id))
summary(model)

# Make boxplots for a few of the most discriminant features from PLS-DA

# Acetylcarnitine (C2:0)
milk_carnitine_C2 <- VIPs_milk_water_Load %>%
  dplyr::filter(ID == "968") 

data_milk_carnitine_C2 <- data_milk %>% 
  dplyr::select("SampleID", milk_carnitine_C2$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_carnitine_C2$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  dplyr::filter(infant_age_days < 15) %>% # Colostrum + transitional milk
  group_by(host_subject_id, water_treatment) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() 

data_milk_carnitine_C2$water_treatment <- factor(data_milk_carnitine_C2$water_treatment,
                                                 levels = c("untreated", "treated"))

plot_milk_carnitine_C2 <- ggboxplot(data = data_milk_carnitine_C2, x = "water_treatment", y = "Log_feat", add = "jitter", 
                                  add.params = list(color = "water_treatment", alpha = 0.5), legend = "none",
                                  palette = palette_wat_rev, xlab = "water_treatment", ylab = "Log2(Peak area)",
                                  title = "carnitine") + 
  stat_compare_means(method = "wilcox.test") + ylim(0, 25) +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Propionylcarnitine (C3:0) 
milk_carnitine_C3 <- VIPs_milk_water_Load %>%
  dplyr::filter(ID == "998") 

data_milk_carnitine_C3 <- data_milk %>% 
  dplyr::select("SampleID", milk_carnitine_C3$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_carnitine_C3$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  dplyr::filter(infant_age_days < 15) %>% # Colostrum + transitional milk
  group_by(host_subject_id, water_treatment) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() 

data_milk_carnitine_C3$water_treatment <- factor(data_milk_carnitine_C3$water_treatment,
                                                 levels = c("untreated", "treated"))

plot_milk_carnitine_C3 <- ggboxplot(data = data_milk_carnitine_C3, x = "water_treatment", y = "Log_feat", add = "jitter", 
                                  add.params = list(color = "water_treatment", alpha = 0.5), legend = "none",
                                  palette = palette_wat_rev, xlab = "water_treatment", ylab = "Log2(Peak area)",
                                  title = "carnitine") + 
  stat_compare_means(method = "wilcox.test") + ylim(0, 25) +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

comb_milk_carnitine <- ggarrange(plot_milk_carnitine_C2, plot_milk_carnitine_C3, nrow = 1)

p1 <- compare_means(Log_feat ~ water_treatment, data = data_milk_carnitine_C2)$p
p2 <- compare_means(Log_feat ~ water_treatment, data = data_milk_carnitine_C3)$p
p_raw <- c(p1, p2)
p_adj <- p.adjust(p_raw, method = "BH") 
#ggsave(plot = comb_milk_carnitine, filename = "milk_carntine.svg", device = "svg", dpi = "retina", width = 2.2, height = 2.5)



# PLSDA milk - Lactation stages
# Colostrum (0-4 days), transitional (5-15), mature (>15)
PLSDA_milk_stage <- mixOmics::plsda(data_milk_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                                    PCA_milk_scores$Stage, ncomp = 2, scale = TRUE)

PLSDA_milk_stage_scores <- data.frame(PLSDA_milk_stage$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_milk_stage_plot <- PLSDA_milk_stage_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Stage", alpha = 0.6, title = "PLSDA milk stage",
            xlab = paste("Component 1 (", round(PLSDA_milk_stage$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_milk_stage$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = palette_stage) +
  geom_point(data = PLSDA_milk_stage_scores %>% group_by(Stage) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Stage), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(plot = PLSDA_milk_stage_plot, filename = "PLSDA_milk_stage_plot.svg", device = "svg", dpi = "retina", width = 4, height = 3)

Loadings_milk_stage <- plotLoadings(PLSDA_milk_stage, plot = FALSE, contrib = "max")$X %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_milk_stage <- perf(PLSDA_milk_stage, validation = "Mfold", folds = 5, nrepeat = 99, progressBar = TRUE) 
#plot(perf_plsda_milk_stage, legend = FALSE)

VIPs_milk_stage <- as.data.frame(mixOmics::vip(PLSDA_milk_stage))
VIPs_milk_stage_filter <- dplyr::filter(VIPs_milk_stage, VIPs_milk_stage$comp1 > 1)
VIPs_milk_stage_filter$ID <- rownames(VIPs_milk_stage_filter)
VIPs_milk_stage_select <- VIPs_milk_stage_filter %>% dplyr::select(ID, comp1)
VIPs_milk_stage_Load <- VIPs_milk_stage_select %>% 
  left_join(Loadings_milk_stage, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete_filter, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#write_csv(x = VIPs_milk_stage_Load, file = "Human_milk_Lactation_stage_VIP.csv")


# Plot some of the key discriminant features with boxplots

### Isoleucine ###
milk_isoleucine <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "296") 

data_milk_isoleucine <- data_milk %>% 
  dplyr::select("SampleID", milk_isoleucine$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_isoleucine$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

my_comparisons <- list(c("Colostrum", "Transitional"),
                       c("Colostrum", "Mature"),
                       c("Transitional", "Mature"))

pairwise_stats <- pairwise.wilcox.test(data_milk_isoleucine$Log_feat, data_milk_isoleucine$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_isoleucine <- ggboxplot(data = data_milk_isoleucine, x = "Stage", y = "Log_feat", add = "jitter", 
                      add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                      palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                      title = "(Iso)leucine") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_isoleucine, filename = "plot_milk_isoleucine.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)

### Methionine ###
milk_methionine <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "680") 

data_milk_methionine <- data_milk %>% 
  dplyr::select("SampleID", milk_methionine$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_methionine$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_methionine$Log_feat, data_milk_methionine$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_methionine <- ggboxplot(data = data_milk_methionine, x = "Stage", y = "Log_feat", add = "jitter", 
                             add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                             palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                             title = "Methionine") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_methionine, filename = "plot_milk_methionine.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Tryptophan ###
milk_trp <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "1768") 

data_milk_trp <- data_milk %>% 
  dplyr::select("SampleID", milk_trp$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_trp$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_trp$Log_feat, data_milk_trp$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_trp <- ggboxplot(data = data_milk_trp, x = "Stage", y = "Log_feat", add = "jitter", 
                                  add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                                  palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                                  title = "Tryptophan") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_trp, filename = "plot_milk_trp.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Glutamine ###
milk_glutamine <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "627") 

data_milk_glutamine <- data_milk %>% 
  dplyr::select("SampleID", milk_glutamine$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_glutamine$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_glutamine$Log_feat, data_milk_glutamine$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_glutamine <- ggboxplot(data = data_milk_glutamine, x = "Stage", y = "Log_feat", add = "jitter", 
                           add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                           palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                           title = "Glutamine") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_glutamine, filename = "plot_milk_glutamine.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Hexose m/z 163.0602 (Spectral match to Galactose) ###
milk_hexose <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "572") 

data_milk_hexose <- data_milk %>% 
  dplyr::select("SampleID", milk_hexose$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_hexose$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_hexose$Log_feat, data_milk_hexose$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_hexose <- ggboxplot(data = data_milk_hexose, x = "Stage", y = "Log_feat", add = "jitter", 
                                 add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                                 palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                                 title = "Hexose ()") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_hexose, filename = "plot_milk_hexose.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Disaccharide m/z 325.1131 (Spectral match to Lactose) ###
milk_lactose <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "504") 

data_milk_lactose <- data_milk %>% 
  dplyr::select("SampleID", milk_lactose$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_lactose$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_lactose$Log_feat, data_milk_lactose$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_lactose <- ggboxplot(data = data_milk_lactose, x = "Stage", y = "Log_feat", add = "jitter", 
                              add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                              palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                              title = "Disaccharide") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_lactose, filename = "plot_milk_lactose.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)

### Ile-Pro ###
milk_ile_pro <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "751") 

data_milk_ile_pro <- data_milk %>% 
  dplyr::select("SampleID", milk_ile_pro$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_ile_pro$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_ile_pro$Log_feat, data_milk_ile_pro$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_ile_pro <- ggboxplot(data = data_milk_ile_pro, x = "Stage", y = "Log_feat", add = "jitter", 
                               add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                               palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                               title = "Ile-Pro") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_ile_pro, filename = "plot_milk_ile_pro.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Val-Pro ###
milk_val_pro <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "985") 

data_milk_val_pro <- data_milk %>% 
  dplyr::select("SampleID", milk_val_pro$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_val_pro$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_val_pro$Log_feat, data_milk_val_pro$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_val_pro <- ggboxplot(data = data_milk_val_pro, x = "Stage", y = "Log_feat", add = "jitter", 
                               add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                               palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                               title = "Val-Pro") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_val_pro, filename = "plot_milk_val_pro.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Arginine ###
milk_arg <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "524") 

data_milk_arg <- data_milk %>% 
  dplyr::select("SampleID", milk_arg$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_arg$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_arg$Log_feat, data_milk_arg$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_arg <- ggboxplot(data = data_milk_arg, x = "Stage", y = "Log_feat", add = "jitter", 
                               add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                               palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                               title = "Arginine") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_arg, filename = "plot_milk_arg.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Tyrosine ###
milk_tyr <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "792") 

data_milk_tyr <- data_milk %>% 
  dplyr::select("SampleID", milk_tyr$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_tyr$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_tyr$Log_feat, data_milk_tyr$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_tyr <- ggboxplot(data = data_milk_tyr, x = "Stage", y = "Log_feat", add = "jitter", 
                           add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                           palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                           title = "Tyrosine") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_tyr, filename = "plot_milk_tyr.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Propionyl-carnitine (C3:0) ###
milk_propionyl_carnitine <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "998") 

data_milk_propionyl_carnitine <- data_milk %>% 
  dplyr::select("SampleID", milk_propionyl_carnitine$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_propionyl_carnitine$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_propionyl_carnitine$Log_feat, data_milk_propionyl_carnitine$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_C3_carnitine <- ggboxplot(data = data_milk_propionyl_carnitine, x = "Stage", y = "Log_feat", add = "jitter", 
                           add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                           palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                           title = "Propionylcarnitine") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_C3_carnitine, filename = "plot_milk_C3_carnitine.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Acetyl-carnitine (C2:0) ###
milk_C2_carnitine <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "968") 

data_milk_C2_carnitine <- data_milk %>% 
  dplyr::select("SampleID", milk_C2_carnitine$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_C2_carnitine$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_C2_carnitine$Log_feat, data_milk_C2_carnitine$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_C2_carnitine <- ggboxplot(data = data_milk_C2_carnitine, x = "Stage", y = "Log_feat", add = "jitter", 
                                    add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                                    palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                                    title = "Acetylcarnitine") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_C2_carnitine, filename = "plot_milk_C2_carnitine.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Glutamate ###
milk_glutamate <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "428") 

data_milk_glutamate <- data_milk %>% 
  dplyr::select("SampleID", milk_glutamate$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_glutamate$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_glutamate$Log_feat, data_milk_glutamate$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_milk_glutamate <- ggboxplot(data = data_milk_glutamate, x = "Stage", y = "Log_feat", add = "jitter", 
                                    add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                                    palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                                    title = "Glutamate") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_milk_glutamate, filename = "plot_milk_glutamate.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Butyryl-carnitine (C4:0) ###
milk_C4_carnitine <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "1425") 

data_milk_C4_carnitine <- data_milk %>% 
  dplyr::select("SampleID", milk_C4_carnitine$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_C4_carnitine$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_C4_carnitine$Log_feat, data_milk_C4_carnitine$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_data_milk_C4_carnitine <- ggboxplot(data = data_milk_C4_carnitine, x = "Stage", y = "Log_feat", add = "jitter", 
                                 add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                                 palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                                 title = "Butyrylcarnitine") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_data_milk_C4_carnitine, filename = "plot_data_milk_C4_carnitine.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Tetrasaccharide (C24H42O21) - Maltotetraose ###
milk_maltotetraose <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "496") 

data_milk_maltotetraose <- data_milk %>% 
  dplyr::select("SampleID", milk_maltotetraose$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_maltotetraose$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_maltotetraose$Log_feat, data_milk_maltotetraose$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_data_milk_maltotetraose <- ggboxplot(data = data_milk_maltotetraose, x = "Stage", y = "Log_feat", add = "jitter", 
                                     add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                                     palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                                     title = "Maltotetraose") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_data_milk_maltotetraose, filename = "plot_data_milk_maltotetraose.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)



### Disaccharide (C12H22O11) - Melbiose ###
milk_melbiose <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "514") 

data_milk_melbiose <- data_milk %>% 
  dplyr::select("SampleID", milk_melbiose$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_melbiose$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_melbiose$Log_feat, data_milk_melbiose$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_data_milk_melbiose <- ggboxplot(data = data_milk_melbiose, x = "Stage", y = "Log_feat", add = "jitter", 
                                         add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                                         palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                                         title = "Melbiose") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_data_milk_melbiose, filename = "plot_data_milk_melbiose.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)


### Monosaccharide (C16H12O6) - D-Fructose ###
milk_fructose <- VIPs_milk_stage_Load %>%
  dplyr::filter(ID == "1503") 

data_milk_fructose <- data_milk %>% 
  dplyr::select("SampleID", milk_fructose$ID) %>%
  dplyr::mutate(feat_sum = rowSums(dplyr::select(., milk_fructose$ID))) %>% 
  dplyr::mutate(Log_feat = log2(feat_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID") %>%
  group_by(host_subject_id, Stage) %>%
  slice_max(feat_sum, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Stage = factor(Stage, levels = c("Colostrum", "Transitional", "Mature")))

pairwise_stats <- pairwise.wilcox.test(data_milk_fructose$Log_feat, data_milk_fructose$Stage, p.adjust.method = "BH")
pairwise_stats$p.value

plot_data_milk_fructose <- ggboxplot(data = data_milk_fructose, x = "Stage", y = "Log_feat", add = "jitter", 
                                     add.params = list(color = "Stage", alpha = 0.5), legend = "none",
                                     palette = palette_stage, xlab = FALSE, ylab = "Log2(Peak area)",
                                     title = "D-Fructose") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  p.adjust.method = "BH") +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = plot_data_milk_fructose, filename = "plot_data_milk_fructose.svg", device = "svg", dpi = "retina", width = 1.3, height = 2.1)







#########
# STOOL #
#########
data_stool <- data_sample_filter %>% 
  dplyr::filter(SampleID %in% sample_stool$SampleID) %>% 
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")
  
# RCLR transformation
data_stool_clr <- decostand(data_stool %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_stool <- mixOmics::pca(data_stool_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_stool_scores <- data.frame(PCA_stool$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_stool_plots <- list()

for (i in c("mother_age_years", "gestational_age_birth", "mode_delivery", "water_filter",
            "sex", "hmo_Secretor")) {
  
  PCA_plot <- PCA_stool_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - stool", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_stool$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_stool$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_stool_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_stool_plots[[i]] <- PCA_plot
  
}

PCA_stool_plots_final <- wrap_plots(PCA_stool_plots, nrow = 3) 

i <- "LogAge"

PCA_plot_age <- PCA_stool_scores %>%
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - stool", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_stool$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_stool$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) + scale_color_viridis_c() +
  theme(plot.title = element_text(size = 9), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6), legend.position = "none")
#ggsave(plot = PCA_plot_age, filename = "PCA_plot_age.svg", device = "svg", dpi = "retina", width = 2.5, height = 2.5)

# PERMANOVA
dist_metabolites <- vegdist(data_stool_clr, method = "euclidean")
disper_delivery <- betadisper(dist_metabolites, PCA_stool_scores$infant_age_days)
anova(disper_delivery)
permanova <- adonis2(dist_metabolites ~ infant_age_days + mode_delivery + water_treatment + hmo_Secretor + AssetIndex2 +
                       sex +  birthweight + im_ever_abxs_adhoc + fp_long_bin + maternal_antibiotics + 
                       gestational_age_birth + host_subject_id,
                     PCA_stool_scores, na.action = na.omit, by = "terms")

#permanova_df <- as.data.frame(permanova)
#write.csv(permanova_df, "permanova_output.csv", row.names = TRUE)

i <- "mode_delivery"

PCA_plot_delivery <- PCA_stool_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - stool", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_stool$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_stool$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) + scale_color_manual(values = palette_del) +
  theme(plot.title = element_text(size = 9), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6), legend.position = "none")
#ggsave(plot = PCA_plot_delivery, filename = "PCA_plot_delivery.svg", device = "svg", dpi = "retina", width = 2.5, height = 2.5)

i <- "water_treatment"

PCA_plot_water <- PCA_stool_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - stool", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_stool$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_stool$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) + scale_color_manual(values = palette_wat) +
  theme(plot.title = element_text(size = 9), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = PCA_plot_water, filename = "PCA_plot_water.svg", device = "svg", dpi = "retina", width = 2.5, height = 2.5)

i <- "AssetIndex2"

PCA_plot_assetindex <- PCA_stool_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - stool", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_stool$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_stool$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) + scale_color_manual(values = palette_wat) + 
  theme(plot.title = element_text(size = 9), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = PCA_plot_assetindex, filename = "PCA_plot_assetindex.svg", device = "svg", dpi = "retina", width = 2.5, height = 2.5)

i <- "hmo_Secretor"

PCA_plot_secretor <- PCA_stool_scores %>%
  mutate(hmo_Secretor = as.factor(hmo_Secretor)) %>% 
  dplyr::filter(!is.na(hmo_Secretor)) %>% 
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - stool", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_stool$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_stool$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) + scale_color_manual(values = palette_hmo_rev) + 
  theme(plot.title = element_text(size = 9), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
#ggsave(plot = PCA_plot_secretor, filename = "PCA_plot_secretor.svg", device = "svg", dpi = "retina", width = 2.5, height = 2.5)


# PLS-DA stool - Delivery mode
PLSDA_stool_delivery <- mixOmics::plsda(data_stool_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                                        PCA_stool_scores$mode_delivery, ncomp = 2, scale = TRUE)

PLSDA_stool_delivery_scores <- data.frame(PLSDA_stool_delivery$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_stool_delivery_plot <- PLSDA_stool_delivery_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "mode_delivery", alpha = 0.6, title = "PLSDA stool Delivery",
            xlab = paste("Component 1 (", round(PLSDA_stool_delivery$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_stool_delivery$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(),
            palette = palette_del) +
  geom_point(data = PLSDA_stool_delivery_scores %>% group_by(mode_delivery) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mode_delivery), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(filename = "VPLSDA_stool_delivery_plot.svg", plot = PLSDA_stool_delivery_plot, device = "svg", width = 5, height = 2.8, dpi = "retina")

Loadings_stool_delivery <- plotLoadings(PLSDA_stool_delivery, plot = FALSE, contrib = "max")$X %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_stool_delivery <- perf(PLSDA_stool_delivery, validation = "Mfold", folds = 5, nrepeat = 99, progressBar = TRUE) 
#plot(perf_plsda_stool_delivery, legend = FALSE)
#perf_plsda_stool_delivery$error.rate

VIPs_stool_delivery <- as.data.frame(mixOmics::vip(PLSDA_stool_delivery))
VIPs_stool_delivery_filter <- dplyr::filter(VIPs_stool_delivery, VIPs_stool_delivery$comp1 > 1)
VIPs_stool_delivery_filter$ID <- rownames(VIPs_stool_delivery_filter)
VIPs_stool_delivery_select <- VIPs_stool_delivery_filter %>% dplyr::select(ID, comp1)
VIPs_stool_delivery_Load <- VIPs_stool_delivery_select %>% 
  left_join(Loadings_stool_delivery, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete_filter, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#write_csv(x = VIPs_stool_delivery_Load, file = "Infant_feces_Delivery_mode_VIP.csv")

# Check top most relevant features
VIPs_stool_vaginal <- VIPs_stool_delivery_Load %>% dplyr::filter(GroupContrib == "vaginal") %>% head(100)
VIPs_stool_csection <- VIPs_stool_delivery_Load %>% dplyr::filter(GroupContrib == "csection") %>% head(100)

data_stool_ratio_delivery <- data_stool %>%
  dplyr::select("SampleID", VIPs_stool_vaginal$ID, VIPs_stool_csection$ID) %>%
  dplyr::mutate(va = rowSums(select(., VIPs_stool_vaginal$ID))) %>%
  dplyr::mutate(cs = rowSums(select(., VIPs_stool_csection$ID))) %>%
  dplyr::mutate(Ratio = log(cs/va + 1)) %>%
  left_join(metadata_metabolomics)

# Plot ratio
plot_ratio_stool_delivery <- data_stool_ratio_delivery %>% 
  ggboxplot(x = "mode_delivery", y = "Ratio", add = "jitter", ylab = "Ln(cs/va)",
            add.params = list(color = "mode_delivery", alpha = 0.6),
            palette = palette_del, xlab = "stool",
            title = "Differential features from PLS-DA model") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

plot_ratio_stool_delivery_time <- data_stool_ratio_delivery %>% 
  ggscatter(x = "infant_age_days", y = "Ratio", add = "loess", ylab = "Ln(cs/va)",
            add.params = list(color = "mode_delivery", alpha = 0.6), legend = "none",
            palette = palette_del, xlab = "Infant age (days)", color = "mode_delivery", alpha = 0.5,
            title = "Stool Ratio - Delivery") +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Linear mixed effect model
model <- data_stool_ratio_delivery %>% 
  lmer(formula = Ratio ~ mode_delivery + infant_age_days + (1|host_subject_id))
summary(model)


# Features enriched between 10 and 30 days of life for mode of delivery
filtered_ids_stool <- metadata_metabolomics %>%
  filter(sample_type == "feces", infant_age_days > 10, infant_age_days < 30) %>%
  select(SampleID, infant_age_days, host_subject_id, mode_delivery) %>%
  left_join(data_stool_ratio_delivery %>% select(SampleID, Ratio), by = "SampleID") %>%
  group_by(host_subject_id) %>%
  slice_max(order_by = Ratio, n = 1, with_ties = FALSE) %>% # pick the sample with the highest ratio
  ungroup()

mean(filtered_ids_stool$infant_age_days)
sd(filtered_ids_stool$infant_age_days)

feat_stool_delivery_sum <- data_stool %>%
  semi_join(filtered_ids_stool, by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  dplyr::select(where(~ sum(. != 0, na.rm = TRUE) >= 13)) %>% # Want to only keep features that are present in at least 20% of samples
  #dplyr::select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% 
  rownames_to_column("SampleID") %>%
  left_join(metadata_metabolomics %>% dplyr::select(SampleID, mode_delivery), by = "SampleID") %>% 
  dplyr::select(-SampleID) %>% group_by(mode_delivery) %>%
  summarise(across(everything(), mean))

feat_stool_delivery_fold <- feat_stool_delivery_sum %>%
  pivot_longer(cols = -mode_delivery, names_to = "Feature", values_to = "Value") %>%
  pivot_wider(names_from = mode_delivery, values_from = Value) %>%
  mutate(across(c("csection", "vaginal"), ~ifelse(. == 0, 1e-9, .))) %>%
  mutate(Fold_Change = csection/vaginal) %>%
  mutate(Log2FC = log2(Fold_Change)) %>%
  arrange(desc(Log2FC)) %>% 
  left_join(info_feature_complete_filter) %>% 
  filter(csection > 1e-9 & vaginal > 1e-9) %>% 
  filter(!is.na(Compound_Name)) %>% 
  left_join(VIPs_stool_delivery_Load %>% select(ID, GroupContrib),
            by = c("Feature" = "ID")) 

data_delivery_1030 <- data_stool %>% 
  semi_join(filtered_ids_stool, by = "SampleID") %>%
  column_to_rownames("SampleID") %>% 
  dplyr::select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>%
  dplyr::select(intersect(names(.), feat_stool_delivery_fold$Feature)) %>% 
  dplyr::mutate(across(everything(), ~ log(. + 1))) %>%
  rownames_to_column("SampleID") %>%
  left_join(metadata_metabolomics %>% dplyr::select("SampleID", "mode_delivery"), by = "SampleID") %>%
  dplyr::select(-SampleID) %>%
  rename_with(~ paste0("Feature_", .), -mode_delivery) 

data_data_delivery_1030_long <- data_delivery_1030 %>%
  pivot_longer(cols = -mode_delivery, names_to = "Feature", values_to = "Value") 

# Wilcoxon test for each feature
wilcox_results <- data_data_delivery_1030_long %>%
  group_by(Feature) %>%
  summarise(p_value = wilcox.test(Value[mode_delivery == "csection"], Value[mode_delivery == "vaginal"], exact = FALSE)$p.value,
            .groups = "drop") %>%
  dplyr::mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  dplyr::mutate(Feature = gsub("Feature_", "", Feature)) %>%
  left_join(feat_stool_delivery_fold %>% dplyr::select(Feature, Log2FC, Compound_Name))

wilcox_results_group <- wilcox_results %>%
  left_join(feat_stool_delivery_fold %>% select(Feature, GroupContrib), by = "Feature") %>% 
  dplyr::filter(adj_p_value < 0.07) %>% 
  mutate(Compound_Name = substr(Compound_Name, 1, 60)) %>% 
  mutate(GroupContrib = if_else(Log2FC > 0, "csection", "vaginal")) %>% arrange(Log2FC) %>% 
  dplyr::filter(!Feature %in% c("1838", "13587", "13547", "17126", "13705", "1768", "8586", "9227", "9357"))

mode_delivery_barplot_stool <- ggplot(wilcox_results_group, aes(x = reorder(Compound_Name, Log2FC), y = Log2FC, fill = GroupContrib)) +
  geom_bar(stat = "identity", width = 0.7) + 
  scale_fill_manual(values = c("csection" = "#6A66A3", "vaginal" = "#84A9C0")) +
  theme_minimal() +
  coord_flip() +
  labs(
    title = "Log2FC - Stool Mode Delivery",
    x = "Feature",
    y = "Log2 Fold Change"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), axis.title = element_text(size = 14),
    axis.text = element_text(size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"))

#ggsave(filename = "mode_delivery_barplot.svg", plot = mode_delivery_barplot_stool, device = "svg", width = 8.9, height = 2.8, dpi = "retina")

# Boxplot of carnitines with VIP > 1 delivery mode
# FEAT 968 - ACETYLCARNITINE
data_stool_delivery_C2_carnitine <- data_stool %>% 
  dplyr::select(SampleID, `968`) %>%
  dplyr::mutate(Log_C2_carnitine = log2(`968` + 1)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

data_delivery_C2_carnitine_30 <- data_stool_delivery_C2_carnitine %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_C2_carnitine, mode_delivery) %>%
  dplyr::filter(infant_age_days > 13 & infant_age_days < 29) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_C2_carnitine == max(Log_C2_carnitine, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE)

data_delivery_C2_carnitine_30_plot <- data_delivery_C2_carnitine_30 %>% arrange(mode_delivery) %>%
  ggboxplot(x = "mode_delivery", y = "Log_C2_carnitine", add = "jitter", legend = "none",
            add.params = list(color = "mode_delivery", alpha = 0.5), palette = palette_del,
            xlab = "Days (14-28)", ylab = "log2(peak area)") + stat_compare_means() + ylim (-0.1, 23) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# FEAT 998 - PROPIONYLCARNITINE
data_stool_delivery_C3_carnitine <- data_stool %>% 
  dplyr::select(SampleID, `998`) %>%
  dplyr::mutate(Log_C3_carnitine = log2(`998` + 1)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

data_delivery_C3_carnitine_30 <- data_stool_delivery_C3_carnitine %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_C3_carnitine, mode_delivery) %>%
  dplyr::filter(infant_age_days > 13 & infant_age_days < 29) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_C3_carnitine == max(Log_C3_carnitine, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE)

data_delivery_C3_carnitine_30_plot <- data_delivery_C3_carnitine_30 %>% arrange(mode_delivery) %>%
  ggboxplot(x = "mode_delivery", y = "Log_C3_carnitine", add = "jitter", legend = "none",
            add.params = list(color = "mode_delivery", alpha = 0.5), palette = palette_del,
            xlab = "Days (14-28)", ylab = "log2(peak area)") + stat_compare_means() + ylim (-0.1, 23) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# FEAT 3969 - HEXANOYLCARNITINE
data_stool_delivery_C6_carnitine <- data_stool %>% 
  dplyr::select(SampleID, `3969`) %>%
  dplyr::mutate(Log_C6_carnitine = log2(`3969` + 1)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

data_delivery_C6_carnitine_30 <- data_stool_delivery_C6_carnitine %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_C6_carnitine, mode_delivery) %>%
  dplyr::filter(infant_age_days > 13 & infant_age_days < 29) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_C6_carnitine == max(Log_C6_carnitine, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE)

data_delivery_C6_carnitine_30_plot <- data_delivery_C6_carnitine_30 %>% arrange(mode_delivery) %>%
  ggboxplot(x = "mode_delivery", y = "Log_C6_carnitine", add = "jitter", legend = "none",
            add.params = list(color = "mode_delivery", alpha = 0.5), palette = palette_del,
            xlab = "Days (14-28)", ylab = "log2(peak area)") + stat_compare_means() + ylim (-0.1, 23) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))


ratio_carnitines_delivery <- ggarrange(data_delivery_C2_carnitine_30_plot, data_delivery_C3_carnitine_30_plot, data_delivery_C6_carnitine_30_plot, nrow = 1)

p1 <- compare_means(Log_C2_carnitine ~ mode_delivery, data = data_delivery_C2_carnitine_30)$p
p2 <- compare_means(Log_C3_carnitine ~ mode_delivery, data = data_delivery_C3_carnitine_30)$p
p3 <- compare_means(Log_C6_carnitine ~ mode_delivery, data = data_delivery_C6_carnitine_30)$p
p_raw <- c(p1, p2,p3)
p_adj <- p.adjust(p_raw, method = "BH") 
#ggsave(plot = ratio_time_data_stool_feat1, filename = "ratio_time_data_stool_feat1.svg", device = "svg", dpi = "retina", width = 3, height = 2)


# PLS-DA stool - maternal secretor status
PCA_stool_scores_filter <- PCA_stool_scores  %>% 
  dplyr::filter(hmo_Secretor != "NA")

PLSDA_stool_secretor <- mixOmics::plsda(data_stool_clr %>% rownames_to_column("SampleID") %>%
                                          dplyr::filter(SampleID %in% PCA_stool_scores_filter$SampleID)%>%
                                          column_to_rownames("SampleID") %>%
                                          select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                                        PCA_stool_scores_filter$hmo_Secretor, ncomp = 2, scale = TRUE)

PLSDA_stool_secretor_scores <- data.frame(PLSDA_stool_secretor$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_stool_secretor_plot <- PLSDA_stool_secretor_scores %>%
  mutate(hmo_Secretor = as.factor(hmo_Secretor)) %>%
  ggscatter(x = "comp1", y = "comp2", color = "hmo_Secretor", alpha = 0.6,
            title = "PLSDA stool secretor",
            xlab = paste("Component 1 (", round(PLSDA_stool_secretor$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_stool_secretor$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(),
            palette = palette_hmo_rev) +
  geom_point(data = PLSDA_stool_secretor_scores %>%
               mutate(hmo_Secretor = as.factor(hmo_Secretor)) %>%
               group_by(hmo_Secretor) %>% 
               summarise(across(matches("comp"), mean)),
             aes(comp1, comp2, color = hmo_Secretor), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(filename = "PLSDA_stool_secretor_plot.svg", plot = PLSDA_stool_secretor_plot, device = "svg", width = 5, height = 2.8, dpi = "retina")

Loadings_stool_secretor <- plotLoadings(PLSDA_stool_secretor, plot = FALSE, contrib = "max")$X %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_stool_secretor <- perf(PLSDA_stool_secretor, validation = "Mfold", folds = 5, nrepeat = 99, progressBar = TRUE) 
#plot(perf_plsda_stool_secretor, legend = FALSE)

VIPs_stool_secretor <- as.data.frame(mixOmics::vip(PLSDA_stool_secretor))
VIPs_stool_secretor_filter <- dplyr::filter(VIPs_stool_secretor, VIPs_stool_secretor$comp1 > 1)
VIPs_stool_secretor_filter$ID <- rownames(VIPs_stool_secretor_filter)
VIPs_stool_secretor_select <- VIPs_stool_secretor_filter %>% dplyr::select(ID, comp1)
VIPs_stool_secretor_Load <- VIPs_stool_secretor_select %>% 
  left_join(Loadings_stool_secretor, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete_filter, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#VIPs_stool_secretor_Load2 <- VIPs_stool_secretor_Load
#VIPs_stool_secretor_Load2$GroupContrib <- ifelse(VIPs_stool_secretor_Load2$GroupContrib == 1, "secretor", "non-secretor")
#write_csv(x = VIPs_stool_secretor_Load2, file = "Infant_feces_Secretor_status_VIP.csv")


# Check top most relevant features
VIPs_stool_seretoryes <- VIPs_stool_secretor_Load %>% dplyr::filter(GroupContrib == "1") %>% head(100)
VIPs_stool_seretorno <- VIPs_stool_secretor_Load %>% dplyr::filter(GroupContrib == "0") %>% head(100)

data_stool_ratio_secretor <- data_stool %>%
  dplyr::select("SampleID", VIPs_stool_seretoryes$ID, VIPs_stool_seretorno$ID) %>%
  dplyr::mutate(yes = rowSums(select(., VIPs_stool_seretoryes$ID))) %>%
  dplyr::mutate(no = rowSums(select(., VIPs_stool_seretorno$ID))) %>%
  dplyr::mutate(Ratio = log(yes/no + 1)) %>%
  left_join(metadata_metabolomics)

# Plot ratio
plot_ratio_stool_secretor <- data_stool_ratio_secretor %>% 
  dplyr::filter(hmo_Secretor %in% c(0,1)) %>%
  ggboxplot(x = "hmo_Secretor", y = "Ratio", add = "jitter", ylab = "Ln(cs/va)",
            add.params = list(color = "hmo_Secretor", alpha = 0.6), legend = "none",
            xlab = "stool", palette = palette_hmo_rev,
            title = "Differential features from PLS-DA model") +
  stat_compare_means() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

plot_ratio_stool_secretor_time <- data_stool_ratio_secretor %>% 
  dplyr::filter(hmo_Secretor %in% c(0,1)) %>%
  dplyr::mutate(hmo_Secretor = factor(hmo_Secretor, levels = c(1, 0))) %>%
  dplyr::filter(Ratio < 1) %>%
  ggscatter(x = "infant_age_days", y = "Ratio", add = "loess", ylab = "Ln(yes/no)",
            add.params = list(color = "hmo_Secretor", alpha = 0.6), legend = "none",
            palette = palette_hmo, xlab = "Infant age (days)", color = "hmo_Secretor", alpha = 0.5,
            title = "Stool Ratio - Secretor") +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Linear mixed effect model
model <- data_stool_ratio_secretor %>% 
  dplyr::filter(hmo_Secretor != "unknown") %>%
  lmer(formula = Ratio ~ infant_age_days + hmo_Secretor + (1|host_subject_id))
summary(model)


# Bile acids stool - maternal secretor status

# 12 oxo BA
data_stool_bas_oxo_10001 <- data_stool %>% 
  dplyr::select("SampleID", "10001") %>%
  dplyr::mutate(Log_oxo = log2(`10001` + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID")

data_ba_10 <- data_stool_bas_oxo_10001 %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_oxo, hmo_Secretor) %>%
  dplyr::filter(infant_age_days < 8) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_oxo == max(Log_oxo, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_ba_10_plot <- data_ba_10 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_oxo", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev, 
            xlab = "Days (0-7)", ylab = "Log2(peak area)") + stat_compare_means() + ylim (-0.1, 20) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_ba_30 <- data_stool_bas_oxo_10001 %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_oxo, hmo_Secretor) %>%
  dplyr::filter(infant_age_days > 13 & infant_age_days < 29) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_oxo == max(Log_oxo, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_ba_30_plot <- data_ba_30 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_oxo", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Days (14-28)", ylab = FALSE) + stat_compare_means() + ylim (-0.1, 20) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_ba_60 <- data_stool_bas_oxo_10001 %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_oxo, hmo_Secretor) %>%
  dplyr::filter(infant_age_days > 29 & infant_age_days < 46) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_oxo == max(Log_oxo, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_ba_60_plot <- data_ba_60 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_oxo", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Days (30-45)", ylab = FALSE) + stat_compare_means() +  ylim (-0.1, 20) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

ratio_time_oxoba_comb <- ggarrange(data_ba_10_plot, data_ba_30_plot, data_ba_60_plot, nrow = 1)

p1 <- compare_means(Log_oxo ~ hmo_Secretor, data = data_ba_10)$p
p2 <- compare_means(Log_oxo ~ hmo_Secretor, data = data_ba_30)$p
p3 <- compare_means(Log_oxo ~ hmo_Secretor, data = data_ba_60)$p
p_raw <- c(p1, p2, p3)
p_adj <- p.adjust(p_raw, method = "BH") 
#ggsave(plot = ratio_time_oxoba_comb, filename = "Ratio_box_oxoBA_HMO.svg", device = "svg", dpi = "retina", width = 3, height = 2.1)


# Trihydroxylated BA HMO secretor status
data_stool_bas_10329 <- data_stool %>% 
  dplyr::select("SampleID", "10329") %>%
  dplyr::mutate(Log_unknown_tri = log2(`10329` + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID")

data_ba_10 <- data_stool_bas_10329 %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_unknown_tri, hmo_Secretor) %>%
  dplyr::filter(infant_age_days < 8) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_unknown_tri == max(Log_unknown_tri, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_ba_10_plot <- data_ba_10 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_unknown_tri", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev, 
            xlab = "Days (0-7)", ylab = "Log2(peak area)") + stat_compare_means() + ylim (-0.1, 20) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_ba_30 <- data_stool_bas_10329 %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_unknown_tri, hmo_Secretor) %>%
  dplyr::filter(infant_age_days > 13 & infant_age_days < 29) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_unknown_tri == max(Log_unknown_tri, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_ba_30_plot <- data_ba_30 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_unknown_tri", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Days (14-28)", ylab = FALSE) + stat_compare_means() + ylim (-0.1, 20) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_ba_60 <- data_stool_bas_10329 %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_unknown_tri, hmo_Secretor) %>%
  dplyr::filter(infant_age_days > 29 & infant_age_days < 46) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_unknown_tri == max(Log_unknown_tri, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_ba_60_plot <- data_ba_60 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_unknown_tri", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Days (30-45)", ylab = FALSE) + stat_compare_means() +  ylim (-0.1, 20) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_ba_100 <- data_stool_bas_10329 %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_unknown_tri, hmo_Secretor) %>%
  dplyr::filter(infant_age_days > 179 & infant_age_days < 206) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_unknown_tri == max(Log_unknown_tri, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_ba_100_plot <- data_ba_100 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_unknown_tri", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Days (180-205)", ylab = FALSE) + stat_compare_means() + ylim (-0.1, 20) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

ratio_time_unknown_tri_comb <- ggarrange(data_ba_10_plot, data_ba_30_plot, data_ba_60_plot, data_ba_100_plot, nrow = 1)

p1 <- compare_means(Log_unknown_tri ~ hmo_Secretor, data = data_ba_10)$p
p2 <- compare_means(Log_unknown_tri ~ hmo_Secretor, data = data_ba_30)$p
p3 <- compare_means(Log_unknown_tri ~ hmo_Secretor, data = data_ba_60)$p
p4 <- compare_means(Log_unknown_tri ~ hmo_Secretor, data = data_ba_100)$p
p_raw <- c(p1, p2, p3, p4)
p_adj <- p.adjust(p_raw, method = "BH") 
#ggsave(plot = ratio_time_unknown_tri_comb, filename = "Ratio_box_unknown tri_HMO.svg", device = "svg", dpi = "retina", width = 4, height = 2.1)

# Boxplot 2FL stool

data_stool_2FL <- data_stool %>%  
  dplyr::select(SampleID, `583`) %>%
  dplyr::mutate(Log_2FL = log2(`583` + 1)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

data_stool_2FL_10 <- data_stool_2FL %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days < 8) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_2FL == max(Log_2FL, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_stool_2FL_10_plot <- data_stool_2FL_10 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_2FL", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev, 
            xlab = "Days (0-7)", ylab = "Log2(peak area)") + stat_compare_means() + ylim (-0.1, 22) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_stool_2FL_30 <- data_stool_2FL %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days > 13 & infant_age_days < 29) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_2FL == max(Log_2FL, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_stool_2FL_30_plot <- data_stool_2FL_30 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_2FL", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Days (14-28)", ylab = FALSE) + stat_compare_means() + ylim (-0.1, 22) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_stool_2FL_60 <- data_stool_2FL %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_2FL, hmo_Secretor) %>%
  dplyr::filter(infant_age_days > 29 & infant_age_days < 46) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_2FL == max(Log_2FL, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_stool_2FL_60_plot <- data_stool_2FL_60 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_2FL", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Days (30-45)", ylab = FALSE) + stat_compare_means() +  ylim (-0.1, 22) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

ratio_time_2FL_comb <- ggarrange(data_stool_2FL_10_plot, data_stool_2FL_30_plot, data_stool_2FL_60_plot, nrow = 1)

p1 <- compare_means(Log_2FL ~ hmo_Secretor, data = data_stool_2FL_10)$p
p2 <- compare_means(Log_2FL ~ hmo_Secretor, data = data_stool_2FL_30)$p
p3 <- compare_means(Log_2FL ~ hmo_Secretor, data = data_stool_2FL_60)$p
p_raw <- c(p1, p2, p3)
p_adj <- p.adjust(p_raw, method = "BH")
#ggsave(plot = ratio_time_2FL_comb, filename = "Ratio_time_2FL_comb.svg", device = "svg", dpi = "retina", width = 4, height = 2.1)


# Boxplot DFLac stool
data_stool_dflac <- data_stool %>% 
  dplyr::select(SampleID, `769`) %>%
  dplyr::mutate(Log_dflac = log2(`769` + 1)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

data_stool_dflac_10 <- data_stool_dflac %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_dflac, hmo_Secretor) %>%
  dplyr::filter(infant_age_days < 8) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_dflac == max(Log_dflac, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_stool_dflac_10_plot <- data_stool_dflac_10 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_dflac", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev, 
            xlab = "Days (0-7)", ylab = "Log2(peak area)") + stat_compare_means() + ylim (-0.1, 22) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_stool_dflac_30 <- data_stool_dflac %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_dflac, hmo_Secretor) %>%
  dplyr::filter(infant_age_days > 13 & infant_age_days < 29) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_dflac == max(Log_dflac, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_stool_dflac_30_plot <- data_stool_dflac_30 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_dflac", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Days (14-28)", ylab = FALSE) + stat_compare_means() + ylim (-0.1, 22) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_stool_dflac_60 <- data_stool_dflac %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_dflac, hmo_Secretor) %>%
  dplyr::filter(infant_age_days > 29 & infant_age_days < 46) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_dflac == max(Log_dflac, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_stool_dflac_60_plot <- data_stool_dflac_60 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_dflac", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Days (30-45)", ylab = FALSE) + stat_compare_means() +  ylim (-0.1, 22) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_stool_dflac_100 <- data_stool_dflac %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_dflac, hmo_Secretor) %>%
  dplyr::filter(infant_age_days > 179 & infant_age_days < 206) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_dflac == max(Log_dflac, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE) %>% 
  dplyr::filter(!hmo_Secretor == "NA")

data_stool_dflac_100_plot <- data_stool_dflac_100 %>% arrange(hmo_Secretor) %>%
  ggboxplot(x = "hmo_Secretor", y = "Log_dflac", add = "jitter", legend = "none",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev,
            xlab = "Days (180-205)", ylab = FALSE) + stat_compare_means() +  ylim (-0.1, 22) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

ratio_time_dflac_comb <- ggarrange(data_stool_dflac_10_plot, data_stool_dflac_30_plot, data_stool_dflac_60_plot, data_stool_dflac_100_plot, nrow = 1)

p1 <- compare_means(Log_dflac ~ hmo_Secretor, data = data_stool_dflac_10)$p
p2 <- compare_means(Log_dflac ~ hmo_Secretor, data = data_stool_dflac_30)$p
p3 <- compare_means(Log_dflac ~ hmo_Secretor, data = data_stool_dflac_60)$p
p4 <- compare_means(Log_dflac ~ hmo_Secretor, data = data_stool_dflac_100)$p
p_raw <- c(p1, p2, p3, p4)
p_adj <- p.adjust(p_raw, method = "BH")
#ggsave(plot = ratio_time_dflac_comb, filename = "Ratio_time_DFLac_comb.svg", device = "svg", dpi = "retina", width = 4, height = 2.5)


##################
# PIECHART DISCRIMINANT FEATURES BY MATERNAL SECRETOR STATUS - STOOL
##################
canopus <- read.delim("data/canopus_formula_summary.tsv")
canopus$mappingFeatureId <- as.character(canopus$mappingFeatureId)

VIPs_stool_secretor_Load_canopus <- VIPs_stool_secretor_Load %>% 
  left_join(canopus, by = c("ID" = "mappingFeatureId")) %>% 
  select(ID, comp1, GroupContrib, mz, RT, Compound_Name, Adduct, molecularFormula, adduct, NPC.pathway, NPC.pathway.Probability,
         NPC.superclass, NPC.superclass.Probability, NPC.class, NPC.class.Probability) %>% 
  dplyr::mutate(NPC.combined = dplyr::case_when(
    NPC.class.Probability >= 0.7 ~ NPC.class,
    NPC.class.Probability < 0.7 & NPC.superclass.Probability >= 0.7 ~ NPC.superclass,
    NPC.class.Probability < 0.7 & NPC.superclass.Probability < 0.7 & NPC.pathway.Probability >= 0.7 ~ NPC.pathway,
    TRUE ~ NA_character_)) %>% 
  dplyr::filter(comp1 >2.5) %>% 
  mutate(NPC.superclass = ifelse(is.na(NPC.superclass), "NA", NPC.superclass))

class_order <- c("Aminosugars and aminoglycosides",
                 "Saccharides",
                 "Fatty acyl glycosides",
                 "Oligopeptides",
                 "Small peptides",
                 "Fatty amides",
                 "Fatty Acids and Conjugates",
                 "Fatty acyls",
                 "Fatty esters",
                 "Glycerolipids",
                 "Steroids",
                 "Triterpenoids",
                 "Spingolipids",
                 "Eicosanoids",
                 "Histidine alkaloids",
                 "Proline alkaloids",
                 "Ornithine alkaloids",
                 "Tyrosine alkaloids",
                 "Pseudoalkaloids (transamidation)",
                 "Macrolides",
                 "Nucleosides",
                 "NA")

class_colors <- c(
  # Carbohydrate-related (pink-purple)
  "Aminosugars and aminoglycosides" = "#fda9c1",
  "Saccharides" = "#fe90ad",
  "Fatty acyl glycosides" = "#ff82ab",
  
  # Peptides (blue-teal)
  "Oligopeptides" = "#fe558c",
  "Small peptides" = "#ff6b9b",
  
  # Lipids (green gradient)
  "Fatty amides" = "#d3ecbc",
  "Fatty Acids and Conjugates" = "#b1d697",
  "Fatty acyls" = "#97ce7d",
  "Fatty esters" = "#8db670",
  "Glycerolipids" = "#6db474",
  
  # Steroids and terpenoids (olive to brown)
  "Steroids" = "#6ab26c",
  "Triterpenoids" = "#467750",
  
  # Complex lipids
  "Spingolipids" = "#2a5e33",
  
  # Bioactive lipid class
  "Eicosanoids" = "#034007",
  
  # Alkaloids (different hues of blue)
  "Histidine alkaloids" = "#668397",
  "Proline alkaloids" = "#4c6e86",
  "Ornithine alkaloids" = "#325a75",
  "Tyrosine alkaloids" = "#194564",
  "Pseudoalkaloids (transamidation)" = "#003153",
  
  # Others
  "Macrolides" = "#ccdaee",
  "Nucleosides" = "#ccdaee",
  "NA" = "#F5F5F5")


pie_data <- VIPs_stool_secretor_Load_canopus %>%
  count(NPC.superclass) %>%
  mutate(NPC.superclass = factor(NPC.superclass, levels = class_order),
    percent = n / sum(n) * 100,
    label = paste0(NPC.superclass, " (", round(percent, 1), "%)"))

piechart_secretor <- ggplot(pie_data, aes(x = "", y = percent, fill = NPC.superclass)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = class_colors) +
  labs(title = "NPC Superclass Distribution", fill = "Superclass") +
  theme_void() +
  theme(legend.position = "right")

#ggsave(plot = piechart_secretor, filename = "piechart_secretor.svg", device = "svg", dpi = "retina", width = 9, height = 4)


# PLS-DA stool - AssetIndex
PLSDA_stool_assetindex <- mixOmics::plsda(data_stool_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                                        PCA_stool_scores$AssetIndex2, ncomp = 2, scale = TRUE)

PLSDA_stool_assetindex_scores <- data.frame(PLSDA_stool_assetindex$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_stool_assetindex_plot <- PLSDA_stool_assetindex_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "AssetIndex2", alpha = 0.6, title = "PLSDA stool AssetIndex",
            xlab = paste("Component 1 (", round(PLSDA_stool_assetindex$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_stool_assetindex$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(),
            palette = palette_wat) +
  geom_point(data = PLSDA_stool_assetindex_scores %>% group_by(AssetIndex2) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = AssetIndex2), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(filename = "PLSDA_stool_assetindex_plot.svg", plot = PLSDA_stool_assetindex_plot, device = "svg", width = 5, height = 2.8, dpi = "retina")

Loadings_stool_assetindex <- plotLoadings(PLSDA_stool_assetindex, plot = FALSE, contrib = "max")$X %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_stool_assetindex <- perf(PLSDA_stool_assetindex, validation = "Mfold", folds = 5, nrepeat = 99, progressBar = TRUE) 
#plot(perf_plsda_stool_assetindex, legend = FALSE)

VIPs_stool_assetindex <- as.data.frame(mixOmics::vip(PLSDA_stool_assetindex))
VIPs_stool_assetindex_filter <- dplyr::filter(VIPs_stool_assetindex, VIPs_stool_assetindex$comp1 > 1)
VIPs_stool_assetindex_filter$ID <- rownames(VIPs_stool_assetindex_filter)
VIPs_stool_assetindex_select <- VIPs_stool_assetindex_filter %>% dplyr::select(ID, comp1)
VIPs_stool_assetindex_Load <- VIPs_stool_assetindex_select %>% 
  left_join(Loadings_stool_assetindex, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete_filter, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#write_csv(x = VIPs_stool_assetindex_Load, file = "Infant_feces_AssetIndex_VIP.csv")

VIPs_stool_Q12 <- VIPs_stool_assetindex_Load %>% dplyr::filter(GroupContrib == "Q12") %>% head(100)
VIPs_stool_Q35 <- VIPs_stool_assetindex_Load %>% dplyr::filter(GroupContrib == "Q35") %>% head(100)

data_stool_ratio_assetindex <- data_stool %>%
  dplyr::select("SampleID", VIPs_stool_Q12$ID, VIPs_stool_Q35$ID) %>%
  dplyr::mutate(Q12 = rowSums(select(., VIPs_stool_Q12$ID))) %>%
  dplyr::mutate(Q35 = rowSums(select(., VIPs_stool_Q35$ID))) %>%
  dplyr::mutate(Ratio = log(Q35/Q12 + 1)) %>%
  left_join(metadata_metabolomics)

# Plot ratio
plot_ratio_stool_assetindex_time <- data_stool_ratio_assetindex %>% 
  ggscatter(x = "infant_age_days", y = "Ratio", add = "loess", ylab = "Ln(Q35/Q12)",
            add.params = list(color = "AssetIndex2", alpha = 0.6),
            palette = palette_wat, xlab = "Infant age (days)", color = "AssetIndex2", alpha = 0.5,
            title = "Stool Ratio - AssetIndex2") +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Linear mixed effect model
model <- data_stool_ratio_assetindex %>% 
  lmer(formula = Ratio ~ AssetIndex2 + infant_age_days + (1|host_subject_id))
summary(model)


# Look closer at the drink_water_safe_variable (water treatment)
# Count and calculate percentages to plot
drink_dist <- baseline %>%
  count(drink_water_safe) %>%
  mutate(percent = 100 * n / sum(n))

ggplot(drink_dist, aes(x = reorder(drink_water_safe, -percent), y = percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(round(percent, 1), "%")), vjust = -0.3, size = 3) +
  labs(x = "Water Treatment Method", y = "Percentage", title = "Water Safety Practices at Baseline") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# PLS-DA stool - water treatment (any treatment vs. no treatment)
PLSDA_stool_water <- mixOmics::plsda(data_stool_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                                          PCA_stool_scores$water_treatment, ncomp = 2, scale = TRUE)

PLSDA_stool_water_scores <- data.frame(PLSDA_stool_water$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_stool_water_plot <- PLSDA_stool_water_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "water_treatment", alpha = 0.6, title = "PLSDA stool water",
            xlab = paste("Component 1 (", round(PLSDA_stool_water$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_stool_water$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(),
            palette = palette_wat) +
  geom_point(data = PLSDA_stool_water_scores %>% group_by(water_treatment) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = water_treatment), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(filename = "PLSDA_stool_water_plot.svg", plot = PLSDA_stool_water_plot, device = "svg", width = 5, height = 2.8, dpi = "retina")

Loadings_stool_water <- plotLoadings(PLSDA_stool_water, plot = FALSE, contrib = "max")$X %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_stool_water <- perf(PLSDA_stool_water, validation = "Mfold", folds = 5, nrepeat = 99, progressBar = TRUE) 
#plot(perf_plsda_stool_water, legend = FALSE)

VIPs_stool_water <- as.data.frame(mixOmics::vip(PLSDA_stool_water))
VIPs_stool_water_treatment <- dplyr::filter(VIPs_stool_water, VIPs_stool_water$comp1 > 1)
VIPs_stool_water_treatment$ID <- rownames(VIPs_stool_water_treatment)
VIPs_stool_water_select <- VIPs_stool_water_treatment %>% dplyr::select(ID, comp1)
VIPs_stool_water_Load <- VIPs_stool_water_select %>% 
  left_join(Loadings_stool_water, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete_filter, by = c("ID" = "Feature")) %>% arrange(desc(comp1)) 

#write_csv(x = VIPs_stool_water_Load, file = "Infant_feces_Water_treatment_VIP.csv")

VIPs_stool_treated <- VIPs_stool_water_Load %>% dplyr::filter(GroupContrib == "treated") %>% head(100)
VIPs_stool_untreated <- VIPs_stool_water_Load %>% dplyr::filter(GroupContrib == "untreated") %>% head(100)

data_stool_ratio_water <- data_stool %>%
  dplyr::select("SampleID", VIPs_stool_treated$ID, VIPs_stool_untreated$ID) %>%
  dplyr::mutate(treated = rowSums(select(., VIPs_stool_treated$ID))) %>%
  dplyr::mutate(untreated = rowSums(select(., VIPs_stool_untreated$ID))) %>%
  dplyr::mutate(Ratio = log(treated/untreated + 1)) %>%
  left_join(metadata_metabolomics)

# Plot ratio
plot_ratio_stool_water_time <- data_stool_ratio_water %>% 
  ggscatter(x = "infant_age_days", y = "Ratio", add = "loess", ylab = "Ln(treated/untreated)",
            add.params = list(color = "water_treatment", alpha = 0.6), 
            palette = palette_wat, xlab = "Infant age (days)", color = "water_treatment", alpha = 0.5,
            title = "Stool Ratio - water_treatment") +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Linear mixed effect model
model <- data_stool_ratio_water %>% 
  lmer(formula = Ratio ~ water_treatment + infant_age_days + (1|host_subject_id))
summary(model)



# PLS-DA stool - boil water vs. not
PLSDA_stool_boil <- mixOmics::plsda(data_stool_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                                     PCA_stool_scores$water_boil, ncomp = 2, scale = TRUE)

PLSDA_stool_boil_scores <- data.frame(PLSDA_stool_boil$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_stool_boil_plot <- PLSDA_stool_boil_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "water_boil", alpha = 0.6, title = "PLSDA stool boil",
            xlab = paste("Component 1 (", round(PLSDA_stool_boil$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_stool_boil$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(),
            palette = palette_wat_rev) +
  geom_point(data = PLSDA_stool_boil_scores %>% group_by(water_boil) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = water_boil), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(filename = "PLSDA_stool_boil_plot.svg", plot = PLSDA_stool_boil_plot, device = "svg", width = 5, height = 2.8, dpi = "retina")

Loadings_stool_boil <- plotLoadings(PLSDA_stool_boil, plot = FALSE, contrib = "max")$X %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_stool_boil <- perf(PLSDA_stool_boil, validation = "Mfold", folds = 5, nrepeat = 99, progressBar = TRUE) 
#plot(perf_plsda_stool_boil, legend = FALSE)

VIPs_stool_boil <- as.data.frame(mixOmics::vip(PLSDA_stool_boil))
VIPs_stool_boil_filter <- dplyr::filter(VIPs_stool_boil, VIPs_stool_boil$comp1 > 1)
VIPs_stool_boil_filter$ID <- rownames(VIPs_stool_boil_filter)
VIPs_stool_boil_select <- VIPs_stool_boil_filter %>% dplyr::select(ID, comp1)
VIPs_stool_boil_Load <- VIPs_stool_boil_select %>% 
  left_join(Loadings_stool_boil, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete_filter, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#write_csv(x = VIPs_stool_boil_Load, file = "Feces_boil_purified_VIP.csv")

VIPs_stool_boil <- VIPs_stool_boil_Load %>% dplyr::filter(GroupContrib == "Yes_boil") %>% head(100)
VIPs_stool_no_boil <- VIPs_stool_boil_Load %>% dplyr::filter(GroupContrib == "No") %>% head(100)

data_stool_ratio_boil <- data_stool %>%
  dplyr::select("SampleID", VIPs_stool_boil$ID, VIPs_stool_no_boil$ID) %>%
  dplyr::mutate(Yes_boil= rowSums(select(., VIPs_stool_boil$ID))) %>%
  dplyr::mutate(No = rowSums(select(., VIPs_stool_no_boil$ID))) %>%
  dplyr::mutate(Ratio = log(No/Yes_boil + 1)) %>%
  left_join(metadata_metabolomics)

# Plot ratio
plot_ratio_stool_boil_time <- data_stool_ratio_boil %>% 
  ggscatter(x = "infant_age_days", y = "Ratio", add = "loess", ylab = "Ln(Yes_boil/No)",
            add.params = list(color = "water_boil", alpha = 0.6), 
            palette = palette_wat_rev, xlab = "Age (days)", color = "water_boil", alpha = 0.5,
            title = "Stool Ratio - water_boil") +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Linear mixed effect model
model <- data_stool_ratio_boil %>% 
  lmer(formula = Ratio ~ water_boil + infant_age_days + (1|host_subject_id))
summary(model)


# PLSDA stool - filtered water vs. not
PLSDA_stool_filter <- mixOmics::plsda(data_stool_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                                    PCA_stool_scores$water_filter, ncomp = 2, scale = TRUE)
PLSDA_stool_filter_scores <- data.frame(PLSDA_stool_filter$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_stool_filter_plot <- PLSDA_stool_filter_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "water_filter", alpha = 0.6, title = "PLSDA stool filter",
            xlab = paste("Component 1 (", round(PLSDA_stool_filter$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_stool_filter$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(),
            palette = palette_wat_rev) +
  geom_point(data = PLSDA_stool_filter_scores %>% group_by(water_filter) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = water_filter), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(filename = "PLSDA_stool_filter_plot.svg", plot = PLSDA_stool_filter_plot, device = "svg", width = 5, height = 2.8, dpi = "retina")


Loadings_stool_filter <- plotLoadings(PLSDA_stool_filter, plot = FALSE, contrib = "max")$X %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_stool_filter <- perf(PLSDA_stool_filter, validation = "Mfold", folds = 5, nrepeat = 99, progressBar = TRUE) 
#plot(perf_plsda_stool_filter, legend = FALSE)


VIPs_stool_filter <- as.data.frame(mixOmics::vip(PLSDA_stool_filter))
VIPs_stool_filter_filter <- dplyr::filter(VIPs_stool_filter, VIPs_stool_filter$comp1 > 1)
VIPs_stool_filter_filter$ID <- rownames(VIPs_stool_filter_filter)
VIPs_stool_filter_select <- VIPs_stool_filter_filter %>% dplyr::select(ID, comp1)
VIPs_stool_filter_Load <- VIPs_stool_filter_select %>% 
  left_join(Loadings_stool_filter, by = c("ID" = "rowname")) %>%
  left_join(info_feature_complete_filter, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#write_csv(x = VIPs_stool_filter_Load, file = "Feces_filter_purified_VIP.csv")

VIPs_stool_filter <- VIPs_stool_filter_Load %>% dplyr::filter(GroupContrib == "Yes_filter") %>% head(100)
VIPs_stool_no_filter <- VIPs_stool_filter_Load %>% dplyr::filter(GroupContrib == "No") %>% head(100)

data_stool_ratio_filter <- data_stool %>%
  dplyr::select("SampleID", VIPs_stool_filter$ID, VIPs_stool_no_filter$ID) %>%
  dplyr::mutate(Yes_filter= rowSums(select(., VIPs_stool_filter$ID))) %>%
  dplyr::mutate(No = rowSums(select(., VIPs_stool_no_filter$ID))) %>%
  dplyr::mutate(Ratio = log(No/Yes_filter + 1)) %>%
  left_join(metadata_metabolomics)

# Plot ratio
plot_ratio_stool_filter_time <- data_stool_ratio_filter %>% 
  ggscatter(x = "infant_age_days", y = "Ratio", add = "loess", ylab = "Ln(Yes_filter/No)",
            add.params = list(color = "water_filter", alpha = 0.6), 
            palette = palette_wat_rev, xlab = "Age (days)", color = "water_filter", alpha = 0.5,
            title = "Stool Ratio - water_filter") +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Linear mixed effect model
model <- data_stool_ratio_filter %>% 
  lmer(formula = Ratio ~ water_filter + infant_age_days + (1|host_subject_id))
summary(model)


# Consistent signature
combined_vip <- VIPs_stool_boil_Load %>% 
  dplyr::select(ID, GroupContrib, comp1) %>%
  inner_join(VIPs_stool_filter_Load %>% dplyr::select(ID, GroupContrib, comp1), by = "ID") %>%
  inner_join(VIPs_stool_water_Load %>% dplyr::select(ID, GroupContrib, comp1), by = "ID") %>%
  dplyr::mutate(concordance = case_when(GroupContrib.x == "No" & GroupContrib.y == "No" & GroupContrib == "untreated" ~ "Yes",
                                        GroupContrib.x == "Yes_boil" & GroupContrib.y == "Yes_filter" & GroupContrib == "treated" ~ "Yes",
                                        TRUE ~ "No")) %>%
  dplyr::filter(concordance == "Yes") %>% left_join(info_feature_complete_filter, by = c("ID" = "Feature")) %>% 
  mutate(mean_comp1 = rowMeans(across(c(comp1, comp1.x, comp1.y)), na.rm = TRUE)) %>%
  arrange(desc(mean_comp1))

#write_csv(x = combined_vip, file = "Infant_feces_combined_vip_water.csv")
### Export in cytoscape to create molecular network

# Combined plots
comb_stool_ratio <- ggarrange(plot_ratio_stool_delivery_time, plot_ratio_stool_secretor_time, plot_ratio_stool_water_time, nrow = 1)
#ggsave(plot = comb_stool_ratio, filename = "Stool_ratios.svg", device = "svg", dpi = "retina", height = 2.5, width = 7)

# Boxplot of carnitines with VIP > 1 water treatment

# FEAT 5927 - ACETYLCARNITINE
data_stool_acetylcarnitine <- data_stool %>% 
  dplyr::select(SampleID, `968`) %>%
  dplyr::mutate(Log_acetylcarnitine = log2(`968` + 1)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

data_acetylcarnitine_10 <- data_stool_acetylcarnitine %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_acetylcarnitine, water_treatment) %>%
  dplyr::filter(infant_age_days < 8) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_acetylcarnitine == max(Log_acetylcarnitine, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE)

data_acetylcarnitine_10_plot <- data_acetylcarnitine_10 %>% arrange(water_treatment) %>%
  mutate(water_treatment = factor(water_treatment, levels = c("untreated", "treated"))) %>%
  ggboxplot(x = "water_treatment", y = "Log_acetylcarnitine", add = "jitter", legend = "none",
            add.params = list(color = "water_treatment", alpha = 0.5), palette = palette_wat_rev, 
            xlab = "Days (0-7)", ylab = "Log2(peak area)") + stat_compare_means() + ylim (-0.1, 25) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# FEAT 1425 BUTYRYLCARNITINE 
data_stool_butyrylcarnitine <- data_stool %>% 
  dplyr::select(SampleID, `1425`) %>% 
  dplyr::mutate(Log_butyrylcarnitine = log2(`1425` + 1)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

butyrylcarnitine_10 <- data_stool_butyrylcarnitine %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_butyrylcarnitine, water_treatment) %>%
  dplyr::filter(infant_age_days < 8) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_butyrylcarnitine == max(Log_butyrylcarnitine, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE)

butyrylcarnitine_10_plot <- butyrylcarnitine_10 %>% arrange(water_treatment) %>%
  mutate(water_treatment = factor(water_treatment, levels = c("untreated", "treated"))) %>%
  ggboxplot(x = "water_treatment", y = "Log_butyrylcarnitine", add = "jitter", legend = "none",
            add.params = list(color = "water_treatment", alpha = 0.5), palette = palette_wat_rev, 
            xlab = "Days (0-7)", ylab = "Log2(peak area)") + stat_compare_means() + ylim (-0.1, 25) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# FEAT 2792 - VALERYLCARNITINE 
data_stool_valerylcarnitine <- data_stool %>% 
  dplyr::select(SampleID, `2792`) %>%
  dplyr::mutate(Log_valerylcarnitine = log2(`2792` + 1)) %>%
  left_join(metadata_metabolomics, by = "SampleID")

data_stool_valerylcarnitine_10 <- data_stool_valerylcarnitine %>% 
  dplyr::select(host_subject_id, infant_age_days, Log_valerylcarnitine, water_treatment) %>%
  dplyr::filter(infant_age_days < 8) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(Log_valerylcarnitine == max(Log_valerylcarnitine, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE)

data_stool_valerylcarnitine_10_plot <- data_stool_valerylcarnitine_10 %>% arrange(water_treatment) %>%
  mutate(water_treatment = factor(water_treatment, levels = c("untreated", "treated"))) %>%
  ggboxplot(x = "water_treatment", y = "Log_valerylcarnitine", add = "jitter", legend = "none",
            add.params = list(color = "water_treatment", alpha = 0.5), palette = palette_wat_rev, 
            xlab = "Days (0-7)", ylab = "Log2(peak area)") + stat_compare_means() + ylim (-0.1, 25) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

ratio_time_data_stool_carnitine <- ggarrange(data_acetylcarnitine_10_plot, butyrylcarnitine_10_plot, data_stool_valerylcarnitine_10_plot, nrow = 1)

p1 <- compare_means(Log_acetylcarnitine ~ water_treatment, data = data_acetylcarnitine_10)$p
p2 <- compare_means(Log_butyrylcarnitine ~ water_treatment, data = butyrylcarnitine_10)$p
p3 <- compare_means(Log_valerylcarnitine ~ water_treatment, data = data_stool_valerylcarnitine_10)$p
p_raw <- c(p1, p2, p3)
p_adj <- p.adjust(p_raw, method = "BH")  
#ggsave(plot = ratio_time_data_stool_carnitine, filename = "ratio_time_data_stool_carnitine.svg", device = "svg", dpi = "retina", width = 3.5, height = 2.3)


##################################################################################################################################################################
##################
### BILE ACIDS ###
##################
bile_acids_annotations <- annotations %>% 
  dplyr::filter(str_detect(pattern = regex("bile|cholic|-CA|-DCA|CHOLATE", ignore_case = TRUE), Compound_Name)) %>%
  dplyr::select(1:2, 15) %>%
  dplyr::filter(!(str_detect(pattern = regex("carboxy|Carnitine", ignore_case = TRUE), Compound_Name))) %>%
  dplyr::mutate(Bile_Type = case_when(str_detect(Compound_Name, regex("tauro", ignore_case = TRUE)) ~ "taurine",
                                      str_detect(Compound_Name, regex("gly", ignore_case = TRUE)) ~ "glycine",
                                      str_detect(Compound_Name, "-CA") ~ "Amino acid",
                                      str_detect(Compound_Name, "oxo") ~ "oxo",
                                      TRUE ~ "Other")) %>% dplyr::mutate_at("X.Scan.", as.character)
colnames(BA_query_filter)[1] <- "X.Scan."

# Remove features annotated as BAs that did not pass validation
ba_interst_filter <- bile_acids_annotations %>% 
  dplyr::filter(!X.Scan. %in% BA_query_filter$X.Scan.)

taurine_conjugated <- ba_interst_filter %>%
  dplyr::filter(Bile_Type == "taurine")

glycine_conjugated <- ba_interst_filter %>%
  dplyr::filter(Bile_Type == "glycine")

aa_conjugated <- ba_interst_filter %>%
  dplyr::filter(Bile_Type == "Amino acid")

oxo_conjugated <- ba_interst_filter %>%
  dplyr::filter(Bile_Type == "oxo")

free_di <- ba_interst_filter %>%
  dplyr::filter(Bile_Type == "Other") %>%
  dplyr::filter(str_detect(pattern = regex("deoxy", ignore_case = TRUE), Compound_Name))

free_tri <- ba_interst_filter %>%
  dplyr::filter(Bile_Type == "Other") %>%
  dplyr::filter(!(str_detect(pattern = regex("deoxy", ignore_case = TRUE), Compound_Name))) %>%
  dplyr::filter(str_detect(pattern = regex("chol", ignore_case = TRUE), Compound_Name)) %>%
  dplyr::filter(!(str_detect(pattern = regex("penta|tetra|sulfated|acetylated|putrescine", ignore_case = TRUE), Compound_Name))) %>%
  dplyr::filter(`X.Scan.` != 10567)

data_stool_bas <- data_stool %>% 
  dplyr::select("SampleID", taurine_conjugated$X.Scan., glycine_conjugated$X.Scan., 
                aa_conjugated$X.Scan., oxo_conjugated$X.Scan., free_di$X.Scan., free_tri$X.Scan.) %>%
  dplyr::mutate(tau_sum = rowSums(dplyr::select(., taurine_conjugated$X.Scan.))) %>%
  dplyr::mutate(gly_sum = rowSums(dplyr::select(., glycine_conjugated$X.Scan.))) %>%
  dplyr::mutate(aa_sum = rowSums(dplyr::select(., aa_conjugated$X.Scan.))) %>% 
  dplyr::mutate(oxo_sum = rowSums(dplyr::select(., oxo_conjugated$X.Scan.))) %>%
  dplyr::mutate(di_sum = rowSums(dplyr::select(., free_di$X.Scan.))) %>%
  dplyr::mutate(tri_sum = rowSums(dplyr::select(., free_tri$X.Scan.))) %>%
  dplyr::select(SampleID, tau_sum, gly_sum, aa_sum, oxo_sum, di_sum, tri_sum) %>%
  dplyr::filter(!(tau_sum == 0 & gly_sum == 0)) %>%
  dplyr::mutate(ratio_tg = log(tau_sum / gly_sum)) %>%
  dplyr::mutate(ratio_at = log(aa_sum / tau_sum+ 1)) %>%
  dplyr::mutate(ratio_ga = log(gly_sum / aa_sum)) %>%
  dplyr::mutate(Log_tau = log2(tau_sum + 1)) %>% 
  dplyr::mutate(Log_gly = log2(gly_sum + 1)) %>% 
  dplyr::mutate(Log_aa = log2(aa_sum + 1)) %>% 
  dplyr::mutate(Log_oxo = log2(oxo_sum + 1)) %>% 
  dplyr::mutate(Log_di = log2(di_sum + 1)) %>% 
  dplyr::mutate(Log_tri = log2(tri_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID")

data_stool_bas %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_tau", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_bas %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_gly", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_bas %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_aa", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_bas %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_oxo", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_bas %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_di", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_bas %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_tri", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_bas %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "ratio_at", add = "loess", alpha = 0.1) #+ ylim(0, 28)

# Combine in two plots
free_ba_plot <- data_stool_bas %>% pivot_longer(cols = c("Log_di", "Log_tri"), names_to = "BA", values_to = "Value") %>%
  ggscatter(x = "infant_age_days", y = "Value", color = "BA",  legend = "none", title = "Free Bile Acids",
            add = "loess", alpha = 0.1, palette = c("#1fa187", "#a0da39"), xlab = "Age (days)", ylab = "Log2(Peak Area)") + ylim(0, 28) +   
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

modified_ba_plot <- data_stool_bas %>% pivot_longer(cols = c("Log_tau", "Log_gly", "Log_aa", "Log_oxo"), names_to = "BA", values_to = "Value") %>%
  ggscatter(x = "infant_age_days", y = "Value", color = "BA",  legend = "none", title = "Conjugated Bile Acids",
            add = "loess", alpha = 0.1, palette = c("#46327e", "#277f8e", "#4ac16d", "#FDE725FF"), xlab = "Age (days)", ylab = "Log2(Peak Area)") + ylim(0, 28) +   
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

comb_ba_time_plot <- ggarrange(free_ba_plot, modified_ba_plot, nrow = 1)
#ggsave(plot = comb_ba_time_plot, filename = "BA_time.svg", device = "svg", dpi = "retina", height = 2, width = 5)

# Mixed effect models
model <- data_stool_bas %>% 
  lmer(formula = Log_di ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_bas %>% 
  lmer(formula = Log_tri ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_bas %>% 
  lmer(formula = Log_tau ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_bas %>% 
  lmer(formula = Log_gly ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_bas %>% 
  lmer(formula = Log_aa ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_bas %>% 
  lmer(formula = Log_oxo ~ infant_age_days + (1|host_subject_id))
summary(model)


# Check ratio AA/Tau
data_stool_bas %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "ratio_at", add = "loess", alpha = 0.1, color = "mode_delivery") # small diff

# Ratio AA/Tau for mode of delivery
# Divide into 3 bins, 0 - 7, 14 - 28, 30 - 45 and keep the highest recorded value per subject
data_ba_10 <- data_stool_bas %>% 
  dplyr::select(host_subject_id, infant_age_days, ratio_at, mode_delivery) %>%
  dplyr::filter(infant_age_days < 8) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(ratio_at == max(ratio_at, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE)

data_ba_10_plot <- data_ba_10 %>% arrange(mode_delivery) %>%
  ggboxplot(x = "mode_delivery", y = "ratio_at", add = "jitter", legend = "none",
            add.params = list(color = "mode_delivery", alpha = 0.5), palette = palette_del, 
            xlab = "Days (0-7)", ylab = "Log(AA/Tau)") + stat_compare_means() + ylim (-0.1, 2.5) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_ba_30 <- data_stool_bas %>% 
  dplyr::select(host_subject_id, infant_age_days, ratio_at, mode_delivery, water_filter) %>%
  dplyr::filter(infant_age_days > 13 & infant_age_days < 29) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(ratio_at == max(ratio_at, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE)

data_ba_30_plot <- data_ba_30 %>% arrange(mode_delivery) %>%
  ggboxplot(x = "mode_delivery", y = "ratio_at", add = "jitter", legend = "none",
            add.params = list(color = "mode_delivery", alpha = 0.5), palette = palette_del,
            xlab = "Days (14-28)", ylab = FALSE) + stat_compare_means() + ylim (-0.1, 2.5) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

data_ba_60 <- data_stool_bas %>% 
  dplyr::select(host_subject_id, infant_age_days, ratio_at, mode_delivery) %>%
  dplyr::filter(infant_age_days > 29 & infant_age_days < 46) %>%
  group_by(host_subject_id) %>%
  dplyr::filter(ratio_at == max(ratio_at, na.rm = TRUE)) %>%
  ungroup() %>% distinct(host_subject_id, .keep_all = TRUE)

data_ba_60_plot <- data_ba_60 %>% arrange(mode_delivery) %>%
  ggboxplot(x = "mode_delivery", y = "ratio_at", add = "jitter", legend = "none",
            add.params = list(color = "mode_delivery", alpha = 0.5), palette = palette_del,
            xlab = "Days (30-45)", ylab = FALSE) + stat_compare_means() + ylim (-0.1, 2.5) +
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

ratio_time_ba_comb <- ggarrange(data_ba_10_plot, data_ba_30_plot, data_ba_60_plot, nrow = 1)

p1 <- compare_means(ratio_at ~ mode_delivery, data = data_ba_10)$p
p2 <- compare_means(ratio_at ~ mode_delivery, data = data_ba_30)$p
p3 <- compare_means(ratio_at ~ mode_delivery, data = data_ba_60)$p
p_raw <- c(p1, p2, p3)
p_adj <- p.adjust(p_raw, method = "BH")  

#ggsave(plot = ratio_time_ba_comb, filename = "Ratio_box_delivery.svg", device = "svg", dpi = "retina", width = 3, height = 2)



##################
### N ACYL LIPIDS ###
##################
nacyl_annotations <- annotations %>%
  dplyr::mutate(N_acyl_type = case_when(
    str_detect(Compound_Name, "C2") ~ "C2",
    str_detect(Compound_Name, "C3") ~ "C3",
    str_detect(Compound_Name, "C4") ~ "C4",
    str_detect(Compound_Name, "C5") ~ "C5",
    str_detect(Compound_Name, "C6") ~ "C6",
    str_detect(Compound_Name, "C7") ~ "C7",
    str_detect(Compound_Name, "C8") ~ "C8",
    str_detect(Compound_Name, "C9") ~ "C9",
    str_detect(Compound_Name, "C12") ~ "C12",
    str_detect(Compound_Name, "C14") ~ "C14",
    str_detect(Compound_Name, "C16") ~ "C16",
    str_detect(Compound_Name, "C18") ~ "C18",
    TRUE ~ "Other"))

nacyl_lipids <- nacyl_annotations %>% 
  dplyr::filter(str_detect(LibraryName, 
                           regex("N-ACYL-LIPIDS|ECG-ACYL-AMIDES-C4-C24-LIBRARY", ignore_case = TRUE))) 

### N-acyl lipids C2:0
putrescine_C2 <- nacyl_lipids %>%
  dplyr::filter(str_detect(Compound_Name, "N-acetylputrescine-C2:0"))

valine_C2 <- nacyl_lipids %>%
  dplyr::filter(str_detect(Compound_Name, "Valine-C2:0"))

aminovaleric_c2 <- nacyl_lipids %>%
  dplyr::filter(str_detect(Compound_Name, "2-aminovaleric acid-C2:0"))

cadaverine_c2 <- nacyl_lipids %>%
  dplyr::filter(str_detect(Compound_Name, "Candidate Cadaverine-C2:0"))

acetylcadaverine_c2 <- nacyl_lipids %>%
  dplyr::filter(str_detect(Compound_Name, "N-acetylcadaverine-C2:0"))

data_stool_nacyl_c2 <- data_stool %>% 
  dplyr::select("SampleID", putrescine_C2$X.Scan., valine_C2$X.Scan., 
                aminovaleric_c2$X.Scan., cadaverine_c2$X.Scan., acetylcadaverine_c2$X.Scan) %>%
  dplyr::mutate(putrescine_C2_sum = rowSums(dplyr::select(., putrescine_C2$X.Scan.))) %>%
  dplyr::mutate(valine_C2_sum = rowSums(dplyr::select(., valine_C2$X.Scan.))) %>%
  dplyr::mutate(aminovaleric_c2_sum = rowSums(dplyr::select(., aminovaleric_c2$X.Scan.))) %>% 
  dplyr::mutate(cadaverine_c2_sum = rowSums(dplyr::select(., cadaverine_c2$X.Scan.))) %>%
  dplyr::mutate(acetylcadaverine_c2_sum = rowSums(dplyr::select(., acetylcadaverine_c2$X.Scan.))) %>%
  dplyr::select(SampleID, putrescine_C2_sum, valine_C2_sum, aminovaleric_c2_sum, cadaverine_c2_sum, acetylcadaverine_c2_sum) %>%
  dplyr::mutate(Log_putrescine_C2 = log2(putrescine_C2_sum + 1)) %>% 
  dplyr::mutate(Log_valine_C2 = log2(valine_C2_sum + 1)) %>% 
  dplyr::mutate(Log_aminovaleric_c2 = log2(aminovaleric_c2_sum + 1)) %>% 
  dplyr::mutate(Log_cadaverine_c2 = log2(cadaverine_c2_sum + 1)) %>% 
  dplyr::mutate(Log_acetylcadaverine_c2 = log2(acetylcadaverine_c2_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID")

data_stool_nacyl_c2 %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_putrescine_C2", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_nacyl_c2 %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_valine_C2", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_nacyl_c2 %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_aminovaleric_c2", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_nacyl_c2 %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_cadaverine_c2", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_nacyl_c2 %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>% 
  ggscatter(x = "infant_age_days", y = "Log_acetylcadaverine_c2", add = "loess", alpha = 0.1)  + ylim(0, 28)


nacyl_cad_plot_C2 <- data_stool_nacyl_c2 %>% pivot_longer(cols = c("Log_putrescine_C2", "Log_valine_C2", 
                                                                   "Log_aminovaleric_c2", "Log_cadaverine_c2", 
                                                                   "Log_acetylcadaverine_c2"), 
                                                          names_to = "NACYL", values_to = "Value") %>%
  ggscatter(x = "infant_age_days", y = "Value", color = "NACYL",   title = "N-acyl lipids C2:0", 
            add = "loess", alpha = 0.1, palette = c("#FDE725FF", "#4ac16d", "#a0da39", "#46327e", "#277f8e"), 
            xlab = "Infant age (days)", ylab = "Log2(Peak Area)") + ylim(0, 28) +   
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Mixed effect models
model <- data_stool_nacyl_c2 %>% 
  lmer(formula = Log_putrescine_C2 ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_nacyl_c2 %>% 
  lmer(formula = Log_valine_C2 ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_nacyl_c2 %>% 
  lmer(formula = Log_aminovaleric_c2 ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_nacyl_c2 %>% 
  lmer(formula = Log_cadaverine_c2 ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_nacyl_c2 %>% 
  lmer(formula = Log_acetylcadaverine_c2 ~ infant_age_days + (1|host_subject_id))
summary(model)

### N-acyl lipids C3:0
pheala_c3 <- nacyl_lipids %>%
  dplyr::filter(str_detect(Compound_Name, "Phenylalanine-C3:0")) %>% 
  dplyr::filter(!X.Scan. %in% c("5308", "1480"))

cadaverine_c3 <- nacyl_lipids %>%
  dplyr::filter(str_detect(Compound_Name, "Candidate Cadaverine-C3:0"))

acetylcadaverine_c3 <- nacyl_lipids %>%
  dplyr::filter(str_detect(Compound_Name, "N-acetylcadaverine-C3:0"))

data_stool_nacyl_c3 <- data_stool %>% 
  dplyr::select("SampleID", pheala_c3$X.Scan., cadaverine_c3$X.Scan., acetylcadaverine_c3$X.Scan.) %>%
  dplyr::mutate(pheala_c3_sum = rowSums(dplyr::select(., pheala_c3$X.Scan.))) %>%
  dplyr::mutate(cadaverine_c3_sum = rowSums(dplyr::select(., cadaverine_c3$X.Scan.))) %>%
  dplyr::mutate(acetylcadaverine_c3_sum = rowSums(dplyr::select(., acetylcadaverine_c3$X.Scan.))) %>% 
  dplyr::select(SampleID, pheala_c3_sum, cadaverine_c3_sum, acetylcadaverine_c3_sum) %>%
  dplyr::mutate(Log_pheala_c3 = log2(pheala_c3_sum + 1)) %>% 
  dplyr::mutate(Log_cadaverine_c3= log2(cadaverine_c3_sum + 1)) %>% 
  dplyr::mutate(Log_acetylcadaverine_c3 = log2(acetylcadaverine_c3_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID")

data_stool_nacyl_c3 %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_pheala_c3", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_nacyl_c3 %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_cadaverine_c3", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_nacyl_c3 %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_acetylcadaverine_c3", add = "loess", alpha = 0.1) + ylim(0, 28)

nacyl_cad_plot_C3 <- data_stool_nacyl_c3 %>% pivot_longer(cols = c("Log_pheala_c3", "Log_cadaverine_c3", "Log_acetylcadaverine_c3"), 
                                                          names_to = "NACYL", values_to = "Value") %>%
  ggscatter(x = "infant_age_days", y = "Value", color = "NACYL",   title = "N-acyl lipids C3:0", 
            add = "loess", alpha = 0.1, palette = c("#FDE725FF", "#a0da39", "#4ac16d"), 
            xlab = "Infant age (days)", ylab = "Log2(Peak Area)") + ylim(0, 28) +   
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Mixed effect models
model <- data_stool_nacyl_c3 %>% 
  lmer(formula = Log_pheala_c3 ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_nacyl_c3 %>% 
  lmer(formula = Log_cadaverine_c3 ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_nacyl_c3 %>% 
  lmer(formula = Log_acetylcadaverine_c3 ~ infant_age_days + (1|host_subject_id))
summary(model)

comb_nacyl_plot <- ggarrange(nacyl_cad_plot_C2, nacyl_cad_plot_C3, nrow = 1)
#ggsave(plot = comb_nacyl_plot, filename = "N_acyl_lipid_time.svg", device = "svg", dpi = "retina", height = 2, width = 5)

# N-ACYL LIPIDS BY DELIVERY MODE
data_stool_nacyl_c2 %>%
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>% 
  ggscatter(x = "infant_age_days", y = "Log_valine_C2", add = "loess", alpha = 0.1, color = "mode_delivery") + ylim(0, 28)

data_stool_nacyl_c2 %>%
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>% 
  ggscatter(x = "infant_age_days", y = "Log_aminovaleric_c2", add = "loess", alpha = 0.1, color = "mode_delivery") + ylim(0, 28)

nacyl_delivery_plot <- data_stool_nacyl_c2 %>%
  pivot_longer(cols = c("Log_valine_C2", "Log_aminovaleric_c2"),
               names_to = "NACYL", values_to = "Value") %>%
  ggplot(aes(x = infant_age_days, y = Value, color = mode_delivery, linetype = NACYL)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +
  ylim(0, 20) +
  labs(title = "N-acyl Lipids C2:0 by Delivery Mode",
       x = "Infant age (days)", y = "Log2(Peak Area)", color = "Delivery Mode", linetype = "Feature") +
  scale_color_manual(values = palette_del) +
  theme_classic() +
  theme(plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.position = "none")

#ggsave(plot = nacyl_delivery_plot, filename = "nacyl_delivery_plot.svg", device = "svg", dpi = "retina", height = 2.4, width = 2.5)

model <- data_stool_nacyl_c2 %>% 
  lmer(formula = Log_valine_C2 ~ infant_age_days * mode_delivery + (1|host_subject_id))
summary(model)

model <- data_stool_nacyl_c2 %>% 
  lmer(formula = Log_aminovaleric_c2 ~ infant_age_days * mode_delivery + (1|host_subject_id))
summary(model)

##################################################################################################################################################################
##################
### ACYL CARNITINES ###
##################

# Short chain acylcarnitines
carnitine_c2 <- annotations %>%
  dplyr::filter(str_detect(Compound_Name, "Acetylcarnitine|Acetyl-DL-carnitine"))

carnitine_c3 <- annotations %>%
  dplyr::filter(str_detect(Compound_Name, "Propionylcarnitine"))

carnitine_c4 <- annotations %>%
  dplyr::filter(str_detect(Compound_Name, "Butyrylcarnitine")) 

carnitine_c5 <- annotations %>%
  dplyr::filter(str_detect(Compound_Name, "Isovalerylcarnitine|pentanoylcarnitine"))

data_stool_carnitine_short <- data_stool %>% 
  dplyr::select("SampleID", carnitine_c2$X.Scan., carnitine_c3$X.Scan., 
                carnitine_c4$X.Scan., carnitine_c5$X.Scan.) %>%
  dplyr::mutate(carnitine_c2_sum = rowSums(dplyr::select(., carnitine_c2$X.Scan.))) %>%
  dplyr::mutate(carnitine_c3_sum = rowSums(dplyr::select(., carnitine_c3$X.Scan.))) %>%
  dplyr::mutate(carnitine_c4_sum = rowSums(dplyr::select(., carnitine_c4$X.Scan.))) %>% 
  dplyr::mutate(carnitine_c5_sum = rowSums(dplyr::select(., carnitine_c5$X.Scan.))) %>%
  dplyr::select(SampleID, carnitine_c2_sum, carnitine_c3_sum, carnitine_c4_sum, carnitine_c5_sum) %>%
  dplyr::mutate(Log_carnitine_c2 = log2(carnitine_c2_sum + 1)) %>% 
  dplyr::mutate(Log_carnitine_c3 = log2(carnitine_c3_sum + 1)) %>% 
  dplyr::mutate(Log_carnitine_c4 = log2(carnitine_c4_sum + 1)) %>% 
  dplyr::mutate(Log_carnitine_c5 = log2(carnitine_c5_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID")

data_stool_carnitine_short %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_carnitine_c2", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_carnitine_short %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_carnitine_c3", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_carnitine_short %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_carnitine_c4", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_carnitine_short %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_carnitine_c5", add = "loess", alpha = 0.1) + ylim(0, 28)

carnitine_short_plot <- data_stool_carnitine_short %>% 
  pivot_longer(cols = c("Log_carnitine_c2", "Log_carnitine_c3", "Log_carnitine_c4", "Log_carnitine_c5"), 
               names_to = "Carnitine", values_to = "Value") %>%
  ggscatter(x = "infant_age_days", y = "Value", color = "Carnitine",   title = "Short chain acylcarnitines", 
            add = "loess", alpha = 0.1, palette = c("#4ac16d", "#a0da39", "#46327e", "#277f8e"), 
            xlab = "Infant age (days)", ylab = "Log2(peak area)") + ylim(0, 28) +   
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6), legend.position = "none")

model <- data_stool_carnitine_short %>% 
  lmer(formula = Log_carnitine_c2 ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_carnitine_short %>% 
  lmer(formula = Log_carnitine_c3 ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_carnitine_short %>% 
  lmer(formula = Log_carnitine_c4 ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_carnitine_short %>% 
  lmer(formula = Log_carnitine_c5 ~ infant_age_days + (1|host_subject_id))
summary(model)


# Long chain acylcarnitines
carnitine_c16 <- annotations %>%
  dplyr::filter(str_detect(Compound_Name, "Palmitoylcarnitine"))

carnitine_c18_0 <- annotations %>%
  dplyr::filter(str_detect(Compound_Name, "Stearoyl-L-Carnitine"))

carnitine_c18_1 <- annotations %>%
  dplyr::filter(str_detect(Compound_Name, "Oleoyl L-carnitine")) 

data_stool_carnitine_long <- data_stool %>% 
  dplyr::select("SampleID", carnitine_c16$X.Scan., carnitine_c18_0$X.Scan., carnitine_c18_1$X.Scan.) %>%
  dplyr::mutate(carnitine_c16_sum = rowSums(dplyr::select(., carnitine_c16$X.Scan.))) %>%
  dplyr::mutate(carnitine_c18_0_sum = rowSums(dplyr::select(., carnitine_c18_0$X.Scan.))) %>%
  dplyr::mutate(carnitine_c18_1_sum = rowSums(dplyr::select(., carnitine_c18_1$X.Scan.))) %>% 
  dplyr::select(SampleID, carnitine_c16_sum, carnitine_c18_0_sum, carnitine_c18_1_sum) %>%
  dplyr::mutate(Log_carnitine_c16 = log2(carnitine_c16_sum + 1)) %>% 
  dplyr::mutate(Log_carnitine_c18_0 = log2(carnitine_c18_0_sum + 1)) %>% 
  dplyr::mutate(Log_carnitine_c18_1 = log2(carnitine_c18_1_sum + 1)) %>% 
  left_join(metadata_metabolomics, by = "SampleID")

data_stool_carnitine_long %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_carnitine_c16", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_carnitine_long %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_carnitine_c18_0", add = "loess", alpha = 0.1) + ylim(0, 28)
data_stool_carnitine_long %>% 
  dplyr::mutate(LogAge = log2(infant_age_days + 1)) %>%
  ggscatter(x = "infant_age_days", y = "Log_carnitine_c18_1", add = "loess", alpha = 0.1) + ylim(0, 28)

carnitine_long_plot <- data_stool_carnitine_long %>% 
  pivot_longer(cols = c("Log_carnitine_c16", "Log_carnitine_c18_0", "Log_carnitine_c18_1"), 
               names_to = "Carnitine", values_to = "Value") %>%
  ggscatter(x = "infant_age_days", y = "Value", color = "Carnitine",   title = "Long chain acylcarnitines", 
            add = "loess", alpha = 0.1, palette = c("#46327e", "#277f8e", "#a0da39"), 
            xlab = "Infant age (days)", ylab = "Log2(peak area)") + ylim(0, 28) +   
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6), legend.position = "none")

# LME
model <- data_stool_carnitine_long %>% 
  lmer(formula = Log_carnitine_c16 ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_carnitine_long %>% 
  lmer(formula = Log_carnitine_c18_0 ~ infant_age_days + (1|host_subject_id))
summary(model)

model <- data_stool_carnitine_long %>% 
  lmer(formula = Log_carnitine_c18_1 ~ infant_age_days + (1|host_subject_id))
summary(model)

comb_carnitine_plot <- ggarrange(carnitine_short_plot, carnitine_long_plot, nrow = 1)
#ggsave(plot = comb_carnitine_plot, filename = "comb_carnitine_plot.svg", device = "svg", dpi = "retina", height = 2, width = 5)




##################################################################################################################################################################
##########
# PLASMA #
##########
data_plasma <- data_sample_filter %>% 
  dplyr::filter(SampleID %in% sample_plasma$SampleID) %>% 
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")


##################
### SIRIUS/CANOPUS OUTPUT ###
##################

class_predictions_filtered <- canopus %>%
  dplyr::mutate(
    NPC.combined = dplyr::case_when(
      NPC.class.Probability >= 0.5 ~ NPC.class,
      NPC.class.Probability < 0.5 & NPC.superclass.Probability >= 0.5 ~ NPC.superclass,
      NPC.class.Probability < 0.5 & NPC.superclass.Probability < 0.5 & NPC.pathway.Probability >= 0.5 ~ NPC.pathway,
      TRUE ~ NA_character_)) %>%
  dplyr::select(mappingFeatureId, NPC.combined) %>%
  distinct(mappingFeatureId, .keep_all = TRUE)

### UPSET PLOT ALL BIOFLUIDS
plasma_t <- data_plasma %>% 
  pivot_longer(cols = -c(SampleID), names_to = "Feature", values_to = "Value") %>% 
  dplyr::filter(Value != 0) %>% distinct(Feature) %>% 
  left_join(class_predictions_filtered, by = c("Feature" = "mappingFeatureId"))

milk_t <- data_milk %>% 
  pivot_longer(cols = -c(SampleID), names_to = "Feature", values_to = "Value") %>% 
  dplyr::filter(Value != 0) %>% distinct(Feature) %>% 
  left_join(class_predictions_filtered, by = c("Feature" = "mappingFeatureId"))

stool_t <- data_stool %>% 
  pivot_longer(cols = -c(SampleID), names_to = "Feature", values_to = "Value") %>% 
  dplyr::filter(Value != 0) %>% distinct(Feature) %>% 
  left_join(class_predictions_filtered, by = c("Feature" = "mappingFeatureId"))

all_features_class <- bind_rows(
  plasma_t %>% mutate(matrix = "plasma"),
  milk_t   %>% mutate(matrix = "milk"),
  stool_t  %>% mutate(matrix = "stool")
)

# Create wide presence/absence table
feature_matrix <- all_features_class %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = matrix, values_from = present, values_fill = 0) %>%
  distinct(Feature, NPC.combined, .keep_all = TRUE)

# Calculate % per class
top_classes <- feature_matrix %>%
  count(NPC.combined) %>%
  mutate(pct = n / sum(n)) %>%
  mutate(NPC.simplified = ifelse(pct < 0.008, "z_Other", NPC.combined))

# Join simplified class back to main data
feature_matrix <- feature_matrix %>%
  left_join(top_classes %>% select(NPC.combined, NPC.simplified), by = "NPC.combined")

class_colors <- c(
  # Amino-related (pink/purple tones)
  "Aminoacids" = "#fec1d3",
  "Amino acids and Peptides" = "#fcb3c7",
  "Aminosugars" = "#fda9c1",
  "Dipeptides" = "#fe90ad",
  "Tripeptides" = "#ff6b9b",
  "Small peptides" = "#f75284",
  "Cyclic peptides" = "#fe558c",
  "Oligopeptides" = "#fb2d73",
  
  # Fatty/lipid-related (greens & neutrals)
  "Fatty acids" = "#ffe192",
  "Fatty amides" = "#ffd670",
  "Fatty acyl carnitines" = "#fed376",
  "N-acyl amines" = "#ffb515",
  "Fatty alcohols" = "#b4f1a0",
  "Ceramides" = "#77dd77",
  "Diacylglycerols" = "#32cd32",
  "Triacylglycerols" = "#228b22",
  "Cholane steroids" = "#a4c195",
  "Steroids" = "#77a361",
  "Terpenoids" = "#2f4f4f",
  
  # Polar lipids (turquoise/blue)
  "Glycerophosphocholines" = "#bec6f9",
  "Glycerophosphoethanolamines" = "#aeb8f8",
  "Phosphosphingolipids" = "#9ba7f3",
  "Sphingoid bases" = "#8492f1",
  "Spingolipids" = "#7084ee",
  
  "Alkaloids" = "#c890e5",
  "Pyridine alkaloids" = "#b06dd1",
  
  # Other & signaling (blue-neutral)
  "z_Other" = "#C9E3F8",
  
  # NA
  "NA" = "#F5F5F5"
)

class_order <- c(
  "Aminoacids",
  "Amino acids and Peptides",
  "Aminosugars",
  "Dipeptides",
  "Tripeptides",
  "Small peptides",
  "Cyclic peptides",
  "Oligopeptides",
  "Fatty acids",
  "Fatty amides",
  "Fatty acyl carnitines",
  "N-acyl amines",
  "Fatty alcohols",
  "Ceramides",
  "Diacylglycerols",
  "Triacylglycerols",
  "Cholane steroids",
  "Steroids",
  "Terpenoids",
  "Glycerophosphocholines",
  "Glycerophosphoethanolamines",
  "Phosphosphingolipids",
  "Sphingoid bases",
  "Spingolipids",
  "Alkaloids",
  "Pyridine alkaloids",
  "z_Other",
  "NA"
)

feature_matrix$NPC.simplified <- factor(
  feature_matrix$NPC.simplified,
  levels = class_order
)


# Plot with custom fill
#svg("upset_plot_all_matrices.svg", width = 7, height = 7.5)
upset(
  feature_matrix,
  intersect = c("plasma", "milk", "stool"),
  name = "Matrix Overlap",
  annotations = list(
    'Class %' = (
      ggplot(mapping = aes(x = intersection, fill = NPC.simplified)) +
        geom_bar(position = "fill") +
        ylab("Composition") +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_fill_manual(values = class_colors) +
        theme_minimal() +
        theme()
    )
  ),
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 3),
      bar_number_threshold = 0
    )
  )
)
#dev.off()



##################################################################################################################################################################
##################
### JOINT ANALYSIS ###
##################
joint_table <- read_csv("data/correlation_table_30day.csv")
taxonomy <- read_tsv("data/taxonomy.tsv")
filtered_taxa <- read_csv("data/micro_shortlist_tempted-jointRPCA-ancombc2_overlap.csv")
filtered_taxa <- filtered_taxa %>%
  dplyr::rename(featureid = "Feature ID")

info_feature_complete_filter_2 <- info_feature_complete_filter %>%
  dplyr::mutate(Compound_Name = substr(Compound_Name, 1, 60))

joint_table_interest_filtered <- joint_table %>% 
  dplyr::select(featureid,  '312', '1634', '9267', '1962', '274', '2275', '4488', '4671', '10545', '495',
                '998', '1425', '7905', '7573', '10001', '10329', '8718', '9204', '7577', '7602', '9171',
                '518', '10011', '9939', '9099', '9227', '10466', '9357', '10484', '10023', '10461', '3969', 
                '524', '786', '583', '561', '769', '2792') %>% 
  dplyr::filter(str_detect(pattern = "G", featureid)) %>% 
  dplyr::filter(featureid %in% filtered_taxa$featureid) %>% 
  dplyr::left_join(taxonomy %>% dplyr::select(1:2), by = c("featureid" = "Feature ID")) %>%
  dplyr::select(-featureid) %>% 
  dplyr::mutate(Taxon = gsub(".*; g__", "", Taxon)) %>%
  dplyr::mutate(Taxon = paste(Taxon, seq_along(Taxon), sep = "_")) %>%
  column_to_rownames("Taxon") %>% t() %>% as.data.frame() %>% rownames_to_column("Feature") %>% 
  left_join(info_feature_complete_filter_2 %>% dplyr::select(1,5)) %>%
  dplyr::mutate(Compound_Name = case_when(is.na(Compound_Name) ~ Feature,
                                          TRUE ~ Compound_Name)) %>%
  dplyr::select(-1) %>% column_to_rownames("Compound_Name") 

filtered_heatmap_joint <- pheatmap::pheatmap(joint_table_interest_filtered, color = viridis(100),   fontsize = 6,
                                             fontsize_row = 6,
                                             fontsize_col = 6)
#ggsave(plot = filtered_heatmap_joint, filename = "filtered_heatmap_joint.svg", device = "svg", dpi = "retina", width = 8, height = 8)

joint_table_interest_complete <- joint_table %>% 
  dplyr::select(featureid,  '312', '1634', '9267', '1962', '274', '2275', '4488', '4671', '10545', '495',
                '998', '1425', '7905', '7573', '10001', '10329', '8718', '9204', '7577', '7602', '9171',
                '518', '10011', '9939', '9099', '9227', '10466', '9357', '10484', '10023', '10461', '3969', 
                '524', '786', '583', '561', '769', '2792') %>% 
  dplyr::filter(str_detect(pattern = "G", featureid)) %>% 
  dplyr::left_join(taxonomy %>% dplyr::select(1:2), by = c("featureid" = "Feature ID")) %>%
  dplyr::select(-featureid) %>% 
  dplyr::mutate(Taxon = gsub(".*; g__", "", Taxon)) %>%
  dplyr::mutate(Taxon = paste(Taxon, seq_along(Taxon), sep = "_")) %>%
  column_to_rownames("Taxon") %>% t() %>% as.data.frame() %>% rownames_to_column("Feature") %>% 
  left_join(info_feature_complete_filter_2 %>% dplyr::select(1,5)) %>%
  dplyr::mutate(Compound_Name = case_when(is.na(Compound_Name) ~ Feature,
                                          TRUE ~ Compound_Name)) %>%
  dplyr::select(-1) %>% column_to_rownames("Compound_Name") 


complete_heatmap_joint <- pheatmap::pheatmap(joint_table_interest_complete, color = viridis(100),   fontsize = 6,
                                             fontsize_row = 6,
                                             fontsize_col = 6)

#ggsave(plot = complete_heatmap_joint, filename = "complete_heatmap_joint.svg", device = "svg", dpi = "retina", width = 18, height = 8)

# Export features of interest for microbeMASST analysis
joint_feat_interest <- data_stool %>%
  select(SampleID, '312', '1634', '9267', '1962', '274', '2275', '4488', '4671', '10545', '495',
         '998', '1425', '7905', '7573', '10001', '10329', '8718', '9204', '7577', '7602', '9171',
         '518', '10011', '9939', '9099', '9227', '10466', '9357', '10484', '10023', '10461', '3969', 
         '524', '786', '583', '561', '769', '2792') 

# Import full mgf and filter it
# mgf for microbeMASST 
# mgf for Sirius
#dda <- Spectra("mzmine/obs_others_iimn_fbmn.mgf", source = MsBackendMgf())
#dda_ids <- data.frame(ID = dda@backend@spectraData@listData$FEATURE_ID) %>%
#  dplyr::mutate(Interest = ID %in% colnames(joint_feat_interest))
#dda_filtered <- dda[dda_ids$Interest]
#export(dda_filtered, MsBackendMgf(), file = "joint_features_interest.mgf", exportTitle = FALSE)




##################################################################################################################################################################
##################
### CORRELATION TARGETED AND UNTARGETED HMO DATA ###
##################

# Scatterplot 2FL milk
data_milk_feat_2FL <- data_milk %>% 
  dplyr::select(SampleID, `583`) %>%
  dplyr::mutate(Log2_2FL_untargeted = log2(`583` + 1)) %>%
  left_join(hmo_data, by = "SampleID") %>% 
  dplyr::mutate(Log2_2FL_targeted = log2(`x2FL` + 1)) %>% 
  dplyr::filter(hmo_Secretor=="1") %>% 
  dplyr::filter(!Log2_2FL_untargeted=="0")

scatter_2FL <- ggscatter(data_milk_feat_2FL,
                         x = "Log2_2FL_targeted", y = "Log2_2FL_untargeted",
                         color = "#80B6A3", fill = "#80B6A3",
                         add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson") +  
  theme_classic()

# Scatterplot DFLac milk
data_milk_feat_DFLac <- data_milk %>% 
  dplyr::select(SampleID, `769`) %>%
  dplyr::mutate(Log2_dflac_untargeted = log2(`769` + 1)) %>%
  left_join(hmo_data, by = "SampleID") %>% 
  dplyr::mutate(Log2_DFLac_targeted = log2(`DFLac` + 1)) %>% 
  dplyr::filter(hmo_Secretor=="1") %>% 
  dplyr::filter(!Log2_dflac_untargeted=="0")

scatter_dflac <- ggscatter(data_milk_feat_DFLac,
                           x = "Log2_DFLac_targeted", y = "Log2_dflac_untargeted",
                           color = "#80B6A3", fill = "#80B6A3",
                           add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson") +  
  theme_classic()

combined_scatter_2FL_DFLac <- ggarrange(scatter_2FL, scatter_dflac, nrow = 1)
#ggsave(plot = combined_scatter_2FL_DFLac, filename = "combined_scatter_HMO_targeted_untargeted.svg", device = "svg", dpi = "retina", width = 4.5, height = 2.3)



### BOXPLOT HMOs - TARGETED ANALYSIS ###

hmo_stage_data <- hmo_data %>% 
  dplyr::select(host_subject_id, hmo_Secretor, infant_age_days, 5:23) %>%
  dplyr::mutate(Stage = case_when(infant_age_days < 5 ~ "Colostrum",
                                  infant_age_days >= 5 & infant_age_days < 16 ~ "Translational",
                                  infant_age_days >= 16 & infant_age_days < 90 ~ "Mature",
                                  TRUE ~ "None")) %>% 
  dplyr::filter(!hmo_Secretor == "NA") %>%
  group_by(host_subject_id, Stage) %>%
  dplyr::filter(x2FL == max(x2FL, na.rm = TRUE)) %>% ungroup() 


hmo_plot_all <- hmo_stage_data %>%
  pivot_longer(cols = 4:22, names_to = "hmo", values_to = "value") %>%
  ggboxplot(x = "Stage", y = "value", add = "jitter", legend = "none", fill = "hmo_Secretor",
            add.params = list(color = "hmo_Secretor", alpha = 0.5), palette = palette_hmo_rev, 
            xlab = "Stage", ylab = "nmol/mL", facet.by = "hmo", scales = "free_y") + 
  theme(plot.title = element_text(size = 9),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))


#ggsave(plot = hmo_plot_all, filename = "hmo_plot_all.svg", device = "svg", dpi = "retina", width = 8, height = 8)

pvals_3FL <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(x3FL ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_LNFPI <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(LNFPI ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_LNFPII <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(LNFPII ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_DFLNT <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(DFLNT ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_LNFPIII <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(LNFPIII ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_FLNH <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(FLNH ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_DFLNH <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(DFLNH ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_x3SL <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(x3SL ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_x6SL <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(x6SL ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_LSTb <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(LSTb ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_LSTc <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(LSTc ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_DSLNT <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(DSLNT ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_DSLNH <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(DSLNH ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_FDSLNH <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(FDSLNH ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_LNT <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(LNT ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_LNnT <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(LNnT ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

pvals_LNH <- hmo_stage_data %>%
  group_by(Stage) %>%
  group_modify(~ compare_means(LNH ~ hmo_Secretor, data = .x, method = "wilcox.test")) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))


####################################################################################
### Table to import into cytoscape ###
### Maternal secretor status - milk and feces - network

VIPs_milk_secretor_Load_cyt <- VIPs_milk_secretor_Load %>%
  dplyr::mutate(Milk_secretor = if_else(GroupContrib == 1, "Yes", "No")) %>% 
  dplyr::select(-GroupContrib)

VIPs_stool_secretor_Load_cyt <- VIPs_stool_secretor_Load %>%
  dplyr::mutate(Feces_secretor = if_else(GroupContrib == 1, "Yes", "No")) %>% 
  dplyr::select(-GroupContrib)

VIP_milk_feces_hmo <- VIPs_milk_secretor_Load_cyt %>% 
  full_join(VIPs_stool_secretor_Load_cyt, by = c("ID", "mz", "RT", "Compound_Name")) %>% 
  mutate(
    Se_status = case_when(
      Feces_secretor == "Yes" & Milk_secretor == "Yes" ~ "Yes_both",
      Feces_secretor == "Yes" & Milk_secretor == "No"  ~ "Yes_Feces",
      Feces_secretor == "No"  & Milk_secretor == "Yes" ~ "Yes_Milk",
      Feces_secretor == "No"  & Milk_secretor == "No"  ~ "No_both",
      is.na(Feces_secretor) & Milk_secretor == "Yes"   ~ "Yes_Milk",
      Feces_secretor == "Yes" & is.na(Milk_secretor)   ~ "Yes_Feces",
      is.na(Feces_secretor) & Milk_secretor == "No"   ~ "No_Milk",
      Feces_secretor == "No" & is.na(Milk_secretor)   ~ "No_Feces",
      TRUE ~ NA_character_)) %>% 
  mutate(mean_VIP = rowMeans(across(c(comp1.x, comp1.y)), na.rm = TRUE)) %>% 
  left_join(canopus, by =c ("ID" = "mappingFeatureId"))

#write_csv(x = VIP_milk_feces_hmo, file = "VIP_Milk_Feces_HMO_Canopus.csv")


