
############################################################
# Clinical data plot
# Project: Breast Cancer Immune Heterogeneity
# Author: Saranya Narayana
# Purpose:
#   - Plot different features of clinical data
############################################################

############################################################
setwd("/home/sara/BioInfo_projects/breast_cancer_immune_rnaseq/")
############################################################
##---------------------------------------------------------------------
## Final plots for clinical data
##---------------------------------------------------------------------
colnames(clinical_filtered_4)
dir.create("results/figures/clinical_overview", recursive = TRUE, showWarnings = FALSE)

##---------------------------------------------------------
## Plot 1: DEMOGRAPHIC OVERVIEW-Age distribution
##---------------------------------------------------------

p1 <- ggplot(clinical_filtered_4, aes(x = age_at_diagnosis_years)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(aes(xintercept = median(age_at_diagnosis_years)), 
             color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = median(clinical_filtered_4$age_at_diagnosis_years) + 5, 
           y = Inf, vjust = 2,
           label = paste("Median:", round(median(clinical_filtered_4$age_at_diagnosis_years), 1), "years"),
           color = "red", size = 4) +
  labs(title = "Age Distribution at Diagnosis",
       x = "Age (years)", y = "Number of Patients") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/01_age_distribution.png", 
       p1, width = 8, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 2: DEMOGRAPHIC OVERVIEW- Age by PAM50 subtype
##---------------------------------------------------------

p2 <- ggplot(clinical_filtered_4, aes(x = PAM50_Subtype, y = age_at_diagnosis_years, 
                                        fill = PAM50_Subtype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.8) +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 1)),
               vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Age Distribution by PAM50 Subtype",
       x = "PAM50 Subtype", y = "Age at Diagnosis (years)") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/02_age_by_pam50.png", 
       p2, width = 10, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 3: DEMOGRAPHIC OVERVIEW- Age groups distribution
##---------------------------------------------------------
p3 <- ggplot(clinical_filtered_4, aes(x = age_group, fill = age_group)) +
  geom_bar(alpha = 0.8) +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(title = "Distribution by Age Group",
       x = "Age Group", y = "Number of Patients") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/03_age_groups.png", 
       p3, width = 8, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 4: PAM50 SUBTYPE ANALYSIS:PAM50 subtype distribution
##---------------------------------------------------------
pam50_counts <- clinical_filtered_4 %>%
  count(PAM50_Subtype) %>%
  mutate(percentage = n / sum(n) * 100)

p4 <- ggplot(pam50_counts, aes(x = reorder(PAM50_Subtype, -n), y = n, 
                                 fill = PAM50_Subtype)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")), 
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "PAM50 Molecular Subtype Distribution",
       subtitle = paste0("Total patients: ", nrow(clinical_filtered_4)),
       x = "PAM50 Subtype", y = "Number of Patients") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/04_pam50_distribution.png", 
       p4, width = 10, height = 6, dpi = 300)

##---------------------------------------------------------
## Plot 5: PAM50 SUBTYPE ANALYSIS:PAM50 pie chart
##---------------------------------------------------------  

p5 <- ggplot(pam50_counts, aes(x = "", y = n, fill = PAM50_Subtype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(PAM50_Subtype, "\n", round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "PAM50 Subtype Distribution (Pie Chart)") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        legend.position = "right")

ggsave("results/figures/clinical_overview/05_pam50_pie.png", 
       p5, width = 10, height = 8, dpi = 300)


##---------------------------------------------------------
## Plot 6: STAGE ANALYSIS:Stage distribution
##---------------------------------------------------------

stage_counts <- clinical_filtered_4 %>%
  count(paper_pathologic_stage) %>%
  mutate(percentage = n / sum(n) * 100)

p6 <- ggplot(stage_counts, aes(x = reorder(paper_pathologic_stage, -n), 
                                 y = n, fill = paper_pathologic_stage)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")), 
            vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(title = "Pathologic Stage Distribution",
       x = "Pathologic Stage", y = "Number of Patients") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/figures/clinical_overview/06_stage_distribution.png", 
       p6, width = 10, height = 6, dpi = 300)

##---------------------------------------------------------
## Plot 7: STAGE ANALYSIS:Stage by PAM50 subtype
##---------------------------------------------------------
p7 <- ggplot(clinical_filtered_4, aes(x = PAM50_Subtype, fill = paper_pathologic_stage)) +
  geom_bar(position = "fill", alpha = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(title = "Stage Distribution by PAM50 Subtype",
       x = "PAM50 Subtype", y = "Proportion of Patients",
       fill = "Stage") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/07_stage_by_pam50.png", 
       p7, width = 10, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 8: STAGE ANALYSIS:Stage by PAM50 subtype
##---------------------------------------------------------

p8 <- ggplot(clinical_filtered_4, aes(x = PAM50_Subtype, fill = paper_pathologic_stage)) +
  geom_bar(position = "dodge", alpha = 0.8) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(title = "Stage Distribution by PAM50 Subtype (Counts)",
       x = "PAM50 Subtype", y = "Number of Patients",
       fill = "Stage") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/08_stage_by_pam50_counts.png", 
       p8, width = 10, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 9: Violin plot: PAM50 subtype vs age at diagnosis
##---------------------------------------------------------

p9 <-ggplot(clinical_filtered_4, 
       aes(x = PAM50_Subtype, y = age_at_diagnosis_years, fill = PAM50_Subtype)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8) +
  geom_jitter(alpha = 0.2, width = 0.1, size = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Age Distribution by PAM50 Subtype",
       x = "PAM50 Subtype", y = "Age at Diagnosis (years)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("results/figures/clinical_overview/09_age_violin_pam50.png", p9, width = 10, height = 6, dpi = 300)

##---------------------------------------------------------
## Plot 10: Violin plot: stage vs age at diagnosis
##---------------------------------------------------------

p10 <- ggplot(clinical_filtered_4 %>% filter(paper_pathologic_stage != "Unknown"), 
       aes(x = paper_pathologic_stage, y = age_at_diagnosis_years, 
           fill = paper_pathologic_stage)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8) +
  geom_jitter(alpha = 0.2, width = 0.1, size = 0.5) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(title = "Age Distribution by Stage",
       x = "Stage", y = "Age at Diagnosis (years)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("results/figures/clinical_overview/10_age_violin_stage.png", p10, width = 10, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 11: DEMOGRAPHIC CHARACTERISTICS - Gender distribution
##---------------------------------------------------------

p11 <- ggplot(clinical_filtered_4, aes(x = gender, fill = gender)) +
  geom_bar(alpha = 0.8) +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Female" = "#FF69B4", "Male" = "#4169E1")) +
  labs(title = "Gender Distribution",
       x = "Gender", y = "Number of Patients") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/11_gender_distribution.png", 
       p11, width = 8, height = 6, dpi = 300)


##---------------------------------------------------------
## Plot 11: DEMOGRAPHIC CHARACTERISTICS - Race distribution
##---------------------------------------------------------

race_counts <- clinical_filtered_4 %>%
  count(race) %>%
  arrange(desc(n))

p12 <- ggplot(race_counts, aes(x = reorder(race, n), y = n, fill = race)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = n), hjust = -0.3, size = 4, fontface = "bold") +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Race/Ethnicity Distribution",
       x = "Race", y = "Number of Patients") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/clinical_overview/12_race_distribution.png", 
       p12, width = 10, height = 6, dpi = 300)
