
########################################################################
#clean.workplace
########################################################################
rm(list = ls())
dev.off()
cat("\014")
########################################################################
#Library and Function
########################################################################

library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)
library(parallel)
library(Hmisc)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(cowplot)
library(FactoMineR)

`%ni%` <- Negate(`%in%`)


########################################################################
#Figure 2B
########################################################################

data.plot #the phenotypic data in UKBB
data.plot$age_stage              <- ifelse(data.plot$age_first_visit <= 50, 1, ifelse(data.plot$age_first_visit > 60, 3, 2))

lt        <- list(
  
  Age50   =  data.plot[data.plot$age_stage ==1,],
  Age60   =  data.plot[data.plot$age_stage ==2,],
  Age70   =  data.plot[data.plot$age_stage ==3,]
  
)

p_MS     <- mapply(data = lt, age = names(lt), function(data, age){ #data = lt[[3]] age = names(lt)[1]

  data$match                 <- paste0(data$metabolic_syndrome_i0, "_", data$p31)
  data$match                 <- as.factor(data$match)
  tmp_male                   <- lm(MHS ~ metabolic_syndrome_i0  + body_mass_index_BMI + age_first_visit + initial_50k_release + pc_1 + pc_2 + pc_3 + pc_4 + pc_5 + pc_6 + pc_7 + pc_8 + pc_9 + pc_10, data = data[data$p31 == "Male",])
  summary(tmp_male)
  tmp_female                 <- lm(MHS ~ metabolic_syndrome_i0  + body_mass_index_BMI + age_first_visit + initial_50k_release + pc_1 + pc_2 + pc_3 + pc_4 + pc_5 + pc_6 + pc_7 + pc_8 + pc_9 + pc_10, data = data[data$p31 == "Female",])
  summary(tmp_female)
  #data$metabolic_syndrome_i0 <- as.factor(data$metabolic_syndrome_i0)
  pdensity <- ggdensity(
    data, x = "MHS", 
    color= "match", palette = c("#a1d99b", "#1b7837", "#bcbddc","#af8dc3"),
    size = 1.5,
    alpha = 0
  ) + ggtitle(age) +
    annotate("text", x=-6, y=0.3, label= paste0("Male: ", nrow(data[data$p31 == "Male",]) , "\n", "Female: ",  nrow(data[data$p31 == "Female",]) )) +
    annotate("text", x=4, y=0.3, label= paste0("Coefficient:\nMetS (Male): ", round(tmp_male$coefficients[2], digits = 4), "***\n", "MetS (Female): ", round(tmp_female$coefficients[2], digits = 4), "***")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "left")  +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)), limits = c(-8, 8), breaks = c(-8, -4, 0, 4, 8)) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text.x = element_text(size = 8, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 8),
      strip.background = element_rect( fill = "white"),
      strip.text = element_text(size = 8),
      legend.position = "bottom"
      #panel.border = element_blank()
    ) 
}, SIMPLIFY = F)

p_MS_all  <- plot_grid(plotlist = p_MS, align = "hv", ncol = 1)


########################################################################
#Figure 2C
#check the correlation between metabolic health score and 
#other general traits
########################################################################

library(jtools)
library(broom)
library(cli)
library(forestmangr)

general_traits_data <- lapply(general_traits, function(x){
  ols               <- lm(healthscore_v2 ~ x + body_mass_index_BMI + genetic_sex_male + age + age_squared + age_x_sex_male + age_squared_x_sex_male + initial_50k_release + pc_1 + pc_2 + pc_3 + pc_4 + pc_5 + pc_6 + pc_7 + pc_8 + pc_9 + pc_10, data = data.plot)
  model_output      <- tidy(ols)
  out_conf          <- tidy(ols, conf.int = TRUE)
  lm_model_out      <- round_df(out_conf, digits=2)
  lm_model_out      <- lm_model_out[2,] #remove the intercept 
  lm_model_out$phenotype <- x
  lm_model_out
})

general_traits_data            <- do.call(rbind, general_traits_data)
general_traits_data            <- general_traits_data[!is.na(general_traits_data$estimate), ]
general_traits_data$adjP       <- p.adjust(general_traits_data$p.value, method = "BH")
general_traits_data_select     <- general_traits_data[general_traits_data$adjP < 0.05,]
general_traits_data_select     <- general_traits_data_select[general_traits_data_select$estimate !=0,]


phenotype_choose               <- c("whole_body_fat_percentage (n =217547)", "a_body_shape_index_ABSI (n =216648)", "waist_hip_ratio (n =217522)", "albumin_creatinine_ratio_in_urine (n =57089)",
                                    "Apolipoprotein_B (n =217057)", "Triglycerides (n =217548)", "FSI (n =215954)", "Glucose (n =217548)", "ALAT_ASAT_ratio (n =216796)", 
                                    "Cystatin_C (n =217425)", "LDL_direct (n =217259)", "Cholesterol (n =217525)","eGFR (n =217214)", "HDL_cholesterol (n =217548)",
                                    "Arm_fat_free_mass_right (n =217480)", "Arm_fat_free_mass_left (n =217435)", "Apolipoprotein_A (n =217312)", "Trunk_fat_mass (n =217395)", "Haemoglobin_concentration (n =212318)"
)
phenotype_choose_v2           <- paste0(gsub(" \\(.*$", "", phenotype_choose), collapse = "|")
general_traits_data_select    <- general_traits_data_select[grepl(phenotype_choose_v2, general_traits_data_select$phenotype), ]
p <- ggplot(result.plot, aes(x=reorder(phenotype, estimate), y= estimate)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_point(size = 4, alpha =1) + coord_flip() +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    strip.background = element_rect( fill = "white"),
    strip.text = element_text(size = 8),
    legend.position = "bottom"

  ) 

########################################################################
#Figure 2D
#Cox model
########################################################################
cox_model     <- coxph(Surv(disease_time, 
                                disease_event) ~ MHS_H + body_mass_index_BMI + genetic_sex_male + age + age_squared + age_x_sex_male + age_squared_x_sex_male + initial_50k_release + pc_1 + pc_2 + pc_3 + pc_4 + pc_5 + pc_6 + pc_7 + pc_8 + pc_9 + pc_10, 
                           data = data_use, 
                           singular.ok=T) 
    




