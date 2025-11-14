library(tidyverse)
library(viridisLite)

results <- list()
rev_results <- list()



for(i in 1:25){
  results[[i]]<- readRDS(paste0("results_phq",i,".rds"))
  
  rev_results[[i]]<- readRDS(paste0("results_phq_rev",i,".rds"))
  
}


methods <- names(rev_results[[1]])

get_runtime <- function(i){
  (sapply(1:length(results[[i]]), 
          function(x) results[[i]][[x]]$runtime )+
     sapply(1:length(rev_results[[i]]), 
            function(x) results[[i]][[x]]$runtime ))/2
}

get_badprops <- function(i){
  (sapply(1:length(results[[i]]), 
          function(x) results[[i]][[x]]$bad_props )+
     sapply(1:length(rev_results[[i]]), 
            function(x) results[[i]][[x]]$bad_props ))/2
}

get_values <- function(i){
  sapply(1:length(results[[i]]), 
         function(x) results[[i]][[x]]$value )
}

get_values_rev <- function(i){
  sapply(1:length(rev_results[[i]]), 
         function(x) rev_results[[i]][[x]]$value )
}

get_rules <- function(i){
  sapply(1:length(results[[i]]), 
         function(x) results[[i]][[x]]$rule )
}

get_rules_rev <- function(i){
  sapply(1:length(rev_results[[i]]), 
         function(x) rev_results[[i]][[x]]$rule )
}




values <- data.frame(-1*t(sapply(1:25,get_values)))

values_rev <- data.frame(-1*t(sapply(1:25,get_values_rev)))

rules <- data.frame(t(sapply(1:25,get_rules)))

rules_rev <- data.frame(t(sapply(1:25,get_rules_rev)))

runtimes <- data.frame(t(sapply(1:25,get_runtime)))

selected <- data.frame(
  OAL = (apply(sapply(1:25,function(x) results[[x]]$oal$selected),1,mean)+
           apply(sapply(1:25,function(x) rev_results[[x]]$oal$selected),1,mean))/2,
  GLiDeR = (apply(sapply(1:25,function(x) results[[x]]$glider$selected),1,mean)+
              apply(sapply(1:25,function(x) rev_results[[x]]$glider$selected),1,mean))/2,
  DiPS = (apply(sapply(1:25,function(x) results[[x]]$dip$selected),1,mean)+
            apply(sapply(1:25,function(x) rev_results[[x]]$dip$selected),1,mean))/2
  
)

names(values) <- methods
names(values_rev) <- methods
names(rules) <- methods
names(rules_rev) <- methods
names(runtimes) <- methods




covars = c("AgeAtIndex",
                      #information on patient sex (male or female).
                      "gend",
                      #Health records data on race and ethnicity
                      "race",
                      #insurance type (commercial, Medicaid, Medicare, or private)
                      "ins",
                      #neighborhood educational attainment (less than 25\% college degrees)
                      "educ",
                      #income (median  lower than 40,000 USD)
                      "income",
                      #level of poverty (20\% of households below federal poverty level)
                      "pov",
                      #urban or rural area (1 to 6, with 1 the most urban and 6 the most rural)
                      "urban_rural",
                      #Charlson score at initiation
                      "Charlson",
                      #tobacco use in the year prior to init
                      "tobacco_priorYr",
                      #BMI
                      "baseline_bmi",
                      #past year anxiety
                      "anx_dx_priorYr",
                      #past year non-anxiety mental health or substance use disorder
                      "non_anx_dx_priorYr",
                      #number of suicide attempts 6 months prior
                      "num_SH_prior6m",
                      #number of psych hospitalizations 6 months prior
                      "num_MHIP_prior6m",
                      #number different antidepressants taken 5 years prior
                      "num_AD_Prior5yr",
                      #received psychotherapy in the 5 years prior
                      "psychtherapy_prior5yr",
                      #number of PHQ measurements recorded year prior
                      "PriorPhqCount",
                      #PHQ recorded closest to treatment initiation
                      "baseline_phq",
                      #calendar year of treatment initiation
                      "epiYear")

string <- c("G:/CTRHS/IMATS/Data/MHRN data/Analytic/Imputations/Final datasets/ADstartImp1.rdata")
load(string)

full_dat <- ADstartImp

full_dat <- full_dat %>% mutate(non_anx_dx_priorYr =
                                  aud_dx_priorYr +
                                  asd_dx_priorYr +
                                  ptsd_dx_priorYr +
                                  ocd_dx_priorYr +
                                  oud_dx_priorYr +
                                  per_dx_priorYr +
                                  sed_dx_priorYr > 0,
                                #any psychthreapy in prior 5 years
                                psychtherapy_prior5yr = numMH_prior5yr > 0,
                                #basline phq score
                                baseline_phq = IndexTrtPHQ8_score +
                                  as.numeric(phq9_cat),
                                baseline_bmi = WeightInitMed*0.453592/(HeightInitMed*0.0254)^2)


X <- full_dat[,covars]

f <- paste0(c("~",names(X)),
            c("",rep(" + ",length(names(X))-1),""),
            collapse = "")

#creating indicator variables for categorical covars
X <- model.matrix(as.formula(f), data = X)[,-1]

write.csv(colnames(X),file = "variable_names.csv")

variable_names <- read.csv("variable_names.csv")

selected$Variable <- variable_names$names


values <- pivot_longer(values,1:5, names_to = "Method", values_to = "Value")
values_rev <- pivot_longer(values_rev,1:5, names_to = "Method", values_to = "Value")
rules <- pivot_longer(rules,1:5, names_to = "Method", values_to = "Threshold")
rules_rev <- pivot_longer(rules_rev,1:5, names_to = "Method", values_to = "Threshold")
runtimes <- pivot_longer(runtimes,1:5, names_to = "Method", values_to = "Runtime (s)")
selected <- pivot_longer(selected,1:3, names_to = "Method", values_to = "Proportion")


values$Method <- recode(values$Method,"no_selection" = "A",
                        "oal" = "OAL",
                        "glider" = "GLiDeR",
                        "dip" = "DiPS",
                        "casual_ball" = "CB",
                        "hd_balancing" = "HDCB")

values$Method <- factor(values$Method, levels = c("A","OAL","GLiDeR","DiPS",
                                                  "HDCB","CB"))


values_rev$Method <- recode(values_rev$Method,"no_selection" = "A",
                            "oal" = "OAL",
                            "glider" = "GLiDeR",
                            "dip" = "DiPS",
                            "casual_ball" = "CB",
                            "hd_balancing" = "HDCB")

values_rev$Method <- factor(values_rev$Method, levels = c("A","OAL","GLiDeR","DiPS",
                                                          "HDCB","CB"))

rules$Method <- recode(rules$Method,"no_selection" = "A",
                       "oal" = "OAL",
                       "glider" = "GLiDeR",
                       "dip" = "DiPS",
                       "casual_ball" = "CB",
                       "hd_balancing" = "HDCB")

rules$Method <- factor(rules$Method, levels = c("A","OAL","GLiDeR","DiPS",
                                                "HDCB","CB"))


rules_rev$Method <- recode(rules_rev$Method,"no_selection" = "A",
                           "oal" = "OAL",
                           "glider" = "GLiDeR",
                           "dip" = "DiPS",
                           "casual_ball" = "CB",
                           "hd_balancing" = "HDCB")

rules_rev$Method <- factor(rules_rev$Method, levels = c("A","OAL","GLiDeR","DiPS",
                                                        "HDCB","CB"))

selected$Method <- recode(selected$Method,"no_selection" = "A",
                          "oal" = "OAL",
                          "glider" = "GLiDeR",
                          "dip" = "DiPS",
                          "casual_ball" = "CB",
                          "hd_balancing" = "HDCB")

selected$Method <- factor(selected$Method, levels = c("A","OAL","GLiDeR","DiPS",
                                                      "HDCB","CB"))




pdf("rule_val_no_CB.pdf", height = 10, width = 12)
print(ggplot(values, aes(y = Value, fill = Method, x = Method)) + geom_boxplot() +
        scale_fill_viridis_d(option = "H") + labs(y = "Estimated Average\n6 Month Change in PHQ9") +
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
              text=element_text(size =20)))
dev.off()

pdf("rule_val_rev_no_CB.pdf", height = 10, width = 12)
print(ggplot(values_rev, aes(y = Value, fill = Method, x = Method)) + geom_boxplot() +
        scale_fill_viridis_d(option = "H") + labs(y = "Estimated Average\n6 Month Change in PHQ9") +
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
              text=element_text(size =20)))
dev.off()

pdf("rule_threshold_no_CB.pdf", height = 10, width = 12)
print(ggplot(rules, aes(x = Threshold, fill = Method)) + geom_histogram() +
        facet_grid(rows = ggplot2::vars(Method)) + scale_fill_viridis_d(option = "H") +
        labs(y = "Count", x = "Treat if Baseline PhQ9\nLess Than Threshold") +
        theme(legend.position = "none",
              text=element_text(size =20)))
dev.off()

pdf("rule_threshold_rev_no_CB.pdf", height = 10, width = 12)
print(ggplot(rules_rev, aes(x = Threshold, fill = Method)) + geom_histogram() +
        facet_grid(rows = ggplot2::vars(Method)) + scale_fill_viridis_d(option = "H") +
        labs(y = "Count", x = "Treat if Baseline PhQ9\nGreater Than Threshold") +
        theme(legend.position = "none",
              text=element_text(size =20)))
dev.off()


pdf("selected_no_CB.pdf", height = 10, width = 12)
print(ggplot(selected, aes(x = Proportion, color = Method, y = Variable, shape = Method)) + geom_point(size =4) +
        # facet_grid(cols = ggplot2::vars(Method)) 
        #scale_colour_viridis_d(option = "H") +
        labs(y = "Proportion Selected") +
        theme(
          text=element_text(size =20),
          panel.spacing.x = unit(2, "lines")))
dev.off()
