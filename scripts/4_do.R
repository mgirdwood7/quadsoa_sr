# Do

quadsoa_cleaned <- read.csv("data/processed/Strength OA Review Cleaned.csv")

# Takagi paper presents two populations  - with OA and at risk of OA. These subgroups will be used in a later sensitivity analysis. 
# For now need to combine together for the main analysis - using a 'pre-meta-analysis
takagi <- quadsoa_cleaned %>%
  filter(author == "Takagi") %>% 
  group_by(author, sex, analysis_group) %>% # group by each publication and subgroups of interest in this review, for later 'pre-metaanalyses'
  nest(data = -c(studyname:year, outcome:muscle_test, reference_group, explanatory_group)) %>%# nest data, remove data that is consistent across subgroups
  mutate(rma = map(data, ~tidy(rma(yi, vi, data = .x)) %>% # for each group, perform a rma
                     select(estimate, std.error) %>% # take only the estimated log rr and log se
                     rename(yi = estimate, # rename for consistency with final dataframe
                            sei = std.error) %>%
                     mutate(vi = sei^2)
  )) %>%
  select(-data) %>% # remove the original data as no longer needed
  unnest(cols = c(rma)) # un-nest back to same format

# now need to join the 'pre-meta' data back to the original dataset
quadsoa_analysis <- quadsoa_cleaned %>% 
  filter(author != "Takagi") %>% # take all data not split before
  bind_rows(., takagi) # bind the pre-meta data 

# Order and arrange
quadsoa_analysis <- quadsoa_analysis %>%
  mutate(sex = factor(sex, levels = c("mixed", "Men", "Women", "Overall"))) %>%
  arrange(outcome, analysis_group, muscle, sex) %>% # order and arrange data frame
  mutate(studyname = case_when( # addition of asterisk for 3 studies with notes for plots
    studyname %in% c("Davis 2019", "Chi 2019", "Woolard 2011") ~ paste0(studyname, "*"),
    studyname == "Takagi 2018" ~ paste0(studyname, "**"),
    TRUE ~ studyname
  ))

# 1. Primary Analysis - Grouped by compartment and muscle tested
# Pools all studies together regardless of outcome type used
# sex specific meta-analyses (groups are: mixed, male, female)

primary_meta_sex <- quadsoa_analysis %>%
  filter(analysis_group %in% c("medial tf", "whole tf", "lateral pf", "whole pf") &
         studyname != "Culvenor 2019" & # Exclude Culvenor as overlaps with Segal 2010
         !(studyname == "Dannhauer 2014" & muscle == "quad")) %>% # exclude Dannhauer Quads data as overlaps with DellaIsola
  mutate(analysis_group = case_when(
    str_detect(analysis_group, "tf") ~ "tf", # medial tf and whole tf combined as tf
    str_detect(analysis_group, "pf") ~ "pf" # lateral pf and whole pf combined as pf
    )) %>%
  arrange(sex) %>%
  group_by(sex, analysis_group, muscle) %>%
  nest() %>%
  mutate(rma = map(data, ~rma(yi, vi, data = .x)),
         output = map(rma, ~tidy(.x, conf.int = TRUE, exponentiate = TRUE)),
         model = map(rma, glance), # get full statistics
         predict = map(rma, ~predict(.x, transf = exp) %>% as.data.frame %>% select(starts_with("pi")))) %>% # prediction interval
  unnest(cols = c(output, model, predict)) # unnest

primary_meta_overall <- quadsoa_analysis %>%
  filter(analysis_group %in% c("medial tf", "whole tf", "lateral pf", "whole pf") &
           studyname != "Culvenor 2019" & # Exclude Culvenor as overlaps with Segal 2010
           !(studyname == "Dannhauer 2014" & muscle == "quad")) %>% # Dannhauer and Culvenor removed due to data duplication
  mutate(analysis_group = case_when(
    str_detect(analysis_group, "tf") ~ "tf",
    str_detect(analysis_group, "pf") ~ "pf"
  )) %>%
  arrange(sex) %>%
  group_by(analysis_group, muscle) %>%
  nest() %>%
  mutate(rma = map(data, ~rma(yi, vi, data = .x)),
         output = map(rma, ~tidy(.x, conf.int = TRUE, exponentiate = TRUE)),
         model = map(rma, glance), # get full statistics
         predict = map(rma, ~predict(.x, transf = exp) %>% as.data.frame %>% select(starts_with("pi"))),
         sex = "Overall") %>% # prediction interval
  unnest(cols = c(output, model, predict)) # unnest

primary_meta_combined <- bind_rows(primary_meta_sex, primary_meta_overall) %>%
  mutate(sex = factor(sex, levels = c("mixed", "Women", "Men", "Overall"))) %>%
  select(analysis_group:type, nobs, everything()) %>%
  arrange(analysis_group, muscle, sex) %>% # order and arrange data frame
  ungroup 

## Write results to file
write_csv(primary_meta_combined %>% select(-c(data,rma)), "output/tables/Primary Meta Analysis Full Results.csv")


## Forest Plots for Primary Analysis
png("output/plots/primary tf quad.png", height = 1250, pointsize = 25, width = 800)
forest_plotr(data = primary_meta_combined, analysis_group = "tf", muscle ="quad", extraspace = TRUE, sumtext = TRUE)
mtext(text = "Tibiofemoral OA Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Low Knee Extensor Strength", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/primary tf hs.png", 800, pointsize = 25, width = 800)
forest_plotr(data = primary_meta_combined, analysis_group = "tf", muscle ="hs", supress = TRUE)
mtext(text = "Tibiofemoral OA Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Low Knee Flexor Strength", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/primary pf quad.png", height = 800, pointsize = 25, width = 800)
forest_plotr(data = primary_meta_combined, analysis_group = "pf", muscle ="quad", supress = TRUE)
mtext(text = "Patellofemoral OA Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Low Knee Extensor Strength", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

## Funnel plot for quad <-> tf 
png("output/plots/funnel.png", pointsize = 25, height = 700, width = 700)
par(family = "Karla")
primary_meta_combined %>%
  filter(analysis_group == "tf",
         muscle == "quad", 
         sex == "Overall") %>%
  pluck('rma') %>% pluck(1) %>%
  funnel(xlab = "Risk Ratio",
         main = "Funnel Plot \nTibiofemoral OA and \n Low Knee Extensor Strength")
dev.off()

# Eggers Test
primary_meta_combined %>%
  filter(analysis_group == "tf",
         muscle == "quad", 
         sex == "Overall") %>%
  pluck('rma') %>% pluck(1) %>% regtest()


# 2. Secondary Analyses - Split further by outcome type (e.g. quantitative JSN/radiographic OA etc)
# sex specific meta-analyses (groups are: mixed, male, female)
secondary_meta_sex <- quadsoa_analysis %>%
  group_by(outcome, sex, analysis_group, muscle) %>% # group by each publication and subgroups of interest in this review, for later 'pre-meta-analyses'
  nest(data = -c(outcome:muscle)) %>%# nest data, remove data that is consistent across subgroups
  mutate(rma = map(data, ~rma(yi, vi, data = .x)), # run meta-analysis
         output = map(rma, ~tidy(.x, conf.int = TRUE, exponentiate = TRUE)), # tidy model outputs, including cis and exponentiate for interpretability 
         model = map(rma, glance), # get full statistics
         predict = map(rma, ~predict(.x, transf = exp) %>% as.data.frame %>% select(starts_with("pi")))) %>% # prediction interval
  unnest(cols = c(output, model, predict)) # unnest

# meta-analyses split by outcome, compartment and muscle tested
secondary_meta_overall <- quadsoa_analysis %>%
  group_by(outcome, analysis_group, muscle) %>% # group by each publication and subgroups of interest in this review, for later 'pre-meta-analyses'
  nest(data = -c(outcome, analysis_group:muscle)) %>%# nest data, remove data that is consistent across subgroups
  mutate(rma = map(data, ~rma(yi, vi, data = .x)), # run meta-analysis
         output = map(rma, ~tidy(.x, conf.int = TRUE, exponentiate = TRUE)), # tidy model outputs, including cis and exponentiate for interpretability 
         model = map(rma, glance), # get full statistics
         predict = map(rma, ~predict(.x, transf = exp) %>% as.data.frame %>% select(starts_with("pi")))) %>% # prediction interval
  unnest(cols = c(output, model, predict)) %>% # unnest
  mutate(sex = "Overall")

# combine results from sex specific and overall meta-analyses together
secondary_meta_combined <- bind_rows(secondary_meta_sex, secondary_meta_overall) %>%
  mutate(sex = factor(sex, levels = c("mixed", "Women", "Men", "Overall"))) %>%
  select(outcome:type, nobs, everything()) %>%
  arrange(outcome, analysis_group, muscle, sex) %>% # order and arrange data frame
  ungroup 

# Write results to file
write_csv(secondary_meta_combined %>% select(-c(data,rma)), "output/tables/Secondary Meta Analysis by Outcome Full Results.csv")

# Individual plots with titles added also

png("output/plots/worsening wholetf quad.png", height = 1110, pointsize = 25, width = 800)
forest_plotr(outcome = "worsening JSN/OA", analysis_group = "whole tf", muscle ="quad", extraspace = TRUE)
mtext(text = "(A) Radiographic OA Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Whole Tibiofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/worsening whole pf quad.png", height = 700, pointsize = 25, width = 800)
forest_plotr(outcome = "worsening JSN/OA", analysis_group = "whole pf", muscle ="quad", supress = TRUE)
mtext(text = "(B) Radiographic OA Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Whole Patellofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/quantJSN medialtf quad.png", height = 700, pointsize = 25, width = 800)
forest_plotr(outcome = "quantitative JSN", analysis_group = "medial tf", muscle ="quad", supress = TRUE)
mtext(text = "Quantitative Joint Space Narrowing", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Medial Tibiofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/quant medialtf quad.png", height = 700, pointsize = 25, width = 800)
forest_plotr(outcome = "quantitative cartilage progression", analysis_group = "medial tf", muscle ="quad", supress = TRUE)
mtext(text = "Quantitative Cartilage Thinning", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Medial Tibiofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/semiquant medialtf quad.png", height = 800, pointsize = 25, width = 800)
forest_plotr(outcome = "semi quantitative cartilage progression", analysis_group = "medial tf", muscle ="quad", supress = TRUE)
mtext(text = "(A) Cartilage Lesion Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Medial Tibiofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/semiquant lateraltf quad.png", height = 750, pointsize = 25, width = 800)
forest_plotr(outcome = "semi quantitative cartilage progression", analysis_group = "lateral tf", muscle ="quad", supress = TRUE)
mtext(text = "(B) Cartilage Lesion Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Lateral Tibiofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/semiquant medial pf quad.png", height = 700, pointsize = 25, width = 800)
forest_plotr(outcome = "semi quantitative cartilage progression", analysis_group = "medial pf", muscle ="quad", supress = TRUE)
mtext(text = "(C) Cartilage Lesion Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Medial Patellofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/semiquant lateral pf quad.png", height = 750, pointsize = 25, width = 800)
forest_plotr(outcome = "semi quantitative cartilage progression", analysis_group = "lateral pf", muscle ="quad", supress = TRUE)
mtext(text = "(D) Cartilage Lesion Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Lateral Patellofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/worsening wholetf hs.png", height = 800, pointsize = 25, width = 800)
forest_plotr(outcome = "worsening JSN/OA", analysis_group = "whole tf", muscle ="hs", supress = TRUE)
mtext(text = "Radiographic OA Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Whole Tibiofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

# Risk of Bias Plot
## Adapted from code from dmetar package
rob <- read_csv("data/raw/rob.csv")

rob <- gather(rob,key = "Domain", value = "Value",-Study) %>%
  mutate(val2=if_else(Value == "Low", "+", if_else(Value == "High", "-", "?"))) %>%
  mutate(Domain=factor(Domain,levels=unique(Domain))) %>%
  mutate(Value = factor(Value, levels = c("Low", "Unclear", "High")))

#Plot with traffic light table
rob.table <- ggplot(data = rob, aes(y = Study, x = Domain)) +
  geom_tile(color="black", fill="white", size = 0.8) +
  geom_point(aes(color=as.factor(val2)), size=8) +
  geom_text(aes(label = val2), size = 8) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits=rev(levels(as.factor(rob$Study)))) +
  scale_color_manual(values = c("-" = "#BF0000",
                                "+" = "#02C100",
                                "?" = "#E2DF07")) +
  theme_minimal() +
  coord_equal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 13, color = "black", angle = 60, hjust=0),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave("output/plots/rob.png", plot = rob.table, device = "png", width = 1500, height = 2700, units = "px", dpi = 300)
ggsave("output/plots/rob.tiff", plot = rob.table, device = "tiff", width = 1500, height = 2700, units = "px", dpi = 300)

# Grade Summary Plot
# Combining grade assessments with summary meta-analyses for better presentation of results
grade_ass <- read_csv("data/raw/grade_assessments.csv")

# Get m-a data to use in summary grade plot
grade_table <- secondary_meta_combined %>%
 
  
grade_table <- bind_rows(primary_meta_combined %>% mutate(order = 1), secondary_meta_combined %>% mutate(order = 2)) %>%
  filter(sex == "Overall",
         nobs > 1) %>% # only those with >1 ob performed grade on
  select(outcome, analysis_group, muscle, estimate, conf.low, conf.high, order) %>%
  mutate(log_estimate = log(estimate), # for metafor forest plot need log val
         log_conf.low = log(conf.low),
         log_conf.high = log(conf.high)) %>% 
  left_join(., grade_ass, by = c("outcome", "analysis_group", "muscle")) %>%
  filter(!is.na(GRADE)) %>%
  mutate(outcome_name = factor(outcome_name, levels = c("Cartilage lesion worsening", "Cartilage thinning", "Joint space narrowing", "Radiographic OA worsening", "Any structural worsening"))) %>%
  arrange(desc(order), muscle, outcome_name)

png("output/plots/grade summary.png", height = 1200, pointsize = 25, width = 1400)
forest(x = grade_table$log_estimate,
       ci.lb = grade_table$log_conf.low,
       ci.ub = grade_table$log_conf.high,
       xlim = c(-15,4),
       ylim = c(0,22),
       cex = 0.8,
       slab = grade_table$outcome_name,
       alim = c(log(1/10),log(10)), # limit of data clipping
       at = c(log(1/5), log(1/2), 0, log(2), log(5)), # axis limit (different to above)
       atransf = exp,
       pch = 20, # use small point so hidden behind polygon
       psize = 1,
       rows = c(1,3,4,5,6,7,8,9,10, 14, 16, 17),
       ilab = cbind(grade_table$location_name,
                    grade_table$`No. studies`, # add matrix of data from grade assessments
                    grade_table$`No. of observations`,
                    grade_table$`Risk of bias`,
                    grade_table$`Inconsistency`,
                    grade_table$Indirectness,
                    grade_table$Imprecision,
                    grade_table$`Publication bias`,
                    grade_table$GRADE),
       ilab.xpos = c(-11, -8.5, -7.75, -7, -6.25, -5.5, -4.75, -4, -3),
       ilab.pos = 4,
       fonts = "Karla",
       xlab = "",
       header = c("Outcome", "Risk Ratio [95%CI]"),
       efac = c(0,1))

# Add diamongs for summary estimates
addpoly(x = grade_table$log_estimate,
        ci.lb = grade_table$log_conf.low,
        ci.ub = grade_table$log_conf.high,
        atransf = exp,
        cex = 0.8,
        annotate = FALSE,
        rows = c(1,3,4,5,6,7,8,9,10, 14, 16, 17))

# expand limits so text not clipped
par(xpd=NA)

# Add risk-factor annotations
text(-15, c(2, 11, 15, 18), pos = 4, c("Low Knee Flexor Strength", "Low Knee Extensor Strength", "Low Knee Flexor Strength", "Low Knee Extensor Strength"), font = 2, cex = 0.8)


# Text for primary/secondary
text(-15, c(12), pos = 4, bquote(paste(underline("Secondary Analysis"), " - Outcome specific")), font = 2, cex = 1)
text(-15, c(19), pos = 4, bquote(paste(underline("Primary Analysis"),  " - All outcomes pooled")), font = 2, cex = 1)


# text for axis descriptors
text(log(c(1/5, 5)), -4, c("Lower strength =\ndeccreased risk","Lower strength =\nincreased risk"), pos=c(3,3), cex = 0.8)

# Add annotations for grade matrix
text(x = c(-8.5, -7.75, -7, -6.25, -5.5, -4.75, -4),
     y = 21,
     pos = 4,
     srt = 45, # angle text to help readability
     cex = 0.8,
     labels = c("No. studies", "No. observations", "Risk of bias", "Inconsistency", "Indirectness", "Imprecision", "Publication bias"))
text(x = c(-3),
     y = 21,
     pos = 4,
     cex = 0.8,
     font = 2,
     labels = c("Level of\nEvidence"))
text(x = c(-11),
     y = 21,
     pos = 4,
     cex = 0.8,
     font = 2,
     labels = c("Compartment"))
text(c(-5.5),
     24,
     pos=3,
     c("GRADE"))

dev.off()

# Sensitivity analysis around impact of population subgrouop
# For worsening JSN/OA grade - impact of at risk of OA vs OA group
# Only possible for medial TF joint

population_sa_data <- quadsoa_cleaned %>%
  filter(analysis_group %in% c("medial tf", "whole tf"),
         !studyname %in% c("Dannhauer 2014", "Culvenor 2019"),
         muscle == "quad") %>% # Dannhauer and Culvenor removed due to data duplication
  mutate(analysis_group = "tf")

# Not considering sex in this analysis - need to combine men + women groups together with a MA before analysing.
sa_combine <- population_sa_data %>%
  filter(sex != "mixed") %>% # select 
  group_by(author, population_subgroup) %>% # group by each publication and subgroups of interest in this review, for later 'pre-metaanalyses'
  nest(data = -c(studyname, author, year, outcome, population_subgroup, analysis_group:explanatory_group)) %>%
  mutate(rma = map(data, ~tidy(rma(yi, vi, data = .x)) %>% # for each group, perform a rma
                     select(estimate, std.error) %>% # take only the estimated log rr and log se
                     rename(yi = estimate, # rename for consistency with final dataframe
                            sei = std.error) %>%
                     mutate(vi = sei^2)
  )) %>%
  select(-data) %>% # remove the original data as no longer needed
  unnest(cols = c(rma)) %>% # un-nest back to same format
  mutate(sex = "mixed")

# now need to join the 'pre-meta' data back to the original dataset
population_sa_data <- population_sa_data %>% 
  filter(sex == "mixed") %>% # take all data not split before
  bind_rows(., sa_combine) # bind the pre-meta data 

population_sa <- population_sa_data %>%
  group_by(population_subgroup) %>% # group by each publication and subgroups of interest in this review, for later 'pre-meta-analyses'
  nest(data = -c(population_subgroup, analysis_group:muscle)) %>%# nest data, remove data that is consistent across subgroups
  mutate(rma = map(data, ~rma(yi, vi, data = .x)), # run meta-analysis
         output = map(rma, ~tidy(.x, conf.int = TRUE, exponentiate = TRUE)), # tidy model outputs, including cis and exponentiate for interpretability 
         model = map(rma, glance), # get full statistics
         predict = map(rma, ~predict(.x, transf = exp) %>% as.data.frame %>% select(starts_with("pi")))) %>% # prediction interval
  unnest(cols = c(output, model, predict)) # unnest

population_sa %>% 
  select(-c(data, rma)) %>%
  write_csv("output/tables/Population Subgroup Sensitivity Analysis.csv")

sa_plot <- population_sa_data %>%
  mutate(population_subgroup = factor(population_subgroup, levels = c("risk of OA", "OA", "OA or risk of OA"))) %>%
  arrange(desc(population_subgroup))


png("output/plots/sensitivity population.png", height = 1000, pointsize = 25, width = 800)
forest(sa_plot$yi,
       sa_plot$vi,
       slab = sa_plot$studyname, 
       xlim = c(-5,5),
       ylim = c(-1,22),
       alim = c(log(1/10),log(10)), # limit of data clipping
       at = c(log(1/5), log(1/2), 0, log(2), log(5)), # axis limit (different to above)
       cex = 0.8,
       rows = c(1, 4.5:11.5, 15:18),
       atransf = exp,
       fonts = "Karla",
       xlab = "",
       header = c("Study", "Risk Ratio [95%CI]"), 
       efac = c(0,1)) # first element is the cis, second is the arrow. This eliminates the tick on the ci
# add subgroup labels
text(-5, c(2, 12.5, 19), pos = 4, c("Combined", "OA", "Risk of OA"), font = 2)
par(xpd=NA, font = 1) # remove plot clipping limits
text(log(c(1/5, 5)), -6 , c("Lower strength =\ndecreased risk","Lower strength =\nincreased risk"), pos=c(3,3), cex = 0.8) 
par(font = 2)
addpoly(population_sa %>% filter(population_subgroup == "OA or risk of OA") %>% pluck(5) %>% pluck(1)
        , atransf = exp, cex = 0.8, row = 0, mlab = "    RE Subgroup Model")
addpoly(population_sa %>% filter(population_subgroup == "OA") %>% pluck(5) %>% pluck(1)
        , atransf = exp, cex = 0.8, row = 3.5, mlab = "    RE Subgroup Model")
addpoly(population_sa %>% filter(population_subgroup == "risk of OA") %>% pluck(5) %>% pluck(1)
        , atransf = exp, cex = 0.8, row = 14, mlab = "    RE Subgroup Model")
par(font = 1)
mtext(text = "Sensitivity Analysis - Population Subgroup", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "OA Worsening - in presence of low knee extensor strength", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()


# Leave one out analysis
## First run the influence analysis with metafor package - gives all indepth statistics and measures
## Join these to original data
## Meta package provides nicer influence plot, but use the data from the metafor output in the plot.

### For Primary Analyses

primary_meta_loo <- quadsoa_analysis %>%
  filter(analysis_group %in% c("medial tf", "whole tf", "lateral pf", "whole pf") &
           studyname != "Culvenor 2019" & # Exclude Culvenor as overlaps with Segal 2010
           !(studyname == "Dannhauer 2014" & muscle == "quad")) %>% # Dannhauer and Culvenor removed due to data duplication
  mutate(analysis_group = case_when(
    str_detect(analysis_group, "tf") ~ "tf",
    str_detect(analysis_group, "pf") ~ "pf"
  )) %>%
  group_by(analysis_group, muscle) %>% # group by each publication and subgroups of interest in this review, for later 'pre-meta-analyses'
  mutate(k = length(studyname)) %>%
  nest(data = -c(analysis_group:muscle, k)) %>%# nest data, remove data that is consistent across subgroups
  ungroup() %>%
  mutate(rma = map(data, ~rma(yi, vi, data = .x)), # runs metafor meta-analysis
         metaforinf = map_if(rma,  # ruin a influence analysis from metafor package if k>1
                             k > 1, 
                             ~influence(.x) %>% pluck(1) %>% data.frame, # take the output
                             .else = ~data.frame(inf = NA))) %>% # for k <=1 return df with NAs
  mutate(new = map2(data, metaforinf, ~bind_cols(.x, .y))) %>% # combine influence data with original data
  mutate(new = map(new, ~.x %>%
                     mutate(studyname = case_when(
                       sex != "mixed" ~ paste(studyname, sex),
                       TRUE ~ studyname)))) %>%
  select(-c(data, rma, metaforinf)) %>% # remove all nested columns now not needed
  # now run meta with metagen (inorder to get nice influence plots)
  mutate(rma = map(new, ~metagen(TE = yi, seTE = sei, studlab = studyname, data = .x, sm = "RR")),
         influence = map_if(rma, k > 1, ~metainf(.x, pooled = "random"))) # run influence analysis with meta package

  
           

### For Secondary Analyses
secondary__meta_loo <- quadsoa_analysis %>%
  group_by(outcome, analysis_group, muscle) %>% # group by each publication and subgroups of interest in this review, for later 'pre-meta-analyses'
  mutate(k = length(studyname)) %>%
  nest(data = -c(outcome, analysis_group:muscle, k)) %>%# nest data, remove data that is consistent across subgroups
  ungroup() %>%
  mutate(rma = map(data, ~rma(yi, vi, data = .x)), # runs metafor meta-analysis
         metaforinf = map_if(rma,  # ruin a influence analysis from metafor package if k>1
                             k > 1, 
                             ~influence(.x) %>% pluck(1) %>% data.frame, # take the output
                             .else = ~data.frame(inf = NA))) %>% # for k <=1 return df with NAs
  mutate(new = map2(data, metaforinf, ~bind_cols(.x, .y))) %>% # combine influence data with original data
  mutate(new = map(new, ~.x %>%
                     mutate(studyname = case_when(
                       sex != "mixed" ~ paste(studyname, sex),
                       TRUE ~ studyname)))) %>%
  select(-c(data, rma, metaforinf)) %>% # remove all nested columns now not needed
  # now run meta with metagen (inorder to get nice influence plots)
  mutate(rma = map(new, ~metagen(TE = yi, seTE = sei, studlab = studyname, data = .x, sm = "RR")),
         influence = map_if(rma, k > 1, ~metainf(.x, pooled = "random"))) # run influence analysis with meta package

# Leave-one-out analysis plot for each (change position of pluck)
pluck(secondary_meta_loo$influence) %>% pluck(2) %>%
forest.meta(., 
            lower.equi = exp(last(.$lower)), 
            upper.equi = exp(last(.$upper)), 
            col.equi = "green", 
            fill.equi = "green",
            leftcols = c("studlab", "inf"),
            leftlabs = c("Analysis", "Sig Influence"),
            text.random = "Original random effects model\n(All studies included)",
            rightcols = c("effect.ci", "I2"),
            backtransf = TRUE,
            fontfamily = "Karla",
            xlim = c(0.33, 3))

