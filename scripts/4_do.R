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

# Primary Analyses
# sex specific meta-analyses (groups are: mixed, male, female)
quadsoa_meta_sex <- quadsoa_analysis %>%
  group_by(outcome, sex, analysis_group, muscle) %>% # group by each publication and subgroups of interest in this review, for later 'pre-meta-analyses'
  nest(data = -c(outcome:muscle)) %>%# nest data, remove data that is consistent across subgroups
  mutate(rma = map(data, ~rma(yi, vi, data = .x)), # run meta-analysis
         output = map(rma, ~tidy(.x, conf.int = TRUE, exponentiate = TRUE)), # tidy model outputs, including cis and exponentiate for interpretability 
         model = map(rma, glance), # get full statistics
         predict = map(rma, ~predict(.x, transf = exp) %>% as.data.frame %>% select(starts_with("pi")))) %>% # prediction interval
  unnest(cols = c(output, model, predict)) # unnest

# meta-analyses split by outcome, compartment and muscle tested
quadsoa_meta_overall <- quadsoa_analysis %>%
  group_by(outcome, analysis_group, muscle) %>% # group by each publication and subgroups of interest in this review, for later 'pre-meta-analyses'
  nest(data = -c(outcome, analysis_group:muscle)) %>%# nest data, remove data that is consistent across subgroups
  mutate(rma = map(data, ~rma(yi, vi, data = .x)), # run meta-analysis
         output = map(rma, ~tidy(.x, conf.int = TRUE, exponentiate = TRUE)), # tidy model outputs, including cis and exponentiate for interpretability 
         model = map(rma, glance), # get full statistics
         predict = map(rma, ~predict(.x, transf = exp) %>% as.data.frame %>% select(starts_with("pi")))) %>% # prediction interval
  unnest(cols = c(output, model, predict)) %>% # unnest
  mutate(sex = "Overall")

# combine results from sex specific and overall meta-analyses together
quadsoa_meta_combined <- bind_rows(quadsoa_meta_sex, quadsoa_meta_overall) %>%
  mutate(sex = factor(sex, levels = c("mixed", "Women", "Men", "Overall"))) %>%
  select(outcome:type, nobs, everything()) %>%
  arrange(outcome, analysis_group, muscle, sex) %>% # order and arrange data frame
  ungroup 

# Write results to file
write_csv(quadsoa_meta_combined %>% select(-c(data,rma)), "output/tables/Meta Analysis Full Results.csv")


# Sensitivity analysis around impact of population subgrouop
# For worsening JSN/OA grade - impact of at risk of OA vs OA group
# Only possible for medial TF joint
population_sa <- quadsoa_cleaned %>%
  filter(outcome == "worsening JSN/OA",
         analysis_group == "whole tf",
         muscle == "quad") %>%
  group_by(outcome, population_subgroup) %>% # group by each publication and subgroups of interest in this review, for later 'pre-meta-analyses'
  nest(data = -c(population_subgroup, outcome, analysis_group:muscle)) %>%# nest data, remove data that is consistent across subgroups
  mutate(rma = map(data, ~rma(yi, vi, data = .x)), # run meta-analysis
         output = map(rma, ~tidy(.x, conf.int = TRUE, exponentiate = TRUE)), # tidy model outputs, including cis and exponentiate for interpretability 
         model = map(rma, glance), # get full statistics
         predict = map(rma, ~predict(.x, transf = exp) %>% as.data.frame %>% select(starts_with("pi")))) %>% # prediction interval
  unnest(cols = c(output, model, predict)) # unnest

population_sa %>% 
  select(-c(data, rma)) %>%
  write_csv("output/tables/Population Subgroup Sensitivity Analysis.csv")

sa_plot <- quadsoa_cleaned %>%
  filter(outcome == "worsening JSN/OA",
         analysis_group == "whole tf",
         muscle == "quad") %>%
  mutate(population_subgroup = factor(population_subgroup, levels = c("risk of OA", "OA", "OA or risk of OA"))) %>%
  arrange(desc(population_subgroup))


png("output/plots/sensitivity population.png", height = 1000, pointsize = 25, width = 800)
forest(sa_plot$yi,
       sa_plot$vi,
       slab = sa_plot$studyname, 
       xlim = c(-5,5),
       ylim = c(-1,20),
       alim = c(log(1/10),log(10)), # limit of data clipping
       at = c(log(1/5), log(1/2), 0, log(2), log(5)), # axis limit (different to above)
       cex = 0.8,
       rows = c(1:2, 5.5:8.5, 12:16),
       atransf = exp,
       fonts = "Karla",
       xlab = "",
       header = c("Study", "Risk Ratio [95%CI]"), 
       efac = c(0,1)) # first element is the cis, second is the arrow. This eliminates the tick on the ci
# add subgroup labels
text(-5, c(3, 9.5, 17), pos = 4, c("Combined", "OA", "Risk of OA"), font = 2)
par(xpd=NA, font = 1) # remove plot clipping limits
text(log(c(1/5, 5)), -6 , c("Lower strength =\ndecreased risk","Lower strength =\nincreased risk"), pos=c(3,3), cex = 0.8) 
par(font = 2)
addpoly(population_sa %>% filter(population_subgroup == "OA or risk of OA") %>% pluck(6) %>% pluck(1)
, atransf = exp, cex = 0.8, row = 0, mlab = "    Random Effects Subgroup Model")
addpoly(population_sa %>% filter(population_subgroup == "OA") %>% pluck(6) %>% pluck(1)
        , atransf = exp, cex = 0.8, row = 4.5, mlab = "    Random Effects Subgroup Model")
addpoly(population_sa %>% filter(population_subgroup == "risk of OA") %>% pluck(6) %>% pluck(1)
        , atransf = exp, cex = 0.8, row = 11, mlab = "    Random Effects Subgroup Model")
par(font = 1)
mtext(text = "Sensitivity Analysis - Population Subgroup", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Radiographic OA - in presence of low knee extensor strength", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

# Function for forest plots
# Customises (heavily) based on the forest.default from metafor
# Most ideas for this customisation taken from metafor webpage as well as forest plots from "meta" package
forest_plotr <- function(outcome_a, analysis_group_a, muscle_a, supress = FALSE, extraspace = FALSE){
  
  # get outcome/predictor 
  data_filt <- quadsoa_meta_combined %>% filter(outcome == outcome_a, 
                               analysis_group == analysis_group_a,                 
                               muscle == muscle_a)
  
  # data for forest
  forestdata <- data_filt %>%
    filter(sex == "Overall") %>%
    select(rma) %>%
    pluck(1) %>% pluck(1)
  
  # labels for forest plot
  labels <- data_filt %>%
    filter(sex == "Overall") %>%
    select(data) %>%
    pluck(1) %>% 
    pluck(1) %>%
    select(studyname) %>%
      pluck(1)
  
  # original data for forest plot
  dat <- data_filt %>%
    filter(sex == "Overall") %>%
    select(data) %>%
    pluck(1) %>% 
    pluck(1) %>%
    mutate(ci.lb = yi - 1.96*sei,
           ci.ub = yi + 1.96*sei)
  
  # size of points
  psize <- weights(forestdata)
  psize <- 1.2 + (psize - min(psize)) / (max(psize) - min(psize))
  
  # prediction interval data
  pi_low <- data_filt$pi.lb[data_filt$sex == "Overall"]
  pi_upper <- data_filt$pi.ub[data_filt$sex == "Overall"]
  pi_est <- pi_upper - ((pi_upper - pi_low) / 2)
  pi_se <- (pi_upper - pi_low) / 3.92
  pi_text <- paste0("[", round(pi_low,2), ", ", round(pi_upper,2), "]")
  
  # get nobs for each
  femaleobs <- data_filt$nobs[data_filt$sex == "Women"] -1
  maleobs <- data_filt$nobs[data_filt$sex == "Men"] - 1
  mixedobs <- data_filt$nobs[data_filt$sex == "mixed"] - 1
  
  # get total obs
  nobs <- femaleobs + maleobs + mixedobs + 3
  
  # adjust rows
  mixedspace <- if(mixedobs == 0 | supress == TRUE) 2 else 3
  malespace <- if(maleobs == 0 | supress == TRUE) 3 else 5
  femalespace <- if(femaleobs == 0 | supress == TRUE) 3 else 5
  
  # starting row of each 
  mixedstart <- mixedspace
  malestart <- mixedspace + malespace + mixedobs
  femalestart <- mixedspace + malespace + femalespace + mixedobs + maleobs
  
  # rows to add studies on forest plot
  rows <- c(mixedspace:(mixedspace + mixedobs), 
            (mixedspace + malespace + mixedobs):(mixedspace + malespace + mixedobs + maleobs), 
            (mixedspace + malespace + femalespace + mixedobs + maleobs):(mixedspace + malespace + femalespace + mixedobs + maleobs + femaleobs))
  
  #subgroup label height
  mixedheight <- if(mixedobs == 0) 3.5 else (3 + 1.5 + mixedobs) 
  maleheight <- if(maleobs == 0) ((mixedspace + mixedobs) + 3.5 + 1) else ((mixedspace + mixedobs) + 5 + 1.5 + maleobs) 
  femaleheight <- if(femaleobs == 0) ((mixedspace + malespace + mixedobs + maleobs) + 3.5 + 1) else ((mixedspace + malespace  + mixedobs + maleobs) + 5 + 1.5 + femaleobs) 
  
  mixedheight <-if(supress == TRUE & mixedobs > 0) mixedheight - 1 else mixedheight
  maleheight <-if(supress == TRUE & maleobs > 0) maleheight - 1 else maleheight
  femaleheight <-if(supress == TRUE & femaleobs > 0) femaleheight - 1 else femaleheight

  grouplabs <- c(mixedheight, maleheight, femaleheight)
  
  # individual data for each subgroup
  mixeddata <- data_filt %>%
    filter(sex == "mixed") %>%
    select(rma) %>%
    pluck(1) %>% pluck(1)
  
  maledata <- data_filt %>%
    filter(sex == "Men") %>%
    select(rma) %>%
    pluck(1) %>% pluck(1)
  
  femaledata <- data_filt %>%
    filter(sex == "Women") %>%
    select(rma) %>%
    pluck(1) %>% pluck(1)
  
  # labeling function for heterogeneity
  mlabfun2 <- bquote(paste(
                      #" (Q = ", .(formatC(forestdata$QE, digits=2, format="f")),
                      #", df = ", .(forestdata$k - forestdata$p),
                      #", p ", .(metafor:::.pval(forestdata$QEp, digits=2, showeq=TRUE, sep=" ")), 
                      I^2, " = ", .(formatC(forestdata$I2, digits=1, format="f")), "%, ",
                      tau^2, " = ", .(formatC(forestdata$tau2, digits=2, format="f"))))
  
  par(xpd = FALSE)
  
  # plot - use raw data points instead of from rma object to allow customisation
  plot <- forest(dat$yi,
         dat$vi,
         slab = labels, 
         xlim = c(-5,5),
         ylim = c(-2,femaleheight + 3.5), # adjust y axis based on number of studies
         alim = c(log(1/10),log(10)), # limit of data clipping
         at = c(log(1/5), log(1/2), 0, log(2), log(5)), # axis limit (different to above)
         cex = 0.8,
         psize = psize,
         rows = rows,
         atransf = exp,
         fonts = "Karla",
         xlab = "",
         header = c("Study", "Risk Ratio [95%CI]"), 
         efac = c(0,1)) # first element is the cis, second is the arrow. This eliminates the tick on the ci
  
  # manually adding points for each study back in so can colour
  points(dat$yi, rows, pch = 15, cex=psize, col = "grey")
  
  #segments(dat$ci.lb, rows, dat$ci.ub, rows, lwd=1.5)
  
  # add subgroup labels
  text(-5, grouplabs, pos = 4, c("Mixed", "Men", "Women"), font = 2)
  
  par(xpd=NA, font = 1) # remove plot clipping limits
  
  # text for axis descriptors
  if(extraspace == TRUE) text(log(c(1/5, 5)), -7.5, c("Lower strength =\ndecreased risk","Lower strength =\nincreased risk"), pos=c(3,3), cex = 0.8) else
    text(log(c(1/5, 5)), -7, c("Lower strength =\ndecreased risk","Lower strength =\nincreased risk"), pos=c(3,3), cex = 0.8)
  
  #Add text for summary descriptors
  text(-5, -1, "Overall RE Model", pos = 4, cex = 0.8, font = 2)
  text(-5, -2, "Prediction Interval", pos = 4, cex = 0.8)
  text(-5, -3, mlabfun2, pos = 4, cex = 0.8)
  text(5, -2, pi_text, pos = 2, cex = 0.8) # add prediction interval text
  
  # Conditionally add summary polygons depending on whether > 1 study and if elect to suppres in function call
  if(mixedobs > 0 & supress != TRUE) addpoly(mixeddata, atransf = exp, cex = 0.8, row = mixedstart - 1.5, mlab = "")
  if(maleobs > 0) addpoly(maledata, atransf = exp, cex = 0.8, row = malestart - 1.5, mlab = "")
  if(femaleobs > 0) addpoly(femaledata, atransf = exp, cex = 0.8, row = femalestart - 1.5, mlab = "")
  
  # Add prediction interval
  segments(x0 = log(pi_low), x1 = log(pi_upper), lwd = 5,  cex = 0.8, y0 = -2, col = "dark grey")
  
  # Add overall summary effet
  addpoly(forestdata, row = -1, atransf = exp, cex = 0.8,  mlab = "", efac=2)
}


# Individual plots with titles added also

png("output/plots/worsening wholetf quad.png", height = 1110, pointsize = 25, width = 800)
forest_plotr("worsening JSN/OA", "whole tf", "quad", extraspace = TRUE)
mtext(text = "(A) Radiographic OA Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Whole Tibiofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/worsening whole pf quad.png", height = 740, pointsize = 25, width = 800)
forest_plotr("worsening JSN/OA", "whole pf", "quad", supress = TRUE)
mtext(text = "(B) Radiographic OA Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Whole Patellofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()


png("output/plots/quantJSN medialtf quad.png", height = 740, pointsize = 25, width = 800)
forest_plotr("quantitative JSN", "medial tf", "quad", supress = TRUE)
mtext(text = "Quantitative Joint Space Narrowing", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Medial Tibiofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()


png("output/plots/quant medialtf quad.png", height = 740, pointsize = 25, width = 800)
forest_plotr("quantitative cartilage progression", "medial tf", "quad", supress = TRUE)
mtext(text = "Quantitative Cartilage Thinning", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Medial Tibiofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/semiquant medialtf quad.png", height = 850, pointsize = 25, width = 800)
forest_plotr("semi quantitative cartilage progression", "medial tf", "quad", supress = TRUE)
mtext(text = "(A) Cartilage Lesion Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Medial Tibiofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/semiquant lateraltf quad.png", height = 790, pointsize = 25, width = 800)
forest_plotr("semi quantitative cartilage progression", "lateral tf", "quad", supress = TRUE)
mtext(text = "(B) Cartilage Lesion Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Lateral Tibiofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/semiquant medial pf quad.png", height = 740, pointsize = 25, width = 800)
forest_plotr("semi quantitative cartilage progression", "medial pf", "quad", supress = TRUE)
mtext(text = "(C) Cartilage Lesion Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Medial Patellofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/semiquant lateral pf quad.png", height = 790, pointsize = 25, width = 800)
forest_plotr("semi quantitative cartilage progression", "lateral pf", "quad", supress = TRUE)
mtext(text = "(D) Cartilage Lesion Worsening", side = 3, adj = 0.01, line = 2, font = 2, cex = 1.3)
mtext(text = "Lateral Patellofemoral", side = 3, adj = 0.01, line = 1, font = 1, cex = 1)
dev.off()

png("output/plots/worsening wholetf hs.png", height = 850, pointsize = 25, width = 800)
forest_plotr("worsening JSN/OA", "whole tf", "hs", supress = TRUE)
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
3

# Leave one out analysis
## First run the influence analysis with metafor package - gives all indepth statistics and measures
## Join these to original data
## Meta package provides nicer influence plot, but use the data from the metafor output in the plot.

quadsoa_meta_loo <- quadsoa_analysis %>%
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
pluck(quadsoa_meta_loo$influence) %>% pluck(2) %>%
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

