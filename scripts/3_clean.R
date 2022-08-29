# Clean

# remove whitespace around entered text
quadsoa <- quadsoa %>%
  mutate(across(c(outcome:alignment_subgroup), str_trim)) %>%
  mutate(baseline_risk = case_when(
    is.na(baseline_risk) ~ mean(baseline_risk, na.rm = TRUE),
    TRUE ~ baseline_risk
  ))

# invert any measures that need inverting 
# first need to convert any RR back to OR as required
quadsoa <- quadsoa %>%
  mutate(or = case_when(
    invert_required == 1 & is.na(or) ~ rr_to_or(rr, baseline_risk),
    TRUE ~ or),
    or_ci_low = case_when(
      invert_required == 1 & is.na(or_ci_low) ~ rr_to_or(rr_ci_low, baseline_risk),
      TRUE ~ or_ci_low),
    or_ci_high = case_when(
      invert_required == 1 & is.na(or_ci_high) ~ rr_to_or(rr_ci_high, baseline_risk),
      TRUE ~ or_ci_high)
  ) %>%
  # now can invert the odds ratios
  mutate(across(c(or, or_ci_low, or_ci_high), 
                ~ case_when(
                  invert_required == 1 ~ 1/.x,
                  TRUE ~ .x))) %>%
  # need to switch high/low cis for inverted rows
  transform(., 
            or_ci_low = case_when(
              invert_required == 1 ~ or_ci_high,
              TRUE ~ or_ci_low),
            or_ci_high = case_when(
              invert_required == 1 ~ or_ci_low,
              TRUE ~ or_ci_high)) %>%
  # need to switch reference groups for any rows where invert was required
  transform(., 
            reference_group = case_when(
              invert_required == 1 ~ explanatory_group,
              TRUE ~ reference_group), 
            explanatory_group = case_when(
              invert_required == 1 ~ reference_group,
              TRUE ~ explanatory_group))

# calculate SMD for continuous data
quadsoa <- quadsoa %>%
  filter(!is.na(prog_mean)) %>% # do this only for rows that have continuous mean/sd etc
  escalc(measure = "SMD", n2i = prog_n, m2i = prog_mean, sd2i = prog_sd,
         n1i = nonprog_n, m1i = nonprog_mean, sd1i = nonprog_sd, data = ., var.names = c("smd", "smd_variance")) %>%
  summary(out.names=c("sei","zi","pval","smd_ci_low","smd_ci_high")) %>% # get full summary including cis
  select(-c(sei, zi, pval)) %>% # remove un-needed data
  bind_rows(., quadsoa %>% filter(is.na(prog_mean))) # join back to original data

# calculate OR from SMD
# ln(OR) = pi/sqrt(3)*SMD
quadsoa <- quadsoa %>%
  mutate(or = case_when(
    !is.na(smd) ~ exp(pi/sqrt(3)*smd),
    TRUE ~ or),
    or_ci_low = case_when(
      !is.na(smd_ci_low) ~ exp(pi/sqrt(3)*smd_ci_low),
      TRUE ~ or_ci_low),
    or_ci_high = case_when(
      !is.na(smd_ci_high) ~ exp(pi/sqrt(3)*smd_ci_high),
      TRUE ~ or_ci_high))

# calcualte RR from OR
quadsoa <- quadsoa %>%
  mutate(rr = case_when(
    is.na(rr) & !is.na(or) ~ or_to_rr(or, baseline_risk),
    TRUE ~ rr),
    rr_ci_low = case_when(
      is.na(rr_ci_low) & !is.na(or_ci_low) ~ or_to_rr(or_ci_low, baseline_risk),
      TRUE ~ rr_ci_low),
    rr_ci_high = case_when(
      is.na(rr_ci_high) & !is.na(or_ci_high) ~ or_to_rr(or_ci_high, baseline_risk),
      TRUE ~ rr_ci_high))

# calculate RR from contingency tables
quadsoa <- quadsoa %>%
  mutate(rr = case_when(
    is.na(rr) & !is.na(ref_progress) ~ (exp_progress/exp_total)/(ref_progress/ref_total),
    TRUE ~ rr),
  rr_ci_low = case_when(
    is.na(rr_ci_low) & !is.na(ref_progress) ~ 
      exp(log(rr)-1.96*sqrt(1/ref_progress + 1/exp_progress - 1/ref_total - 1/exp_total)), 
    TRUE ~ rr_ci_low),
  rr_ci_high = case_when(
    is.na(rr_ci_high) & !is.na(ref_progress) ~ 
      exp(log(rr)+1.96*sqrt(1/ref_progress + 1/exp_progress - 1/ref_total - 1/exp_total)), 
    TRUE ~ rr_ci_high))

# Calculate log rr and log variance, log se

quadsoa <- quadsoa %>%
  mutate(yi = log(rr), # log risk ratio
         sei = (log(rr_ci_high) - log(rr_ci_low))/3.92, # log standard error
         vi = sei^2) # log variance

# Filter out any excluded pieces of data
quadsoa <- quadsoa %>%
  filter(is.na(exclude))

# MA of studies with alignment subgroups
# Some studies include different subgroups not of interest in this review
# Need to combine together with a 'pre' meta-analysis

premeta <- quadsoa %>%
  filter(alignment_subgroup != "All") %>% # find studies that include alignments subgroups
  group_by(author, sex, analysis_group) %>% # group by each publication and subgroups of interest in this review, for later 'pre-metaanalyses'
  nest(data = -c(outcome:muscle_test, reference_group, explanatory_group)) %>%# nest data, remove data that is consistent across subgroups
  mutate(rma = map(data, ~tidy(rma(yi, vi, data = .x)) %>% # for each group, perform a rma
                     select(estimate, std.error) %>% # take only the estimated log rr and log se
                     rename(yi = estimate, # rename for consistency with final dataframe
                            sei = std.error) %>%
                     mutate(vi = sei^2)
                     )) %>%
  select(-data) %>% # remove the original data as no longer needed
  unnest(cols = c(rma)) # un-nest back to same format

# now need to join the 'pre-meta' data back to the original dataset
quadsoa_cleaned <- quadsoa %>% 
  filter(alignment_subgroup == "All") %>% # take all data not split before
  bind_rows(., premeta) %>% # bind the pre-meta data 
  mutate(studyname = paste(author, year, sep = " ")) %>%
  select(c(studyname, author, year, population_subgroup, outcome, sex:reference_group, explanatory_group, or:vi))
  
write_csv(quadsoa_cleaned, "data/processed/Strength OA Review Cleaned.csv")
  
