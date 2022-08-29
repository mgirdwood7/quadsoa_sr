# Functions

# Convert odds ratio to risk ratio
or_to_rr <- function(or, baselinerisk){
  # or = extracted odds ratio
  # baselinerisk = baseline risk probability)
  or/(1 - baselinerisk + (baselinerisk * or))
}

# Convert risk ratio to odds ratio
rr_to_or <- function(rr, baselinerisk){
  # rr = extracted risk ratio
  # baselinerisk = baseline risk probablity
  ((1 - baselinerisk) * rr) / (1 - rr * baselinerisk)
}

smd_cell <- function(n1, mean1, sd1, n2, mean2, sd2, variable){
  escalc("SMD", n1i = n1, m1i = mean1, sd1i = sd1, n2i = n2, m2i = mean2, sd2i = sd2, data = quadsoa) %>%
    summary() %>%
    select()
}


foresthelper <- function(data, outcome_a, muscle_a){
  
  # get outcome/predictor 
  data_filt <- data %>% filter(outcome == outcome_a, 
                          muscle == muscle_a)
  # get nobs for each
  femaleobs <- data_filt$nobs[data_filt$sex == "Women"] -1
  maleobs <- data_filt$nobs[data_filt$sex == "Men"] - 1
  mixedobs <- data_filt$nobs[data_filt$sex == "mixed"] - 1
  
  # adjust rows
  mixedspace <- if(mixedobs == 0) 2 else 3
  malespace <- if(maleobs == 0) 3 else 5
  femalespace <- if(femaleobs == 0) 3 else 5
  
  rows <- c(mixedspace:(mixedspace + mixedobs), 
            (mixedspace + malespace + mixedobs):(mixedspace + malespace + mixedobs + maleobs), 
            (mixedspace + malespace + femalespace + mixedobs + maleobs):(mixedspace + malespace + femalespace + mixedobs + maleobs + femaleobs))
  return(rows)
}
