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



# Function for forest plots
# Customises (heavily) based on the forest.default from metafor
# Most ideas for this customisation taken from metafor webpage as well as forest plots from "meta" package
forest_plotr <- function(data = NULL, outcome = NULL, analysis_group = NULL, muscle = NULL, supress = FALSE, extraspace = FALSE, sumtext = FALSE){
 
  outcome_a <- outcome
  analysis_group_a <- analysis_group
  muscle_a <- muscle
  
  # define filter call, if no parameter specified then NULL to filter. rlang::exprs is passed to dplyr filter later.
  outcomefilt <- if(is.null(outcome_a)) NULL else rlang::exprs(outcome == !!outcome_a)
  analysisfilt <- if(is.null(analysis_group_a)) NULL else rlang::exprs(analysis_group == !!analysis_group_a)
  musclefilt <- if(is.null(muscle_a)) NULL else rlang::exprs(muscle ==  !!muscle_a)
  
  # get data by outcome/predictor 
  data_filt <- if(is.null(data)) secondary_meta_combined else data
     
  
  # get data by outcome/predictor 
  data_filt <- data_filt %>% filter(!!!outcomefilt, 
                                     !!!analysisfilt,                 
                                     !!!musclefilt)
  
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
  
  # labeling function for combo plots
  mlabfun_mixed <- bquote(paste("RE Model - ",
                                #" (Q = ", .(formatC(forestdata$QE, digits=2, format="f")),
                                #", df = ", .(forestdata$k - forestdata$p),
                                #", p ", .(metafor:::.pval(forestdata$QEp, digits=2, showeq=TRUE, sep=" ")), 
                                I^2, " = ", .(formatC(mixeddata$I2, digits=1, format="f")), "%, ",
                                tau^2, " = ", .(formatC(mixeddata$tau2, digits=2, format="f"))))
  
  mlabfun_male <- bquote(paste("RE Model - ",
                               #" (Q = ", .(formatC(forestdata$QE, digits=2, format="f")),
                               #", df = ", .(forestdata$k - forestdata$p),
                               #", p ", .(metafor:::.pval(forestdata$QEp, digits=2, showeq=TRUE, sep=" ")), 
                               I^2, " = ", .(formatC(maledata$I2, digits=1, format="f")), "%, ",
                               tau^2, " = ", .(formatC(maledata$tau2, digits=2, format="f"))))
  
  mlabfun_female <- bquote(paste("RE Model - ",
                                 #" (Q = ", .(formatC(forestdata$QE, digits=2, format="f")),
                                 #", df = ", .(forestdata$k - forestdata$p),
                                 #", p ", .(metafor:::.pval(forestdata$QEp, digits=2, showeq=TRUE, sep=" ")), 
                                 I^2, " = ", .(formatC(femaledata$I2, digits=1, format="f")), "%, ",
                                 tau^2, " = ", .(formatC(femaledata$tau2, digits=2, format="f"))))
  
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
  
  if(sumtext == TRUE) text(-5, mixedstart -1.5, pos = 4, mlabfun_mixed, cex = 0.8)
  if(sumtext == TRUE) text(-5, malestart -1.5, pos = 4, mlabfun_male, cex = 0.8)
  if(sumtext == TRUE) text(-5, femalestart -1.5, pos = 4, mlabfun_female, cex = 0.8)
  
  # Add overall summary effet
  addpoly(forestdata, row = -1, atransf = exp, cex = 0.8,  mlab = "", efac=2)
}
