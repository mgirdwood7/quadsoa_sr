# Load

library(tidyverse)
library(metafor)
library(broom)
library(meta)

quadsoa <- read.csv("data/raw/dataextraction_20220822.csv", header = TRUE, na = c("NA", ""))