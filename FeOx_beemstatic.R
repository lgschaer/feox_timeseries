#packages used for this script
library(devtools)
#install.packages("foreach", repos="http://R-Forge.R-project.org")
library(foreach)
#install.packages("doMC")   #only works on mac or linux os
library(doMC)
#install.packages("lokern")
library(lokern)
#install.packages("pspline")
library(pspline)
#install.packages("monomvn")
library(monomvn)
library(beemStatic)
library(tidyverse)
library(csv)
#trying beemstatic instead
library(beemStatic)
data("beemDemo")
attach(beemDemo)
#?beemDemo
#install.packages("ggraph")
library(ggraph)

#-----------"beemstatic" package -------------#
#https://github.com/lch14forever/BEEM-static

#start with rarefied phyloseq object
pso <- readRDS("/home/lgschaer/old/FeOx_demux/feox-phyloseq.rds")
pso

#filter out eukaryotes and mitochondria
justbacteria <- pso %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
justbacteria

pso_subset <- justbacteria %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance >= 0.0001) %>%                      # Filter out low abundance taxa - CHANGED THIS THROUGHOUT ANALYSIS
  mutate(
    Count1 = Abundance * 10000,                        # Multiply all abundances by number large enough so rel. abundance is not a decimal
    Count2 = as.integer(Count1),
    Count3 = as.numeric(Count2),                       # Convert counts first to integer, then to numeric characters
    Taxa = as.character(Family)                        # Change "Family" to character (instead of factor), rename to "Taxa"
  ) %>%
  group_by(Taxa) %>%                                   # Group by "Taxa"
  add_count(Taxa) %>%                                  # Count how many times each taxa appears in the data set
  filter(n >= 17) %>%                                  # Filter out the taxa that appear less than x number of times, this parameter can be changed depending on taxa level and abundance chosen
  ungroup() %>%
  unite(TaxID, Class:Family) %>%                       # Merge desired taxonomic levels together to create labels for network graph
  select(TaxID, Sample, Count3) %>%                    # Select desired columns
  spread(key = Sample, value = Count3, fill = 0) %>%   # Spread data frame to get desired OTU table format
  column_to_rownames(var = "TaxID")                    # Convert "TaxID" column to rownames
dim(pso_subset)
#head(pso_subset)

counts <- data.matrix(pso_subset, rownames.force = TRUE)   # Convert to matrix format
#dim(counts)
#head(counts)
#View(pso_subset)

#test data set
#data("beemDemo")
#attach(beemDemo)
#res <- func.EM(dat.w.noise, ncpu=4, scaling=median(biomass.true), max.iter=200, epsilon = 1e-4)
#res
#class(res)
#showInteraction(res, dat.w.noise)
#dim(dat.w.noise)
#dim(counts)

#run beemstatic with my data
res <- func.EM(counts, dev = 1000, ncpu=4, scaling=median(biomass.true), max.iter=200, epsilon = 1e-4)

showInteraction(res, counts)