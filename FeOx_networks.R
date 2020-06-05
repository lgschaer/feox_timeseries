# Required Packages

#install.packages("Hmisc")
library(Hmisc)
library(phyloseq)
library(tidyverse)
library(vegan)
library(ape)

# Start with rarefied phyloseq object
pso <- readRDS("/home/lgschaer/old/FeOx_demux/feox-phyloseq.rds")
pso

# Filter out eukaryotes and mitochondria
justbacteria <- pso %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
justbacteria

# Modify data to get matrix of sampleIDs as rownames and Taxa classification as colnames
pso_subset <- justbacteria %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  mutate(
    Abundance = Abundance*100,
    Abundance = as.numeric(as.integer(Abundance))
  ) %>%
  filter(Location == "GLRC" & Metal != "Swab") %>%
  unite(TaxID, Class:Family, sep = "_") %>%
  select(TaxID, Sample, Abundance) %>%
  pivot_wider(names_from = "TaxID", values_from = "Abundance", values_fn = list(Abundance = sum)) %>%
  arrange(Sample)%>%
  column_to_rownames(var = "Sample")
View(pso_subset)
class(pso_subset)

# Use Hmisc to do spearman correlation on matrix made in last step
res <- rcorr(as.matrix(pso_subset), type = "spearman")
summary(res)

# Pull out Spearman coefficients and mutate to get in usable format
corr_co <- res$r %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "Spearman_Co") %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
View(corr_co)

# Pull out P-values and mutate to get in usable format
p_val <- res$P %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "P_Value", values_drop_na = TRUE) %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
View(p_val)

# Combine spearman coefficients and p-values and filter to keep only significant observations
spearman <- p_val %>%
  left_join(corr_co, by = "Interaction")%>%
  filter(!is.na(Spearman_Co)) %>%
  filter(Spearman_Co != 1) %>%
  filter(P_Value <= 0.01) %>%
  separate(Interaction, into = c("Taxa1", "Taxa2"), sep = "_and_") %>%
  mutate(
    Relationship = ifelse(Spearman_Co >= 0, "Positive", "Negative")
  )
head(spearman)

# Make heatmap to show correlations between taxa
ggplot(data = spearman)+
  geom_tile(aes(x = Taxa1, y = Taxa2, fill = Relationship))+
  theme(axis.text.y.left = element_blank(),
        axis.text.x = element_blank())


# Make a matrix of spearman correlation between taxa
spearman2 <- spearman %>%
  select(-c("P_Value", "Relationship")) %>%
  pivot_wider(names_from = "Taxa2", values_from = "Spearman_Co", values_fill = list("Spearman_Co" = 0)) %>%
  group_by(Taxa1) %>%
  summarise_if(is.numeric, abs) %>%
  column_to_rownames(var = "Taxa1")
head(spearman2)

# Make ordination, bray curtis
dist <- vegdist(spearman2,  method = "bray")

# Calculate Axis values
PCOA <- pcoa(dist)

# Pull axis vectors out of PCoA list
spearman_dist <- PCOA$vectors
head(spearman_dist)

#select axis 1 and axis 2
spearman_dist2 <- spearman_dist %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa1") %>%
  mutate(Taxa1 = as.character(Taxa1),
         Taxa1b = Taxa1) %>%
  as_tibble() %>%
  select(c("Taxa1", "Taxa1b", "Axis.1", "Axis.2"))
head(spearman_dist2)

# Match axis 1 and axis 2 up with Taxa classification for both Taxa 1 and Taxa 2
spearman_net <- spearman_dist2 %>%
  right_join(spearman, by = "Taxa1") %>%
  mutate(
    Origin_Taxa = Taxa1,
    Origin_X = Axis.1, 
    Origin_Y = Axis.2,
    Taxa1b = Taxa2
  ) %>%
  select(c("Origin_Taxa", "Origin_X", "Origin_Y", "Taxa1b", "Spearman_Co", "Relationship")) %>%
  left_join(spearman_dist2, by = "Taxa1b") %>%
  mutate(
    Dest_Taxa = Taxa1b,
    Dest_X = Axis.1, 
    Dest_Y = Axis.2) %>%
  select(c("Origin_Taxa", "Origin_X", "Origin_Y", "Spearman_Co", "Relationship", "Dest_Taxa", "Dest_X", "Dest_Y"))
View(spearman_net)

# Make a list of points we want to label (since we only want to label some of the points)
spearman_labels_df <- spearman_net %>% filter(abs(Spearman_Co) >= 0.6)
spearman_labels <- unique(spearman_labels_df$Origin_Taxa)
spearman_labels
length(spearman_labels)

# Label using geom_text
ggplot(data = spearman_net, aes(x = Origin_X, y = Origin_Y))+
  geom_segment(aes(x = Origin_X, y = Origin_Y, xend = Dest_X, yend = Dest_Y, color = Relationship))+
  geom_point(shape = 21, fill = "black", aes(size = abs(Spearman_Co)))+
  geom_text(aes(label = Origin_Taxa), check_overlap = TRUE, vjust = -1, 
            color = "gray15", size = 3,
            data = filter(spearman_net, Origin_Taxa %in% spearman_labels))+
  theme_classic()


