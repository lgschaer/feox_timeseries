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
  select(TaxID, Sample, Time, Abundance) %>%
  pivot_wider(names_from = "TaxID", values_from = "Abundance", values_fn = list(Abundance = sum)) %>%
  arrange(Sample, Time)%>%
  column_to_rownames(var = "Sample")
#View(pso_subset)
class(pso_subset)

pso1 <- pso_subset %>%
  filter(Time == 1)

pso1 <- pso_subset %>%
  filter(Time == 1)

pso2 <- pso_subset %>%
  filter(Time == 2)

pso3 <- pso_subset %>%
  filter(Time == 3)

pso4 <- pso_subset %>%
  filter(Time == 4)

pso5 <- pso_subset %>%
  filter(Time == 5)

pso6 <- pso_subset %>%
  filter(Time == 6)

# Use Hmisc to do spearman correlation on matrix made in last step
res1 <- rcorr(as.matrix(pso1), type = "spearman")
summary(res1)

res2 <- rcorr(as.matrix(pso2), type = "spearman")
summary(res2)

res3 <- rcorr(as.matrix(pso3), type = "spearman")
summary(res3)

res4 <- rcorr(as.matrix(pso4), type = "spearman")
summary(res4)

res5 <- rcorr(as.matrix(pso5), type = "spearman")
summary(res5)

res6 <- rcorr(as.matrix(pso6), type = "spearman")
summary(res6)

# Pull out Spearman coefficients and mutate to get in usable format
corr_co1 <- res1$r %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "Spearman_Co") %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(corr_co1)

corr_co2 <- res2$r %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "Spearman_Co") %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(corr_co2)

corr_co3 <- res3$r %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "Spearman_Co") %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(corr_co3)

corr_co4 <- res4$r %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "Spearman_Co") %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(corr_co4)

corr_co5 <- res5$r %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "Spearman_Co") %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(corr_co5)

corr_co6 <- res6$r %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "Spearman_Co") %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(corr_co6)

# Pull out P-values and mutate to get in usable format
p_val1 <- res1$P %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "P_Value", values_drop_na = TRUE) %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(p_val1)

p_val2 <- res2$P %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "P_Value", values_drop_na = TRUE) %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(p_val2)

p_val3 <- res3$P %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "P_Value", values_drop_na = TRUE) %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(p_val3)

p_val4 <- res4$P %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "P_Value", values_drop_na = TRUE) %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(p_val4)

p_val5 <- res5$P %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "P_Value", values_drop_na = TRUE) %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(p_val5)

p_val6 <- res6$P %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Taxa2", values_to = "P_Value", values_drop_na = TRUE) %>%
  rownames_to_column(var = "Taxa1") %>%
  separate(Taxa1, into = c("Taxa1", "RowNum"), sep = "\\...") %>%
  select(-c("RowNum")) %>%
  unite(Interaction, Taxa1, Taxa2, sep = "_and_")
#View(p_val6)

# Combine spearman coefficients and p-values and filter to keep only significant observations
spearman1 <- p_val1 %>%
  left_join(corr_co1, by = "Interaction")%>%
  filter(!is.na(Spearman_Co)) %>%
  filter(Spearman_Co != 1) %>%
  filter(P_Value <= 0.01) %>%
  separate(Interaction, into = c("Taxa1", "Taxa2"), sep = "_and_") %>%
  mutate(
    Relationship = ifelse(Spearman_Co >= 0, "Positive", "Negative")
  )
head(spearman1)

spearman2 <- p_val2 %>%
  left_join(corr_co2, by = "Interaction")%>%
  filter(!is.na(Spearman_Co)) %>%
  filter(Spearman_Co != 1) %>%
  filter(P_Value <= 0.01) %>%
  separate(Interaction, into = c("Taxa1", "Taxa2"), sep = "_and_") %>%
  mutate(
    Relationship = ifelse(Spearman_Co >= 0, "Positive", "Negative")
  )
head(spearman2)

spearman3 <- p_val3 %>%
  left_join(corr_co3, by = "Interaction")%>%
  filter(!is.na(Spearman_Co)) %>%
  filter(Spearman_Co != 1) %>%
  filter(P_Value <= 0.01) %>%
  separate(Interaction, into = c("Taxa1", "Taxa2"), sep = "_and_") %>%
  mutate(
    Relationship = ifelse(Spearman_Co >= 0, "Positive", "Negative")
  )
head(spearman3)

spearman4 <- p_val4 %>%
  left_join(corr_co4, by = "Interaction")%>%
  filter(!is.na(Spearman_Co)) %>%
  filter(Spearman_Co != 1) %>%
  filter(P_Value <= 0.01) %>%
  separate(Interaction, into = c("Taxa1", "Taxa2"), sep = "_and_") %>%
  mutate(
    Relationship = ifelse(Spearman_Co >= 0, "Positive", "Negative")
  )
head(spearman4)

spearman5 <- p_val5 %>%
  left_join(corr_co5, by = "Interaction")%>%
  filter(!is.na(Spearman_Co)) %>%
  filter(Spearman_Co != 1) %>%
  filter(P_Value <= 0.01) %>%
  separate(Interaction, into = c("Taxa1", "Taxa2"), sep = "_and_") %>%
  mutate(
    Relationship = ifelse(Spearman_Co >= 0, "Positive", "Negative")
  )
head(spearman5)

spearman6 <- p_val6 %>%
  left_join(corr_co6, by = "Interaction")%>%
  filter(!is.na(Spearman_Co)) %>%
  filter(Spearman_Co != 1) %>%
  filter(P_Value <= 0.01) %>%
  separate(Interaction, into = c("Taxa1", "Taxa2"), sep = "_and_") %>%
  mutate(
    Relationship = ifelse(Spearman_Co >= 0, "Positive", "Negative")
  )
head(spearman6)

# Make heatmap to show correlations between taxa
ggplot(data = spearman1)+
  geom_tile(aes(x = Taxa1, y = Taxa2, fill = Relationship))+
  theme(axis.text.y.left = element_blank(),
        axis.text.x = element_blank())


# Make a matrix of spearman correlation between taxa
spearmanB1 <- spearman1 %>%
  select(-c("P_Value", "Relationship")) %>%
  pivot_wider(names_from = "Taxa2", values_from = "Spearman_Co", values_fill = list("Spearman_Co" = 0)) %>%
  group_by(Taxa1) %>%
  summarise_if(is.numeric, abs) %>%
  column_to_rownames(var = "Taxa1")
head(spearmanB1)

spearmanB2 <- spearman2 %>%
  select(-c("P_Value", "Relationship")) %>%
  pivot_wider(names_from = "Taxa2", values_from = "Spearman_Co", values_fill = list("Spearman_Co" = 0)) %>%
  group_by(Taxa1) %>%
  summarise_if(is.numeric, abs) %>%
  column_to_rownames(var = "Taxa1")
head(spearmanB2)

spearmanB3 <- spearman3 %>%
  select(-c("P_Value", "Relationship")) %>%
  pivot_wider(names_from = "Taxa2", values_from = "Spearman_Co", values_fill = list("Spearman_Co" = 0)) %>%
  group_by(Taxa1) %>%
  summarise_if(is.numeric, abs) %>%
  column_to_rownames(var = "Taxa1")
head(spearmanB3)

spearmanB4 <- spearman4 %>%
  select(-c("P_Value", "Relationship")) %>%
  pivot_wider(names_from = "Taxa2", values_from = "Spearman_Co", values_fill = list("Spearman_Co" = 0)) %>%
  group_by(Taxa1) %>%
  summarise_if(is.numeric, abs) %>%
  column_to_rownames(var = "Taxa1")
head(spearmanB4)

spearmanB5 <- spearman5 %>%
  select(-c("P_Value", "Relationship")) %>%
  pivot_wider(names_from = "Taxa2", values_from = "Spearman_Co", values_fill = list("Spearman_Co" = 0)) %>%
  group_by(Taxa1) %>%
  summarise_if(is.numeric, abs) %>%
  column_to_rownames(var = "Taxa1")
head(spearmanB5)

spearmanB6 <- spearman6 %>%
  select(-c("P_Value", "Relationship")) %>%
  pivot_wider(names_from = "Taxa2", values_from = "Spearman_Co", values_fill = list("Spearman_Co" = 0)) %>%
  group_by(Taxa1) %>%
  summarise_if(is.numeric, abs) %>%
  column_to_rownames(var = "Taxa1")
head(spearmanB6)

# Make ordination, bray curtis
dist1 <- vegdist(spearmanB1,  method = "bray")
dist2 <- vegdist(spearmanB2,  method = "bray")
dist3 <- vegdist(spearmanB3,  method = "bray")
dist4 <- vegdist(spearmanB4,  method = "bray")
dist5 <- vegdist(spearmanB5,  method = "bray")
dist6 <- vegdist(spearmanB6,  method = "bray")

# Calculate Axis values
PCOA1 <- pcoa(dist1)
PCOA2 <- pcoa(dist2)
PCOA3 <- pcoa(dist3)
PCOA4 <- pcoa(dist4)
PCOA5 <- pcoa(dist5)
PCOA6 <- pcoa(dist6)

# Pull axis vectors out of PCoA list
spearman_dist1 <- PCOA1$vectors
#head(spearman_dist1)

spearman_dist2 <- PCOA2$vectors
#head(spearman_dist2)

spearman_dist3 <- PCOA3$vectors
#head(spearman_dist3)

spearman_dist4 <- PCOA4$vectors
#head(spearman_dist4)

spearman_dist5 <- PCOA5$vectors
#head(spearman_dist5)

spearman_dist6 <- PCOA6$vectors
#head(spearman_dist6)

#select axis 1 and axis 2
spearman_distB1 <- spearman_dist1 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa1") %>%
  mutate(Taxa1 = as.character(Taxa1),
         Taxa1b = Taxa1) %>%
  as_tibble() %>%
  select(c("Taxa1", "Taxa1b", "Axis.1", "Axis.2"))
#head(spearman_distB1)

spearman_distB2 <- spearman_dist2 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa1") %>%
  mutate(Taxa1 = as.character(Taxa1),
         Taxa1b = Taxa1) %>%
  as_tibble() %>%
  select(c("Taxa1", "Taxa1b", "Axis.1", "Axis.2"))
#head(spearman_distB2)

spearman_distB3 <- spearman_dist3 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa1") %>%
  mutate(Taxa1 = as.character(Taxa1),
         Taxa1b = Taxa1) %>%
  as_tibble() %>%
  select(c("Taxa1", "Taxa1b", "Axis.1", "Axis.2"))
#head(spearman_distB3)

spearman_distB4 <- spearman_dist4 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa1") %>%
  mutate(Taxa1 = as.character(Taxa1),
         Taxa1b = Taxa1) %>%
  as_tibble() %>%
  select(c("Taxa1", "Taxa1b", "Axis.1", "Axis.2"))
#head(spearman_distB4)

spearman_distB5 <- spearman_dist5 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa1") %>%
  mutate(Taxa1 = as.character(Taxa1),
         Taxa1b = Taxa1) %>%
  as_tibble() %>%
  select(c("Taxa1", "Taxa1b", "Axis.1", "Axis.2"))
#head(spearman_distB5)

spearman_distB6 <- spearman_dist6 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa1") %>%
  mutate(Taxa1 = as.character(Taxa1),
         Taxa1b = Taxa1) %>%
  as_tibble() %>%
  select(c("Taxa1", "Taxa1b", "Axis.1", "Axis.2"))
#head(spearman_distB6)

# Match axis 1 and axis 2 up with Taxa classification for both Taxa 1 and Taxa 2
spearman_net1 <- spearman_distB1 %>%
  right_join(spearman1, by = "Taxa1") %>%
  mutate(
    Origin_Taxa = Taxa1,
    Origin_X = Axis.1, 
    Origin_Y = Axis.2,
    Taxa1b = Taxa2
  ) %>%
  select(c("Origin_Taxa", "Origin_X", "Origin_Y", "Taxa1b", "Spearman_Co", "Relationship")) %>%
  left_join(spearman_distB1, by = "Taxa1b") %>%
  mutate(
    Dest_Taxa = Taxa1b,
    Dest_X = Axis.1, 
    Dest_Y = Axis.2,
    Time = "1") %>%
  select(c("Origin_Taxa","Time", "Origin_X", "Origin_Y", "Spearman_Co", "Relationship", "Dest_Taxa", "Dest_X", "Dest_Y"))
head(spearman_net1)

spearman_net2 <- spearman_distB2 %>%
  right_join(spearman2, by = "Taxa1") %>%
  mutate(
    Origin_Taxa = Taxa1,
    Origin_X = Axis.1, 
    Origin_Y = Axis.2,
    Taxa1b = Taxa2
  ) %>%
  select(c("Origin_Taxa", "Origin_X", "Origin_Y", "Taxa1b", "Spearman_Co", "Relationship")) %>%
  left_join(spearman_distB2, by = "Taxa1b") %>%
  mutate(
    Dest_Taxa = Taxa1b,
    Dest_X = Axis.1, 
    Dest_Y = Axis.2,
    Time = "2") %>%
  select(c("Origin_Taxa", "Time","Origin_X", "Origin_Y", "Spearman_Co", "Relationship", "Dest_Taxa", "Dest_X", "Dest_Y")) 
head(spearman_net2)

spearman_net3 <- spearman_distB3 %>%
  right_join(spearman3, by = "Taxa1") %>%
  mutate(
    Origin_Taxa = Taxa1,
    Origin_X = Axis.1, 
    Origin_Y = Axis.2,
    Taxa1b = Taxa2
  ) %>%
  select(c("Origin_Taxa", "Origin_X", "Origin_Y", "Taxa1b", "Spearman_Co", "Relationship")) %>%
  left_join(spearman_distB3, by = "Taxa1b") %>%
  mutate(
    Dest_Taxa = Taxa1b,
    Dest_X = Axis.1, 
    Dest_Y = Axis.2,
    Time = "3") %>%
  select(c("Origin_Taxa","Time", "Origin_X", "Origin_Y", "Spearman_Co", "Relationship", "Dest_Taxa", "Dest_X", "Dest_Y"))
head(spearman_net3)

spearman_net4 <- spearman_distB4 %>%
  right_join(spearman4, by = "Taxa1") %>%
  mutate(
    Origin_Taxa = Taxa1,
    Origin_X = Axis.1, 
    Origin_Y = Axis.2,
    Taxa1b = Taxa2
  ) %>%
  select(c("Origin_Taxa", "Origin_X", "Origin_Y", "Taxa1b", "Spearman_Co", "Relationship")) %>%
  left_join(spearman_distB4, by = "Taxa1b") %>%
  mutate(
    Dest_Taxa = Taxa1b,
    Dest_X = Axis.1, 
    Dest_Y = Axis.2,
    Time = "4") %>%
  select(c("Origin_Taxa", "Time", "Origin_X", "Origin_Y", "Spearman_Co", "Relationship", "Dest_Taxa", "Dest_X", "Dest_Y"))
head(spearman_net4)

spearman_net5 <- spearman_distB5 %>%
  right_join(spearman5, by = "Taxa1") %>%
  mutate(
    Origin_Taxa = Taxa1,
    Origin_X = Axis.1, 
    Origin_Y = Axis.2,
    Taxa1b = Taxa2
  ) %>%
  select(c("Origin_Taxa", "Origin_X", "Origin_Y", "Taxa1b", "Spearman_Co", "Relationship")) %>%
  left_join(spearman_distB5, by = "Taxa1b") %>%
  mutate(
    Dest_Taxa = Taxa1b,
    Dest_X = Axis.1, 
    Dest_Y = Axis.2,
    Time = "5") %>%
  select(c("Origin_Taxa", "Time", "Origin_X", "Origin_Y", "Spearman_Co", "Relationship", "Dest_Taxa", "Dest_X", "Dest_Y"))
head(spearman_net5)

spearman_net6 <- spearman_distB6 %>%
  right_join(spearman6, by = "Taxa1") %>%
  mutate(
    Origin_Taxa = Taxa1,
    Origin_X = Axis.1, 
    Origin_Y = Axis.2,
    Taxa1b = Taxa2
  ) %>%
  select(c("Origin_Taxa", "Origin_X", "Origin_Y", "Taxa1b", "Spearman_Co", "Relationship")) %>%
  left_join(spearman_distB6, by = "Taxa1b") %>%
  mutate(
    Dest_Taxa = Taxa1b,
    Dest_X = Axis.1, 
    Dest_Y = Axis.2,
    Time = "6") %>%
  select(c("Origin_Taxa","Time", "Origin_X", "Origin_Y", "Spearman_Co", "Relationship", "Dest_Taxa", "Dest_X", "Dest_Y"))
head(spearman_net6)

master_net <- spearman_net1 %>%
  full_join(spearman_net2) %>%
  full_join(spearman_net3) %>%
  full_join(spearman_net4) %>%
  full_join(spearman_net5) %>%
  full_join(spearman_net6)
#View(master_net)

# Plot network
ggplot(data = master_net, aes(x = Origin_X, y = Origin_Y))+
  facet_wrap(nrow = 2, ncol = 3, facets = vars(Time))+
  geom_segment(aes(x = Origin_X, y = Origin_Y, xend = Dest_X, yend = Dest_Y, color = Relationship))+
  geom_jitter(shape = 21, aes(size = abs(Spearman_Co), fill = Origin_Taxa), show.legend = FALSE)+
  #geom_text(aes(label = Origin_Taxa), check_overlap = TRUE, hjust = -0.05, 
           # color = "gray15", size = 3,
            #data = spearman_net1)+
  theme_classic()


# Another variation
ggplot(data = master_net, aes(x = Origin_X, y = Origin_Y))+
  facet_wrap(nrow = 2, ncol = 3, facets = vars(Time))+
  geom_segment(aes(x = Origin_X, y = Origin_Y, xend = Dest_X, yend = Dest_Y, color = Relationship))+
  geom_jitter(shape = 21, aes(fill = Origin_Taxa), show.legend = FALSE)+
  #geom_text(aes(label = Origin_Taxa), check_overlap = TRUE, vjust = -1, 
  #          color = "gray15", size = 3,
  #          data = filter(spearman_net, Origin_Taxa %in% spearman_labels))+
  theme_classic()


# Make heatmap to show correlations between taxa
ggplot(data = master_net)+
  facet_wrap(nrow = 2, ncol = 3, facets = vars(Time), shrink = TRUE)+
  geom_point(aes(x = Origin_Taxa, y = Dest_Taxa, fill = Relationship, size = abs(Spearman_Co)), shape = 21)+
  theme_bw()+
  theme(axis.text.y.left = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

ggplot(data = master_net)+
  #facet_wrap(nrow = 2, ncol = 3, facets = vars(Time), shrink = TRUE)+
  geom_tile(aes(x = Origin_Taxa, y = Dest_Taxa, fill = Relationship))+
  theme(axis.text.y.left = element_blank(),
        axis.text.x = element_blank())
