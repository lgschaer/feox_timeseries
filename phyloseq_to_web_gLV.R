# Load required packages
library(tidyverse)
library(phyloseq)

# Load rarefied phyloseq object
pso <- readRDS("/home/lgschaer/old/FeOx_demux/feox-phyloseq.rds")
pso                                           # View summary of phyloseq object

# Filter out eukaryotes and mitochondria
justbacteria <- pso %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   # Only include bacteria
      Family  != "mitochondria" &             # Filter out mitochondria
      Class   != "Chloroplast"                # Filter out chloroplasts
  )
justbacteria                                  # View summary of phyloseq object

# Convert to the right format for web-based gLV tool found at: https://web.rniapps.net/webglv/
pso_subset <- justbacteria %>%
  tax_glom(taxrank = "Family") %>%                     # Agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  mutate(
    Abundance = Abundance*100,                         # Multiply relative abundance by 100 to get whole numbers
    Abundance = as.numeric(as.integer(Abundance))      # Round values to nearest whole number
  ) %>%
  filter(Location == "GLRC" & Metal != "Swab") %>%     # Filter, keep only GLRC samples and remove all swabs
  unite(TaxID, Class:Family, sep = "_") %>%            # Join taxonomy columns Class through Family to get a single taxonomic identifier for each ASV
  select(TaxID, Time, Abundance) %>%                   # Select variables that will be in final spreadsheet
  pivot_wider(names_from = "TaxID", values_from = "Abundance", values_fn = list(Abundance = sum)) %>% # Mutate the table so the first column is time and the colnames are ASVs
  arrange(Time)%>%                                     # Arrange rows in chronological order
  column_to_rownames(var = "Time")                     # Convert the first column to row names
head(pso_subset)                                       # Confirm that everything is how it should be

# Save final table as tsv file for input into web-gLV tool
write.table(pso_subset, file = "/home/lgschaer/old/FeOx_demux/gLV_otu_table.tsv", row.names=TRUE, sep="\t")
