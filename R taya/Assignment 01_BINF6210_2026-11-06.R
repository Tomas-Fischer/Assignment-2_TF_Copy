##***************************
## Assignment 2
##
## TF Edits of TC's Assignment 1
##
## 2026-11-12
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
library(ggridges)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())

# Startup ends here

## _ Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# Code organization (## XXXXX -----)
# Justification for a section of code ## XXX
# Dead end analyses because it did not work, or not pursuing this line of inquiry (but leave it in as a trace of it, to potentially solve this issue, or avoid making the same mistake in the future # (>_<) 
# Solutions/results/interpretations (#==> XXX)
# Reference to manuscript pieces, figures, results, tables, ... # (*_*)
# TODO items #TODO
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)

## code for assignment part 1 ----

##### TF Edit rationale ----
#Edit 1: Deleting/merging redundant code (lists of non transposible elements,lists of retrotransposons,importing spreadsheets, removing redundant script from )

#Edit 2: I found figure 1 hard to read due to all of the datapoints, so I wanted to upgrade it by using geom_bin_2d, and have a separate chart for each TE type and species.

#Edit 3: I figured that it might be interesting to include multiple p

df_trich_sp <- read.csv("../Assignment-2_TF_Copy_Good/data taya/trichinella_spiralis_te.csv")
df_trich_sp

plot(sort(table(df_trich_sp$matching_class), decreasing = T)[1:10])

df_trich_t9 <- read.csv("../Assignment-2_TF_Copy_Good/data taya/GCA_001447505.1.csv")
df_trich_t9

plot(sort(table(df_trich_t9$matching_class), decreasing = T)[1:10])

# # filtering for only TEs
# df_trich_sp$matching_class |> unique()
# 
# non_te_sp = c("simple_repeat", "Simple_repeat", "Low_complexity", "Satellite", "rRNA", "snRNA", "ARTEFACT", "Unknown")
# 
# df_trich_s9$matching_class |> unique()
# non_te_s9 = c("Simple_repeat", "Satellite", "Low_complexity", "rRNA", "snRNA", "ARTEFACT", "Unknown")


##### TF Edit #1 ----
#we can simply make one list of non-transposible elements inclusive to both isolates to remove redundant code.

ls_non_te = c("simple_repeat", "Simple_repeat",
           "Low_complexity", "Satellite",
           "rRNA", "snRNA", "ARTEFACT",
           "Unknown") 

### The following section remains unchanged from the original script, besides changing the list of non-transposable elements to the common one i made above.

# visualizing the number of TEs per contig
df_trich_sp |> 
  group_by(qry_id) |> 
  summarise(count = n()) |> # this is the number of TEs per contig
  arrange(desc(count)) |> 
  filter(count < 1260) |> # this may not be correct, s9 has max count of 1259
  pull(count) |> 
  hist()

df_trich_t9 |> 
  group_by(qry_id) |> 
  summarise(count = n()) |> # this is the number of TEs per contig
  arrange(desc(count)) |> 
  pull(count) |> 
  hist()


# histogram of the percentage deletion

df_trich_sp |> 
  filter(matching_class %notin% ls_non_te) |>
  ggplot(aes(x = perc_del)) +
  geom_histogram()

df_trich_t9 |> 
  filter(matching_class %notin% ls_non_te) |>
  ggplot(aes(x = perc_del)) +
  geom_histogram()

# getting the TEs

df_te_sp <- df_trich_sp |> 
  filter(matching_class %notin% ls_non_te)

df_te_t9 <- df_trich_t9 |>
  filter(matching_class %notin% ls_non_te)

# getting DNA- and retro- transposons
df_te_sp$matching_class |> unique()

df_te_t9$matching_class |> unique()

# defining retrotransposons
# 
# retro_sp <- c("PLE/Naiad", "PLE", "LINE/CR1", "LTR/Pao", "LTR/Gypsy", "LINE/RTE-BovB",
#               "SINE/tRNA-RTE", "LTR/Ngaro", "SINE/tRNA-Deu-L2", "LTR/DIRS")
# retro_s9 <- c("LTR/Gypsy", "PLE", "PLE/Naiad", "LINE/CR1", "LTR/Pao", "LTR/Ngaro",
#               "SINE/tRNA-Deu-L2", "SINE/tRNA-RTE", "LINE/RTE-BovB", "LINE")

##### TF Edit #1 cont'd ----
#we can just make one single list of retrotransposons inclusive to both species instead of making 2 nearly identical lists.

ls_retro_all <- c("PLE/Naiad", "PLE", "LINE/CR1", "LTR/Pao", "LTR/Gypsy", "LINE/RTE-BovB", "SINE/tRNA-RTE", "LTR/Ngaro", "SINE/tRNA-Deu-L2", "LTR/DIRS", "LINE")

### TF - The below code was original. Below it I have made my implementations.
# 
# new_sp <- df_trich_sp |>
#   filter(matching_class %notin% non_te_sp) |>
#   mutate(class_te = if_else(matching_class %in% retro_sp, "retrotrans", "DNAtrans"))
# new_sp |>
#   ggplot(aes( x = perc_del, y = perc_div, colour = class_te)) +
#   geom_point() +
#   geom_smooth()
# 
# new_s9 <- df_trich_s9 |>
#   filter(matching_class %notin% non_te_s9) |>
#   mutate(class_te = if_else(matching_class %in% retro_s9, "retrotrans", "DNAtrans"))
# new_s9 |>
#   ggplot(aes( x = perc_del, y = perc_div, colour = class_te)) +
#   geom_point() +
#   geom_smooth()

# do TE families evolve differently
# compare  perc del vs perc div of DNA trans between species
# need to have a column with species

# df_all_te <- bind_rows(
#   new_sp |>
#     select(perc_div, perc_del, matching_class, class_te) |>
#     mutate(species = "spiralis"),
#   new_s9 |>
#     select(perc_div, perc_del, matching_class, class_te) |>
#     mutate(species = "T9"))
# 
# names(df_all_te)
# 
# ggplot(df_all_te, aes(x = perc_del, y = perc_div, colour = species)) +
#   geom_point() +
#   geom_smooth() +
#   facet_wrap(~class_te)

##### TF Edit #2 -----
#I found the original figure 1 hard to read due to the thousands of datapoints. So, i decided to use geom_bin_2d to bin the datapoints that are nearly identical together using a logarithmic scale to make the figure easier to interpret (in my opinion)

#Step 1 - filter out all of the non-transposible elements and define which are retrotransposons and which are DNA transposons.

df_newplot_density_sp <- df_trich_sp |> 
  filter(matching_class %notin% ls_non_te) |> 
  mutate(class_te = if_else(matching_class %in% ls_retro_all, "Retrotransposons", "DNA transposons"))

df_newplot_density_t9 <- df_trich_t9 |> 
  filter(matching_class %notin% ls_non_te) |> 
  mutate(class_te = if_else(matching_class %in% ls_retro_all, "Retrotransposons", "DNA transposons"))

#Step 2 - make plots for each isolate

df_newplot_density_sp |> 
  ggplot(aes(x = perc_del, y = perc_div)) +
  geom_bin_2d(bins=30) +
  geom_smooth(aes(colour = class_te)) +
  scale_fill_viridis_c(trans = "log10") +
  scale_color_viridis_d(option = "plasma") +
  facet_wrap(~ class_te)

df_newplot_density_t9 |> 
  ggplot(aes(x = perc_del, y = perc_div)) +
  geom_bin_2d(bins=30) +
  geom_smooth(aes(colour = class_te)) +
  scale_fill_viridis_c(trans = "log10") +
  scale_color_viridis_d(option = "plasma") +
  facet_wrap(~ class_te)

#Step 3 - merge the two dataframes into one, and define them by the species

df_all_te <- bind_rows(
  df_newplot_density_sp |>
    select(perc_div, perc_del, matching_class, class_te) |>
    mutate(species = "spiralis"),
  df_newplot_density_t9 |>
    select(perc_div, perc_del, matching_class, class_te) |>
    mutate(species = "T9"))

#Step 4 - Make a final plot comparing perc_del and perc_div between t. spirlis and T9, sorting between DNA transposons and retrotransposons.

ggplot(df_all_te, aes(x = perc_del, y = perc_div)) +
  geom_bin_2d(bins=30) +
  geom_smooth(aes(colour = species)) +
  scale_fill_viridis_c(trans = "log10") +
  scale_color_viridis_d(option = "plasma") +
  facet_wrap(species ~ class_te) +
  labs(title = "Comparing Percent Deletion and Percent Divergence in Transposible Elements in two Trichinella isolates")

#Step 5 - Save the figure!
ggsave("Figure1_TF.png", width = 10, height = 8)

# End of Edit 2

###Stats - unchanged

# trying to test it
dna_transposons <- df_all_te |>
  filter(class_te == "DNA transposons")

model <- lm(perc_div ~ 1 + perc_del + species + perc_del*species, data = dna_transposons)
summary(model)

# the species does have a significant effect on the divergence
model2 <- lm(perc_div ~ 1 + species, data = dna_transposons)
summary(model2)

# for retrotransposons
retro_transposons <- df_all_te |>
  filter(class_te == "Retrotransposons")

model3 <- lm(perc_div ~ 1 + perc_del + species + perc_del*species, data = retro_transposons)
summary(model3)
# the species does not have a significant effect on the divergence

model4 <- lm(perc_div ~ 1 + perc_del, data = retro_transposons)
summary(model4)
# the percentage deletion has a significant effect

# df_all_te |>
#   ggplot(aes(x = perc_div, fill = species)) +
#   geom_density(bins = 30, alpha = 0.6, position = "identity") +
#   scale_fill_manual(values = c("spiralis" = "steelblue", "T9" = "tomato")) +
#   labs(x = "Percent Divergence", y = "Count", fill = "Species")

##### TF Edit #3 ----
#Changing the figure from geom_density to geom_density ridges, and comparing retrotransposons and transposons.

#Step 1: Create a new column within our datafrage with merged species and te class to plot in separate ridges
df_all_te <- df_all_te |>  
  mutate(species_teclass = paste(species, class_te, sep = " ")) #|> 
  #view()

#Step 2: Plot a ridge graph, and calculate quartiles to easily compare the differences of TE composition with respect to percent divergence, between T9 and T. spiralis

df_all_te |> 
  ggplot(aes(x = perc_div, y = species_teclass, 
             fill = after_stat(quantile))) + #after_stat means that it will essentially colour the plots after the quartiles are calculated
  stat_density_ridges(quantile_lines =  FALSE, calc_ecdf = TRUE, 
                      geom = "density_ridges_gradient") + 
#we are computing the stats using calc_ecdf to determine the quartiles
#We also dont need quantile lines since we are colouring by the 4 quartiles already.
  scale_fill_viridis_d(name = "Quartile") +
  labs(x = "Percent Divergence", 
       y = "Transposon class per Species", 
       title = "Mapping the differences in Percent Divergence by Quartile in two Trichinella Isolates")

#Step 3: Save the figure

ggsave("Figure2_TF.png", height = 6, width = 10)

### End of Edit #3


common_te <- c("DNA/MULE-MuDR", "DNA/PIF-Harbinger", "DNA/TcMar-Tc1", "LINE/CR1", "LTR/Gypsy", "LTR/Ngaro", "LTR/Pao", "PLE", "PLE/Naiad")

df_common_te <- df_all_te |>
  filter(matching_class %in% common_te)
df_common_te

df_common_te |>
  ggplot(aes(x = matching_class, y = perc_div, fill = species)) +
  geom_boxplot() +
  scale_fill_viridis_d()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "TE Class", y = "Pecrent Divergence", fill = "Species") 

## code for assignment part 2 ----

#TE Composition

# Loading datasets

##### TF Edit 1 Contd -----
#This is already imported, so I replaced these dataframes with the original ones for above. We know they work and this is redundant.
#A lot of the code below was done to generate the first 3 figures, so I've commented out the redundant code and implemented the lists and dataframes that already exist.

# df_trichinella_sp <- read.csv("../Assignment-2_TF_Copy_Good/data taya/trichinella_spiralis_te.csv")
# 
# df_trichinella_t9 <- read.csv("../Assignment-2_TF_Copy_Good/data taya/GCA_001447505.1.csv")


# Getting to know about the datasets

# dim(df_trichinella_sp)
# dim(df_trichinella_t9)
# 
# 
# names(df_trichinella_sp)
# names(df_trichinella_t9)
# 
# head(df_trichinella_sp)
# head(df_trichinella_t9)
# 
# unique(df_trichinella_sp$matching_class)
# unique(df_trichinella_t9$matching_class)


# Filtering for TE annotations

# non_te <- c("Simple_repeat",
#             "simple_repeat",
#             "Low_complexity",
#             "Satellite",
#             "Unknown",
#             "rRNA",
#             "snRNA",
#             "ARTEFACT")


# df_sp_te <- df_trichinella_sp |>
#   filter(!(matching_class %in% non_te))
# 
# 
# df_t9_te <- df_trichinella_t9 |>
#   filter(!(matching_class %in% non_te))

# Adding dataset names
# 
# df_sp_te <- df_sp_te |>
#   mutate(dataset = "T. spiralis")
# 
# 
# df_t9_te <- df_t9_te |>
#   mutate(dataset = "T9")
# 
# 
# # Combining datasets
# 
# df_all_tes <- bind_rows(df_sp_te, df_t9_te)


#TF Edit - We already have df_all_te from above, which with all of the transformations, is functionally identical to this code used above. To simplify, I simply use df_all_te, which was used extensively above.
#Note that I have changed dataset to species for the following lines of code.

# Checking TE annotation counts

### I have made small changes to the code below. The scructure is from the original script, but the lists and dataframes have been replaced by onces already used earlier.

dim(df_trich_sp)
dim(df_trich_t9)


names(df_trich_sp)
names(df_trich_t9)

head(df_trich_sp)
head(df_trich_t9)

unique(df_trich_sp$matching_class)
unique(df_trich_t9$matching_class)
table(df_all_te$species)

# Counting TE families in each dataset

te_composition <- df_all_te |>
  count(species, matching_class)

te_composition

# Calculating TE family proportions

te_composition <- te_composition |>
  group_by(species) |>
  mutate(proportion = n / sum(n),
         percentage = proportion * 100) |>
  ungroup()

te_composition

# Visualizing TE family composition

te_composition |>
  ggplot(aes(x = matching_class,
             y = proportion,
             colour = species)) +
  geom_point() +
  theme_dark() + #theme_dark was chosen to make the light viridis T9 points stand out - TF
  scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "TE Family/Class",
       y = "Proportion of TE Annotations",
       colour = "Species")

# Ranking TE families by abundance

te_composition |>
  arrange(species, n)
# Gives us from the lowest abundance to the highest

# Top five TE families in T. spiralis

te_composition |>
  filter(species == "spiralis") |>
  arrange(n) |>
  tail(5) |>
  summarise(top_five_percentage = sum(percentage))
# Gives us the combined percentage of top five, %88.5

# Top five TE families in T9

te_composition |>
  filter(species == "T9") |>
  arrange(n) |>
  tail(5) |>
  summarise(top_five_percentage = sum(percentage))
# Gives us the combined percentage of top five, %91.7
# test
