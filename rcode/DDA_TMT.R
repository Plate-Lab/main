############################## READ ME #########################################

# MS Proteomics Pipeline - DDA TMT-Labeled

#Version 1.1 - Released August 13th, 2024

#Written by Sarah Garcia, Lea Barny, Jake Hermanson - 2024

#Set Working Directory - unique to each system
setwd("path/to/your/data.csv OR data.tsv/file")

#Goal: Analyze Mass Spectrometry Data (TMT-Labeled DDA)
#       and perform statistical testing

#LOAD ALL FUNCTIONS AND LIBRARIES THEN START AT STEP ONE

############################ Load Tools ########################################

#BEFORE BEGINNING: Load ALL listed libraries and functions below:
#                    (Collapse functions for easier viewing)

#Load libraries
library(tidyverse)        #used in almost every function
library(data.table)       #used in almost every function
library(broom)            #Used for t.tests
library(ggfortify)        #Used for the PCA plot

### (OPTIONAL LIBRARIES) ### - Load for each optional step if desired
library(EnhancedVolcano)  #used for EnhancedVolcano plot
library(VennDiagram)      #used for making VennDiagrams
library(NbClust)          ####
library(clusterCrit)      # ^
library(factoextra)       #Used for k-means clustering and heatmap
library(cluster)          # V
library(fpc)              # V
library(pheatmap)         ####

#Load functions
load_file <- function(file){
  as.data.frame(fread(file, stringsAsFactors = FALSE))
}                         #works with .csv and .tsv; not .xlsx
tmt.dda.qc <- function(df){
  
  raw.data <- df %>% mutate(`Gene Symbol` = ifelse(`Gene Symbol` == "", NA, `Gene Symbol`))
  
  #If missing a gene name, fill in name with UniProt Accession number
  raw.data$"Gene Symbol" <- ifelse(is.na(raw.data$"Gene Symbol"),
                                   raw.data$Accession, raw.data$"Gene Symbol")
  
  #If there is a gene name listed in the description column
  if(any(str_detect(raw.data$Description, pattern = "GN=\\S+"))){  #if detected pattern 'GN='
    data.sorted <- raw.data %>%                                    #take our raw file
      dplyr::rename(Gene_Symbol = `Gene Symbol`) %>%
      select(Accession, Description, Gene_Symbol,                #Select relevant columns
             starts_with("Abundance:")) %>%
      mutate(Gene_Symbol = (str_extract(.$Description,             #replace name with pattern found in column
                                        pattern = "GN=\\S+"))) %>%
      mutate(Gene_Symbol = str_replace(string = .$Gene_Symbol,     #get rid of 'GN='
                                       pattern = "GN=", 
                                       replacement = "")) %>%
      pivot_longer(cols = starts_with("Abundance:"),               #pivot longer
                   names_to = "Sample", values_to = "Intensity") %>%
      mutate(Sample = str_replace_all(                             #remove 'Abundance:'
        Sample, pattern = "Abundance: ", replacement = "")) %>%
      mutate(Sample = str_replace_all(                             #remove 'Sample'
        Sample, pattern = "Sample, ", replacement = "")) %>%     
      mutate(Sample = str_replace_all(                             #Change ':' to ','
        Sample, pattern = ":", replacement = ","))
    
    max_parts <- max(str_count(data.sorted$Sample, ",")) + 1       #get max size of separated string @ ','
    
    data.sorted <- data.sorted %>%                                 #take our sorted data
      mutate(Run = Sample) %>%                                     #overwrite 'Run' column with 'Sample' col
      separate_wider_delim(cols = Sample, delim = ",",             #split column at every ','
                           names = paste0("X_", 1:(max_parts)),    ### make names dynamic to fit varying size data
                           too_few = "align_start") %>%
      dplyr::rename(MS_Run = "X_1", TMT = "X_2", Condition = "X_3") %>%   #Replace column names
      mutate(Protein.Info = paste0(Gene_Symbol,'*',Accession)) %>% #Combine columns and separate by '*'
      select(-Gene_Symbol, -Accession)                             #remove irrelevant columns
  }else{                                                            #if no pattern 'GN=' detected
    data.sorted <- raw.data %>%                                    #take our raw file
      dplyr::rename(Gene_Symbol = `Gene Symbol`) %>%
      select(Accession, Description, Gene_Symbol,                #select relevant columns
             starts_with("Abundance:")) %>%
      pivot_longer(                                                #pivot longer
        cols = starts_with("Abundance:"), 
        names_to = "Sample", values_to = "Intensity") %>%
      mutate(Sample = str_replace_all(                             #remove 'Abundance'
        Sample, pattern = "Abundance: ", replacement = "")) %>%
      mutate(Sample = str_replace_all(                             #remove 'Sample'
        Sample, pattern = "Sample, ", replacement = "")) %>%
      mutate(Sample = str_replace_all(                             #Change ':' to ','
        Sample, pattern = ":", replacement = ","))            
    
    max_parts <- max(str_count(data.sorted$Sample, ",")) + 1       #get max size of separate string @ ','
    
    data.sorted <- data.sorted %>%                                 #take our sorted data
      mutate(Run = Sample) %>%                                     #overwrite 'Run' column with 'Sample col
      separate_wider_delim(cols = Sample, delim = ",",             #split column at every ','  
                           names = paste0("X_", 1:(max_parts)),    ### makes names dynamic to fit varying size data
                           too_few = "align_start") %>% 
      dplyr::rename(MS_Run = "X_1", TMT = "X_2", Condition = "X_3") %>%   #Replace column names
      mutate(Protein.Info = paste0(Gene_Symbol,'*',Accession)) %>% #Combine columns and separate by '*'
      select(-Gene_Symbol, -Accession)                             #remove irrelevant columns
  }
  
  before <- length(unique(data.sorted$Protein.Info))               #store size for comparison
  
  raw.data.sorted <- data.sorted %>%                               #Take our almost-finished file
    unite("Condition", Condition, starts_with("X_"),               #Combine condition columns into single
          sep = ",", remove = TRUE) %>%
    filter(!str_detect(Protein.Info, 
                       regex("contaminant|cont|Contaminant|Cont")))#remove any proteins that are listed as contaminant
  
  after <- length(unique(raw.data.sorted$Protein.Info))            #store size for comparison
  
  removed <- (1-(after/before))*100                                #get size comparison
  
  cat("\nRemoved", removed, "% of the data due to:\n identified as contaminant.\n")
  
  raw.data.sorted <- raw.data.sorted %>% mutate_at(vars(2:ncol(.)), ~ str_replace_all(., "\\s+", "")) %>%
    mutate(Intensity = as.numeric(Intensity)) %>% filter(!is.na(Intensity))
  
}                    #import PD / TMT file and apply quality control filters
medNorm <- function(df, samples = Run, 
                    values = Precursor.Normalised) { 
  df %>% data.frame() %>%                                #take your specified data.frame
    mutate(global_median = median({{ values }},          #calculate the global median
                                  na.rm = TRUE)) %>% 
    group_by({{ samples }}) %>%                          #group by sample name
    mutate(x_median = median({{ values }},               #calculate each sample's median
                             na.rm=TRUE)) %>% 
    mutate(norm_factor = global_median / x_median) %>%   #calculate the norm_factor
    mutate(Norm_abundance = norm_factor * {{ values }}   #multiply values by the norm_factor
    ) %>%
    select(-global_median, -x_median, -norm_factor)      #remove irrelevant columns
} #medianNormalize data
info.group <- function(df){
  group.size <- data.frame(df) %>%                       #get group size (rows)
    group_by(Run) %>% group_size()
  group.name <- data.frame(df) %>%                       #get group labels
    group_by(Run) %>% group_keys()
  group.info <- bind_cols(group.name, group.size)        #create group info matrix
  names(group.info) <- c("Run", "Count")                 #rename column names
  group.info <- group.info                               #export file to global environment 
}                          #function to obtain sample sizes
filter.NA <- function(df, points = 3){
  
  test <- list()                                         #initialize empty list
  
  df_subset <- df %>% dplyr::select(Run, Protein.Info, Intensity)
  
  df.wide <- df_subset %>% data.frame() %>%
    pivot_wider(                         #take the data.frame, pivot wider
    names_from = 'Run', 
    values_from = 'Intensity',
    values_fn = list(Intensity = median) #if there are any list-cols present, the median will be taken
  )
  
  for(i in treatments){                                  #for each treatment specified
    sample.number <- df.wide %>%                         #take our wide data.frame
      dplyr::select(Protein.Info,                        #select treatment columns
                    ends_with(i)
      ) %>% ncol(.)-1                                    #store number of treatment columns
    NA.amount <- sample.number - points                  #number of treatment columns minus minimum number of points
    test[[i]] <- df.wide %>%                             #take our wide data.frame
      dplyr::select(Protein.Info,                        #select treatment columns  
                    ends_with(i)) %>% 
      filter(rowSums(is.na(.)) <= NA.amount)             #filter for columns that meet the criteria
  }
  
  combined_df <- Reduce(function(x, y) full_join(x, y, by = c("Protein.Info")), test)
  sample.amount <- ncol(combined_df) - 1
  combined_df_filter <- combined_df %>% unique() %>% filter(rowSums(is.na(.)) < sample.amount) %>%
    pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Intensity") %>%
    filter(Intensity != 0) %>%
    inner_join(df, by = c("Protein.Info", "Run", "Intensity"))
}               #filter for minimum # of points per observation
average.FC <- function(data){
  test <- list()
  data.wide <- data %>% pivot_wider(names_from = "Run", values_from = "Log2") %>%
    column_to_rownames(var = "Protein.Info")
  for(i in treatments){
    test[[i]] <- data.wide %>% mutate(AverageFC = rowMeans(dplyr::select(., ends_with(i)), na.rm = TRUE)) %>% 
      dplyr::select(AverageFC) %>% rownames_to_column() %>% setNames(c("Protein.Info", paste0("Avg.", i)))
  }
  Avg.FC.log2 <- test %>% purrr::reduce(full_join, by = "Protein.Info")
}                        #function to calculate treatment avg.Log2FC

### (OPTIONAL FUNCTIONS) ### - Load for each optional step if desired
volcano.curve <- function(result.volcano,                #function to plot curvature volcano
                          c.95 = 0.4,             #Define curvature
                          c.99 = 0.8,             #Define high-confidence curvature
                          x0 = 0.5,               #Define minimum fold change
                          adjpCutoff = 0.01){
  
  #adding curve cutoff
  curve_cutoff <- function(lfc, c, x0, sigma) {
    return(c * sigma / (lfc - x0))
  }
  
  mirrored_function <- function(x, c, x0, sigma) {
    ifelse(abs(x) > x0, curve_cutoff(abs(x), c, x0, sigma), NA)
  }
  
  # Calculate standard deviation of log2 fold changes
  sigma <- sd(result.volcano$estimate, na.rm = TRUE)
  
  # Calculate curve cutoff values
  result.volcano$curve_cutoff_95 <- curve_cutoff(abs(result.volcano$estimate), c.95, x0, sigma)
  result.volcano$curve_cutoff_99 <- curve_cutoff(abs(result.volcano$estimate), c.99, x0, sigma)
  result.volcano$significant_med <- -log10(result.volcano$p.adj) > result.volcano$curve_cutoff_95 &
    abs(result.volcano$estimate) >= x0
  result.volcano$significant_high <- -log10(result.volcano$p.adj) > result.volcano$curve_cutoff_99 &
    abs(result.volcano$estimate) >= x0
  
  #Calculate graph sizes
  min.x <- as.integer(min(result.volcano$estimate, na.rm = TRUE)) - 1
  max.x <- as.integer(max(result.volcano$estimate, na.rm = TRUE)) + 1
  
  min.y <- 0
  max.y <- as.integer(max(-log10(result.volcano$p.adj), na.rm = TRUE)) + 1
  
  # ggplot volcano plot with specific labels
  volcano_plot <<- ggplot(result.volcano, aes(x = estimate, y = -log10(p.adj)#, label = rownames(result.volcano)
  )) +
    geom_point(alpha = 0.2, colour = 'grey30') +
    stat_function(fun = function(x) mirrored_function(x, c = c.95, x0 = x0, sigma = sigma), color = "black") +
    stat_function(fun = function(x) mirrored_function(x, c = c.99, x0 = x0, sigma = sigma), color = "blue") +
    theme_bw() +
    labs(title = "Volcano Plot with Curve Cutoff Based on Standard Deviation",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted p-value") +
    scale_x_continuous(breaks = seq(min.x, max.x, by = 1), limits = c(min.x, max.x))+
    scale_y_continuous(breaks = seq(min.y, max.y, by = 1), limits = c(min.y, max.y)) +
    geom_point(data = subset(result.volcano, significant_med), aes(color = "med-conf"), size = 2) +
    geom_point(data = subset(result.volcano, significant_high), aes(color = "high-conf"), size = 2) +
    geom_text_repel(data = subset(result.volcano, significant_med), aes(label = Protein), size = 3)
  
  return(volcano_plot)
} #Define adj.p-value cutoff

############### Step One: Input and Reformat Data ##############################

#Load data files - .csv or .tsv (doesn't like .xlsx)
raw.data <- load_file("your.data.file.csv OR .tsv")

#Clean up data files
raw.data.qc <- tmt.dda.qc(raw.data)

#Store Sample Treatments
treatments <- unique(raw.data.qc$Condition)

#Optional Filtering (points = 1 for removal of only NA rows)
raw.data.qc.filtered <- filter.NA(raw.data.qc, points = 1)

### Visualization ###

#Boxplot 1 - Raw abundances across TMT labels grouped by treatment or TMT Label
ggplot(raw.data.qc.filtered, aes(x= TMT, y=Intensity, fill= Condition)) +
  geom_hline(yintercept = median(raw.data.qc.filtered$Intensity, na.rm = TRUE),  #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot() +
  labs(title= "Distribution of Values Grouped by Treatment", 
       x= "Sample",
       y= "Intensity") +
  scale_y_log10() +
  theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
  facet_grid(cols = vars(Condition), scales = "free_x")

#Same plot, grouped instead by TMT Label 
ggplot(raw.data.qc.filtered, aes(x= Condition, y=Intensity, fill= TMT)) +
  geom_hline(yintercept = median(raw.data.qc.filtered$Intensity, na.rm = TRUE), #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot()+ 
  labs(title= "Distribution of Values Grouped by TMT Label", 
       x= "Sample",
       y= "Intensity") +
  scale_y_log10() +
  theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
  facet_grid(cols = vars(TMT), scales = "free_x")

############## Step Two: Normalization #########################################

#Global Median Normalize (per Sample/TMT Label)
# Specify: samples = column_samplenames    values = column_values  
medNorm.data <- medNorm(raw.data.qc.filtered, samples = Run, values = Intensity)

#Log2 Transformation
medNorm.log2.data <- medNorm.data %>% mutate(Log2 = log2(Norm_abundance))

### Visualization ###

#Boxplot 2 - Normalized Log2 Abundances grouped by Treatment or TMT Label
ggplot(medNorm.log2.data, aes(x= TMT, y=Log2, fill= Condition)) +
  geom_hline(yintercept = median(medNorm.log2.data$Log2, na.rm = TRUE),  #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot() +
  labs(title= "Distribution of Values Grouped by Treatment", 
       x= "Sample",
       y= "Intensity") +
  theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
  facet_grid(cols = vars(Condition), scales = "free_x")

#Same plot, grouped instead by TMT Label 
ggplot(medNorm.log2.data, aes(x= Condition, y=Log2, fill= TMT)) +
  geom_hline(yintercept = median(medNorm.log2.data$Log2, na.rm = TRUE), #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot()+ 
  labs(title= "Distribution of Values Grouped by TMT Label", 
       x= "Sample",
       y= "Intensity") +
  theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
  facet_grid(cols = vars(TMT), scales = "free_x")  

#Barplot 1 - Unique Protein IDs Identified per Sample
protein.count <- medNorm.log2.data %>% na.omit() %>% info.group()
protein.count.info <-medNorm.log2.data %>% dplyr::select(Run, MS_Run) %>% 
  distinct() %>% inner_join(protein.count, by = "Run")

ggplot(protein.count.info, aes(x= Run, y=Count, fill = MS_Run)) +
  geom_col()+
  labs(title= "Unique Protein IDs Identified per Sample", 
       x= "Sample",
       y= "Count") +
  theme(axis.text.x= element_text(angle = 45, hjust =1))

#PCA Plot
data <- medNorm.log2.data %>% filter(is.finite(Log2)) %>%
  select(-Intensity, -Norm_abundance, -starts_with("X_"), -Description) %>%
  data.frame() %>% filter(grepl(".*\\*.*", Protein.Info)) %>% 
  pivot_wider(names_from = Protein.Info, values_from = Log2) %>%
  select(-MS_Run, -TMT, -Condition) %>% column_to_rownames(var = "Run")

#Replace any missing values with 1
data[is.na(data)] <- 1

#Store our sample names
sample_names <- data %>% rownames()

#Construct sample design matrix
sample_groups <- medNorm.log2.data %>% 
  dplyr::select(-Protein.Info, -Intensity, -Description, -Norm_abundance, -Log2) %>%
  unique()

#Run the PCA analysis
data.prcomp <- data %>% prcomp(scale = TRUE)

#Label our samples with their treatment groups
data.groupInfo <- data %>%
  rownames_to_column(var = "Run") %>%
  inner_join(sample_groups, by = "Run")

#Plot the PCA analysis - change 'colour' to grouping variable (MS_Run, Condition, etc.)
autoplot(data.prcomp, data=data.groupInfo, label= TRUE, size=4, colour = 'MS_Run') +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  labs(title = "A. Principal Component Analysis Plot")

############## Step Three: Statistics ##########################################

### Options: ###
# Multiple t.tests w/ p.adj(method = "fdr")  -  For comparison of 2 groups
# ANOVA / Mixed Model  -  For comparison of >= 3 groups

  # #If needed, use these lines to separate the "Condition" column into individual pieces for statistical testing
  # medNorm.log2.data.sep <- medNorm.log2.data %>% 
  #   separate_wider_delim(cols = Condition, names = c("Condition", "Condition_1"), delim = ",") #add more columns if needed

#Perform Multiple unpaired t.tests with equal variance - Visualize via Volcano Plot
comparison <- c("treatment", "control") #define in order of [treatment - control]

#Filter for specified comparison and perform t.tests
proteins.log2.long.stat <- medNorm.log2.data %>% data.frame() %>%
  dplyr::select(Run, Protein.Info, Condition, Log2) %>%
  filter(Condition %in% comparison) %>%
  mutate(Condition = factor(Condition, levels = comparison)) %>%
  group_by(Protein.Info, Condition) %>% filter(n() >= 2) %>% 
  
  #Note: Default is minimum 3 points present in both groups for t.test. [filter(n() >= 3)]
  #However, R allows a minimum of 2 points to output a p.value. [filter(n() >= 2)]
  #Change at your own discretion.
  
  ungroup() %>%
  group_by(Protein.Info) %>% filter(n_distinct(Condition) == 2) %>%
  do(tidy(t.test(Log2 ~ Condition, data = ., var.equal = TRUE)))

#Adjust p-values
result.volcano <- proteins.log2.long.stat %>% data.frame() %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  select(Protein.Info, estimate, p.adj) %>% 
  separate_wider_delim(cols = "Protein.Info", names = c("Protein", "Accession"), delim = "*") %>%
  data.frame()

    ###   Simple Volcano Plot - EnhancedVolcano package   ###
    EnhancedVolcano(result.volcano, lab = result.volcano$Protein, 
                    x = "estimate", y = "p.adj", pointSize = 3.0,
                    FCcutoff = 0.5,   #set fold change cutoff
                    pCutoff = 0.01,   #set p.adj cutoff
                    xlim = c((min(result.volcano$estimate)-0.5), abs(min(result.volcano$estimate)-0.5)),
                    ylim = c(0, (max(-log10(result.volcano$p.adj))+0.5))
    )
    
    #Filter for significant proteins - Enhanced Volcano
    sig.prt.EV <- result.volcano %>% filter(abs(estimate) > 0.5 & p.adj < 0.01) %>%
      rownames_to_column(var = "Protein.Info")
    
    ###                                                  ###
    
    ###   Volcano plot using curvature function          ###
    volcano.curve(result.volcano)
    
    #Filter for significant proteins - Curvature Volcano
    volcano_plot.data <- volcano_plot$data %>% 
      filter(significant_med == TRUE | significant_high == TRUE) %>%
      rownames_to_column(var = "Protein.Info")
    
    ###                                                 ###
    
    ###   Venn Diagram Comparing Both Volcano Plots ###
    venn.plot <- venn.diagram(
      x = list(Set1 = sig.prt.EV$Protein.Info, Set2 = volcano_plot.data$Protein.Info),
      category.names = c("EV", "Curve"),
      filename = NULL,
      output = FALSE,
      fill = c("blue", "yellow"),
      alpha = 0.5,
      cex = 1.5,
      cat.cex = 1.5,
      cat.col = c("black", "black")
    )
    
    grid.newpage()
    grid.draw(venn.plot)
    
    ###                                            ###

#Perform One-way ANOVA
  temp <- medNorm.log2.data %>% filter(is.finite(Log2)) %>%
    dplyr::select(Run, Protein.Info, Condition, Log2)
  
  anova_result <- aov(Log2 ~ Condition, data = temp)  
  
  summary(anova_result)
  
  tukey_result <- TukeyHSD(anova_result)
  
  tukey_df <- as.data.frame(tukey_result$Condition) 
  
  tukey_df_sig <- filter(tukey_df, `p adj` < 0.01 )

############### Additional Analysis ############################################

#Calculate Fold Changes normalized to user-defined control

#Define the control group
control <- treatments[1]  

#Get Control information
data.control <- medNorm.log2.data %>% 
  dplyr::select(Run, Protein.Info, Log2) %>%
  pivot_wider(names_from = "Run", values_from = "Log2") %>%
  dplyr::select(Protein.Info, ends_with(control)) %>% 
  column_to_rownames(var = "Protein.Info") %>%
  mutate(sd = apply(., 1, sd, na.rm = TRUE)) %>% 
  mutate(mean = rowMeans(dplyr::select(., -sd), na.rm = TRUE)) %>%
  dplyr::select(sd, mean) %>%
  rownames_to_column(var = "Protein.Info")

#Calculate Fold Changes for Original Data Frame
FC.log2 <- medNorm.log2.data %>% 
  dplyr::select(Run, Protein.Info, Log2) %>%
  inner_join(data.control, by = "Protein.Info") %>%
  mutate_at(vars(Log2), ~ (.-mean)) %>% 
  dplyr::select(-mean, -sd)

#Plot Fold Changes per sample
FC.log2 %>%
  inner_join(sample_groups, by = "Run") %>%
  ggplot(aes(x = Run, y = Log2, fill = Condition)) + 
  geom_boxplot(outliers = TRUE) + theme_bw() +
  labs(title = "A. Boxplot of Average Log2 Fold Changes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  theme(axis.text.x= element_text(angle = 45, hjust =1)) +
  facet_grid(cols = vars(Condition), scales = "free_x")

#This block calculates the average fold change of each protein per treatment group
Avg.FC.log2 <- FC.log2 %>% average.FC() %>% filter(rowSums(!is.na(.)) > 2) %>%
  separate_wider_delim(cols = Protein.Info, delim = '*', names = c("Protein", "UniProt.Id"))

#This block reformats the data into long format, then plots a boxplot of the average fold changes per treatment
Avg.FC.log2 %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "Treatment", values_to = "AvgFC") %>%
  ggplot(aes(x = Treatment, y = AvgFC, fill = Treatment)) + 
  geom_boxplot(outliers = TRUE) + theme_bw() +
  labs(title = "A. Boxplot of Average Log2 Fold Changes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value")

#Clustering (Using K-Means Clustering)

  #Copy the data
  clusterDF <- Avg.FC.log2
  
  #Replace NAs with zero
  clusterDF[is.na(clusterDF)] <- 0
  
  #Reformat
  clusterDF.mean <- clusterDF %>%
    unite(Protein.Info, Protein, UniProt.Id, sep = "*") %>% column_to_rownames(var = "Protein.Info")
  
  ######## Run this when deciding the number of kmeans clusters ###############
  #Change this to change the limit to amount of clusters you calculate for
  clusLength <- 20
  
  result <- fviz_nbclust(clusterDF.mean, kmeans, method = "wss", k.max = clusLength)
  print(result)
  elbow <- result$data
  
  #start for the silhouette score calculation
  silPoints <- 1
  df.silhouette <- data.frame(k = 2:clusLength)
  
  for(i in 1:20){  for(k in 2:clusLength){
    # Perform k-means clustering
    kmeans_result <- kmeans(clusterDF.mean, centers = k)
    # Compute silhouette information using silhouette function
    sil_info <- silhouette(kmeans_result$cluster, dist(as.matrix(clusterDF.mean)))
    silPoints[k-1] <- summary(sil_info)$avg.width
  }
    silhouetteDf <- data.frame(k = 2:clusLength, Score = silPoints)
    df.silhouette <- left_join(df.silhouette, silhouetteDf, by = 'k')
  }
  
  row_sds <- apply(df.silhouette, 1, sd)
  row_mean <- apply(df.silhouette, 1, median)
  
  df.silMean <- data.frame(k = df.silhouette$k, means = row_mean, stdev = row_sds)
  
  silhouttePlot <- ggplot(silhouetteDf, aes(x = k, y = Score)) +
    geom_line() +
    theme_classic()
  silhouttePlot
  
  silhouttePlot <- ggplot(df.silMean, aes(x = k, y = means)) +
    geom_line() +
    theme_classic() + labs(title = "Average of 100 Silhouette Scores")
  silhouttePlot
  #############################################################################
  
  #Change clusAmount to be the amount of clusters you decide on, based on the above plots
  clusAmount = 10 

  # Create the heatmap
  set.seed(52)   
  pheatmap.kmean<-pheatmap(clusterDF.mean, kmeans_k = clusAmount, cluster_rows = TRUE, cluster_cols = FALSE,
                           clustering_method = "complete",  # You can change the method as needed
                           cutree_rows = 1,  # Number of clusters for rows, estimated from dendrogram above
                           main = "Heatmap of Proteins with Hierarchical Clustering using K-means", 
                           labels_row = clusterDF.mean$Protein.Group,
                           color = hcl.colors(50, "Temps"), fontsize_row = 10, display_numbers = TRUE)
  
  cluster.values <- pheatmap.kmean$kmeans$cluster %>% data.frame() %>% rownames_to_column() %>%
    setNames(c("Protein.Info","Cluster"))
  
  clusterDF.mean.fixed <- clusterDF.mean %>% 
    rownames_to_column(var = "Protein.Info") %>% full_join(cluster.values, by = "Protein.Info") %>%
    separate_wider_delim(col = Protein.Info, names = c("Protein", "UniProt.AC"), delim = "*")
  
  #Export proteins with cluster assignment
  #write.csv(clusterDF.mean.fixed, "Avg.FC.log2.A.Kmean.csv")

