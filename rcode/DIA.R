############################### READ ME #########################################

# MS Proteomics Pipeline - DIA Unlabeled - Written for DIA-NN report.tsv files ONLY
  #DIA-NN:  Demichev, V., et al. (2020). https://doi.org/10.1038/s41592-019-0638-x

#Written by Sarah Garcia, Lea Barny, Jake Hermanson - 2024

#Set Working Directory - unique to each system
setwd("path/to/your/report.tsv/file")

#Goal: Analyze Mass Spectrometry Data (Un-labeled DIA)
#       and perform statistical testing

#LOAD ALL FUNCTIONS AND LIBRARIES THEN START AT STEP ONE

############################ Load Tools ########################################

#BEFORE BEGINNING: Load ALL listed libraries and functions below:
#                    (Collapse functions for easier viewing)

#Load libraries
library(tidyverse)        #used in almost every function
library(data.table)       #used in almost every function
library(diann)            #DIA only (peptide to protein roll-up)
library(broom)            #Used for t.tests
library(ggfortify)        #Used for the PCA plot

#Load functions
load_file <- function(file){
  as.data.frame(fread(file, stringsAsFactors = FALSE))
}                         #works with .csv and .tsv; not .xlsx
diann.dia.qc <- function(raw.data, cont.name){
  
  chomp.contaminants<-function(file){  
    
    fasta.file <- fread(file, stringsAsFactors = FALSE, header = FALSE)
    
    i=0                                             #Initialize empty objects
    flag=0
    y<-data.frame()
    protein.size<-list()
    sequence.all<-list()
    alist<-list()
    b<-''
    c<-''
    
    for(i in fasta.file$V1){                      #for each line of inputted proteomeFASTA
      if(str_detect(i,"^>")){                     #if row starts with >
        if(flag==1){                              #flag check to concatenate sequence lines
          sequence<-str_c(alist,collapse = "")    #collapse
          aasize<-nchar(sequence)                 #get amino acid size
          protein.size<-rbind.data.frame(protein.size,aasize)  #store size
          sequence.all<-rbind.data.frame(sequence.all,sequence)  #store sequence
          alist<-list()                           #empty temp variable
          flag=0                                  #reset flag
        }
        flag=0                                    #reset
        pattern.species<-paste0("OS=\\S+")        #make species name pattern
        pattern.gn<-paste0("GN=\\S+")             #make gene name pattery
        pattern.type<-paste0("^>\\S+")            #make >tr or >sp pattern 
        a<-str_extract(i,pattern = pattern.species) #find and store species name (user supplies species name)
        b<-str_extract(i,pattern = pattern.gn)    #find and store gene name
        c<-str_extract(i,pattern = pattern.type)  #find and store pattern type (>tr or >sp, should be only two options)
        abc<-c(a,b,c)                             #make above identifiers into easy to bind object
        y<-rbind(y,abc)                           #rbind to fill rows of new df
      }
      else{
        flag=1                                    #Raise the flag!
        b<-str_trim(i,side = c("right"))          #Cut off new line character
        alist[[i]]<-b                             #Storing all sequence lines of current protein
      }
    }
    if(flag==1){                                  #this catches the last protein's sequence VVV
      sequence<-str_c(alist,collapse = "")        #
      aasize<-nchar(sequence)                     #
      protein.size<-rbind.data.frame(protein.size,aasize) #
      sequence.all<-rbind.data.frame(sequence.all,sequence) #
      alist<-list()                               #
      flag=0                                      #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    }
    
    y[,4]<-protein.size                             #add sizes to df
    y[,5]<-sequence.all                             #add sequences to df
    colnames(y)<-c("Genus","Gene.Name","UniProt Type","AA Length","Sequence") #add colnames to object 'y'
    y$Genus <- str_replace_all(y$Genus, "OS=", "")  #remove OS= from species column
    y$`Gene.Name` <- str_replace_all(y$`Gene.Name`, "GN=", "")  #remove GN= from gene name column
    input.fasta<-y
    
    fasta.info <- input.fasta %>%
      separate(`UniProt Type`, into = c("discard", "UniProt.AC","discard2"), sep = "\\|") %>%
      dplyr::select(-discard,-discard2) %>% filter(startsWith(UniProt.AC, "Cont")) %>%
      separate(UniProt.AC, into = c("cat", "UniProt.AC"), sep = "_") %>%
      dplyr::select(-cat)
    
  }  #function written to reformat .fasta files
  contaminants.info <- chomp.contaminants(file = cont.name) #import DIA contaminant file
  
  #Add list of contaminants to be removed from data before protein quantification: 
  contaminants <- contaminants.info$UniProt.AC
  
  cat("Cleaning up raw data...\n")                  #print message to console
  
  before <- length(unique(raw.data$Protein.Group))  #store size for comparison
  
  # #Data Clean-up - Parameters set at: Q.value <= 0.01; 1+ Peptide per Protein
  # raw.data.filtered <- raw.data %>%                                       #Takes our imported raw data
  #   group_by(Protein.Group) %>%                                           #Sort peptides into groups based on sequence
  #   filter(sum(grepl(Protein.Group[1], contaminants)) == 0) %>%           #filter out groups that match a given contaminant
  #   filter(Global.Q.Value <= 0.01) %>% filter(Q.Value <= 0.01) %>%        #Filter for qvalues
  #   group_by(Run, Protein.Group) %>%                                      #group by Run and Protein.Group
  #   filter(length(unique(paste0(Modified.Sequence, Precursor.Charge)))>1) #Filters for lengths > 1 of pasted sequence 
  # 
  
  #Data Clean-up - Parameters set at: Q.value <= 0.01; 1+ Peptide per Protein
  raw.data.filtered <- raw.data %>%                                       #Takes our imported raw data
    group_by(Protein.Group) %>%                                           #Sort peptides into groups based on sequence
    filter(!str_detect(Protein.Group, "^Cont_")) %>%                      #filter out groups that match a given contaminant
    filter(Global.Q.Value <= 0.01) %>% filter(Q.Value <= 0.01) %>%        #Filter for qvalues
    group_by(Run, Protein.Group) %>%                                      #group by Run and Protein.Group
    filter(length(unique(paste0(Modified.Sequence, Precursor.Charge)))>1) #Filters for lengths > 1 of pasted sequence 
  
  after <- length(unique(raw.data.filtered$Protein.Group))                #store size for comparison
  
  removed <- (1-(after/before))*100                                       #get size comparison
  
  cat("\nRemoved", removed, "% of the data due to q.value > 0.01 \n and/or only 1 peptide per protein observed \n and/or determined to be a contaminant.\n")
  
  export.file <- raw.data.filtered %>%                                    #take our qc filtered data                   
    dplyr::select(Run, Protein.Group, Genes, Stripped.Sequence,           #select relevant columns
                  Precursor.Id, Precursor.Normalised) %>%                 
    tidyr::unite(Protein.Info, Genes, Protein.Group, sep = "*")                  #combine two columns, separated by '*'
  
}       #import diann file and contaminants file; apply quality control filters
assign.group <- function(df, groups = treatments){
  
  sample.names <- unique(df$Run)                                     #store unique sample names
  
  temp <- sample.names %>% data.frame(Run = .) %>%                   #make names a data.frame
    arrange(Run) %>%                                                 #arrange names alphabetically
    mutate(group = NA, sample = NA) %>% data.frame()                 #make two empty columns
  
  cat("Groups:", paste0(groups[1:length(groups)]))                   #print message to console
  
  for(i in 1:nrow(temp)){                                            #for each row in 'temp'
    user_input <- readline(prompt = paste0(                          #prompt the user to give
      "Assign ", temp[i,"Run"], " to group [1 - ",                   ## each sample a group number
      length(groups), "] : "))
    temp[i,"group"] <- groups[as.numeric(user_input)]                #store the sample's group assignment
  }
  
  name.fix <- function(df){                                          #function definition
    df <- df %>% data.frame()                                        #make df a data.frame
    max_parts <- max(str_count(df[ ,1], "_"))                        #count max pieces after separation
    separated_df <- df %>%                                           #take data.frame
      separate_wider_delim(col = 1,                                  #separate wider
                           names = paste0("X_", 1:(max_parts + 1)),  ### make names dynamic for varying sizes
                           delim = "_", too_few = "align_start") %>% 
      dplyr::select(where(~length(unique(.)) > 1)) %>%               #select columns with at least 2 different entries
      unite(Run, starts_with("X_"), sep = "_", na.rm = TRUE) %>%     #combine Run column with columns starting with 'X_'
      arrange(Run)                                                   #arrange by the 'Run' column alphabetically
  } 
  
  temp[ ,"sample"] <- name.fix(temp$Run)                             #run previously define function on sample names
  
  sample.info <- temp %>% mutate(                                    #take our almost finished data.frame
    sample = paste0(temp$sample, '_', group)                         #combine sample name with associated treatment group
  )
  
}   #Assign treatment groups
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
max_LFQ <- function(df, values = "Norm_abundance") {
  diann_maxlfq(df, sample.header = "Run", group.header = "Protein.Info",
               id.header = "Precursor.Id", quantity.header = values,
  )
} #simpler wrapper for diann_maxLFQ
filter.NA <- function(df, points = 3){
  
  test <- list()                                         #initialize empty list
  
  df <- df %>% rownames_to_column(var = "Protein.Info")
  
  for(i in treatments){                                  #for each treatment specified
    sample.number <- df %>%                              #take our wide data.frame
      dplyr::select(Protein.Info, ends_with(i)) %>%      #select treatment columns
      ncol(.) - 1                                       #store number of treatment columns
    NA.amount <- sample.number - points                  #number of treatment columns minus minimum number of points
    test[[i]] <- df %>%                                  #take our wide data.frame
      dplyr::select(Protein.Info, ends_with(i)) %>%      #select treatment columns  
      filter(rowSums(is.na(.)) <= NA.amount)             #filter for columns that meet the criteria
  }
  
  combined_df <- Reduce(function(x, y) full_join(x, y, by = c("Protein.Info")), test)
  sample.amount <- ncol(combined_df) - 1
  combined_df_filter <- combined_df %>% unique() %>% filter(rowSums(is.na(.)) < sample.amount)
}               #filter for minimum # of points per observation
info.group <- function(df){
  group.size <- data.frame(df) %>%                       #get group size (rows)
    group_by(Run) %>% group_size()
  group.name <- data.frame(df) %>%                       #get group labels
    group_by(Run) %>% group_keys()
  group.info <- bind_cols(group.name, group.size)        #create group info matrix
  names(group.info) <- c("Run", "Count")                 #rename column names
  group.info <- group.info                               #export file to global environment 
}                          #function to obtain sample sizes
average.FC <- function(data){
  test <- list()
  for(i in treatments){
    test[[i]] <- data %>% mutate(AverageFC = rowMeans(dplyr::select(., ends_with(i)), na.rm = TRUE)) %>% 
      dplyr::select(AverageFC) %>% rownames_to_column() %>% setNames(c("Protein.Info", paste0("Avg.", i)))
  }
  Avg.FC.log2 <- test %>% purrr::reduce(full_join, by = "Protein.Info")
}                        #function to calculate treatment avg.Log2FC

############### Step One: Input and Reformat Data ##############################

#Load data files - .csv or .tsv (doesn't like .xlsx)
raw.data <- load_file("your.diann.report.tsv")

    #For contaminants, recommend fasta file provided by this paper:
    # Frankenfield, A., et al. (2022). https://doi.org/10.1021/acs.jproteome.2c00145

#Clean up data files
raw.data.qc <- diann.dia.qc(raw.data, cont.name = "your.contaminants.file.fasta")

#Define treatment groups
treatments <- c("treatmentA", "treatmentB", "treatmentC", "etc.")

#Generate treatment groups - since A and B are from one data file, only need to run this once for both species
sample.groups <- assign.group(raw.data.qc, groups = treatments)

#Replace sample names with reformatted ones
raw.data.qc.f <- raw.data.qc %>% inner_join(sample.groups, by = "Run") %>%
  mutate(Run = sample) %>% select(-sample)

    ### (OPTIONAL) Export peptide groups ###
    peptides.raw.qc.f <- raw.data.qc.f %>% dplyr::select(Run, Protein.Info, Precursor.Id, Precursor.Normalised) %>%
      pivot_wider(names_from = Run, values_from = Precursor.Normalised)
    
    #write.csv(peptides.raw.qc, "peptides.raw.qc.csv")

### Visualization ###

#Boxplot 1 - A. Raw abundances across samples grouped by treatment
ggplot(raw.data.qc.f, aes(x= Run, y=Precursor.Normalised, fill = group)) +
  geom_hline(yintercept = median(raw.data.qc.f$Precursor.Normalised, na.rm = TRUE),  #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot()+
  labs(title= "A. Distribution of Values Grouped by Treatment", 
       x= "Sample",
       y= "Abundance") +
  scale_y_log10() +
  theme(axis.text.x= element_text(angle = 45, hjust =1)) +
  facet_grid(cols = vars(group), scales = "free_x")

############## Step Two: Normalization #########################################

#Global Median Normalize (per Sample/TMT Label)
# Specify: samples = column_samplenames    values = column_values
peptides.medNorm <- medNorm(raw.data.qc.f, samples = Run, values = Precursor.Normalised)

#Protein Roll-up using maxLFQ algorithm
proteins.medNorm <- max_LFQ(peptides.medNorm, values = "Norm_abundance")

#Log2 Transformation
proteins.medNorm.log2 <- proteins.medNorm %>% data.frame() %>% log2()

##Prepare the data for plotting - Run only one of the two lines below
  #Run if sample names start with a number in sample.groups
  sample.groups.v2 <- sample.groups %>% mutate(Run = paste0("X", sample)) %>% select(-sample)
  
  #Run if sample names start with a letter in sample.groups
  sample.groups.v2 <- sample.groups %>% mutate(Run = paste0(sample)) %>% select(-sample)
  ###

      ### (OPTIONAL) Bait Normalization ###
      
      #Define bait using UniProt Accession ID (Ex. Human GAPDH)
      bait <- "your.bait.UniProt.accession.ID"
      
      #Retrieve bait protein information
      bait.quant <- proteins.medNorm.log2 %>% rownames_to_column(var = "Protein.Info") %>%
        separate_wider_delim(cols = "Protein.Info", delim = '*', names = c("Protein", "Accession")) %>%
        filter(Accession %in% bait) %>% 
        pivot_longer(cols = 3:ncol(.), names_to = "Run", values_to = "bait") %>% 
        dplyr::select(Run, bait)
      
      #Barplot - Visualize Bait Levels Across Samples
      bait.quant %>% inner_join(sample.groups.v2, by = "Run") %>%
        ggplot(aes(x= Run, y=bait, fill = group)) +
        geom_col()+
        labs(title= "A. Bait Protein Levels", 
             x= "Sample",
             y= "Abundance") +
        theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
        facet_grid(cols = vars(group), scales = "free_x")
      
      #Normalize to bait protein levels
      proteins.log2.baitNorm <- proteins.medNorm.log2 %>% rownames_to_column(var = "Protein.Info") %>%
        pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2" ) %>%
        inner_join(bait.quant, by = "Run") %>% mutate(Log2norm = Log2 - bait) %>%
        dplyr::select(Protein.Info, Run, Log2norm) %>%
        pivot_wider(names_from = "Run", values_from = "Log2norm") %>%
        column_to_rownames(var = "Protein.Info")
      
      boxplot.proteins.log2 <- proteins.log2.baitNorm %>%
        rownames_to_column(var = "Protein.Info") %>%
        pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2_medNorm") %>%
        inner_join(sample.groups.v2, by = "Run") %>% na.omit()

      #Boxplot 2 - A.Normalized Log2 Abundances grouped by Treatment
      ggplot(boxplot.proteins.log2, aes(x= Run, y=Log2_medNorm, fill = group)) +
        geom_hline(yintercept = median(boxplot.proteins.log2$Log2_medNorm, na.rm = TRUE), #this adds a line at the global median
                   color = "red4", linetype = "solid", size = 1) +
        geom_boxplot()+
        labs(title= "A. Distribution of Values Grouped by Treatment", 
             x= "Sample",
             y= "Abundance") +
        theme(axis.text.x= element_text(angle = 45, hjust =1)) +
        facet_grid(cols = vars(group), scales = "free_x")

### Visualization ### - If bait norm was used, change 'proteins.medNorm.log2' to 'proteins.log2.baitNorm'
      
boxplot.proteins.log2 <- proteins.medNorm.log2 %>% 
  rownames_to_column(var = "Protein.Info") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2_medNorm") %>%
  inner_join(sample.groups.v2, by = "Run") %>% na.omit()

#Boxplot 2 - A.Normalized Log2 Abundances grouped by Treatment
ggplot(boxplot.proteins.log2, aes(x= Run, y=Log2_medNorm, fill = group)) +
  geom_hline(yintercept = median(boxplot.proteins.log2$Log2_medNorm, na.rm = TRUE), #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot()+
  labs(title= "A. Distribution of Values Grouped by Treatment", 
       x= "Sample",
       y= "Abundance") +
  theme(axis.text.x= element_text(angle = 45, hjust =1)) +
  facet_grid(cols = vars(group), scales = "free_x")

#Barplot 1 - Unique Protein IDs Identified per Sample
protein.count <- proteins.medNorm.log2 %>% 
  rownames_to_column(var = "Protein.Info") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>% 
  na.omit() %>% info.group() %>% inner_join(sample.groups.v2, by = "Run") %>% na.omit()

ggplot(protein.count, aes(x= Run, y=Count, fill = group)) +
  geom_col()+
  labs(title= "A. Unique Protein IDs Identified per Sample", 
       x= "Sample",
       y= "Count") +
  theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
  facet_grid(cols = vars(group), scales = "free_x")

#PCA Plot
  #Copy our data
  data <- proteins.medNorm.log2
  
  #Replace all missing values with 1
  data[is.na(data)] <- 1
  
  #Transpose the data
  data.t <- data %>% t() %>% data.frame()
  
  #Store our sample names
  sample_names <- data.t %>% rownames()
  
  #Run the PCA analysis
  data.prcomp <- data.t %>% prcomp(scale = TRUE)
  
  #Label our samples with their treatment groups
  data.groupInfo <- data.t %>%
    rownames_to_column(var = "Run") %>%
    inner_join(sample.groups.v2, by = "Run")
  
  #Plot the PCA analysis
  autoplot(data.prcomp, data=data.groupInfo, label= FALSE, size=4, colour = 'group') +
    theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
    labs(title = "A. Principal Component Analysis Plot")

############## Step Three: Statistics ##########################################

#Protein Filtering - If no NA filtering is needed, set points = 1
#Filters for proteins with at least x observations per treatment group
proteins.filtered <- filter.NA(proteins.medNorm.log2, points = 3)

### Options: ###
# Multiple t.tests w/ p.adj(method = "fdr")  -  For comparison of 2 groups
# ANOVA / Mixed Model  -  For comparison of >= 3 groups

#Perform Multiple unpaired t.tests with equal variance - Visualize via Volcano Plot [Ex. DIA]
comparison <- c("treatment", "control") #define in order of [treatment - control]

proteins.log2.long.stat <- proteins.filtered %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>%
  na.omit() %>% inner_join(sample.groups.v2, by = "Run") %>%
  filter(group %in% comparison) %>%
  mutate(group = factor(group, levels = comparison)) %>%
  group_by(Protein.Info, group) %>% filter(n() > 2) %>% ungroup() %>%
  group_by(Protein.Info) %>% filter(n_distinct(group) == 2) %>%
  do(tidy(t.test(Log2 ~ group, data = ., paired = FALSE, var.equal = TRUE)))

#Adjust p-values
result.volcano <- proteins.log2.long.stat %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr", n = length(p.value))) %>%
  select(Protein.Info, estimate, p.adj) %>% 
  separate_wider_delim(cols = "Protein.Info", delim = '*', names = c("Protein", "Accession")) %>%
  select(-Accession) %>% mutate(Protein = paste0(1:nrow(.), '-', Protein)) %>%
  column_to_rownames(var = 'Protein') %>% data.frame()

#Simple Volcano Plot - EnhancedVolcano package
library(EnhancedVolcano)

proteins.log2.long.stat %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr", n = length(p.value))) %>%
  select(Protein.Info, estimate, p.adj) %>% 
  separate_wider_delim(cols = "Protein.Info", delim = '*', names = c("Protein", "Accession")) %>%
  select(-Accession) %>% mutate(Protein = paste0(1:nrow(.), '-', Protein)) %>%
  column_to_rownames(var = 'Protein') %>% data.frame() %>%
  EnhancedVolcano(lab = rownames(.), x = "estimate", y = "p.adj", pointSize = 3.0,
                  FCcutoff = 0.5,  #set fold change cutoff
                  pCutoff = 0.01   #set p.adj cutoff
                  )

  #Volcano plot using curvature function
  volcano.curve <- function(result.volcano, 
                            c.95 = 1.96,        #Define curvature
                            c.99 = 2.576,       #Define high-confidence curvature
                            x0 = 0.5,           #Define minimum fold change
                            adjpCutoff = 0.01   #Define adj.p-value cutoff
                            ){
    
                            #adding curve cutoff
                            curve_cutoff <- function(lfc, c, x0, sigma) {
                              return(c * sigma / (lfc - x0))
                            }
                            
                            mirrored_function <- function(x, c, x0, sigma) {
                              ifelse(abs(x) > x0, curve_cutoff(abs(x), c, x0, sigma), NA)
                            }
                            
                            # Calculate standard deviation of log2 fold changes
                            sigma <- sd(result.volcano$estimate, na.rm = TRUE)
                            
                            # # Define constants for the curve
                            # c.95 <- 1.96 # Adjust based on your confidence interval (e.g., 1.96 for 95%)
                            # c.99 <- 2.576 # Adjust based on your confidence interval (e.g., 2.576 for 99%)
                            # x0 <- 0.5 # Minimum fold change for enriched proteins
                            # adjpCutoff <- 0.01
                            
                            # Calculate curve cutoff values
                            result.volcano$curve_cutoff_95 <- curve_cutoff(abs(result.volcano$estimate), c.95, x0, sigma)
                            result.volcano$curve_cutoff_99 <- curve_cutoff(abs(result.volcano$estimate), c.99, x0, sigma)
                            result.volcano$significant_med <- -log10(result.volcano$p.adj) > result.volcano$curve_cutoff_95 &
                              abs(result.volcano$estimate) >= x0
                            result.volcano$significant_high <- -log10(result.volcano$p.adj) > result.volcano$curve_cutoff_99 &
                              abs(result.volcano$estimate) >= x0
                            
                            #Calculate graph sizes
                            min.x <- as.integer(min(result.volcano$estimate)) - 1
                            max.x <- as.integer(max(result.volcano$estimate)) + 1
                            
                            min.y <- 0
                            max.y <- as.integer(max(-log10(result.volcano$p.adj))) + 1
                            
                            # ggplot volcano plot with specific labels
                            volcano_plot <- ggplot(result.volcano, aes(x = estimate, y = -log10(p.adj)#, label = rownames(result.volcano)
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
                              geom_point(data = subset(result.volcano, significant_med), aes(color = "95% C.I."), size = 2) +
                              geom_point(data = subset(result.volcano, significant_high), aes(color = "99% C.I."), size = 2)
                            
                            return(volcano_plot)
  }
  
  volcano.curve(result.volcano)


#Perform One-way ANOVA - [Ex. DIA]

temp <- proteins.filtered %>% 
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>%
  inner_join(sample.groups.v2, by = "Run") %>%
  filter(is.finite(Log2))

anova_result <- aov(Log2 ~ group, data = temp)  

summary(anova_result)

tukey_result <- TukeyHSD(anova_result)

tukey_df <- as.data.frame(tukey_result$group)

tukey_df_sig <- filter(tukey_df, `p adj` < 0.01 )

############### Additional Analysis ############################################

#Step Seven: Calculate Fold Changes normalized to user-defined control

#Define the control group
control <- treatments[1]  

#Get Control information
data.control <- proteins.medNorm.log2 %>% dplyr::select(ends_with(control)) %>% 
  mutate(sd = apply(., 1, sd, na.rm = TRUE)) %>% mutate(mean = rowMeans(dplyr::select(., -sd), na.rm = TRUE))

#Calculate Fold Changes for Original Data Frame
FC.log2 <- proteins.medNorm.log2 %>% mutate(Mean = data.control$mean, StDev = data.control$sd) %>%
  mutate_at(vars(-Mean, -StDev), ~ (.-Mean)) %>% dplyr::select(-Mean, -StDev) %>% filter(rowSums(!is.na(.)) > 0)

#Plot Fold Changes per sample
FC.log2 %>% rownames_to_column(var = "Protein.Info") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2FC") %>%
  inner_join(sample.groups.v2, by = "Run") %>%
  ggplot(aes(x = Run, y = Log2FC, fill = group)) + 
  geom_boxplot(outliers = TRUE) + theme_bw() +
  labs(title = "A. Boxplot of Average Log2 Fold Changes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  theme(axis.text.x= element_text(angle = 45, hjust =1)) +
  facet_grid(cols = vars(group), scales = "free_x")

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

#Step Eight: Clustering (Using K-Means Clustering)
library(NbClust)
library(clusterCrit)
library(factoextra)
library(cluster)
library(fpc)
library(pheatmap)

clusterDF <- Avg.FC.log2
clusterDF[is.na(clusterDF)] <- 0

clusterDF.mean <- clusterDF %>% #filter(Avg.50 != 0 | Avg.60 != 0) %>%
  unite(Protein.Info, Protein, UniProt.Id, sep = "*") %>% column_to_rownames(var = "Protein.Info")

######## Run this when deciding the number of kmeans clusters ###############
#Change this to change the limit to amount of clusters you calculate for
clusLength <- 20

# #remove the normalized condtion column that is all 0
# clusterDF.mean <- clusterDF.mean.all[ ,-1]

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

clusAmount = 10 #Change clusAmount to be the amount of clusters you decide on, based on the above plots

d <- clusterDF.mean

# Perform hierarchical clustering
dist_matrix <- dist(d)  # Compute distance matrix
hc <- hclust(dist_matrix, method = "complete")  # Perform hierarchical clustering

# Perform hierarchical clustering
clustering_rows <- hclust(dist(d)) 

# Create the heatmap
set.seed(52)   
pheatmap.kmean<-pheatmap(d, kmeans_k = clusAmount, cluster_rows = TRUE, cluster_cols = FALSE,
                         clustering_method = "complete",  # You can change the method as needed
                         cutree_rows = 1,  # Number of clusters for rows, estimated from dendrogram above
                         main = "Heatmap of Proteins with Hierarchical Clustering using K-means", 
                         labels_row = d$Protein.Group,
                         color = hcl.colors(50, "Temps"), fontsize_row = 10, display_numbers = TRUE)

cluster.values <- pheatmap.kmean$kmeans$cluster %>% data.frame() %>% rownames_to_column() %>%
  setNames(c("Protein.Info","Cluster"))

clusterDF.mean.fixed <- clusterDF.mean %>% 
  rownames_to_column(var = "Protein.Info") %>% full_join(cluster.values, by = "Protein.Info") %>%
  separate_wider_delim(col = Protein.Info, names = c("Protein", "UniProt.AC"), delim = "*")

#Export proteins with cluster assignment
#write.csv(clusterDF.mean.fixed, "Avg.FC.log2.A.Kmean.csv")

