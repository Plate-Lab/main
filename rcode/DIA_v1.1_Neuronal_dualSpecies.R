############################### READ ME #########################################

# MS Proteomics Pipeline - DIA Unlabeled - Written for DIA-NN report.tsv files ONLY
  #DIA-NN:  Demichev, V., et al. (2020). https://doi.org/10.1038/s41592-019-0638-x

#Written by Sarah Garcia, Lea Barny, Jake Hermanson - 2024

#Set Working Directory - unique to each system
#setwd("path/to/your/report.tsv/file")

setwd("C:/Users/nyanc/OneDrive/Documents/All Files - personal/Work_Vandy/R_misc/Generic_MS_pipeline/Version1.0_Generic")

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
library(cleaver)          #used for species separation

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
diann.dia.qc <- function(raw.data, cont.name){
  
  cat("Cleaning up raw data...\n")                  #print message to console
  
  before <- length(unique(raw.data$Protein.Group))  #store size for comparison
  
  #Data Clean-up - Parameters set at: Q.value <= 0.01; 1+ Peptide per Protein
  raw.data.filtered <- raw.data %>%                                       #Takes our imported raw data
    group_by(Protein.Group) %>%                                           #Sort peptides into groups based on sequence
    filter(!str_detect(Protein.Group, "^Cont_")) %>%                      #filter out groups that match a given contaminant
    filter(Global.Q.Value <= 0.01) %>% filter(Q.Value <= 0.01) %>%        #Filter for qvalues
    group_by(Run, Protein.Group) %>%                                      #group by Run and Protein.Group
    filter(length(unique(Stripped.Sequence))>1)                           #Filters for lengths > 1 of pasted sequence 

  
  after <- length(unique(raw.data.filtered$Protein.Group))                #store size for comparison
  
  removed <- (1-(after/before))*100                                       #get size comparison
  
  cat("\nRemoved", removed, "% of the data due to q.value > 0.01 \n and/or only 1 peptide per protein observed \n and/or determined to be a contaminant.\n")
  
  export.file <- raw.data.filtered %>%                                    #take our qc filtered data                   
    dplyr::select(Run, Protein.Group, Genes, Stripped.Sequence,           #select relevant columns
                  Precursor.Id, Precursor.Normalised) %>%   
    mutate(Precursor.Normalised = na_if(Precursor.Normalised, 0)) %>%     #replace any zeros with NA
    tidyr::unite(Protein.Info, Genes, Protein.Group, sep = "*")           #combine two columns, separated by '*'
  
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
  combined_df_filter <- combined_df %>% distinct() %>% filter(rowSums(is.na(.)) < sample.amount)
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
    geom_point(data = subset(result.volcano, significant_high), aes(color = "high-conf"), size = 2)
  
  return(volcano_plot)
} #Define adj.p-value cutoff

############### Step One: Input and Reformat Data ##############################

#Load data files - .csv or .tsv (doesn't like .xlsx)
#raw.data <- load_file("your.diann.report.tsv")

raw.data <- load_file("report_noRT_newContam_13samples_neuro.tsv")

#Clean up data files
#raw.data.qc <- diann.dia.qc(raw.data, cont.name = "your.contaminants.file.fasta")
raw.data.qc <- diann.dia.qc(raw.data, cont.name = "contaminants.fasta")

    ### Dual-species ###
    
    #Define species list and proteome directory
    species.list<-c("Homo sapiens","Mus musculus")           #define species in order of input file
    dir<-"Proteome_Dual/"                                    #input file folder location (folder within working directory) 
    
    #Load functions
    digest.two.species <- function(species.list, dir) {
      
      chomp.protein<-function(species.list,dir){
        files<-list.files(dir, pattern = '*.fasta')     #inputs .fa filenames in specified folder
        #Read in each file in data folder
        for (file in 1:length(files)){                  #1:length(files) makes exact range to fit data
          assign(files[file],
                 read.csv(paste(dir,files[file],sep=''),header = FALSE,sep = "\n",na.strings = "N/A") #default sep=',' changed it to new line
          )
        }
        proteome.list<-mget(ls(pattern='*.fasta'))      #gets all objects with '.fa', pushes into list
        
        i=0                                             #Initialize empty objects for faster processing time
        flag=0
        y<-data.frame()
        protein.size<-list()
        sequence.all<-list()
        alist<-list()
        species.name<-species.list                      #this will be user supplied, must match how it looks in FASTA
        b<-''
        c<-''
        entry<-0
        
        for(j in proteome.list){                        #for each 'entry' i.e. given proteome
          entry<-entry+1                                #make marker for easy printing/speciation
          cat(paste("Now importing:",names(proteome.list[entry]), "\n"))  #basic progress checker
          for(i in j[,1]){                              #for each line of inputted proteomeFASTA
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
              pattern.species<-paste0("OS=",species.name[entry]) #make species name pattern
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
        }
        y[,4]<-protein.size                             #add sizes to df
        y[,5]<-sequence.all                             #add sequences to df
        colnames(y)<-c("Species","Gene.Name","UniProt Type","AA Length","Sequence") #add colnames to object 'y'
        y$Species <- str_replace_all(y$Species, "OS=", "")  #remove OS= from species column
        y$`Gene.Name` <- str_replace_all(y$`Gene.Name`, "GN=", "")  #remove GN= from gene name column
        proteomes<<-y %>% separate(`UniProt Type`, into = c("discard", "UniProt.AC","discard2"), sep = "\\|") %>%
          dplyr::select(-discard,-discard2)              #Drop the 'discard' columns
      }   #Load 'chomp' function to import proteome files
      chomp.protein(species.list,dir)                  #Generate 'Proteomes" object (both species combined)
      
      species.A <<- species.list[[1]]                  #Store species names
      species.B <<- species.list[[2]]                  #Store species names
      
      proteomes.A <- proteomes %>%                     #Take our Proteomes object
        filter(Species == species.A)
      
      proteomes.B <- proteomes %>%                     #Take our Proteomes ojbect
        filter(Species == species.B)
      
      #If missing a gene name, fill in name with UniProt Accession number
      proteomes.A$Gene.Name <- ifelse(is.na(proteomes.A$Gene.Name), 
                                      proteomes.A$UniProt.AC, proteomes.A$Gene.Name)
      proteomes.B$Gene.Name <- ifelse(is.na(proteomes.B$Gene.Name), 
                                      proteomes.B$UniProt.AC, proteomes.B$Gene.Name)

      cat("Now cleaving peptides...\n")
      
      #Create filtered peptide list
      tryp.digest.A <- proteomes.A$Sequence %>%       #Looks at protein sequences one at a time
        cleave(enzym = "trypsin", missedCleavages = 0:1, unique = TRUE) %>% #Specify trypsin enzyme, 0-1 missed cleavages, and only unique peptides (per protein)
        lapply(data.frame) %>%                        #Makes all list entries into data frames
        setNames(str_c(proteomes.A$Gene.Name, proteomes.A$UniProt.AC, proteomes.A$Species, sep = '*')) %>% #Rename data frames with the associated protein name
        lapply(bind_rows) %>%                         #Collapses the peptides into single protein entries (all peptides are grouped under their associated protein)
        lapply(function(df) {                         #Make a small function to apply across the list
          filtered_df <- filter(df, str_length(X..i..) > 6 & str_length(X..i..) < 31) #filter for peptides with charLength 6 < x < 31
          colnames(filtered_df) <- "Peptide"          #Rename the column names of each protein data frame
          return(filtered_df)                         #Gives the edited filtered_df back to the parent pipe sequence
        }) %>% discard(~ nrow(.) == 0)                #discard proteins with zero peptides (after filtering for length)
      
      #Create filtered peptide list - repeating above code block for the second species
      tryp.digest.B <- proteomes.B$Sequence %>%
        cleave(enzym = "trypsin", missedCleavages = 0:1, unique = TRUE) %>%
        lapply(data.frame) %>%
        setNames(str_c(proteomes.B$Gene.Name, proteomes.B$UniProt.AC, proteomes.B$Species, sep = '*')) %>%
        lapply(bind_rows) %>%
        lapply(function(df) {
          filtered_df <- filter(df, str_length(X..i..) > 6 & str_length(X..i..) < 31)
          colnames(filtered_df) <- "Peptide"
          return(filtered_df)
        }) %>% discard(~ nrow(.) == 0)

      cat("Performing species identification...\n")
      
      #Assemble Species.A Peptides
      species.A.peptides <<- tryp.digest.A %>%         #Takes our tryptic digestion
        bind_rows(.id = "source") %>%                 #Collapses protein data frames into a single data frame
        dplyr::select(source, Peptide) %>%            #Selects for specified columns in specified order
        separate_wider_delim(cols = source, delim = '*', names = c("Gene.Name", "UniProt.Id", "Species")) #Separates 'source' column into three columns (descriptors)
      
      #Assemble Species.B Peptides
      species.B.peptides <<- tryp.digest.B %>%
        bind_rows(.id = "source") %>%
        dplyr::select(source, Peptide) %>%
        separate_wider_delim(cols = source, delim = '*', names = c("Gene.Name", "UniProt.Id", "Species"))
      
      #Group by similar peptides, filter for groups with only one entry (proteotypic)
      proteotypic.A <- species.A.peptides %>% group_by(Peptide) %>% filter(n() == 1) %>% mutate(P=1)
      proteotypic.B <- species.B.peptides %>% group_by(Peptide) %>% filter(n() == 1) %>% mutate(P=1)
      
      #Group by similar peptides, filter for groups with one+ entry (razor)
      multi.A <- species.A.peptides %>% group_by(Peptide) %>% filter(n() > 1) %>% mutate(P=0)
      multi.B <- species.B.peptides %>% group_by(Peptide) %>% filter(n() > 1) %>% mutate (P=0)
      
      #We now have a master list of all peptides with Proteotypic Status
      all.species.all.peptides<-bind_rows(proteotypic.A, proteotypic.B, multi.A, multi.B)
      
      #Group by peptide, filter for Species.A
      A.only<<-all.species.all.peptides %>% group_by(Peptide) %>% 
        filter(any(str_detect(Species, species.A))) %>%
        filter(!(any(str_detect(Species, species.B))))
      
      #Group by peptide, filter for Species.B
      B.only<<-all.species.all.peptides %>% group_by(Peptide) %>% 
        filter(any(str_detect(Species, species.B))) %>%
        filter(!(any(str_detect(Species, species.A))))
      
      #Compile list of only shared peptides
      shared<<- all.species.all.peptides %>% group_by(Peptide) %>%
        filter(any(str_detect(Species, species.B))) %>%
        filter(any(str_detect(Species, species.A))) %>%
        ungroup() %>% distinct() %>% mutate(Species = 'Both')
      
      #Compile our peptide lists
      all.species.all.peptides<<-bind_rows(A.only, B.only, shared)

      cat("Done!\n")
      
}  #outputs tryptic database generated from supplied proteomes (2)
    match.data.sm <- function(raw.data) {
      
      cat("Matching input data to tryptic database...\n")
      
      #Pull out stripped and precursor.id sequence, along with run and precursor.normalised
      raw.data.filtered.subset <- subset(raw.data[, c('Run', 'Precursor.Id', 'Stripped.Sequence', 'Precursor.Normalised')]) %>%
        filter(Precursor.Normalised != 0) %>% distinct()
      
      #Rename for easier matching
      names(raw.data.filtered.subset) <- c("Run", "Precursor.Id", "Peptide", "Precursor.Normalised")
      
      #Match raw data to tryptic database - Specifies only proteotypic peptides
      raw.data.filtered.matched <- data.frame(raw.data.filtered.subset) %>% 
        inner_join(all.species.all.peptides, by = "Peptide", relationship = "many-to-many") %>%
        filter(P == 1) %>% distinct()
      
      cat("Identifying species-specific peptides...\n")
      
      # #Now to compile our species groups
      raw.data.matched.A <<- data.frame(raw.data.filtered.matched) %>%
        filter(Species == species.A) %>% filter(P == 1) %>% dplyr::select(-Peptide, -Species, -P) %>% 
        distinct() %>% mutate(Protein.Info = paste(Gene.Name, UniProt.Id, sep = '*')) %>%
        dplyr::select(Run, Precursor.Id, Precursor.Normalised, Protein.Info) %>% 
        filter(Precursor.Normalised != 0)
      
      raw.data.matched.B <<- data.frame(raw.data.filtered.matched) %>%
        filter(Species == species.B) %>% filter(P == 1) %>% dplyr::select(-Peptide, -Species, -P) %>% 
        distinct() %>% mutate(Protein.Info = paste(Gene.Name, UniProt.Id, sep = '*')) %>%
        dplyr::select(Run, Precursor.Id, Precursor.Normalised, Protein.Info) %>% 
        filter(Precursor.Normalised != 0)
      
      raw.data.matched.shared <<- data.frame(raw.data.filtered.matched) %>%
        filter(Species == "Both") %>% filter(P ==1) %>% dplyr::select(-Peptide, -Species, -P, -UniProt.Id)  %>%
        mutate(Gene.Name = toupper(Gene.Name)) %>% distinct()
      
      cat("Compiling group summaries...\n")
      
      #Need to compile numbers of species-uniqueness per run
      info.group <- function(df){
        group.size <- data.frame(df) %>%  #get group size (rows)
          group_by(Run) %>% group_size()
        group.name <- data.frame(df) %>%  #get group labels
          group_by(Run) %>% group_keys()
        group.info <- bind_cols(group.name, group.size)       #create group info matrix
        names(group.info) <- c("Run", "Peptides")
        group.info <<- group.info
      }
      
      suppressMessages(  {info.group.A <- info.group(raw.data.matched.A)
      info.group.B <- info.group(raw.data.matched.B)
      info.group.shared <- info.group(raw.data.matched.shared)
      info.group.allspecies.unique <- bind_cols(info.group.A, info.group.B, info.group.shared)
      info.group.allspecies.unique <- info.group.allspecies.unique %>% data.frame() %>%
        dplyr::select(-Run...3, -Run...5) %>% setNames(c("Run", species.A, species.B, "Both"))})
      
      group.info <<- info.group.allspecies.unique
      
      cat("Done!\n")
      
}                #Input is DIA-NN file, output is peptides matched to species
    reMatch <- function(raw.data) {
      
      cat("Compiling unmatched peptides...\n")
      
      raw.data.filtered.subset <- raw.data %>% data.frame() %>%
        dplyr::select(Run, Precursor.Id, Stripped.Sequence, Precursor.Normalised) %>%
        filter(Precursor.Normalised != 0) %>% distinct() %>% 
        setNames(c("Run", "Precursor.Id", "Peptide", "Precursor.Normalised"))
      
      unmatched <<- data.frame(raw.data.filtered.subset) %>%
        anti_join(all.species.all.peptides, by = "Peptide") %>% distinct()
      
      unmatched.wide <- unmatched %>% pivot_wider(names_from = 'Run', values_from = 'Precursor.Normalised')
      
      cat("Mapping unmatched peptides. This can take a few minutes...\n")
      
      df <- unmatched.wide %>% mutate(protein_match = map(Peptide, ~ proteomes %>% filter(str_detect(Sequence, .x)))) %>%
        unnest(protein_match) %>% filter(n() > 0) %>% distinct()   #this can take awhile...
      
      cat("Gathering peptide information...\n")
      
      df.edit <- df %>% dplyr::select(-Sequence, -'AA Length')
      
      df.edit.uniqueP <- df.edit %>% group_by(Precursor.Id) %>% filter(n()==1)

      cat("Identifying species-unique and proteotypic status...\n")
      
      uniqueP.A <- df.edit.uniqueP %>% filter(Species == species.A) %>% dplyr::select(-Species, -Peptide) %>% distinct() %>%
        pivot_longer(cols = where(is.numeric), names_to = "Run", values_to = "Precursor.Normalised") %>%
        filter(Precursor.Normalised != 0) %>% mutate(Protein.Info = paste(Gene.Name, UniProt.AC, sep = "*")) %>%
        dplyr::select(-Gene.Name, -UniProt.AC) %>% distinct()
      
      uniqueP.B <- df.edit.uniqueP %>% filter(Species == species.B) %>% dplyr::select(-Species, -Peptide) %>% distinct() %>%
        pivot_longer(cols = where(is.numeric), names_to = "Run", values_to = "Precursor.Normalised") %>%
        filter(Precursor.Normalised != 0) %>% mutate(Protein.Info = paste(Gene.Name, UniProt.AC, sep = "*")) %>%
        dplyr::select(-Gene.Name, -UniProt.AC) %>% distinct()
      
      raw.data.A.reMatched <<- rbind(raw.data.matched.A, uniqueP.A)
      raw.data.B.reMatched <<- rbind(raw.data.matched.B, uniqueP.B)
      
      added.A <- length(unique(uniqueP.A$Precursor.Id))
      added.B <- length(unique(uniqueP.B$Precursor.Id))
      
      cat("Done! \n Rematched", added.A, "peptides to", paste0(species.A, "."),
          "\n Rematched", added.B, "peptides to", paste0(species.B, "."))
      
}                      #Input is DIA-NN file, output is peptides matched and rematched
    
    #Run Functions
    digest.two.species(species.list, dir)
    match.data.sm(raw.data = raw.data.qc)
    reMatch(raw.data = raw.data.qc)
    
    # ### Code to Calculate Database Numbers ###
    # print(nrow(proteomes %>% filter(Species == "Homo sapiens")))
    # print(nrow(proteomes %>% filter(Species == "Mus musculus")))
    # 
    # print(nrow(species.A.peptides))
    # print(nrow(species.B.peptides))
    # 
    # print(length(unique(species.A.peptides$Gene.Name)))
    # print(length(unique(species.B.peptides$Gene.Name)))
    # 
    # print(nrow(all.species.all.peptides %>% filter(Species == "Homo sapiens")))
    # print(nrow(all.species.all.peptides %>% filter(Species == "Mus musculus")))
    # print(nrow(all.species.all.peptides %>% filter(Species == "Both")))
    # 
    # print(nrow(all.species.all.peptides %>% filter(Species == "Homo sapiens") %>% 
    #              ungroup() %>% dplyr::select(Gene.Name) %>% unique()))
    # print(nrow(all.species.all.peptides %>% filter(Species == "Mus musculus") %>% 
    #              ungroup() %>% dplyr::select(Gene.Name) %>% unique()))
    # print(nrow(all.species.all.peptides %>% filter(Species == "Both") %>% 
    #              ungroup() %>% dplyr::select(Gene.Name) %>% unique()))
    # 
    # print(nrow(all.species.all.peptides %>% filter(Species == "Homo sapiens" & P == 1)))
    # print(nrow(all.species.all.peptides %>% filter(Species == "Mus musculus" & P == 1)))
    # print(nrow(all.species.all.peptides %>% filter(Species == "Both" & P == 1)))
    # 
    # print(nrow(all.species.all.peptides %>% filter(Species == "Homo sapiens" & P == 1) %>% 
    #              ungroup() %>% dplyr::select(Gene.Name) %>% unique()))
    # print(nrow(all.species.all.peptides %>% filter(Species == "Mus musculus" & P == 1) %>% 
    #              ungroup() %>% dplyr::select(Gene.Name) %>% unique()))
    # print(nrow(all.species.all.peptides %>% filter(Species == "Both" & P == 1) %>% 
    #              ungroup() %>% dplyr::select(Gene.Name) %>% unique()))
    # ###                                    ###
    # 
    # ### Code to Calculate Database-Rematch Data Numbers ###
    # print(length(unique(raw.data$Precursor.Id)))
    # print(length(unique(raw.data$Protein.Group)))
    # 
    # print(length(unique(raw.data.qc$Precursor.Id)))
    # print(length(unique(raw.data.qc$Protein.Info)))
    # 
    # print(length(unique(raw.data.matched.A$Precursor.Id)))
    # print(length(unique(raw.data.matched.B$Precursor.Id)))
    # print(length(unique(raw.data.matched.shared$Precursor.Id)))
    # 
    # print(length(unique(raw.data.matched.A$Protein.Info)))
    # print(length(unique(raw.data.matched.B$Protein.Info)))
    # print(length(unique(raw.data.matched.shared$Gene.Name)))
    # 
    # print(length(unique(raw.data.A.reMatched$Precursor.Id)))
    # print(length(unique(raw.data.B.reMatched$Precursor.Id)))
    # 
    # print(length(unique(raw.data.A.reMatched$Protein.Info)))
    # print(length(unique(raw.data.B.reMatched$Protein.Info)))
    # 
    # unknown <- unmatched %>% dplyr::select(Precursor.Id) %>% unique() %>%
    #   anti_join(raw.data.A.reMatched, by = "Precursor.Id") %>%
    #   anti_join(raw.data.B.reMatched, by = "Precursor.Id")
    # ###                                                 ###

#Define treatment groups
treatments <- c("D30", "D50", "D60")

#Generate treatment groups - since A and B are from one data file, only need to run this once for both species
sample.groups <- assign.group(raw.data.A.reMatched, groups = treatments)

#Replace sample names with reformatted ones
raw.data.qc.A <- raw.data.A.reMatched %>% inner_join(sample.groups, by = "Run") %>%
  mutate(Run = sample) %>% select(-sample)

raw.data.qc.B <- raw.data.B.reMatched %>% inner_join(sample.groups, by = "Run") %>%
  mutate(Run = sample) %>% select(-sample)

raw.data.qc.shared <- raw.data.matched.shared %>% inner_join(sample.groups, by = "Run") %>%
  mutate(Run = sample) %>% select(-sample)

  ### (OPTIONAL) Export peptide groups ###
  peptides.A <- raw.data.qc.A %>% dplyr::select(Run, Protein.Info, Precursor.Id, Precursor.Normalised) %>%
    pivot_wider(names_from = Run, values_from = Precursor.Normalised)
  
  peptides.B <- raw.data.qc.B %>% dplyr::select(Run, Protein.Info, Precursor.Id, Precursor.Normalised) %>%
    pivot_wider(names_from = Run, values_from = Precursor.Normalised)
  
  peptides.shared <- raw.data.qc.shared %>% dplyr::select(Run, Gene.Name, Precursor.Id, Precursor.Normalised) %>%
    pivot_wider(names_from = Run, values_from = Precursor.Normalised)
  
  #write.csv(peptides.A, "peptides.A.csv")
  #write.csv(peptides.B, "peptides.B.csv")
  #write.csv(peptides.shared, "peptides.shared.csv")

### Visualization ###

#Boxplot 1 - A. Raw abundances across samples grouped by treatment
ggplot(raw.data.qc.A, aes(x= Run, y=Precursor.Normalised, fill = group)) +
  geom_hline(yintercept = median(raw.data.qc.A$Precursor.Normalised, na.rm = TRUE),  #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot()+
  labs(title= "A. Distribution of Values Grouped by Treatment", 
       x= "Sample",
       y= "Abundance") +
  scale_y_log10() +
  theme(axis.text.x= element_text(angle = 45, hjust =1)) +
  facet_grid(cols = vars(group), scales = "free_x")

#Boxplot 1 - B. Raw abundances across samples grouped by treatment
ggplot(raw.data.qc.B, aes(x= Run, y=Precursor.Normalised, fill = group)) +
  geom_hline(yintercept = median(raw.data.qc.B$Precursor.Normalised, na.rm = TRUE),  #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot()+
  labs(title= "B. Distribution of Values Grouped by Treatment", 
       x= "Sample",
       y= "Abundance") +
  scale_y_log10() +
  theme(axis.text.x= element_text(angle = 45, hjust =1)) +
  facet_grid(cols = vars(group), scales = "free_x")

############## Step Two: Normalization #########################################

#Global Median Normalize (per Sample)
# Specify: samples = column_samplenames    values = column_values
peptides.medNorm.A <- medNorm(raw.data.qc.A, samples = Run, values = Precursor.Normalised)
peptides.medNorm.B <- medNorm(raw.data.qc.B, samples = Run, values = Precursor.Normalised)

#Protein Roll-up using maxLFQ algorithm
proteins.medNorm.A <- max_LFQ(peptides.medNorm.A, values = "Norm_abundance")
proteins.medNorm.B <- max_LFQ(peptides.medNorm.B, values = "Norm_abundance")

#Log2 Transformation
proteins.medNorm.log2.A <- proteins.medNorm.A %>% data.frame() %>% log2()
proteins.medNorm.log2.B <- proteins.medNorm.B %>% data.frame() %>% log2()

##Prepare the data for plotting
sample.groups.v2 <- sample.groups %>% 
  mutate(Run = ifelse(grepl("^[0-9]", sample), paste0("X", sample), sample)) %>% 
  select(-sample)

    ### (OPTIONAL) Bait Normalization ###
    
    #Define bait using UniProt Accession ID (Ex. Human GAPDH)
    bait.A <- "P04406"
    bait.B <- "P16858"
    
    #Retrieve bait protein information
    bait.quant.A <- proteins.medNorm.log2.A %>% rownames_to_column(var = "Protein.Info") %>%
      separate_wider_delim(cols = "Protein.Info", delim = '*', names = c("Protein", "Accession")) %>%
      filter(Accession %in% bait.A) %>% 
      pivot_longer(cols = 3:ncol(.), names_to = "Run", values_to = "bait") %>% 
      dplyr::select(Run, bait)
    
    bait.quant.B <- proteins.medNorm.log2.B %>% rownames_to_column(var = "Protein.Info") %>%
      separate_wider_delim(cols = "Protein.Info", delim = '*', names = c("Protein", "Accession")) %>%
      filter(Accession %in% bait.B) %>% 
      pivot_longer(cols = 3:ncol(.), names_to = "Run", values_to = "bait") %>% 
      dplyr::select(Run, bait)
    
    #Barplot - Visualize Bait Levels Across Samples
    bait.quant.A %>% inner_join(sample.groups.v2, by = "Run") %>%
      ggplot(aes(x= Run, y=bait, fill = group)) +
      geom_col()+
      labs(title= "A. Bait Protein Levels", 
           x= "Sample",
           y= "Abundance") +
      theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
      facet_grid(cols = vars(group), scales = "free_x")
    
    bait.quant.B %>% inner_join(sample.groups.v2, by = "Run") %>%
      ggplot(aes(x= Run, y=bait, fill = group)) +
      geom_col()+
      labs(title= "A. Bait Protein Levels", 
           x= "Sample",
           y= "Abundance") +
      theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
      facet_grid(cols = vars(group), scales = "free_x")
    
    #Normalize to bait protein levels
    proteins.log2.baitNorm.A <- proteins.medNorm.log2.A %>% rownames_to_column(var = "Protein.Info") %>%
      pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2" ) %>%
      inner_join(bait.quant.A, by = "Run") %>% mutate(Log2norm = Log2 - bait) %>%
      dplyr::select(Protein.Info, Run, Log2norm) %>%
      pivot_wider(names_from = "Run", values_from = "Log2norm") %>%
      column_to_rownames(var = "Protein.Info")
    
    proteins.log2.baitNorm.B <- proteins.medNorm.log2.B %>% rownames_to_column(var = "Protein.Info") %>%
      pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2" ) %>%
      inner_join(bait.quant.B, by = "Run") %>% mutate(Log2norm = Log2 - bait) %>%
      dplyr::select(Protein.Info, Run, Log2norm) %>%
      pivot_wider(names_from = "Run", values_from = "Log2norm") %>%
      column_to_rownames(var = "Protein.Info")
    
    boxplot.proteins.log2.A <- proteins.log2.baitNorm.A %>%
      rownames_to_column(var = "Protein.Info") %>%
      pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2_medNorm") %>%
      inner_join(sample.groups.v2, by = "Run") %>% na.omit()
    
    boxplot.proteins.log2.B <- proteins.log2.baitNorm.B %>%
      rownames_to_column(var = "Protein.Info") %>%
      pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2_medNorm") %>%
      inner_join(sample.groups.v2, by = "Run") %>% na.omit()
    
    #Boxplot 2 - A.Normalized Log2 Abundances grouped by Treatment
    ggplot(boxplot.proteins.log2.A, aes(x= Run, y=Log2_medNorm, fill = group)) +
      geom_hline(yintercept = median(boxplot.proteins.log2.A$Log2_medNorm, na.rm = TRUE), #this adds a line at the global median
                 color = "red4", linetype = "solid", size = 1) +
      geom_boxplot()+
      labs(title= "A. Distribution of Values Grouped by Treatment", 
           x= "Sample",
           y= "Abundance") +
      theme(axis.text.x= element_text(angle = 45, hjust =1)) +
      facet_grid(cols = vars(group), scales = "free_x")
    
    ggplot(boxplot.proteins.log2.B, aes(x= Run, y=Log2_medNorm, fill = group)) +
      geom_hline(yintercept = median(boxplot.proteins.log2.B$Log2_medNorm, na.rm = TRUE), #this adds a line at the global median
                 color = "red4", linetype = "solid", size = 1) +
      geom_boxplot()+
      labs(title= "B. Distribution of Values Grouped by Treatment", 
           x= "Sample",
           y= "Abundance") +
      theme(axis.text.x= element_text(angle = 45, hjust =1)) +
      facet_grid(cols = vars(group), scales = "free_x")

### Visualization ###
    
#Update proteins.medNorm.log2.A/B to proteins.log2.baitNorm.A/B if desired (or vice-versa)

boxplot.proteins.log2.A <- proteins.log2.baitNorm.A %>% 
  rownames_to_column(var = "Protein.Info") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2_medNorm") %>%
  inner_join(sample.groups.v2, by = "Run") %>% na.omit()

boxplot.proteins.log2.B <- proteins.log2.baitNorm.B %>% 
  rownames_to_column(var = "Protein.Info") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2_medNorm") %>%
  inner_join(sample.groups.v2, by = "Run") %>% na.omit()

#Boxplot 2 - A.Normalized Log2 Abundances grouped by Treatment
ggplot(boxplot.proteins.log2.A, aes(x= Run, y=Log2_medNorm, fill = group)) +
  geom_hline(yintercept = median(boxplot.proteins.log2.A$Log2_medNorm, na.rm = TRUE), #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot()+
  labs(title= "A. Distribution of Values Grouped by Treatment", 
       x= "Sample",
       y= "Abundance") +
  theme(axis.text.x= element_text(angle = 45, hjust =1)) +
  facet_grid(cols = vars(group), scales = "free_x")

#Boxplot 2 - B.Normalized Log2 Abundances grouped by Treatment
ggplot(boxplot.proteins.log2.B, aes(x= Run, y=Log2_medNorm, fill = group)) +
  geom_hline(yintercept = median(boxplot.proteins.log2.B$Log2_medNorm, na.rm = TRUE), #this adds a line at the global median
             color = "red4", linetype = "solid", size = 1) +
  geom_boxplot()+
  labs(title= "B. Distribution of Values Grouped by Treatment", 
       x= "Sample",
       y= "Abundance") +
  theme(axis.text.x= element_text(angle = 45, hjust =1)) +
  facet_grid(cols = vars(group), scales = "free_x")

#Barplot 1 - Unique Protein IDs Identified per Sample
protein.count.A <- proteins.log2.baitNorm.A %>% 
  rownames_to_column(var = "Protein.Info") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>% 
  na.omit() %>% info.group() %>% inner_join(sample.groups.v2, by = "Run") %>% na.omit()

protein.count.B <- proteins.log2.baitNorm.B %>% 
  rownames_to_column(var = "Protein.Info") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>% 
  na.omit() %>% info.group() %>% inner_join(sample.groups.v2, by = "Run") %>% na.omit()

ggplot(protein.count.A, aes(x= Run, y=Count, fill = group)) +
  geom_col()+
  labs(title= "A. Unique Protein IDs Identified per Sample", 
       x= "Sample",
       y= "Count") +
  theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
  facet_grid(cols = vars(group), scales = "free_x")

ggplot(protein.count.B, aes(x= Run, y=Count, fill = group)) +
  geom_col()+
  labs(title= "A. Unique Protein IDs Identified per Sample", 
       x= "Sample",
       y= "Count") +
  theme(axis.text.x= element_text(angle = 45, hjust =1), legend.position = "none") +
  facet_grid(cols = vars(group), scales = "free_x")

#PCA Plot
  #Copy our data
  data.A <- proteins.log2.baitNorm.A
  data.B <- proteins.log2.baitNorm.B
  
  #Replace all missing values with 1
  data.A[is.na(data.A)] <- 1
  data.B[is.na(data.B)] <- 1
  
  #Transpose the data
  data.A.t <- data.A %>% t() %>% data.frame()
  data.B.t <- data.B %>% t() %>% data.frame()
  
  #Store our sample names
  sample_names <- data.A.t %>% rownames()
  
  #Run the PCA analysis
  data.A.prcomp <- data.A.t %>% select(where(~ var(.) != 0)) %>% prcomp(scale = TRUE)
  data.B.prcomp <- data.B.t %>% select(where(~ var(.) != 0)) %>% prcomp(scale = TRUE)
  
  #Label our samples with their treatment groups
  data.A.groupInfo <- data.A.t %>%
    rownames_to_column(var = "Run") %>%
    inner_join(sample.groups.v2, by = "Run")
  
  data.B.groupInfo <- data.B.t %>%
    rownames_to_column(var = "Run") %>%
    inner_join(sample.groups.v2, by = "Run")
  
  #Plot the PCA analysis
  autoplot(data.A.prcomp, data=data.A.groupInfo, label= FALSE, size=4, colour = 'group') +
    theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
    labs(title = "A. Principal Component Analysis Plot")
  
  autoplot(data.B.prcomp, data=data.B.groupInfo, label= FALSE, size=4, colour = 'group') +
    theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
    labs(title = "B. Principal Component Analysis Plot")

############## Step Three: Statistics ##########################################

#Protein Filtering - If no NA filtering is needed, set points = 1
#Filters for proteins with at least x observations per treatment group
proteins.filtered.A <- filter.NA(proteins.log2.baitNorm.A, points = 3)
proteins.filtered.B <- filter.NA(proteins.log2.baitNorm.B, points = 3)

### Options: ###
# Multiple t.tests w/ p.adj(method = "fdr")  -  For comparison of 2 groups
# ANOVA / Mixed Model  -  For comparison of >= 3 groups

#Perform Multiple unpaired t.tests with equal variance - Visualize via Volcano Plot [Ex. DIA]
#comparison <- c("treatment", "control") #define in order of [treatment - control]

comparison <- c("D60", "D30")

proteins.log2.long.stat.A <- proteins.filtered.A %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>%
  na.omit() %>% inner_join(sample.groups.v2, by = "Run") %>%
  filter(group %in% comparison) %>%
  mutate(group = factor(group, levels = comparison)) %>%
  group_by(Protein.Info, group) %>% filter(n() > 2) %>% ungroup() %>%
  group_by(Protein.Info) %>% filter(n_distinct(group) == 2) %>%
  do(tidy(t.test(Log2 ~ group, data = ., var.equal = TRUE)))

proteins.log2.long.stat.B <- proteins.filtered.B %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>%
  na.omit() %>% inner_join(sample.groups.v2, by = "Run") %>%
  filter(group %in% comparison) %>%
  mutate(group = factor(group, levels = comparison)) %>%
  group_by(Protein.Info, group) %>% filter(n() > 2) %>% ungroup() %>%
  group_by(Protein.Info) %>% filter(n_distinct(group) == 2) %>%
  do(tidy(t.test(Log2 ~ group, data = ., var.equal = TRUE)))

#Adjust p-values
result.volcano.A <- proteins.log2.long.stat.A %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr", n = length(p.value))) %>%
  select(Protein.Info, estimate, p.adj) %>% 
  column_to_rownames(var = 'Protein.Info') %>% data.frame()

result.volcano.B <- proteins.log2.long.stat.B %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr", n = length(p.value))) %>%
  select(Protein.Info, estimate, p.adj) %>% 
  column_to_rownames(var = 'Protein.Info') %>% data.frame()

  ###   Simple Volcano Plot - EnhancedVolcano package   ###
  EnhancedVolcano(result.volcano.A, lab = rownames(result.volcano.A), 
                  x = "estimate", y = "p.adj", pointSize = 3.0,
                  FCcutoff = 0.5,  #set fold change cutoff
                  pCutoff = 0.01   #set p.adj cutoff
  )
  
  EnhancedVolcano(result.volcano.B, lab = rownames(result.volcano.B), 
                  x = "estimate", y = "p.adj", pointSize = 3.0,
                  FCcutoff = 0.5,  #set fold change cutoff
                  pCutoff = 0.01   #set p.adj cutoff
  )
    
  #Filter for significant proteins - Enhanced Volcano
  sig.prt.EV.A <- result.volcano.A %>% filter(abs(estimate) > 0.5 & p.adj < 0.01) %>%
    rownames_to_column(var = "Protein.Info")
  
  sig.prt.EV.B <- result.volcano.B %>% filter(abs(estimate) > 0.5 & p.adj < 0.01) %>%
    rownames_to_column(var = "Protein.Info")
  
  ###                                                  ###
  
  ###   Volcano plot using curvature function          ###
  volcano.curve.A <- volcano.curve(result.volcano.A)
  
  volcano.curve.A
  
  volcano.curve.B <- volcano.curve(result.volcano.B)
  
  volcano.curve.B
  
  #Filter for significant proteins - Curvature Volcano
  volcano.curve.A.data <- volcano.curve.A$data %>% 
    filter(significant_med == TRUE | significant_high == TRUE) %>%
    rownames_to_column(var = "Protein.Info")
  
  volcano.curve.B.data <- volcano.curve.B$data %>% 
    filter(significant_med == TRUE | significant_high == TRUE) %>%
    rownames_to_column(var = "Protein.Info")
    
  ###                                                 ###
  
  ###   Venn Diagram Comparing Both Volcano Plots ###
    venn.plot <- venn.diagram(
      x = list(Set1 = sig.prt.EV.A$Protein.Info, Set2 = volcano.curve.A.data$Protein.Info),
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
  prt.A <- proteins.filtered.A %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>%
    inner_join(sample.groups.v2, by = "Run") %>%
    filter(is.finite(Log2))
  
  prt.B <- proteins.filtered.B %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2") %>%
    inner_join(sample.groups.v2, by = "Run") %>%
    filter(is.finite(Log2))
  
  anova_result.A <- aov(Log2 ~ group, data = prt.A)
  anova_result.B <- aov(Log2 ~ group, data = prt.B)
  
  summary(anova_result.A)
  summary(anova_result.B)
  
  tukey_result.A <- TukeyHSD(anova_result.A)
  tukey_result.B <- TukeyHSD(anova_result.B)
  
  tukey_df.A <- as.data.frame(tukey_result.A$group)
  tukey_df.B <- as.data.frame(tukey_result.B$group)
  
  tukey_df_sig.A <- filter(tukey_df.A, `p adj` < 0.01 )
  tukey_df_sig.B <- filter(tukey_df.B, `p adj` < 0.01 )

############### Additional Analysis ############################################

### Calculate Fold Changes normalized to user-defined control ###

  #Define the control group
  control <- treatments[1]  
  
  #Copy our data
  data.A <- proteins.log2.baitNorm.A
  data.B <- proteins.log2.baitNorm.B
  
  #Get Control information
  data.control.A <- data.A %>% dplyr::select(ends_with(control)) %>% 
    mutate(sd = apply(., 1, sd, na.rm = TRUE)) %>% mutate(mean = rowMeans(dplyr::select(., -sd), na.rm = TRUE))
  
  data.control.B <- data.B %>% dplyr::select(ends_with(control)) %>% 
    mutate(sd = apply(., 1, sd, na.rm = TRUE)) %>% mutate(mean = rowMeans(dplyr::select(., -sd), na.rm = TRUE))
  
  #Calculate Fold Changes for Original Data Frame
  FC.log2.A <- data.A %>% mutate(Mean = data.control.A$mean, StDev = data.control.A$sd) %>%
    mutate_at(vars(-Mean, -StDev), ~ (.-Mean)) %>% dplyr::select(-Mean, -StDev) %>% filter(rowSums(!is.na(.)) > 0)
  
  FC.log2.B <- data.B %>% mutate(Mean = data.control.B$mean, StDev = data.control.B$sd) %>%
    mutate_at(vars(-Mean, -StDev), ~ (.-Mean)) %>% dplyr::select(-Mean, -StDev) %>% filter(rowSums(!is.na(.)) > 0)
  
  #Plot Fold Changes per sample
  FC.log2.A %>% rownames_to_column(var = "Protein.Info") %>%
    pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2FC") %>%
    inner_join(sample.groups.v2, by = "Run") %>%
    ggplot(aes(x = Run, y = Log2FC, fill = group)) + 
    geom_boxplot(outliers = TRUE) + theme_bw() +
    labs(title = "A. Boxplot of Log2 Fold Changes",
         x = "Sample",
         y = "Log2 Fold Change") +
    theme(axis.text.x= element_text(angle = 45, hjust =1)) +
    facet_grid(cols = vars(group), scales = "free_x")
  
  FC.log2.B %>% rownames_to_column(var = "Protein.Info") %>%
    pivot_longer(cols = 2:ncol(.), names_to = "Run", values_to = "Log2FC") %>%
    inner_join(sample.groups.v2, by = "Run") %>%
    ggplot(aes(x = Run, y = Log2FC, fill = group)) + 
    geom_boxplot(outliers = TRUE) + theme_bw() +
    labs(title = "B. Boxplot of Log2 Fold Changes",
         x = "Sample",
         y = "Log2 Fold Change") +
    theme(axis.text.x= element_text(angle = 45, hjust =1)) +
    facet_grid(cols = vars(group), scales = "free_x")
  
  #This block calculates the average fold change of each protein per treatment group
  Avg.FC.log2.A <- FC.log2.A %>% average.FC() %>% filter(rowSums(!is.na(.)) > 2) %>%
    separate_wider_delim(cols = Protein.Info, delim = '*', names = c("Protein", "UniProt.Id"))
  
  Avg.FC.log2.B <- FC.log2.B %>% average.FC() %>% filter(rowSums(!is.na(.)) > 2) %>%
    separate_wider_delim(cols = Protein.Info, delim = '*', names = c("Protein", "UniProt.Id"))
  
  #This block reformats the data into long format, then plots a boxplot of the average fold changes per treatment
  Avg.FC.log2.A %>% 
    pivot_longer(cols = 3:ncol(.), names_to = "Treatment", values_to = "AvgFC") %>%
    ggplot(aes(x = Treatment, y = AvgFC, fill = Treatment)) + 
    geom_boxplot(outliers = TRUE) + theme_bw() +
    labs(title = "A. Boxplot of Average Log2 Fold Changes",
         x = "Sample",
         y = "Log2 Fold Change")
  
  Avg.FC.log2.B %>% 
    pivot_longer(cols = 3:ncol(.), names_to = "Treatment", values_to = "AvgFC") %>%
    ggplot(aes(x = Treatment, y = AvgFC, fill = Treatment)) + 
    geom_boxplot(outliers = TRUE) + theme_bw() +
    labs(title = "A. Boxplot of Average Log2 Fold Changes",
         x = "Sample",
         y = "Log2 Fold Change")

### Clustering (Using K-Means Clustering) ###

  #Copy our data
  clusterDF.A <- Avg.FC.log2.A
  clusterDF.B <- Avg.FC.log2.B
  
  #Replace any NAs with zeros
  clusterDF.A[is.na(clusterDF.A)] <- 0
  clusterDF.B[is.na(clusterDF.B)] <- 0
  
  #Remake Protein.Info column
  clusterDF.mean.A <- clusterDF.A %>% #filter(Avg.50 != 0 | Avg.60 != 0) %>%
    unite(Protein.Info, Protein, UniProt.Id, sep = "*") %>% column_to_rownames(var = "Protein.Info")
  
  clusterDF.mean.B <- clusterDF.B %>% #filter(Avg.50 != 0 | Avg.60 != 0) %>%
    unite(Protein.Info, Protein, UniProt.Id, sep = "*") %>% column_to_rownames(var = "Protein.Info")
  
  ######## Run this when deciding the number of kmeans clusters ###############
  #Change this to change the limit to amount of clusters you calculate for
  clusLength <- 20

  result.A <- fviz_nbclust(clusterDF.mean.A, kmeans, method = "wss", k.max = clusLength)
  result.B <- fviz_nbclust(clusterDF.mean.B, kmeans, method = "wss", k.max = clusLength)
  
  print(result.A)
  print(result.B)
  
      #Do for first species
  
        silPoints <- 1
        df.silhouette <- data.frame(k = 2:clusLength)
        
        for(i in 1:20){  
          for(k in 2:clusLength){
            # Perform k-means clustering
            kmeans_result <- kmeans(clusterDF.mean.A, centers = k)
            # Compute silhouette information using silhouette function
            sil_info <- silhouette(kmeans_result$cluster, dist(as.matrix(clusterDF.mean.A)))
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
        
        silPoints <- 1
        df.silhouette <- data.frame(k = 2:clusLength)
        
        #Set number of clusters based on graphs above
        clusAmount.A <- 13
        
      #Repeat for second species
        
        silPoints <- 1
        df.silhouette <- data.frame(k = 2:clusLength)
        
        for(i in 1:20){  
          for(k in 2:clusLength){
            # Perform k-means clustering
            kmeans_result <- kmeans(clusterDF.mean.B, centers = k)
            # Compute silhouette information using silhouette function
            sil_info <- silhouette(kmeans_result$cluster, dist(as.matrix(clusterDF.mean.B)))
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
        
        #Set number of clusters based on graphs above
        clusAmount.B <- 11
        
  #############################################################################
  
  A <- clusterDF.mean.A
  B <- clusterDF.mean.B
  
  # # Perform hierarchical clustering
  # dist_matrix.A <- dist(A)  # Compute distance matrix
  # dist_matrix.B <- dist(B)  # Compute distance matrix
  # 
  # hc.A <- hclust(dist_matrix.A, method = "complete")  # Perform hierarchical clustering
  # hc.B <- hclust(dist_matrix.B, method = "complete")  # Perform hierarchical clustering
  
  # # Perform hierarchical clustering
  # clustering_rows.A <- hclust(dist(A))
  # clustering_rows.B <- hclust(dist(B))
  
  # Create the heatmap
  set.seed(52)   
  pheatmap.kmean.A<-pheatmap(A, kmeans_k = clusAmount.A, cluster_rows = TRUE, cluster_cols = FALSE,
                           clustering_method = "complete",  # You can change the method as needed
                           cutree_rows = 1,  # Number of clusters for rows, estimated from dendrogram above
                           main = "A. Heatmap of Proteins with Hierarchical Clustering using K-means", 
                           labels_row = A$Protein.Group,
                           color = hcl.colors(50, "Temps"), fontsize_row = 10, display_numbers = TRUE)
  
  set.seed(52)   
  pheatmap.kmean.B<-pheatmap(B, kmeans_k = clusAmount.B, cluster_rows = TRUE, cluster_cols = FALSE,
                           clustering_method = "complete",  # You can change the method as needed
                           cutree_rows = 1,  # Number of clusters for rows, estimated from dendrogram above
                           main = "B. Heatmap of Proteins with Hierarchical Clustering using K-means", 
                           labels_row = B$Protein.Group,
                           color = hcl.colors(50, "Temps"), fontsize_row = 10, display_numbers = TRUE)
  
  cluster.values.A <- pheatmap.kmean.A$kmeans$cluster %>% data.frame() %>% rownames_to_column() %>%
    setNames(c("Protein.Info","Cluster"))
  
  cluster.values.B <- pheatmap.kmean.B$kmeans$cluster %>% data.frame() %>% rownames_to_column() %>%
    setNames(c("Protein.Info","Cluster"))
  
  clusterDF.mean.fixed.A <- clusterDF.mean.A %>% 
    rownames_to_column(var = "Protein.Info") %>% full_join(cluster.values.A, by = "Protein.Info") %>%
    separate_wider_delim(col = Protein.Info, names = c("Protein", "UniProt.AC"), delim = "*")
  
  clusterDF.mean.fixed.B <- clusterDF.mean.B %>% 
    rownames_to_column(var = "Protein.Info") %>% full_join(cluster.values.B, by = "Protein.Info") %>%
    separate_wider_delim(col = Protein.Info, names = c("Protein", "UniProt.AC"), delim = "*")
  
  #Export proteins with cluster assignment
  #write.csv(clusterDF.mean.fixed.A, "Avg.FC.log2.A.Kmean.csv")
  #write.csv(clusterDF.mean.fixed.B, "Avg.FC.log2.B.Kmean.csv")

