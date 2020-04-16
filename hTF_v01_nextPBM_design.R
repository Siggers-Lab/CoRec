#!/usr/bin/env Rscript

### DESCRIPTION ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Script to generate the PBM design, associated annotations, and session info for the expanded human TF (hTF) CoRec array



### INSTRUCTIONS ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# USAGE:
# Rscript hTF_v01_nextPBM_design.R args[1]
#
# ARGUMENTS:
#     args[1]     prefix to use for naming



### INITIALIZE ENVIRONMENT AND FUNCTIONS ---------------------------------------------------------------------------------------------------------------------------------------------

# set up environment
rm(list=ls())
options(scipen=999)
options(digits=22)

# load libraries
library(TFBSTools)
library(JASPAR2018)
library(plyr)

# function to fetch consensus sequence from a formatted pwm
get_consensus_seq <- function(pwm_mat) {
  
  # collect the number representation of the top-scoring base at each position
  top_scores <- apply(pwm_mat, 2, which.max)
  
  # bind all of the top-scoring positions together into a single vector
  top_scores <- as.character(top_scores)
  top_scores <- paste(top_scores, collapse="")
  
  # translate the numbers to their corresponding nucleotide
  consensus_seq <- chartr("1234", "ACGT", top_scores)
  
  # return consensus sequence
  return(consensus_seq)
  
}

# function to obtain reverse complement of a DNA sequence
rev_compl <- function(DNA_seq){
  # reverses and complements the input DNA sequence
  temp_seq <- unlist(strsplit(DNA_seq, split=''))
  temp_seq <- rev(temp_seq)
  temp_seq <- paste(temp_seq, collapse='')
  temp_seq <- chartr("ACGTN", "TGCAN", temp_seq)
  return(temp_seq)
}

# save the folder that contains files relevant to background processing/generation
bg_dir <- "C:/Users/bray/Google Drive/Boston University/Spring 2016/Siggers Rotation/hTF_array_project/data/PBM_design/hTF_v01/bg_probe_generation"

# TEST: local folder change (uncomment if using as a script on the cluster)
setwd("C:/Users/bray/Google Drive/Boston University/Spring 2016/Siggers Rotation/hTF_array_project/data/PBM_design")
args <- c("hTF_v01")

# save the prefix argument as a named variable
prefix <- as.character(args[1])

# create a directory to house all of the output files generated below
dir.create(paste(getwd(), prefix, sep='/'))
setwd(paste(getwd(), prefix, sep='/'))



### OBTAIN ALL CORE MOTIF SEEDS --------------------------------------------------------------------------------

# initialize options and obtain the human PWMs
opts <- list()
opts[["species"]] <- 9606 # Homo sapiens
JASPAR_list <- getMatrixSet(JASPAR2018, opts)
JASPAR_list <- toPWM(JASPAR_list, type="prob", pseudocounts = 0.1, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))

# initialize a list to hold the matrices
PWM_list <- list()

# initialize a dataframe to hold the metadata of interest
PWM_df <- data.frame(ID=character(),
                     TF_name=character(),
                     cons_seq=character(),
                     TF_class=character(),
                     alias=character(),
                     descr=character(),
                     family=character(),
                     symbol=character(),
                     exp_type=character(),
                     species=character(),
                     stringsAsFactors=F)
PWM_df_names <- names(PWM_df)

# iterate through the JASPAR_list to isolate matrices and metadata
for (i in 1:length(JASPAR_list)) {
  
  # isolate the name of the current entry
  curr_ID <- names(JASPAR_list[i])
  
  # use the name to isolate other info
  curr_TF_name <- JASPAR_list[i][[curr_ID]]@name
  curr_cons_seq <- get_consensus_seq(as.matrix(JASPAR_list[i][[curr_ID]]@profileMatrix))
  curr_TF_class <- JASPAR_list[i][[curr_ID]]@matrixClass
  curr_alias <- JASPAR_list[i][[curr_ID]]@tags$alias
  curr_descr <- JASPAR_list[i][[curr_ID]]@tags$description
  curr_family <- JASPAR_list[i][[curr_ID]]@tags$family
  curr_symbol <- JASPAR_list[i][[curr_ID]]@tags$symbol
  curr_exp_type <- JASPAR_list[i][[curr_ID]]@tags$type
  curr_species <- JASPAR_list[i][[curr_ID]]@tags$species
  
  # clean up the entries if necessary (replace NULL entries)
  if (is.null(curr_TF_name)) { curr_TF_name <- "-"}
  if (is.null(curr_TF_class)) { curr_TF_class <- "-"}
  if (is.null(curr_alias)) { curr_alias <- "-"}
  if (is.null(curr_descr)) { curr_descr <- "-"}
  if (is.null(curr_family)) { curr_family <- "-"}
  if (is.null(curr_symbol)) { curr_symbol <- "-"}
  if (is.null(curr_exp_type)) { curr_exp_type <- "-"}
  if (is.null(curr_species)) { curr_species <- "-"}
  
  # push the matrix to the list
  PWM_list[[i]] <- as.matrix(JASPAR_list[i][[curr_ID]]@profileMatrix)
  
  # push the metadata to the dataframe
  PWM_df <- rbind(PWM_df,
                  list(curr_ID,
                       curr_TF_name,
                       curr_cons_seq,
                       paste(curr_TF_class, collapse=","),
                       curr_alias,
                       curr_descr,
                       paste(curr_family, collapse=","),
                       curr_symbol,
                       curr_exp_type,
                       paste(curr_species, collapse=",")),
                  stringsAsFactors=F)
  names(PWM_df) <- PWM_df_names
  
  # clean up temp variables
  rm(curr_ID,
     curr_TF_name,
     curr_cons_seq,
     curr_TF_class,
     curr_alias,
     curr_descr,
     curr_family,
     curr_symbol,
     curr_exp_type,
     curr_species)
}

# write all of the motifs to individual files for easy reference if needed (or if the above changes)
dir.create(paste(getwd(), "JASPAR_CORE_motifs", sep='/'))
for (i in 1:nrow(PWM_df)) {
  write.table(PWM_list[[i]], file=paste(getwd(), "/JASPAR_CORE_motifs/", PWM_df[i, "ID"], ".txt", sep=''), row.names=T, col.names=F, quote=F, sep='\t')
}


### FLAG EQUIVALENT MOTIF SEEDS FOR FILTERING ----------------------------------------------------------

# sort the df from the previous step by the size of the consensus sequences (needed for algo below)
PWM_df$cons_seq_size <- nchar(PWM_df$cons_seq)
PWM_df <- PWM_df[order(-PWM_df$cons_seq_size, PWM_df$TF_name),]
PWM_df <- PWM_df[,1:(ncol(PWM_df)-1)]

# clone the initial df to be able to modify it
final_df <- PWM_df

# add a FILTER flag column and equivalent IDs columns
final_df$equiv_core_ID <- rep("-", nrow(final_df))
final_df$equiv_core_TF_name <- rep("-", nrow(final_df))
final_df$filter_flag <- rep("-", nrow(final_df))

# initialize a size threshold
size_thresh <- 0.90

# for each possibly duplicated entry in the original table
for (i in 2:nrow(PWM_df)) {
  
  # create a version of the df that excludes smaller consensus sequences that have yet to be seen
  curr_df <- PWM_df[-(i:nrow(PWM_df)),]
  
  # isolate the current consensus seq in question
  curr_cons_seq <- PWM_df[i, "cons_seq"]
  
  # search for the current cons seq in the mod version of the df
  for_grepl <- grepl(curr_cons_seq, curr_df$cons_seq)
  rev_grepl <- grepl(rev_compl(curr_cons_seq), curr_df$cons_seq)
  seq_grepl <- for_grepl | rev_grepl
  
  # isolate matches if one or more exist
  if (sum(seq_grepl>0)) {
    
    # isolate the matching entries
    curr_matches <- curr_df[which(seq_grepl), ]
    
    # for each of the matches
    for (j in 1:nrow(curr_matches)) {
      
      # determine the size of the cons seq relative to the match
      rel_size <- nchar(curr_cons_seq)/nchar(curr_matches[j, "cons_seq"])
      
      # if the size is comparable and this is a duplicate not yet encountered
      if (rel_size > size_thresh) {
        
        # flag the curr seq motif's annotation for filtering
        final_df[i, "filter_flag"] <- "FILTER"
        
      }
    }
  }
}

# sort resulting df by TF name
final_df <- final_df[order(final_df$TF_name),]

# annotate each entry of the final_df with all equivalent sites
for (i in 1:nrow(final_df)) {
  
  # create a version of the df that excludes only the current entry
  curr_df <- final_df[-i, ]
  
  # isolate the curent consensus seq in question
  curr_cons_seq <- final_df[i, "cons_seq"]
  
  # search for the current cons seq in the mod version of the df
  for_grepl <- grepl(curr_cons_seq, curr_df$cons_seq)
  rev_grepl <- grepl(rev_compl(curr_cons_seq), curr_df$cons_seq)
  seq_grepl <- for_grepl | rev_grepl
  
  # isolate matches if one or more exist
  if (sum(seq_grepl>0)) {
    
    # isolate the matching entries
    curr_matches <- curr_df[which(seq_grepl), ]
    
    # for each of the matches
    for (j in 1:nrow(curr_matches)) {
      
      # determine the final df index that corresponds to the current match
      match_idx <- which(final_df$ID==curr_matches[j, "ID"])
      
      # modify the match's annotation to show the curr cons seq motif is equivalent
      if (final_df[match_idx, "equiv_core_ID"]=="-") {
        final_df[match_idx, "equiv_core_ID"] <- final_df[i, "ID"]
        final_df[match_idx, "equiv_core_TF_name"] <- final_df[i, "TF_name"]
      } else {
        final_df[match_idx, "equiv_core_ID"] <- paste(final_df[match_idx, "equiv_core_ID"], final_df[i, "ID"], sep=",")
        final_df[match_idx, "equiv_core_TF_name"] <- paste(final_df[match_idx, "equiv_core_TF_name"], final_df[i, "TF_name"], sep=",")
      }
      
    }
    
  }
  
}









### BUILD SEED DATAFRAME -------------------------------------------------------

# filter out the equivalent seeds and scrap the flag column (no longer needed after filtering)
df <- final_df[which(final_df$filter_flag!="FILTER"), ]
df <- df[, 1:(ncol(df)-1)]

# start dummy columns needed for the SNV probe annotations in the following section (using CASCADE-esque conventions)
df$seed_names <- paste(df$ID, df$TF_name, sep="_")
df$SNV_pos_offset <- rep(0, nrow(df))
df$seed_nuc <- rep(NA, nrow(df))
df$SNV_nuc <- rep(NA, nrow(df))
df$probe_type <- rep("MOTIF", nrow(df))



### PAD cons_seq TO CREATE A target_seq COLUMN ---------------------------------------------------

# initialize a variable to hold the set of all possible non-repeating flanking dinucleotides
non_rep_dinuc <- c("AC", "AG", "AT",
                   "CA", "CG", "CT",
                   "GA", "GC", "GT",
                   "TA", "TC", "TG")

# obtain a vector of random set.seeds integers (so that the following step could be reproduced if needed)
set.seed(1234)
set_seed_ints <- sample(1:nrow(df), nrow(df), replace=F) 

# generate target_seq by padding cons_seq with flanking AT
df$target_seq <- rep("", nrow(df))
for (i in 1:nrow(df)) {
  
  # sample a set of random pads
  set.seed(set_seed_ints[i])
  rand_pads <- non_rep_dinuc[sample(1:length(non_rep_dinuc), 2, replace=T)]
  
  # replace the target_seq with the padded consensus_seq
  df[i, "target_seq"] <- paste(rand_pads[1], df[i, "cons_seq"], rand_pads[2], sep="") 
  
}

# remove clutter
rm(non_rep_dinuc,
   set_seed_ints)



### GENERATE SNV PROBES FOR MOTIF MODELING -------------------------------------------------------------------------------------------------------------------------------------------

# helper method to perform all the SNV modeling for a set of seed probes
generate.SNVs <- function(seeds) {
  
  # initialize vectors needed for the annotation
  ID <- vector(mode="character")
  TF_name <- vector(mode="character")
  cons_seq <- vector(mode="character")
  TF_class <- vector(mode="character")
  alias <- vector(mode="character")
  descr <- vector(mode="character")
  family <- vector(mode="character")
  symbol <- vector(mode="character")
  exp_type <- vector(mode="character")
  species <- vector(mode="character")
  equiv_core_ID <- vector(mode="character")
  equiv_core_TF_name <- vector(mode="character")
  seed_names <- vector(mode="character")
  SNV_pos_offset <- vector(mode="integer")
  seed_nuc <- vector(mode="character")
  SNV_nuc <- vector(mode="character")
  probe_type <- vector(mode="character")
  target_seq <- vector(mode="character")
  
  # initialize nucleotide alphabet
  nucleotides <- c("A", "C", "G", "T")
  
  # for each seed sequence in the input vector of seeds
  for (i in 1:nrow(seeds)) {
    # for each position in the probe seed sequence of interest
    for (j in 1:nchar(seeds[i,"target_seq"])) {
      # for each possible nucleotide (excluding the one already in the seed sequence)
      for (k in 1:length(nucleotides)) {
        # if the current nucleotide in the cycle is not equal to the seed nucleotide
        if (nucleotides[k]!=substring(seeds[i,"target_seq"], j, j)){
          # create a probe sequence with the nucleotide substition at the given position
          SNV_current <- paste(substring(seeds[i,"target_seq"], 1, j-1), nucleotides[k], substring(seeds[i,"target_seq"], j+1, nchar(seeds[i,"target_seq"])), sep="")
          
          # append the current SNV probe to the growing array of SNV probes
          target_seq <- c(target_seq, SNV_current)
          
          # append the current seed nucleotide (the original nucleotide in the reference seed probe)
          seed_nuc <- c(seed_nuc, substring(seeds[i,"target_seq"], j, j))
          
          # append the current SNV nucleotide (the nucleotide that the seed is mutated to)
          SNV_nuc <- c(SNV_nuc, nucleotides[k])
          
          # append the genomic position corresponding to where the SNV was applied
          SNV_pos_offset <- c(SNV_pos_offset, j)
          
          # generate rest of columns
          ID <- c(ID, as.character(seeds[i, "ID"]))
          TF_name <- c(TF_name, as.character(seeds[i, "TF_name"]))
          cons_seq <- c(cons_seq, as.character(seeds[i, "cons_seq"]))
          TF_class <- c(TF_class, as.character(seeds[i, "TF_class"]))
          alias <- c(alias, as.character(seeds[i, "alias"]))
          descr <- c(descr, as.character(seeds[i, "descr"]))
          family <- c(family, as.character(seeds[i, "family"]))
          symbol <- c(symbol, as.character(seeds[i, "symbol"]))
          exp_type <- c(exp_type, as.character(seeds[i, "exp_type"]))
          species <- c(species, as.character(seeds[i, "species"]))
          equiv_core_ID <- c(equiv_core_ID, as.character(seeds[i, "equiv_core_ID"]))
          equiv_core_TF_name <- c(equiv_core_TF_name, as.character(seeds[i, "equiv_core_TF_name"]))
          seed_names <- c(seed_names, as.character(seeds[i, "seed_names"]))
          probe_type <- c(probe_type, as.character(seeds[i, "probe_type"]))
          
        }
      }
    }
  }
  
  # return the resulting SNV probes
  return(data.frame(ID,
                    TF_name,
                    cons_seq,
                    TF_class,
                    alias,
                    descr,
                    family,
                    symbol,
                    exp_type,
                    species,
                    equiv_core_ID,
                    equiv_core_TF_name,
                    seed_names,
                    SNV_pos_offset,
                    seed_nuc,
                    SNV_nuc,
                    probe_type,
                    target_seq,
                    stringsAsFactors = F))
  
  # remove clutter
  rm(ID,
     TF_name,
     cons_seq,
     TF_class,
     alias,
     descr,
     family,
     symbol,
     exp_type,
     species,
     equiv_core_ID,
     equiv_core_TF_name,
     seed_names,
     SNV_pos_offset,
     seed_nuc,
     SNV_nuc,
     probe_type,
     target_seq)
  
}

# generate SNV probes for the tiling seeds and add to the table
df_SNV <- generate.SNVs(df)
df <- rbind(df, df_SNV)
rm(df_SNV)



### INSERT BACKGROUND PROBES ---------------------------------------------------------------------------------------------------------

# read in the background probe coordinates
x <- read.table(paste(bg_dir, "hg38_random_34_bases.bed", sep="/"), stringsAsFactors = F)

# bind the uppercase sequences to them
y <- read.table(paste(bg_dir, "hg38_random_34_bases_seqs.bed", sep="/"), stringsAsFactors = F)
x <- cbind(x, toupper(y[, ncol(y)]))
rm(y)

# modify the coordinate to match UCSC conventions
names(x) <- c("bg_chr", "bg_start", "bg_end", "target_seq")
x$target_seq <- as.character(x$target_seq)
x$bg_start <- x$bg_start+1

# pre-filter them for uniqueness
x <- x[which(!duplicated(x$target_seq)),]

# collect all important sequences used so far in the array design
all_seq <- final_df[which(final_df$filter_flag!="FILTER"), "cons_seq"]
all_seq <- c(all_seq, df$target_seq)
all_seq <- unique(all_seq)

# sequentially filter out each seed sequence from the background
for (i in 1:length(all_seq)) {
  
  # check to see if the current cons_seq is present in the background probes
  for_grepl <- grepl(all_seq[i], x$target_seq)
  rev_grepl <- grepl(rev_compl(all_seq[i]), x$target_seq)
  seq_grepl <- for_grepl | rev_grepl
  
  # if one or more matches exist
  if (sum(seq_grepl>0)) {
    
    # update the background probe data frame
    x <- x[which(!seq_grepl),]
    
  }
  
}

# keep a sample of 261 of the remaining background probes to include in the final design (number needed to fill rest of design)
x <- x[1:261,] # they were already randomized to begin with, no need to re-sample them
x$probe_order <- 1:261

# generate seed names for the resulting background probes
x$seed_names <- paste("BG", as.character(x$probe_order), x$bg_chr, x$bg_start, x$bg_end, sep="-")

# tag them with a "BACKGROUND" probe_type
x$probe_type <- rep("BACKGROUND", nrow(x))

# save a copy of the background probes with their coordinates (in case this is needed in the future)
write.table(x, file=paste(prefix, "BG_PROBES_ANNOT.txt", sep="_"), row.names=F, col.names=T, sep='\t', quote=F)

# keep only the columns relevant to the PBM design
x <- x[, c("seed_names", "probe_type", "target_seq")]

# bind background probes to the array design
df <- rbind.fill(df, x)

# update the all_seq variable
all_seq <- c(all_seq, x$target_seq)
all_seq <- unique(all_seq)

# clear clutter
rm(x)



### GENERATE REVERSE COMPLEMENT PROBES --------------------------------------------------------------------------------------------------------------------------------

# ensure that the target_seq column is a character vector (not factor)
df$target_seq <- as.character(df$target_seq)

# annotate each probe so far with the baseline orientation
df$probe_or <- rep("o1", nrow(df))

# copy the current design but replace the orientation
df_rev <- df
df_rev$probe_or <- rep("o2", nrow(df))

# obtain the reverse complement of each probe
for (i in 1:nrow(df)) {
  df_rev[i, "target_seq"] <- rev_compl(df[i, "target_seq"])
}

# add the reverse oriented probes into the array design
df <- rbind(df, df_rev)
rm(df_rev)



### GENERATE REPLICATE PROBES -----------------------------------------------------------------------------------------------------------------------------------------

# tag each probe in the design with the replicate 1 tag
df$probe_repl <- rep("r1", nrow(df))

# specify the number of replicates to include in the array
n_repl <- 5

# copy the current design
x <- df

# create the desired number of probe replicates and tag them with their replicate number
for (i in 2:n_repl) {
  # create a temporary copy of the current design
  y <- x
  
  # assign a replicate name to the current set
  y$probe_repl <- rep(paste("r", as.character(i), sep=''), nrow(y))
  
  # rbind the new replicates to the array design
  df <- rbind(df, y)
}

# remove the temporary clutter
rm(x, y, i, n_repl)



### GENERATE A BACKBONE SEQUENCE --------------------------------------------------------------------------------------------------------------------------------------

# helper function to generate a backbone sequence
generate.backbone <- function(seq_len, all_seq) {
  
  # initialize a test flag
  test_flag <- TRUE
  
  # while the backbone 
  while(test_flag) {
    
    # flip the test_flag
    test_flag <- FALSE
    
    # initialize a backbone sequence vector
    backbone_seq <- ""
    
    # fill the first position with a random nucleotide
    backbone_seq <- paste(backbone_seq, sample(1:4, 1), sep="")
    
    # for each remaining position in the desired sequence
    for (i in 2:seq_len) {
      
      # generate a random nucleotide that doesn't match the previous position
      possible_nucs <- 1:4
      prev_nuc <- as.integer(substr(backbone_seq, i-1, i-1))
      possible_nucs <- possible_nucs[-prev_nuc]
      backbone_seq <- paste(backbone_seq, sample(possible_nucs, 1), sep="")
      
    }
    
    # translate the integer sequence to nucleotides
    backbone_seq <- chartr("1234", "ACGT", backbone_seq)
    
    # test if any of the target_seqs are in the backbone_seq
    for (i in 1:length(all_seq)) {
      
      # check to see if the current cons_seq is present in the background probes
      for_grepl <- grepl(all_seq[i], backbone_seq)
      rev_grepl <- grepl(rev_compl(all_seq[i]), backbone_seq)
      seq_grepl <- for_grepl | rev_grepl
      
      # if one or more matches exist
      if (sum(seq_grepl>0)) {
        
        # update the test flag if the backbone needs to be regenerated
        test_flag <- TRUE
        break
      }
      
    }
    
  }
  
  # return the sequence to the parent process
  return(backbone_seq)
}

# helper function to test backbone after it generates
test.backbone <- function(backbone) {
  
  # initialize test variable
  test <- ""
  
  # test if any of the target_seqs are in the backbone_seq
  for (i in 1:length(all_seq)) {
    
    # check to see if the current cons_seq is present in the background probes
    for_grepl <- grepl(all_seq[i], backbone)
    rev_grepl <- grepl(rev_compl(all_seq[i]), backbone)
    seq_grepl <- for_grepl | rev_grepl
    
    test <- c(test, seq_grepl)
    
  }
  
  # test to see if backbone passes the test
  test <- test[2:length(test)]
  test <- as.logical(test)
  return(sum(test))
  
}

# add "TATA" motif to all_seq manually to prevent this motif from appearing in the backbone
all_seq <- c(all_seq, "TATA")

# generate a backbone sequence excluding matches from all_seq
backbone <- generate.backbone(34, all_seq)

# test the backbone sequence
test_int <- test.backbone(backbone)
stopifnot(test_int==0)

# update all_seq
all_seq <- c(all_seq, backbone)
all_seq <- unique(all_seq)



### APPEND CONSTANT PROBE REGIONS AND PRIMER SEQUENCES ----------------------------------------------------------------------------------------------------------------

# initialize a probe_seq column
df$probe_seq <- rep("", nrow(df))

# for each probe in the design
for (i in 1:nrow(df)) {
  
  # contruct the current probe seq using the target seq and backbone
  curr_target_seq <- as.character(df[i, "target_seq"])
  curr_probe_seq <- paste(curr_target_seq, substr(backbone, nchar(curr_target_seq)+1, nchar(backbone)), sep="")
  
  # push to the table
  df[i, "probe_seq"] <- curr_probe_seq
  
  # remove clutter
  rm(curr_target_seq,
     curr_probe_seq)
  
}

# helper function to test the primer (check to see if there are target sequences in it)
test.primer <- function(primer, all_seq) {
  
  # initialize test variable
  test <- ""
  
  # test if any of the target_seqs are in the backbone_seq
  for (i in 1:length(all_seq)) {
    
    # check to see if the current cons_seq is present in the background probes
    for_grepl <- grepl(all_seq[i], primer)
    rev_grepl <- grepl(rev_compl(all_seq[i]), primer)
    seq_grepl <- for_grepl | rev_grepl
    
    test <- c(test, seq_grepl)
    
  }
  
  # test to see if primer passes the test
  test <- test[2:length(test)]
  test <- as.logical(test)
  return(all_seq[which(test)])
  
}
primer_check <- test.primer("GTCTTGATTCGCTTGACGCTGCTG", all_seq)
stopifnot(length(primer_check)==0)

# apply GC cap and double-stranding primer
df$probe_seq <- paste("GC", df$probe_seq, "GTCTTGATTCGCTTGACGCTGCTG", sep="")



### GENERATE PROBE NAMES ----------------------------------------------------------------------------------------------------------------------------------------------

# use the probe category, ID, orientation, and replicate to create a unique name for each probe
df$probe_name <- paste(prefix, gsub("[[:punct:]]", "-", df$seed_names), df$SNV_pos_offset, df$seed_nuc, df$SNV_nuc, df$probe_or, df$probe_repl, sep="_")

# assign a probe ID based on the probe_name
df$probeID <- gsub("_o[1-5]_r[1-5]", "", df$probe_name)



### SAVE PBM DESIGN, ANNOTATIONS, AND SESSION INFO --------------------------------------------------------------------------------------------------------------------

# assign an order to the current configuration so it can be reconstituted if necessary
df$probe_order <- seq.int(1, nrow(df))

# write the complete JASPAR annotations to file
write.table(final_df, file=paste(prefix, "JASPAR_FULL_ANNOT.txt", sep="_"), row.names=F, col.names=T, sep='\t', quote=F)
write.table(final_df[which(final_df$filter_flag!="FILTER"), 1:(ncol(final_df)-1)], file=paste(prefix, "JASPAR_FILT_ANNOT.txt", sep="_"), row.names=F, col.names=T, sep='\t', quote=F)

# write full annotation, probes for Agilent, and session info
write.table(df, file=paste(prefix, "PBM_ANNOT.txt", sep="_"), row.names=F, col.names=T, sep='\t', quote=F)
PBM <- df[,c("probe_name", "probe_seq")]
names(PBM) <- c("ProbeID", "Sequence")
write.table(PBM, file=paste(prefix, "PROBES.txt", sep="_"), row.names=F, col.names=T, sep='\t', quote=F)
writeLines(capture.output(sessionInfo()), paste(prefix, "RsessionInfo.txt", sep="_"))

# write backbone used
write(backbone, file=paste(prefix, "BACKBONE.txt"), sep="_")



# ### DONE: RUN LOCAL TESTS ----------------------------------------------------------------------------
# 
# # read in the final PBM annotation to inspect
# setwd("C:/Users/bray/Google Drive/Boston University/Spring 2016/Siggers Rotation/hTF_array_project/data/PBM_design/hTF_v01")
# df <- read.table("hTF_v01_PBM_ANNOT.txt", header=T, sep='\t', stringsAsFactors = F)
# 
# # are all of the probe names less than 100 characters?
# table(nchar(as.character(df$probe_name)))
# 
# # do the expected number of probes appear in each probe_type category?
# x <- df[which(df$probe_repl=="r1" & df$probe_or=="o1"),]
# table(x$probe_type)
# 
# # what is the number of seeds included in design?
# x <- df[which(df$probe_type=="MOTIF"),]
# length(unique(as.character(x$seed_names)))
# 
# # do the expected number of probes appear in each probe orientation and replicate?
# x <- df
# table(x$probe_repl, x$probe_or)
# 
# # are the target sequences all the expected size?
# x <- df
# table(nchar(as.character(x$target_seq)), x$probe_type)
# 
# # are there only standard uppercase nucleotides in the probe sequences?
# x <- df
# length(grep("[ACGT]", as.character(x$target_seq)))   # 174340 (missing dummy probe)
# length(grep("[^ACGT]", as.character(x$target_seq)))  # 0
# length(grep("[acgt]", as.character(x$target_seq)))   # 0
# length(grep("[N]", as.character(x$target_seq)))      # 0
# length(grep("[ACGT]", as.character(x$probe_seq)))   # 174340
# length(grep("[^ACGT]", as.character(x$probe_seq)))  # 0
# length(grep("[acgt]", as.character(x$probe_seq)))   # 0
# length(grep("[N]", as.character(x$probe_seq)))      # 0
# 
# # are the probe names unique?
# x <- df
# length(unique(as.character(x$probe_name)))
# 
# # are the probes all the same length?
# x <- df
# table(nchar(as.character(x$probe_seq)), nchar(as.character(x$target_seq)))
# x <- df
# table(nchar(as.character(x$probe_seq)))
# 
# # assert that the correct number of motif probes and assert correct SNV sequences are present
# x <- df[which(df$probe_type=="MOTIF" & df$probe_repl=="r1" & df$probe_or=="o1"),]
# seeds <- unique(as.character(x$seed_names))
# for (z in 1:length(seeds)) {
#   
#   # collect probes from x that match the current seed name
#   test_probes <- x[which(x$seed_names==seeds[z]),]
#   target_seq <- test_probes[which(is.na(test_probes$seed_nuc)), "target_seq"]
#   
#   # raise error if number of probes doesn't match expectation
#   len_target_seq <- nchar(target_seq)
#   stopifnot(nrow(test_probes)==((3*len_target_seq)+1))
#   
#   # raise error if expected sequences don't exist
#   nucs <- c("A", "C", "G", "T")
#   
#   # collect z-score values for each position in the matrix
#   for (i in 1:4) {
#     for (j in 1:nchar(target_seq)) {
#       # construct the current sequence of interest using the current SNV
#       SNV_seq <- paste(substr(target_seq, 1, j-1), nucs[i], substr(target_seq, j+1, nchar(target_seq)), sep="")
#       
#       # collect the entry that has the SNV_seq as the target_seq
#       test_SNV <- test_probes[which(test_probes$target_seq==SNV_seq), "target_seq"]
#       stopifnot(length(test_SNV)==1)
#       
#     }
#   }
#   
# }
# 
# # verify background sequences are unique probes
# x <- df[which(df$probe_type=="BACKGROUND" & df$probe_repl=="r1" & df$probe_or=="o1"),]
# length(unique(x$probe_seq))
# 
# # manually inspect seeds
# x <- df[which(df$probe_type=="MOTIF"),]
# seeds <- unique(as.character(x$seed_names))
# test_probes <- x[which(x$seed_names==seeds[143]), c("seed_names", "seed_nuc", "SNV_nuc", "SNV_pos_offset", "target_seq", "probe_seq", "probe_or", "probe_repl")]
# write.table(test_probes, "test_probes.txt", row.names=F, col.names=T, sep='\t', quote=F)
