#!/usr/bin/env Rscript

library("optparse")




# Function for filtering RNA epitopes 
filter_calling <- function(df){
    # Different Filter Steps Specific for RNA 
    #     had to remove steps specific to DNA to make it compatable for RNA
    #######################################
    # Coverage Filters
    #######################################
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Binding score filtering 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    #idx = which(df$Best.MT.Score <= 500)
    idx = which(df$Median.MT.Score <= 500)
    filtered_DF = df[idx, ]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # RNA coverage filtering 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    idx = which(filtered_DF$Tumor.RNA.Depth >= 10)
    filtered_DF = filtered_DF[idx, ]
    # #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # # DNA Normal coverage filtering 
    # #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # # This is only present in the DNA mehtod. has normal values 
    # idx = which(filtered_DF$Normal.Depth >= 10)
    # filtered_DF = filtered_DF[idx, ]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # RNA VAF filtering 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    idx = which(filtered_DF$Tumor.RNA.VAF >= 0.25)
    filtered_DF = filtered_DF[idx, ]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # DNA VAF filtering 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    idx = which(filtered_DF$Tumor.DNA.VAF >= 0.25)
    filtered_DF = filtered_DF[idx, ]
    # #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # # DNA Normal VAF filtering
    # #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # idx = which(filtered_DF$Normal.VAF <= 0.02)
    # filtered_DF = filtered_DF[idx, ]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Gene EXP filtering 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    idx = which(filtered_DF$Gene.Expression >= 1)
    filtered_DF = filtered_DF[idx, ]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Transcript EXP filtering 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    idx = which(filtered_DF$Transcript.Expression >= 1)
    filtered_DF = filtered_DF[idx, ]
    #################################################
    # Transcript Suport Level
    #################################################
    # Transcript Suport Level filtering
    idx <- which(!(is.na(filtered_DF$Transcript.Support.Level)))
    filtered_DF = filtered_DF[idx, ]
    
    idx = which(filtered_DF$Transcript.Support.Level <= 1)
    filtered_DF = filtered_DF[idx, ]
    
    return(filtered_DF)
}


option_list = list(
    make_option(c("-f", "--file"), 
                type="character", default=NA, 
                help="Input File",
                dest = "input"),
    make_option(c("-o", "--out"), type="character", default="out.txt", 
                help="output file name [default= %default]", 
                dest = "output")
) 

if (!is.null(option_list)){
    # Set Constants
    args_parsed = optparse::parse_args(optparse::OptionParser(option_list=option_list))
    input <- args_parsed$input
    output <- args_parsed$output
    
    df_RNA = read.table(file = input, sep = "\t",header = TRUE)
    
    output_df = filter_calling(df_RNA)
    
    write.table(output_df, output, quote=FALSE, append = FALSE, sep = "\t", dec = ".")
}





