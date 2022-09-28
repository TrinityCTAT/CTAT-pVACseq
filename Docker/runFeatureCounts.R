#!/usr/bin/env Rscript


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsubread")


# library(Rsubread)
# install.packages("futile.logger")
# library(futile.logger)
# GX_file <- "/home/mbrown/pVACseq/WDL/TEST/cromwell-executions/pVACseq/0049324c-7405-4a09-b912-7ff5559cf07f/call-RunExpressionData/execution/featureCounts_output.GX.txt"

RunFeatureCount <- function( input_bam,
                             GTF,
                             ref_genome,
                             threads,
                             type   = "transcript_id",
                             output){
    ###########################
    # RUN Subread: FeatureCount
    ###########################
    ## df Holds COUNTS and ANNOTATIONS
    ## counts:           Counts matrix 
    ##                   a data matrix containing read counts for each feature or meta-feature for each library.
    ## annotations:      annotation Matrix 
    ##                   a data frame with six columns including GeneID, Chr, Start, End and Length
    futile.logger::flog.info(paste("Running featureCounts "))
    df <- Rsubread::featureCounts(files = input_bam,
                            annot.inbuilt = "hg38",
                            annot.ext = GTF,
                            isPairedEnd = TRUE,
                            GTF.attrType = type,
                            GTF.featureType = "exon",
                            isGTFAnnotationFile = TRUE,
                            genome = ref_genome,
                            nthreads = threads)
    ################
    # Calculate TPM
    ################
    # Check that we have a length for each count
    if (dim(df$counts)[1] != dim(df$annotation)[1]) {
        error_message <- paste("Somthing went wrong.",
                               "Number of counts does not match the number of lengths.")
        futile.logger::flog.error(error_message)
        stop(error_message)
    }
    
    futile.logger::flog.info(paste("Calculating TPM Values"))
    a = df$counts/(df$annotation$Length/1000) 
    TPM = a/sum(a) * 1e6
    
    # Add these newly created TPM values to the counts DataFrame
    counts = data.frame(df$counts)
    counts$TPM = TPM
    
    # Make gene id column 
    counts <- cbind(Geneid = rownames(counts), counts)
    
    if (is.na(output)){
        output = paste0(type,'_counts.tsv')    
    }
    # Save the output 
    write.table(counts, 
                file=output, 
                quote=FALSE, sep='\t', row.names = FALSE)
    
    return(counts)
}



# input_bam <- "/home/mbrown/pVACseq/OUTPUT/recalibrated.bam"
# GTF <- "/home/bhaas/CTAT_GENOME_LIBS/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf"
# ref_genome <- "/home/bhaas/CTAT_GENOME_LIBS/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"
# 
# counts <- RunFeatureCount( input_bam,
#                            GTF,
#                            ref_genome,
#                            threads = 5)



##########################
# Command Line Arguments #
##########################
option_list = list(
    optparse::make_option(c("-B", "--BAM"), 
                action="store", default=NA, type='character',
                dest = "BAM", 
                help="BAM alignment file of interest."),
    optparse::make_option(c("-G", "--GTF"), 
                action="store", default=NA, type='character',
                dest = "GTF", 
                help="GTF annotation file."),
    optparse::make_option(c("-R", "--Reference"), 
                action="store", default=NA, type='character',
                dest = "Reference",
                help="Reference FASTA file."),
    optparse::make_option(c("-t", "--threads"), 
                action="store", default=2, type='integer',
                dest = "threads", 
                help="Number of threads."),
    optparse::make_option(c("-T", "--Type"), 
                action="store", default=NA, type='character',
                dest = "type", 
                help="transcript or gene expression. has to match the GTF file (transcript_id, gene_id)."),
    optparse::make_option(c("-O", "--output"), 
                action="store", default=NA, type='character',
                dest = "output", 
                help="output file name.")
)

if (!is.null(option_list)){
    # Set Constants
    args_parsed = optparse::parse_args(optparse::OptionParser(option_list=option_list))
    BAM <- args_parsed$BAM
    GTF <- args_parsed$GTF
    Reference <- args_parsed$Reference
    threads <- args_parsed$threads
    Type <- args_parsed$type
    output <- args_parsed$output
    RunFeatureCount(input_bam = BAM,
                    GTF = GTF,
                    ref_genome = Reference,
                    threads = threads,
                    type = Type,
                    output = output
                    )
}

