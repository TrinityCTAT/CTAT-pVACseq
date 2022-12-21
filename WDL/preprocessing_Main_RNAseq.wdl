version 1.0



            ############################################################
        ####################################################################
######################################################################################
#                                 pVACseq Main Workflow                             #
######################################################################################
        #####################################################################
            ############################################################

# pVACseq is a cancer immunotherapy pipeline for the identification of personalized Variant Antigens by Cancer Sequencing (pVACseq) 
# that integrates tumor mutation and expression data (DNA- and RNA-Seq). It enables cancer immunotherapy research by using massively 
# parallel sequence data to predicting tumor-specific mutant peptides (neoantigens) that can elicit anti-tumor T cell immunity



#import "Tasks_wdl/annotateVCF.wdl" as annotateVCF
#import "Tasks_wdl/phasing.wdl" as phasing

import "https://raw.githubusercontent.com/brownmp/CTAT-pVACseq/main/WDL/Tasks_wdl/annotateVCF.wdl"
import "https://raw.githubusercontent.com/brownmp/CTAT-pVACseq/main/WDL/Tasks_wdl/phasing.wdl"



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
workflow pVACseqPreprocess {
    input {

        #~~~~~~~~~~~~
        # Sample ID
        #~~~~~~~~~~~~
        String sample_id
      
        
        #~~~~~~~~~~~~
        # VCF Files
        #~~~~~~~~~~~~
        File HaplotypeCaller_VCF
        File HaplotypeCaller_VCF_index

        String Tumor_ID
        
        #~~~~~~~~~~~~
        # BAM Files
        #~~~~~~~~~~~~
        File BAM
        File BAM_index
        
        File GTF
        File RNA_editing_VCF

        #~~~~~~~~~~~~
        # Resources 
        #~~~~~~~~~~~~
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File VEP_Reference

        File gnomadVCF
        File gnomadVCFindex

        #~~~~~~~~~~~~
        # CPU count 
        #~~~~~~~~~~~~
        Int cpus = 10

        #~~~~~~~~~~~~
        # general runtime settings
        #~~~~~~~~~~~~
        Int preemptible = 2
        # String docker = "brownmp/<>:devel"
    }

    parameter_meta {
        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}
        cpus:{help:"CPU count"}
        #docker:{help:"Docker image"}
    }


    #########################
    # Input File Preparation
    #########################

    #~~~~~~~~~~~~~~~~~~~~
    # VEP Annotations
    #~~~~~~~~~~~~~~~~~~~~
    call annotateVCF.RunAnnotateVEP as RunAnnotateVEP{
        input:
            VCF                 = HaplotypeCaller_VCF,
            VEP_Reference       = VEP_Reference,

            ref_fasta           = ref_fasta,
            ref_fasta_index     = ref_fasta_index,
            ref_dict            = ref_dict,
            
            cpus                = cpus,
            preemptible         = preemptible,
            sample_id           = sample_id
    }

    #~~~~~~~~~~~~~~~~~~~~
    # Decompose
    #~~~~~~~~~~~~~~~~~~~~
    call annotateVCF.RunDecompose as RunDecompose{
        input:
            VCF                 = RunAnnotateVEP.VEP_output,

            preemptible         = preemptible,
            sample_id           = sample_id

    }

    #~~~~~~~~~~~~~~~~~~~~
    # Create Exon BED File
    #~~~~~~~~~~~~~~~~~~~~
    call annotateVCF.RunCreateExonBED as RunCreateExonBED{
        input:
            GTF                 = GTF,
            
            cpus                = 1,
            preemptible         = preemptible
    }

    #~~~~~~~~~~~~~~~~~~~~
    # Filter Non-Exon Regions
    #~~~~~~~~~~~~~~~~~~~~
    call annotateVCF.RunFilterNonExons as RunFilterNonExons{
        input:
            Exon_BED_File       = RunCreateExonBED.exon_bed,
            VCF                 = RunDecompose.decomposed_VCF,
            VCF_index           = RunDecompose.decomposed_VCF_index,
            
            cpus                = 1,
            preemptible         = preemptible,
            sample_id           = sample_id
    }


    #~~~~~~~~~~~~~~~~~~~~
    # Filter RNA-editing
    #~~~~~~~~~~~~~~~~~~~~
    call annotateVCF.RunAnnotateRNAediting as RunAnnotateRNAediting{
        input:
            VCF                 = RunFilterNonExons.exon_filtered_VCF,
            VCF_index           = RunFilterNonExons.exon_filtered_VCF_index,
            RNA_editing_VCF     = RNA_editing_VCF,
            
            cpus                = 1,
            preemptible         = preemptible,
            sample_id           = sample_id
    }

    call annotateVCF.RunFilterRNAediting as RunFilterRNAediting{
        input:
            VCF                 = RunAnnotateRNAediting.annotated_VCF,
            VCF_index           = RunAnnotateRNAediting.annotated_VCF_index,
            
            cpus                = 1,
            preemptible         = preemptible,
            sample_id           = sample_id
    }

    #~~~~~~~~~~~~~~~~~~~~
    # gnomAD Annotations and filtering 
    #~~~~~~~~~~~~~~~~~~~~
    call annotateVCF.RunAnnotateGnomad as RunAnnotateGnomad{
        input:
            VCF                 = RunFilterRNAediting.filtered_VCF,
            VCF_index           = RunFilterRNAediting.filtered_VCF_index,
            gnomad_vcf          = gnomadVCF,
            gnomad_vcf_index    = gnomadVCFindex,
            
            preemptible         = preemptible,
    }




    #~~~~~~~~~~~~~~~~~~~~
    # Read Count Data
    #~~~~~~~~~~~~~~~~~~~~
    call annotateVCF.RunBamReadcount as RunBamReadcount{
        input:
            VCF                 = RunAnnotateGnomad.gnomad_filtered_VCF,
            VCF_index           = RunAnnotateGnomad.gnomad_filtered_VCF_index,
            BAM                 = BAM,
            BAM_index           = BAM_index,
            ref_fasta           = ref_fasta,
            ref_fasta_index     = ref_fasta_index,
            ref_dict            = ref_dict,
            
            cpus                = 1,
            preemptible         = preemptible,
            sample_id           = sample_id
    }
    call annotateVCF.RunAddReadcount as RunAddReadcount{
        input:
            VCF                 = RunAnnotateGnomad.gnomad_filtered_VCF,
            readcount_indel     = RunBamReadcount.sample_bam_readcount_indel,
            readcount_snv       = RunBamReadcount.sample_bam_readcount_snv,
            RNAorDNA            = "RNA",

            cpus                = 1,
            preemptible         = preemptible,
            sample_id           = sample_id
    }

    #~~~~~~~~~~~~~~~~~~~~
    # Expression Data
    #~~~~~~~~~~~~~~~~~~~~
    call annotateVCF.RunExpressionData as RunExpressionData{
        input:
            VCF                 = RunAddReadcount.annotated_vcf,
            BAM                 = BAM,
            BAM_index           = BAM_index,
            ref_fasta           = ref_fasta,
            ref_fasta_index     = ref_fasta_index,
            ref_dict            = ref_dict,
            
            Tumor_ID            = Tumor_ID,

            GTF                 = GTF,

            preemptible         = preemptible,
            sample_id           = sample_id,
            cpus                = cpus
    }



    # Add expression data annotation 
    call annotateVCF.RunAddExpressionData as RunAddExpressionData{
        input:
            VCF                 = RunAddReadcount.annotated_vcf,
            featureCounts_GX    = RunExpressionData.featureCounts_GX,
            featureCounts_TX    = RunExpressionData.featureCounts_TX,

            ref_fasta           = ref_fasta,
            ref_fasta_index     = ref_fasta_index,
            ref_dict            = ref_dict,
            
            Tumor_ID            = Tumor_ID,

            cpus                = cpus,
            preemptible         = preemptible,
            sample_id           = sample_id
    }



    output {

        File annotated_VCF = RunAddExpressionData.annotated_vcf
        File phasedVCF     = RunDecompose.decomposed_VCF
        
    }
    
}
