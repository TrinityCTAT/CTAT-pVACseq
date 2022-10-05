version 1.0



            ############################################################
        ####################################################################
######################################################################################
#                                 pVACseq Maine Workflow                             #
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
        # FASTQ Files
        #~~~~~~~~~~~~
        File left
        File right
        
        #~~~~~~~~~~~~
        # VCF Files
        #~~~~~~~~~~~~
        File HaplotypeCaller_VCF
        File HaplotypeCaller_VCF_index
        File Mutect2_VCF
        File Mutect2_VCF_index

        String Tumor_ID
        String Normal_ID
        
        #~~~~~~~~~~~~
        # BAM Files
        #~~~~~~~~~~~~
        File BAM_Tumor
        File BAM_Tumor_index

        File BAM_RNA
        File BAM_RNA_index

        File BAM_Normal
        File BAM_Normal_index
        
        File GTF

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
            VCF                 = Mutect2_VCF,
            VEP_Reference       = VEP_Reference,

            ref_fasta           = ref_fasta,
            ref_fasta_index     = ref_fasta_index,
            ref_dict            = ref_dict,
            
            cpus                = cpus,
            preemptible         = preemptible,
            sample_id           = sample_id
    }

    #~~~~~~~~~~~~~~~~~~~~
    # Read Count Data
    #~~~~~~~~~~~~~~~~~~~~
    # Get the Bam read count information
    call annotateVCF.RunBamReadcount as TumorRunBamReadcount{
        input:
            VCF                 = RunAnnotateVEP.decomposed_VEP_output,
            VCF_index           = RunAnnotateVEP.decomposed_VEP_output_index,
            BAM                 = BAM_Tumor,
            BAM_index           = BAM_Tumor_index,
            ref_fasta           = ref_fasta,
            ref_fasta_index     = ref_fasta_index,
            ref_dict            = ref_dict,
            
            cpus                = 1,
            preemptible         = preemptible,
            sample_id           = Tumor_ID
    }
    call annotateVCF.RunBamReadcount as NormalRunBamReadcount{
        input:
            VCF                 = RunAnnotateVEP.decomposed_VEP_output,
            VCF_index           = RunAnnotateVEP.decomposed_VEP_output_index,
            BAM                 = BAM_Normal,
            BAM_index           = BAM_Normal_index,
            ref_fasta           = ref_fasta,
            ref_fasta_index     = ref_fasta_index,
            ref_dict            = ref_dict,
            
            cpus                = 1,
            preemptible         = preemptible,
            sample_id           = Normal_ID
    }

    call annotateVCF.RunBamReadcount as RnaRunBamReadcount{
        input:
            VCF                 = RunAnnotateVEP.decomposed_VEP_output,
            VCF_index           = RunAnnotateVEP.decomposed_VEP_output_index,
            BAM                 = BAM_RNA,
            BAM_index           = BAM_RNA_index,
            ref_fasta           = ref_fasta,
            ref_fasta_index     = ref_fasta_index,
            ref_dict            = ref_dict,
            
            cpus                = 1,
            preemptible         = preemptible,
            sample_id           = Tumor_ID
    }

    # Add the bam readcount information
    call annotateVCF.RunAddReadcount as TumorRunAddReadcount{
        input:
            VCF                 = RunAnnotateVEP.decomposed_VEP_output,
            readcount_indel     = TumorRunBamReadcount.sample_bam_readcount_indel,
            readcount_snv       = TumorRunBamReadcount.sample_bam_readcount_snv,
            RNAorDNA            = "DNA",

            cpus                = 1,
            preemptible         = preemptible,
            sample_id           = Tumor_ID
    }

    call annotateVCF.RunAddReadcount as NormalRunAddReadcount{
        input:
            VCF                 = TumorRunAddReadcount.annotated_vcf,
            readcount_indel     = NormalRunBamReadcount.sample_bam_readcount_indel,
            readcount_snv       = NormalRunBamReadcount.sample_bam_readcount_snv,
            RNAorDNA            = "DNA",

            cpus                = 1,
            preemptible         = preemptible,
            sample_id           = Normal_ID
    }

    call annotateVCF.RunAddReadcount as RnaRunAddReadcount{
        input:
            VCF                 = NormalRunAddReadcount.annotated_vcf,
            readcount_indel     = RnaRunBamReadcount.sample_bam_readcount_indel,
            readcount_snv       = RnaRunBamReadcount.sample_bam_readcount_snv,
            RNAorDNA            = "RNA",

            cpus                = 1,
            preemptible         = preemptible,
            sample_id           = Tumor_ID
    }

    #~~~~~~~~~~~~~~~~~~~~
    # Expression Data
    #~~~~~~~~~~~~~~~~~~~~
    call annotateVCF.RunExpressionData as RunExpressionData{
        input:
            left                = left,
            right               = right,
            VCF                 = RnaRunAddReadcount.annotated_vcf,
            BAM                 = BAM_RNA,
            BAM_index           = BAM_RNA_index,
            ref_fasta           = ref_fasta,
            ref_fasta_index     = ref_fasta_index,
            ref_dict            = ref_dict,
            
            Tumor_ID            = Tumor_ID,

            GTF                 = GTF,

            cpus                = cpus,
            preemptible         = preemptible,
            sample_id           = sample_id
    }

    ##~~~~~~~~~~~~~~~~~~~~
    ## Phasing 
    ##~~~~~~~~~~~~~~~~~~~~
    call phasing.RunCombineSomaticGermline as RunCombineSomaticGermline{
        input:
            HaplotypeCaller_VCF = HaplotypeCaller_VCF,
            HaplotypeCaller_VCF_index = HaplotypeCaller_VCF_index,
            Somatic_VCF         = Mutect2_VCF,
            Somatic_VCF_index   = Mutect2_VCF_index,
            ref_fasta           = ref_fasta,
            ref_fasta_index     = ref_fasta_index,
            ref_dict            = ref_dict,
            BAM                 = BAM_Tumor,
            BAM_index           = BAM_Tumor_index,
            Tumor_ID            = Tumor_ID,
    
            cpus                = cpus,
            preemptible         = preemptible,
            sample_id           = sample_id
    }
    
    call phasing.RunAnnotateCombinedVCF as RunAnnotateCombinedVCF{
        input:
            VCF                 = RunCombineSomaticGermline.combined_VCF,
            ref_fasta           = ref_fasta,
            ref_fasta_index     = ref_fasta_index,
            ref_dict            = ref_dict,
            VEP_Reference       = VEP_Reference,
    
            cpus                = cpus,
            preemptible         = preemptible,
            sample_id           = sample_id
    }


    output {

        File annotated_VCF = RunExpressionData.annotated_vcf
        File phasedVCF = RunAnnotateCombinedVCF.combined_annotated_VCF
        
    }
    
}
