version 1.0




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Annotations: Add RNA-Editing to VCF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add the RNA-Editing Annotations to the VCF file 
task RunAnnotateRNAediting {
    input {
        File VCF
        File VCF_index
        File RNA_editing_VCF
        String sample_id


        Int preemptible
        Int cpus = 1
        Int disk = ceil((size(VCF, "GB") * 2) + size(RNA_editing_VCF, "GB") + 50)
    }

    command <<<
        set -ex

        echo "####### Creating Exon Bed File ########"

        # index the rna-editing file 
        tabix -p vcf ~{RNA_editing_VCF}

        bcftools annotate \
            --output-type z \
            --annotations ~{RNA_editing_VCF} \
            --columns "INFO/RNAEDIT" \
            --header-line '##INFO=<ID=RNAEDIT,Number=1,Type=String,Description="A known or predicted RNA-editing site">' \
            --output ~{sample_id}_annotated_RNAediting.vcf.gz \
            ~{VCF}

        tabix -p vcf ~{sample_id}_annotated_RNAediting.vcf.gz

    >>>
    output {
        File annotated_VCF       = "~{sample_id}_annotated_RNAediting.vcf.gz"
        File annotated_VCF_index = "~{sample_id}_annotated_RNAediting.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: "brownmp/pvactools:devel"
        memory: "4G"
        preemptible: preemptible
        cpus : cpus
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Annotations: Filter RNA-Editing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add the RNA-Editing Annotations to the VCF file 
task RunFilterRNAediting {
    input {
        File VCF
        File VCF_index
        String sample_id

        Int preemptible
        Int cpus
        Int disk = ceil((size(VCF, "GB") * 2) + 50)
    }

    command <<<
        set -ex

        echo "####### Filter the RNA-editing variants ########"

        bcftools filter \
            --include 'RNAEDIT="."' \
            --output ~{sample_id}_Filtered_RNAediting.vcf.gz \
            --threads ~{cpus} \
            ~{VCF}

        tabix -p vcf ~{sample_id}_Filtered_RNAediting.vcf.gz

    >>>
    output {
        File filtered_VCF       = "~{sample_id}_Filtered_RNAediting.vcf.gz"
        File filtered_VCF_index = "~{sample_id}_Filtered_RNAediting.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: "brownmp/pvactools:devel"
        memory: "4G"
        preemptible: preemptible
        cpus : cpus
    }
}





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Annotations: Create Exon Bed file 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create the exon bed file using a GTF file 

task RunCreateExonBED {
    input {
        File GTF

        Int preemptible
        Int cpus = 1
        Int disk = ceil((size(GTF, "GB") * 2) + 50)
    }

    command <<<
        set -ex

        echo "####### Creating Exon Bed File ########"

        python <<CODE
        
        import pandas as pd 
        import os, sys 

        # Read in the GTF file to a pandas DF
        df = pd.read_csv("~{GTF}", sep = "\t", comment="#", header = None)

        # Subset the GTF file
        ## Gather the Exons 
        tmp = df[(df[2] == "exon")]
        ## get the CHR:START-STOP
        output_bed = tmp.iloc[:,[0,3,4]]
        ## Filter Dublications 
        output_bed = output_bed.drop_duplicates()
        output_bed

        # Save the Output 
        output_bed.to_csv("Exon_bed.bed", sep = "\t", index = False, header = False)

        CODE
    >>>
    output {
        File exon_bed       = "Exon_bed.bed"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: "brownmp/pvactools:devel"
        memory: "4G"
        preemptible: preemptible
        cpus : cpus
    }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Annotations: Filter to only include exon regions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create the exon bed file using a GTF file 

task RunFilterNonExons {
    input {
        File Exon_BED_File
        File VCF
        File VCF_index
        String sample_id

        Int preemptible
        Int cpus = 1
        Int disk = ceil(size(Exon_BED_File, "GB") + (size(VCF, "GB") * 2) + 50)
    }

    command <<<
        set -ex
        echo "####### Filter non-exon region ########"
        tabix -h -R ~{Exon_BED_File} ~{VCF} > ~{sample_id}_exon_filtered.vcf

        bcftools sort ~{sample_id}_exon_filtered.vcf > ~{sample_id}_exon_filtered_sorted.vcf
        bgzip -c ~{sample_id}_exon_filtered_sorted.vcf > ~{sample_id}_exon_filtered_sorted.vcf.gz
        tabix -p vcf ~{sample_id}_exon_filtered_sorted.vcf.gz
    >>>
    output {
        File exon_filtered_VCF       = "~{sample_id}_exon_filtered_sorted.vcf.gz"
        File exon_filtered_VCF_index = "~{sample_id}_exon_filtered_sorted.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: "brownmp/pvactools:devel"
        memory: "4G"
        preemptible: preemptible
        cpus : cpus
    }
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Annotations: Add gnomAD to VCF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adding gnomAD population Allele frequencies to your VCF
## can determine somatic and germiline variants 
task RunAnnotateGnomad {
    input {
        File VCF
        File VCF_index
        File gnomad_vcf
        File gnomad_vcf_index

        Int preemptible
        Int cpus = 1
        Int disk = ceil((size(VCF, "GB") * 2) + size(gnomad_vcf, "GB") + 50)
    }

    command <<<
        set -ex

        echo "####### Annotate gnomAD ########"

        bcftools annotate \
            --output-type z \
            --annotations ~{gnomad_vcf} \
            --columns "INFO/gnomad_RS,INFO/gnomad_AF" \
            --output gnomAD_annotated.vcf.gz \
            ~{VCF}


        #~~~~~~~~~~~~~~~~~~~~~
        # Filter 
        #~~~~~~~~~~~~~~~~~~~~~
        python <<CODE
        from pysam import VariantFile
        import os, sys 

        vcf_in = VariantFile("gnomAD_annotated.vcf.gz", "rb")
        vcf_out = VariantFile('gnomAD_filtered.vcf', 'w', header=vcf_in.header)

        counter, total = 0, 0
        for rec in vcf_in.fetch():
            total += 1

            vcf_out.write(rec)
            counter+=1

            try:
                if rec.info["gnomad_AF"][0] <= .009:
                    vcf_out.write(rec)
                    counter+=1
            except:
                vcf_out.write(rec)
                counter+=1
        CODE

        bgzip -c gnomAD_filtered.vcf > gnomAD_filtered.vcf.gz
        tabix -p vcf gnomAD_filtered.vcf.gz
    >>>
    output {
        File unfiltered_gnomad_annotated_VCF    = "gnomAD_annotated.vcf.gz"
        File gnomad_annotated_VCF               = "gnomAD_filtered.vcf.gz"
        File gnomad_annotated_VCF_index         = "gnomAD_filtered.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: "brownmp/pvactools:devel"
        memory: "4G"
        preemptible: preemptible
        cpus : cpus
    }

}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Annotations: GT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adding genotype sample information to your VCF
## pVACseq requires that the input VCF contains sample genotype information (GT field), 
## which identifies whether or not a variant was called in a specific sample of interest.

task RunAnnotateGenotypeSample{
    input {
        File VCF

        Int cpus
        Int preemptible
        String sample_id
    }

    command <<<
        set -e

        vcf-genotype-annotator \
            ~{VCF} \
            ~{sample_id} 0/1 \
            -o /usr/local/src/output/gt_annotated_vcf.vcf
    >>>

    output {
        
        File gt_annotated_vcf = "/usr/local/src/output/gt_annotated_vcf.vcf"
        
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(VCF, "GB")*2 + 100) + " HDD"
        docker: "brownmp/pvactools:devel"
        cpu: cpus
        memory: "20GB"
    }
}

###########################################################
# Run Annotations: Annotating your VCF with VEP
###########################################################
# The input to the pVACseq pipeline is a VEP-annotated VCF.
#          Variant Effect Predictor (VEP)
## This will add consequence, transcript, and gene information to your VCF.

task RunAnnotateVEP{
    input {
        File VCF
        File VEP_Reference
        
        #~~~~~~~~~~~~
        # FASTQ Files
        #~~~~~~~~~~~~
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int cpus
        Int preemptible
        String sample_id
    }

    command <<<
        set -e

        #~~~~~~~~~~~~~~~~~~~~~~~
        # untar the references 
        #~~~~~~~~~~~~~~~~~~~~~~~
        tar -xvf ~{VEP_Reference}
        rm ~{VEP_Reference}

        /usr/local/src/annotate_VEP.py\
            --VCF ~{VCF} \
            --sample_id ~{sample_id} \
            --output `pwd` \
            --cpu ~{cpus} \
            --reference_fa ~{ref_fasta}

    >>>

    output {
        File VEP_output = "annotated_VEP.vcf"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(ref_fasta, "GB") + size(VEP_Reference, "GB") + size(VCF, "GB")*2 + 100) + " HDD"
        docker: "brownmp/pvactools:devel"
        cpu: cpus
        memory: "50GB"
    }
}

task RunDecompose{
    input {
        File VCF

        Int preemptible
        String sample_id
    }

    command <<<
        set -e

        #~~~~~~~~~~~~~~~~
        # Decompose 
        #~~~~~~~~~~~~~~~~
        # bam-readcount needs to be run separately for snvs and indels
        # Need to split multi-allelic sites
        /opt/vt/vt decompose \
            -s ~{VCF} \
            -o ~{sample_id}_decomposed_output.vcf


        # Index the output VCF file 
        bgzip -c ~{sample_id}_decomposed_output.vcf > ~{sample_id}_decomposed_output.vcf.gz
        tabix -p vcf ~{sample_id}_decomposed_output.vcf.gz
    >>>

    output {
        File decomposed_VCF = "~{sample_id}_decomposed_output.vcf.gz"
        File decomposed_VCF_index   = "~{sample_id}_decomposed_output.vcf.gz.tbi"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(VCF, "GB")*2 + 100) + " HDD"
        docker: "brownmp/pvactools:devel"
        memory: "50GB"
    }
}





###########################################################
# Run Annotations: bam_readcount
###########################################################

task RunBamReadcount{
    input {
        File VCF
        File VCF_index
        File BAM
        File BAM_index
        
        #~~~~~~~~~~~~
        # FASTQ Files
        #~~~~~~~~~~~~
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int cpus
        Int preemptible
        String sample_id
    }

    command <<<
        set -ea

        cp ~{BAM_index} .

        /usr/bin/python /usr/bin/bam_readcount_helper.py \
            ~{VCF} \
            ~{sample_id} \
            ~{ref_fasta} \
            ~{BAM} \
            NOPREFIX \
            `pwd`

        
    >>>

    output {
        File sample_bam_readcount_indel = "~{sample_id}_bam_readcount_indel.tsv"
        File sample_bam_readcount_snv   = "~{sample_id}_bam_readcount_snv.tsv"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(VCF, "GB")*2 + 100) + " HDD"
        docker: "mgibio/bam_readcount_helper-cwl:1.1.1"
        cpu: cpus
        memory: "50GB"
    }
}

###########################################################
# Run Annotations: Read Counts 
###########################################################
task RunAddReadcount{
    input {
        File VCF
        File readcount_indel
        File readcount_snv
        String RNAorDNA

        Int cpus
        Int preemptible
        String sample_id
    }

    command <<<
        set -e

        #~~~~~~~~~~~~~~~~~~~
        # Read counts 
        #~~~~~~~~~~~~~~~~~~~
        # The readcounts for snvs and indels are then added to your VCF separately, by running the vcf-readcount-annotator twice.
        ## ADd the SNV and INDELS seperaetly 
        # SNV
        vcf-readcount-annotator \
            ~{VCF} \
            ~{readcount_snv} \
            ~{RNAorDNA} \
            -s ~{sample_id} \
            -t snv \
            -o annotated_bam_readcount_snv.vcf

        # INDELS 
        vcf-readcount-annotator \
            annotated_bam_readcount_snv.vcf \
            ~{readcount_indel} \
            ~{RNAorDNA} \
            -s ~{sample_id} \
            -t indel \
            -o annotated_bam_readcount.vcf


        # index 
        bgzip -c annotated_bam_readcount.vcf > annotated_bam_readcount.vcf.gz
        tabix -p vcf annotated_bam_readcount.vcf.gz
        
    >>>

    output {
        File annotated_vcf = "annotated_bam_readcount.vcf.gz"
        File annotated_vcf_index = "annotated_bam_readcount.vcf.gz.tbi"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(VCF, "GB")*2 + size(readcount_indel,"GB") + size(readcount_snv,"GB") + 100) + " HDD"
        docker: "brownmp/pvactools:devel"
        cpu: cpus
        memory: "10GB"
    }
}



###########################################################
# Run Annotations: Expression Data
###########################################################
task RunExpressionData{
    input {
        File VCF
        
        #~~~~~~~~~~~~
        # FASTQ Files
        #~~~~~~~~~~~~
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File BAM
        File BAM_index
        File GTF

        String Tumor_ID

        Int cpus
        Int preemptible
        String sample_id
    }

    command <<<
        set -e

        cp ~{BAM_index} .

        ######################################
        # Expression data 
        ######################################
        # GTF/GFF files define genomic regions covered by different types of genomic features, 
        #  e.g. genes, transcripts, exons, or UTRs.
        #
        ## If the VCF is a single-sample VCF, pVACseq assumes that this sample is the tumor sample. 
        ## If the VCF is a multi-sample VCF, pVACseq will look for the sample using the sample_name 
        ##  parameter and treat that sample as the tumor sample.



        #~~~~~~~~~~~~~~~~~~~~
        # Subread-featureCounts: Transcript expression 
        #~~~~~~~~~~~~~~~~~~~~
        # use Subreads featureCounts to obtain the Trascript level counts 
        Rscript /usr/local/src/runFeatureCounts.R \
            --BAM ~{BAM} \
            --GTF ~{GTF} \
            --Reference ~{ref_fasta} \
            --threads ~{cpus} \
            --Type transcript_id \
            --output featureCounts_output.TX.txt

        #~~~~~~~~~~~~~~~~~~~~
        # Subread-featureCounts: Gene expression 
        #~~~~~~~~~~~~~~~~~~~~
        # use Subreads featureCounts to obtain the gene level counts 
        Rscript /usr/local/src/runFeatureCounts.R \
            --BAM ~{BAM} \
            --GTF ~{GTF} \
            --Reference ~{ref_fasta} \
            --threads ~{cpus} \
            --Type gene_id \
            --output featureCounts_output.GX.txt


    >>>

    output {
        File featureCounts_TX = "featureCounts_output.TX.txt"
        File featureCounts_GX = "featureCounts_output.GX.txt"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(VCF, "GB")*2 + size(BAM, "GB")*2 + size(ref_fasta,"GB") + 100) + " HDD"
        docker: "brownmp/pvactools:devel"
        cpu: cpus
        memory: "50GB"
    }
}


task RunAddExpressionData{
    input {
        
        File VCF
        
        #~~~~~~~~~~~~
        # FASTQ Files
        #~~~~~~~~~~~~
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File featureCounts_GX
        File featureCounts_TX

        String Tumor_ID

        Int cpus
        Int preemptible
        String sample_id
        Int disk = ceil( (size(VCF, "GB") * 2) + size(ref_fasta) + size(featureCounts_GX, "GB")*2 + 50)
    }

    command <<<
        set -e

        ######################################
        # Expression data 
        ######################################
        # GTF/GFF files define genomic regions covered by different types of genomic features, 
        #  e.g. genes, transcripts, exons, or UTRs.
        #
        ## If the VCF is a single-sample VCF, pVACseq assumes that this sample is the tumor sample. 
        ## If the VCF is a multi-sample VCF, pVACseq will look for the sample using the sample_name 
        ##  parameter and treat that sample as the tumor sample.



        #bgzip -c decomposed_VEP_output.vcf > decomposed_VEP_output.vcf.gz
        #tabix -p vcf /usr/local/src/decomposed_VEP_output.vcf.gz

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Annotate Expression Data 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        /usr/local/src/annotate_expression.py \
            --VCF ~{VCF} \
            --sample_id ~{Tumor_ID} \
            --cpu ~{cpus} \
            --reference_fa ~{ref_fasta} \
            --output . \
            --featureCounts_GX ~{featureCounts_GX} \
            --featureCounts_TX ~{featureCounts_TX}
    >>>

    output {
        File annotated_vcf = "annotated_TXGX.vcf"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + disk + " HDD"
        docker: "brownmp/pvactools:devel"
        cpu: cpus
        memory: "50GB"
    }
}





