version 1.0



###########################################################
# Run Annotations: Combine Somatic and Germline
###########################################################
# Combine the Somatic and Germline VCF files for phasing pursposes 
##  Mutect2 (Somatic), HaplotypeCaller (Germline)
##  Combine the Somatic and Germiline variants to create the phasesd vcf. 
##  pVACseq will use this to address variants in the proximity of the variant it is assessing in isolation 
task RunCombineSomaticGermline{
    input {
        File HaplotypeCaller_VCF
        File HaplotypeCaller_VCF_index
        File Somatic_VCF
        File Somatic_VCF_index
        
        #~~~~~~~~~~~~
        # FASTQ Files
        #~~~~~~~~~~~~
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File BAM
        File BAM_index

        String Tumor_ID

        Int cpus
        Int preemptible
        String sample_id
    }

    command <<<
        set -e


        # Select the Tumor variants
        ## Mutect2 is run on Tumor-Normal mode so need to selct Tumor variants
        
        /usr/bin/java -jar /usr/GenomeAnalysisTK.jar \
            -T SelectVariants \
            -R ~{ref_fasta} \
            --variant ~{Somatic_VCF} \
            --sample_name ~{Tumor_ID} \
            -o tumor_only.vcf



        # Combine somatic and germline variants using GATK’s CombineVariants
        /usr/bin/java -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
            -R ~{ref_fasta} \
            --variant ~{HaplotypeCaller_VCF} \
            --variant tumor_only.vcf \
            --assumeIdenticalSamples \
            -o combined_somatic_plus_germline.vcf

        # Sort the combined VCF output using PICARD
        java -jar /usr/local/src/picard.jar SortVcf \
            I=combined_somatic_plus_germline.vcf \
            O=combined_somatic_plus_germline.sorted.vcf \
            SEQUENCE_DICTIONARY=~{ref_dict}


        # Phase variants using GATK’s ReadBackedPhasing
        /usr/bin/java -jar /usr/GenomeAnalysisTK.jar \
            -T ReadBackedPhasing \
            -R ~{ref_fasta} \
            -I ~{BAM} \
            --variant combined_somatic_plus_germline.sorted.vcf \
            -L combined_somatic_plus_germline.sorted.vcf \
            -o phased.vcf
  
    >>>

    output {
        File combined_VCF = "phased.vcf"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(HaplotypeCaller_VCF, "GB")*2 + size(Somatic_VCF, "GB")*2 + size(ref_fasta,"GB") + 100) + " HDD"
        docker: "mbrown/gatk3:devel"
        cpu: cpus
        memory: "50GB"
    }
}

###########################################################
# Run Annotations: Annotate Combine Somatic and Germline VCF
###########################################################
# Annotate the previously combined Somatic and Germline VCF file 
task RunAnnotateCombinedVCF{
    input {
        File VCF
        
        #~~~~~~~~~~~~
        # FASTQ Files
        #~~~~~~~~~~~~
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File VEP_Reference

        Int cpus
        Int preemptible
        String sample_id
    }

    command <<<
        set -e
        
        tar -xvf ~{VEP_Reference}

        #~~~~~~~~~~~~~~~~~~~~~~~
        # Run annotation step 
        #~~~~~~~~~~~~~~~~~~~~~~~
        /opt/ensembl-vep/vep \
            --input_file ~{VCF} \
            --output_file VEP_output.vcf \
            --format vcf \
            --vcf --symbol --terms SO --tsl \
            --hgvs \
            --fasta ~{ref_fasta} \
            --offline \
            --cache \
            --dir_cache `pwd` \
            --plugin Downstream \
            --plugin Wildtype \
            --dir_plugins /opt/VEP_plugins/

    >>>

    output {
        File combined_annotated_VCF = "VEP_output.vcf"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(VCF, "GB")*2 + size(ref_fasta,"GB") + 100) + " HDD"
        docker: "mbrown/pvactools:devel"
        cpu: cpus
        memory: "50GB"
    }
}