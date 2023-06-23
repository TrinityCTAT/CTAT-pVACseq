version 1.0


###########################################################
###########################################################
# TASKs
###########################################################
###########################################################
# 
task RunpVACseq{
    input {
        File VCF
        File phase_VCF

        String HLAs
        String fasta_size = 200
        String epitope_prediction_algorithms
        String? epitope_lengths_I
        String? epitope_lengths_II

        String sample_id
        String? Normal_ID

        Int cpus
        Int preemptible
        Int disk = ceil((size(VCF, "GB") * 3) + (size(phase_VCF, "GB") * 2) + 50)
    }

    command {
        set -e
        set -v

        #~~~~~~~~~~~~~~~~~~~~~~~~
        # bgzip and index the input vcf and phased vcf
        #~~~~~~~~~~~~~~~~~~~~~~~~

        if file ~{VCF} | grep -q "Variant Call Format (VCF)"; then
            bgzip -c ~{VCF} > VEP_output.vcf.gz
        elif file ~{VCF} | grep -q "gzip"; then
            mv ~{VCF} VEP_output.vcf.gz    
        fi
        tabix -p vcf VEP_output.vcf.gz

        if file ~{phase_VCF} | grep -q "Variant Call Format (VCF)"; then
            bgzip -c ~{phase_VCF} > phased.vcf.gz
        elif file ~{phase_VCF} | grep -q "gzip"; then
            mv ~{phase_VCF} phased.vcf.gz
        fi
        tabix -p vcf phased.vcf.gz


        # Need to change the tmp directory 
        # Setting the TMPDIR environment variable to a long path will cause an error in Python mulitprocessing library.



        
        tmpDir="/tmp"
        echo "$tmpDir"
        #chmod 777 "$tmpDir"
        #export _JAVA_OPTIONS=-Djava.io.tmpdir="$tmpDir"
        export TMPDIR="$tmpDir"

        mkdir output

        ## possible HLA typers 
        # MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \

        #~~~~~~~~~~~~~~~~~~~~~~~~
        # RUN pvacseq
        #~~~~~~~~~~~~~~~~~~~~~~~~
        pvacseq run \
            VEP_output.vcf.gz \
            ~{sample_id} \
            ~{HLAs} \
            ~{epitope_prediction_algorithms} \
            `pwd`/output \
            --phased-proximal-variants-vcf phased.vcf.gz \
            --iedb-install-directory /opt/iedb \
            -t ~{cpus} \
            --fasta-size ~{fasta_size} \
            ~{"-e1 " + epitope_lengths_I } ~{"-e2 " + epitope_lengths_II } ~{"--normal-sample-name " + Normal_ID}

        tar cvf output.tar output

    }

    output {
        File pvacseq_output   = "output.tar"
        File pvacseq_epitopes = "output/MHC_Class_I/~{sample_id}.all_epitopes.tsv"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + disk + " HDD"
        docker: "cgeorgescu/pvactools:devel"
        cpu: cpus
        memory: "50GB"
    }
}


task RunBindingFilterpVACseq {
    input {
        File epitopes_tsv
        String sample_id

        Int preemptible
        Int cpus = 1
        Int disk = ceil((size(epitopes_tsv, "GB") * 2) + 50)
    }

    command <<<
        set -ex

        echo "####### Creating Exon Bed File ########"

        pvacseq binding_filter \
            --exclude-NAs \
            ~{epitopes_tsv} \
            ~{sample_id}_binding_filter_noNAs.tsv
    >>>
    output {
        File binding_filter       = "~{sample_id}_binding_filter_noNAs.tsv"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: "cgeorgescu/pvactools:devel"
        memory: "4G"
        preemptible: preemptible
        cpus : cpus
    }
}

task RunCoverageFilterpVACseq {
    input {
        File epitopes_tsv
        String sample_id

        Int preemptible
        Int cpus = 1
        Int disk = ceil((size(epitopes_tsv, "GB") * 2) + 50)
    }

    command <<<
        set -ex

        echo "####### Creating Exon Bed File ########"

        pvacseq coverage_filter \
            --exclude-NAs \
            ~{epitopes_tsv} \
            ~{sample_id}_coverage_binding_filter_noNAs.tsv
    >>>
    output {
        File coverage_filter       = "~{sample_id}_coverage_binding_filter_noNAs.tsv"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: "cgeorgescu/pvactools:devel"
        memory: "4G"
        preemptible: preemptible
        cpus : cpus
    }
}


task RunTranscriptSupportFilterpVACseq {
    input {
        File epitopes_tsv
        String sample_id

        Int preemptible
        Int cpus = 1
        Int disk = ceil((size(epitopes_tsv, "GB") * 2) + 50)
    }

    command <<<
        set -ex

        echo "####### Creating Exon Bed File ########"

        pvacseq transcript_support_level_filter \
            --exclude-NAs \
            ~{epitopes_tsv} \
            ~{sample_id}_transcript_support_coverage_binding_filter_noNAs.tsv
    >>>
    output {
        File TSL_filter       = "~{sample_id}_transcript_support_coverage_binding_filter_noNAs.tsv"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: "cgeorgescu/pvactools:devel"
        memory: "4G"
        preemptible: preemptible
        cpus : cpus
    }
}





task RunFilterRNA {
    input {
        File epitopes_tsv
        String sample_id

        Int preemptible
        Int cpus = 1
        Int disk = ceil((size(epitopes_tsv, "GB") * 2) + 50)
        Int memory_gb = ceil((size(epitopes_tsv, "GB") * 2) + 4)
    }

    command <<<
        set -ex

        echo "####### Filter RNA Epitopes ########"
        
        Rscript /usr/local/src/Filter.R \
            --file ~{epitopes_tsv} \
            --out ~{sample_id}_Filtered.tsv
    >>>

    output {
        File TSL_filter       = "~{sample_id}_Filtered.tsv"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: "cgeorgescu/pvactools:devel"
        memory: memory_gb + "GiB"
        preemptible: preemptible
        cpus : cpus
    }
}


            ############################################################
        ####################################################################
######################################################################################
#                                      pVACseq                                       #
######################################################################################
        #####################################################################
            ############################################################

# pVACseq is a cancer immunotherapy pipeline for the identification of personalized Variant Antigens by Cancer Sequencing (pVACseq) 
# that integrates tumor mutation and expression data (DNA- and RNA-Seq). It enables cancer immunotherapy research by using massively 
# parallel sequence data to predicting tumor-specific mutant peptides (neoantigens) that can elicit anti-tumor T cell immunity

##########################
# Workflow
##########################
workflow pVACseq {
    input {

        #~~~~~~~~~~~~
        # Sample ID
        #~~~~~~~~~~~~
        String sample_id
        String? Normal_ID
        String HLAs
        String? fasta_size
        String epitope_prediction_algorithms
        String? epitope_lengths_I
        String? epitope_lengths_II

        #~~~~~~~~~~~~
        # DNA or RNA
        #~~~~~~~~~~~~
        String Type
      
        #~~~~~~~~~~~~
        # VCF Files
        #~~~~~~~~~~~~
        File VCF
        File phased_VCF
        

        #~~~~~~~~~~~~
        # CPU count 
        #~~~~~~~~~~~~
        Int cpus = 10

        #~~~~~~~~~~~~
        # general runtime settings
        #~~~~~~~~~~~~
        Int preemptible = 0
        
    }

    parameter_meta {
        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}
        cpus:{help:"CPU count"}
        #docker:{help:"Docker image"}
    }


    ########################
    ## CALLS
    ########################
    call RunpVACseq as RunpVACseq{
        input:
            VCF                             = VCF,
            phase_VCF                       = phased_VCF,
            HLAs                            = HLAs,
            fasta_size                      = fasta_size,
            epitope_prediction_algorithms   = epitope_prediction_algorithms,
            epitope_lengths_I               = epitope_lengths_I,
            epitope_lengths_II              = epitope_lengths_II,
    
            cpus                            = cpus,
            preemptible                     = preemptible,
            sample_id                       = sample_id,
            Normal_ID                       = Normal_ID
    }

    #~~~~~~~~~~~~~~~~~~~~~~
    #   DNA Filter
    #~~~~~~~~~~~~~~~~~~~~~~
    if (Type == "DNA"){
        call RunBindingFilterpVACseq as RunBindingFilterpVACseq{
            input:
                epitopes_tsv                    = RunpVACseq.pvacseq_epitopes,
                sample_id                       = sample_id,
        
                cpus                            = cpus,
                preemptible                     = preemptible,
                sample_id                       = sample_id
        }

        call RunCoverageFilterpVACseq as RunCoverageFilterpVACseq{
            input:
                epitopes_tsv                    = RunBindingFilterpVACseq.binding_filter,
                sample_id                       = sample_id,
        
                cpus                            = cpus,
                preemptible                     = preemptible,
                sample_id                       = sample_id
        }

        call RunTranscriptSupportFilterpVACseq as RunTranscriptSupportFilterpVACseq{
            input:
                epitopes_tsv                    = RunCoverageFilterpVACseq.coverage_filter,
                sample_id                       = sample_id,
        
                cpus                            = cpus,
                preemptible                     = preemptible,
                sample_id                       = sample_id
        }
    }

    #~~~~~~~~~~~~~~~~~~~~~~~
    #   RNA Filter
    #~~~~~~~~~~~~~~~~~~~~~~~
    if (Type == "RNA"){

        call RunFilterRNA as RunFilterRNA{
                input:
                    epitopes_tsv                    = RunpVACseq.pvacseq_epitopes,
                    sample_id                       = sample_id,
            
                    cpus                            = cpus,
                    preemptible                     = preemptible,
                    sample_id                       = sample_id
        }
    }

    output {
        File? pvacseq_epitopes = RunpVACseq.pvacseq_epitopes
        File? TSL_filter_RNA = RunFilterRNA.TSL_filter
        File? TSL_filter_DNA = RunTranscriptSupportFilterpVACseq.TSL_filter
    }
}


