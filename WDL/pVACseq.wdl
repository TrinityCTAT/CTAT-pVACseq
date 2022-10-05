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
        String epitope_prediction_algorithms
        String? epitope_lengths_I
        String? epitope_lengths_II

        String sample_id
        String? Normal_ID

        Int cpus
        Int preemptible
        Int disk = ceil((size(VCF, "GB") * 3) + (size(phase_VCF, "GB") * 2) + 50)
    }

    command <<<
        set -e

        #~~~~~~~~~~~~~~~~~~~~~~~~
        # bgzip and index the input vcf and phased vcf
        #~~~~~~~~~~~~~~~~~~~~~~~~
        bgzip -c ~{VCF} > VEP_output.vcf.gz
        tabix -p vcf VEP_output.vcf.gz

        bgzip -c ~{phase_VCF} > phased.vcf.gz
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
            ~{"-e1 " + epitope_lengths_I } ~{"-e2 " + epitope_lengths_II } ~{"--normal-sample-name " + Normal_ID}

        tar cvf output.tar output

    >>>

    output {
        File pvacseq_output = "output.tar"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + disk + " HDD"
        docker: "brownmp/pvactools:devel"
        cpu: cpus
        memory: "50GB"
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
        String epitope_prediction_algorithms
        String? epitope_lengths_I
        String? epitope_lengths_II
      
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
        Int preemptible = 2
        
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
            epitope_prediction_algorithms   = epitope_prediction_algorithms,
            epitope_lengths_I               = epitope_lengths_I,
            epitope_lengths_II              = epitope_lengths_II,
    
            cpus                            = cpus,
            preemptible                     = preemptible,
            sample_id                       = sample_id,
            Normal_ID                       = Normal_ID
    }
}
