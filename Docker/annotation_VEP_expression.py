#!/usr/bin/env python3

import os, sys, re
import logging
import argparse
import subprocess
import gzip
import multiprocessing
import time
import numpy as np

## Set up the logging  
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)



def runVEPannotator(input_vcf, idx, count, header, sample_id, reference):
    '''
    Open the VCF fiel and extract the lines of interest
        write these lines to a new vcf (subsetted)
        run the expression analysis on this subsetted file 
        then remove the temp subsetted files 
    
    '''
    with gzip.open(input_vcf, 'r') as vcf:
        # read in the specific lines 
        output_vcf_lines = vcf.readlines()[idx[0]:(idx[-1]+1)]

        file_name = f"tmp{count}.vcf"
        with open(file_name, "w") as tmp:
            # Write the header to the file 
            for header_line in header:
                tmp.write(header_line)#.encode())
            # Write the output to the temp file 
            for line in output_vcf_lines:
                # print(line)
                tmp.write(line.decode())

    # Bgzip and index the subsetted file 
    cmd = f'''bgzip -c {file_name} > {file_name}.gz'''
    # print(cmd)
    subprocess.run(cmd, shell=True)
    cmd = f'''tabix -p vcf {file_name}.gz'''
    # print(cmd)
    subprocess.run(cmd, shell=True)

    cmd = f'''
        /opt/ensembl-vep/vep \
            --input_file {tmp.name}.gz \
            --output_file {count}_VEP_output.vcf \
            --format vcf \
            --vcf --symbol --terms SO --tsl \
            --hgvs \
            --fasta {reference} \
            --offline \
            --cache \
            --dir_cache `pwd` \
            --plugin Frameshift \
            --plugin Wildtype \
            --dir_plugins /opt/VEP_plugins/
        '''
    # print(cmd)
    subprocess.run(cmd, shell=True)


    cmd = f"rm {file_name}*"
    subprocess.run(cmd, shell=True)




def runExpressionAnnotator(input_vcf, idx, count, header, sample_id, featureCounts_GX, featureCounts_TX):
    '''
    Open the VCF fiel and extract the lines of interest
        write these lines to a new vcf (subsetted)
        run the expression analysis on this subsetted file 
        then remove the temp subsetted files 
    
    '''
    with gzip.open(input_vcf, 'r') as vcf:
        # read in the specific lines 
        output_vcf_lines = vcf.readlines()[idx[0]:(idx[-1]+1)]

        #with tempfile.NamedTemporaryFile() as tmp:
        file_name = f"tmp{count}.vcf"
        with open(file_name, "w") as tmp:
            
            # Write the header to the file 
            for header_line in header:
                tmp.write(header_line)#.encode())

            # Write the output to the temp file 
            for line in output_vcf_lines:
                # print(line)
                tmp.write(line.decode())

    # print(tmp.name)
    # tmp.write(output.encode())
    cmd = f'''bgzip -c {file_name} > {file_name}.gz'''
    print(cmd)
    subprocess.run(cmd, shell=True)
    cmd = f'''tabix -p vcf {file_name}.gz'''
    print(cmd)
    subprocess.run(cmd, shell=True)

    cmd = f'''
    vcf-expression-annotator\n\t{tmp.name}.gz\n\tfeatureCounts_output.TX.txt\n\tcustom\n\ttranscript\n\t--id-column Geneid\n\t--expression-column TPM\n\t--ignore-ensembl-id-version\n\t-o {count}_annotated.tx.vcf\n\t-s {sample_id}
    '''
    print(cmd)

    # Add Transcript Expression
    cmd = f'''
    vcf-expression-annotator \
        {file_name}.gz \
        {featureCounts_TX} \
        custom \
        transcript \
        --id-column Geneid \
        --expression-column TPM \
        --ignore-ensembl-id-version \
        -o {count}_annotated.tx.vcf \
        -s {sample_id}
    '''
    subprocess.run(cmd, shell=True)



    cmd = f'''
    vcf-expression-annotator\n\t{count}_annotated.tx.vcf\n\tfeatureCounts_output.GX.txt\n\tcustom\n\tgene\n\t--id-column Geneid\n\t--expression-column TPM\n\t--ignore-ensembl-id-version\n\t-o {count}_annotated.tx.gx.vcf\n\t-s {sample_id}
    '''
    print(cmd)
    # Add Gene Expression
    cmd = f'''
    vcf-expression-annotator \
        {count}_annotated.tx.vcf \
        {featureCounts_GX} \
        custom \
        gene \
        --id-column Geneid \
        --expression-column TPM \
        --ignore-ensembl-id-version \
        -o {count}_annotated.tx.gx.vcf \
        -s {sample_id}
    '''
    subprocess.run(cmd, shell=True)


    message_str = "\t Removing TMP files..."
    cmd = f'''rm {tmp.name}*'''
    print("\t\t",cmd)
    subprocess.run(cmd, shell=True)



class SplitVCF:
    '''
    Class to annotate a given vcf in a multiprocessing fashion 
    '''
    # initialize object
    def __init__(   self, 
                    input_vcf,
                    reference_fa,
                    sample_id,
                    output,
                    cpu,
                    featureCounts_GX = "",
                    featureCounts_TX = ""
                    ): # arguments to class instantiation 
        
        self.input_vcf    = input_vcf
        self.output       = output
        self.sample_id    = sample_id
        self.cpu          = cpu
        self.reference_fa = reference_fa
        self.featureCounts_GX = featureCounts_GX
        self.featureCounts_TX = featureCounts_TX

        message_str = f"Annotating VCF:\n\tVCF : {input_vcf}\n\tSample ID : {sample_id}\n\tCPU count : {cpu}\n\t reference : {reference_fa}"
        logger.info(message_str)

    def getIDs( self ):
        #####################
        # Get contig IDs
        #####################
        # Get the contig IDs that are present in the VCF file 

        message_str = f"\t\tGetting contig IDs."
        logger.info(message_str)

        # Get contig IDs from the VCF header 
        tmp = "#"
        vcf_head = []
        with gzip.open(self.input_vcf, 'rb') as vcf:
            while tmp == "#":
                line = next(vcf).decode('ASCII')
                tmp = line[0]
                vcf_head.append(line)

        r = re.compile("##contig=<ID=*")
        long_IDS = list(filter(r.match, vcf_head))

        IDs = [i.split(",")[0].split("##contig=<ID=")[1] for i in long_IDS]


        self.IDs  = IDs
        return self 

    def getStats( self ):
        # create the dic to hold the counts 
        # contig_dic = {i:0 for i in self.IDs}
        # test = 0
        # with gzip.open(self.input_vcf, 'rb') as vcf:
        #     for line in vcf:
        #         line = line.decode()
        #         if line[0] != "#":
        #             # Add to the dic
        #             val = contig_dic[line.split("\t")[0]] + 1
        #             contig_dic[line.split("\t")[0]] = val

        test = 0
        with gzip.open(self.input_vcf, 'rb') as vcf:
            for line in vcf:
                test+=1
        self.stats = test
        return self 

    def VEP( self ):
        # Apply annotations to the VCF
        ## Split the VCF file by thier contig 
        message_str = f"\tAdding VEP INFO"
        logger.info(message_str)

        idx_range = list(range(len(self.header), self.stats))
        idx_list = np.array_split(idx_range, self.cpu)

        pool = multiprocessing.Pool(self.cpu)

        # get the start time
        start_time = time.process_time()
        message_str = f"\t\tStart Time: {time.asctime( time.localtime(time.time()) )}"
        logger.info(message_str)

        args = [[self.input_vcf, list(i), j, self.header, self.sample_id, self.reference_fa] for j,i in enumerate(idx_list) ]
        pool.starmap(runVEPannotator, args)

        pool.close()

        # get the end time
        end_time = time.process_time()
        message_str = f"\t\tEnd Time: {time.asctime( time.localtime(time.time()) )}"
        logger.info(message_str)
        # get execution time
        time_total = end_time - start_time
        message_str = f"\t\tProcess time: {time_total}"
        logger.info(message_str)

               # save the arguments 
        file_list = []
        for count,i in enumerate(idx_list):
            file_list.append(f"{count}_VEP_output.vcf")

        # SAVE files created 
        self.VEP_files = file_list

        return self 

    def expression( self ):
        # Apply Expession (gene and transcript) to the VCF
        ## Split the VCF file by thier contig
        message_str = f"\tAdding Expression INFO"
        logger.info(message_str)

        idx_range = list(range(len(self.header), self.stats))
        idx_list = np.array_split(idx_range, self.cpu)

        pool = multiprocessing.Pool(self.cpu)
        
        
    
        # get the start time
        start_time = time.process_time()
        message_str = f"\t\tStart Time: {time.asctime( time.localtime(time.time()) )}"
        logger.info(message_str)

        args = [[self.input_vcf, list(i), j, self.header, self.sample_id, self.featureCounts_GX, self.featureCounts_TX] for j,i in enumerate(idx_list) ]
        
        pool.starmap(runExpressionAnnotator, args)

        pool.close()

        # get the end time
        end_time = time.process_time()
        message_str = f"\t\tEnd Time: {time.asctime( time.localtime(time.time()) )}"
        logger.info(message_str)
        # get execution time
        time_total = end_time - start_time
        message_str = f"\t\tProcess time: {time_total}"
        logger.info(message_str)

        # save the arguments 
        file_list = []
        for count,i in enumerate(idx_list):
            file_list.append(f"{count}_annotated.tx.gx.vcf")

        self.txgx_exp_files = file_list

        all_file_list = []
        for count,i in enumerate(idx_list):
            all_file_list.append(f"{count}_annotated.tx.gx.vcf")
            all_file_list.append(f"{count}_annotated.tx.vcf")

        self.tx_exp_files = all_file_list

        return self 


    def getHeader( self ):
        #--------------------------------------------------------
        # Get the header for the vcf file and hold it as a string
        #-------------------------------------------------------- 
        vcf_head = []
        tmp = "#"
        with gzip.open(self.input_vcf, 'rb') as vcf:
            while tmp == "#":
                line = next(vcf).decode('ASCII')
                tmp = line[0]
                if tmp == "#": # had trouble with my while loop
                    vcf_head.append(line)
        self.header = vcf_head 
        return self 

    

    def mergeVCFs( self, file_list, output):
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # Merge the VCFs created in expression annotaion step 
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        message_str = f"\t\tMerging VCFs."
        logger.info(message_str)

        vcfs = [f"I={i}" for i in file_list ]

        vcfs = " ".join(vcfs)

        cmd = f"""java -jar /usr/local/src/picard.jar MergeVcfs {vcfs} O={output}"""

        print(cmd.replace(" ", " \n"))

        subprocess.run(cmd, shell=True)


    def removeVCFs( self, file_list):
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # Remove the unneeded VCFs created in expression annotaion step 
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        message_str = f"\tRemoving unneeded VCFs."
        logger.info(message_str)

        # remove the intermediate files 
        vcfs_remove = " ".join(file_list)

        cmd = f"""rm {vcfs_remove}"""

        subprocess.run(cmd, shell=True)

