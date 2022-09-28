#!/usr/bin/env python3


import os, sys, re
import logging
import argparse
import subprocess
import gzip
import multiprocessing
import time
import numpy as np


from annotation_VEP_expression import *


def main():

    ####################
    # Parse the use supplied information 
    ####################
    # Set the variables to supply 
    parser = argparse.ArgumentParser(description="Rsplit VCF.", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--VCF", type=str, required=True, help="VCF of interest.")
    parser.add_argument("--sample_id", type=str, required=True, help="the id for the sample.")
    parser.add_argument("--output", type=str, required=False, help="output directory.", default = ".")
    parser.add_argument("--cpu", type=int, required=False, help="Number of CPUs to use.", default = "8")
    parser.add_argument("--reference_fa", type=str, required=False, help="Reference used for VEP annotations.")
    parser.add_argument("--featureCounts_GX", type=str, required=False, help="Gene Expression txt file.")
    parser.add_argument("--featureCounts_TX", type=str, required=False, help="Transcript expression txt file.")

    # Parse the variables given 
    args = parser.parse_args()
    VCF = args.VCF
    output_path = args.output
    sample_id = args.sample_id
    cpu = args.cpu
    reference_fa = args.reference_fa
    featureCounts_GX = args.featureCounts_GX
    featureCounts_TX = args.featureCounts_TX

    if output_path ==  ".":
        output_path = os.getcwd()

    message_str = "\n####################################################################################\n\tAnnotating Expression Information\n####################################################################################"
    print(message_str)

    ##############################
    # Load Data
    ##############################
    # initiate the ViFi object 
    
    VCF = SplitVCF(input_vcf = VCF,  reference_fa = reference_fa, sample_id= sample_id,  output= output_path, cpu = cpu, 
        featureCounts_GX =featureCounts_GX, featureCounts_TX = featureCounts_TX)
    VCF = VCF.getIDs()
    VCF = VCF.getStats()
    VCF = VCF.getHeader()
    

    # #EXPRESSION
    VCF.expression()
    VCF.mergeVCFs(file_list = VCF.txgx_exp_files, output = "annotated_TXGX.vcf")
    VCF.removeVCFs(file_list = VCF.tx_exp_files)




if __name__ == "__main__":
    main()



