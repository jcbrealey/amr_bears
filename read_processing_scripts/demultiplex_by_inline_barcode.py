#!/usr/bin/python
from sys import argv
from itertools import izip, islice
import gzip
"""
author: Tom van der Valk, tom.vandervalk@ebc.uu.se
edited: Jaelle Brealey
"""

# TAKES zipped files!! Outputs unzipped
# Edited from demultiplex_by_index_and_barcode_double_3of4_notrim_181222.py
# Assumes samples have already been demultiplexed by SciLifeLab with index, and barcodes only are required
# Note given large number of pooled samples (200 odd) need both barcodes to assign to sample

def demultiplex_by_index(FORWARD_FILE,REVERSE_FILE,INDEX_FILE,indexID):
    
    bar_dict_F = {}
    bar_dict_R = {}
    bar_dict_F_V1 = {}
    bar_dict_R_V1 = {}
    nucleotides = ["N","T","C","G","A"]
    rev_complement = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
    sample_dict_V1 = {}
    sample_dict = {}
 
    with open(INDEX_FILE) as f:
	next(f)
	for line in f:
	    splitted = line.strip().split(",")
	    sample,barcode_F,barcode_R,P5,P7,barF_ID,barR_ID,P5_ID,P7_ID = splitted[0],splitted[1],splitted[2],splitted[3],\
		splitted[4],splitted[5],splitted[6],splitted[7],splitted[8]
            if P5_ID != P7_ID:
                print("Error: unexpected P5 and P7 index combination (%s and %s)" % (P5_ID, P7_ID))
                break

            elif P5_ID == indexID:
                # add to dictionary
                bar_dict_F_V1[barF_ID] = barcode_F
                bar_dict_R_V1[barR_ID] = barcode_R
                sample_dict_V1[sample] = [barF_ID,barR_ID]

            else:
                # sample not part of this sequencing file, skip
                continue

    for key,value in sample_dict_V1.items():
	sample_dict["".join(value)] = key

    # DEBUGGING
#    print("Sample dict:", sample_dict)
#    print("Sample dict 3of4:", sample_dict_U)

    # Barcode/index dictionaries:
    # input: {ID : sequence}
    # generate all possible 1 substitutions per sequence
    # output: {sequence : ID}

    for key, value in bar_dict_F_V1.items():
        for i in range(len(value)):
            barcode = list(value)
            for j in nucleotides:
                barcode[i] = j
                bar_dict_F["".join(barcode)] = key

    for key, value in bar_dict_R_V1.items():
        for i in range(len(value)):
            barcode = list(value)
            for j in nucleotides:
                barcode[i] = j
                bar_dict_R["".join(barcode)] = key

#    print("BC5 dict:", bar_dict_F)
#    print("BC7 dict:", bar_dict_R)

    sample_reads_F = {}
    sample_reads_R = {}
    sample_reads_F["unidentified"] = []
    sample_reads_R["unidentified"] = []
    
    count = 0
    readcount = 0

    with gzip.open(FORWARD_FILE, 'r') as forward_file, gzip.open(REVERSE_FILE,'r') as reverse_file:
        while True:
            fastq_read_F = list(islice(forward_file, 4))
            fastq_read_R = list(islice(reverse_file, 4))
	    readcount += 1
	    if readcount % 1000 == 0:
		print("Processed %i reads, %i unassigned (%0.3f percent)" % (readcount, count, float(count)/readcount*100)) 
            if not fastq_read_F:
                print("File ended at read %i" % (readcount-1))
                break

            barcode_F_number = bar_dict_F.get(fastq_read_F[1].strip()[0:7],"unknown")
            barcode_R_number = bar_dict_R.get(fastq_read_R[1].strip()[0:7],"unknown")
            merged_numbers = barcode_F_number + barcode_R_number
#            print("Current barcode:", merged_numbers)

            if not fastq_read_F[0].strip().split(" ")[-1][-15:] == fastq_read_R[0].strip().split(" ")[-1][-15:]: #check at same point in both read files
                print("Forward and reverse read index headers do not match for:/n", fastq_read_F)
                print("File ended at read %i" % (counter-1))
                break

            if merged_numbers in sample_dict:
		sample = sample_dict[merged_numbers]
            else:
                sample = "unidentified"
                count += 1
                outputfile_failed = open("log_unidentified_barcode_index.txt","a")
                outputfile_failed.write(merged_numbers+"\n")
                outputfile_failed.close()

            # no hard trimming of barcodes
            outputfile_F = open(sample + "_F.fastq", "a")
            outputfile_F.write(fastq_read_F[0] + fastq_read_F[1] + fastq_read_F[2] + fastq_read_F[3])
            outputfile_F.close()

            outputfile_R = open(sample + "_R.fastq", "a")
            outputfile_R.write(fastq_read_R[0] + fastq_read_R[1] + fastq_read_R[2] + fastq_read_R[3])
            outputfile_R.close()

    print("No. of reads unassigned: %i" % (count))

if __name__ == "__main__":
    forward_file = argv[1]
    reverse_file = argv[2]
    index_file = argv[3]
    indexID = argv[4]
    demultiplex_by_index(forward_file,reverse_file,index_file,indexID)
    print("DONE")
