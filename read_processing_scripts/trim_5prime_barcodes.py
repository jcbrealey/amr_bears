#!/usr/bin/python
from sys import argv
from itertools import islice #izip
import gzip

# python 3 assumes all strings are bytes, need to convert to str with proper encoding
# see https://stackoverflow.com/questions/21689365/python-3-typeerror-must-be-str-not-bytes-with-sys-stdout-write/21689447

def trim_fwd_barcodes(FASTQ_FILE,basename,outname):
    outputfile = open(str(basename)+str(outname), "w")
    #fastq_dict = {}

    with gzip.open(FASTQ_FILE, 'r') as readfile:
        while True:
            fastq_read_F = list(islice(readfile, 4))
            if not fastq_read_F:
                break
#            if fastq_read_F[1].strip() in fastq_dict:
#                continue
#            else:
#                fastq_dict[fastq_read_F[1].strip()] = ""
#            outputfile.write(fastq_read_F[0] + fastq_read_F[1][7:] + fastq_read_F[2] + fastq_read_F[3][7:])

            outputfile.write(fastq_read_F[0].decode('utf-8') + fastq_read_F[1][7:].decode('utf-8') + fastq_read_F[2].decode('utf-8') + fastq_read_F[3][7:].decode('utf-8'))

        outputfile.close()

if __name__ == "__main__":
    FASTQ_FILE = argv[1]
    basename = argv[2]
    outname = argv[3]
    trim_fwd_barcodes(FASTQ_FILE,basename,outname)
