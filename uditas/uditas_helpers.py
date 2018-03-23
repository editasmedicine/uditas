#!/usr/bin/env python
"""
    Function definitions to process UDiTaS data
"""

# imports here

from __future__ import print_function

from _version import __version__

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

import pylab

import pandas as pd

import numpy as np

import os

import sys

import gzip

import itertools

import operator

import subprocess

import twobitreader

from Bio.Alphabet import IUPAC

from Bio import SeqIO

from Bio.Seq import Seq

from Bio.SeqRecord import SeqRecord

import pysam



# Metadata
__author__ = "Eugenio Marco"
__credits__ = ["David Kelly"]
__status__ = "Development"


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class StrandError(Error):
    """Exception raised for errors in the strand information.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class ReactionTypeError(Error):
    """Exception raised for errors in the reaction type to be processed.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


################################################################################
# Open .fastq or .fastq.gz files for reading
################################################################################
def open_fastq_or_gz(filename):
    if filename.endswith(".fastq") and os.access(filename, os.F_OK):
        return open(filename, "rU")
    elif filename.endswith(".fastq.gz") and os.access(filename, os.F_OK):
        return gzip.open(filename, "rb")
    elif filename.endswith(".fastq") and os.access(filename + ".gz", os.F_OK):
        return gzip.open(filename + ".gz", "rb")
    elif filename.endswith(".fastq.gz") and os.access(filename[:-3], os.F_OK):
        return open(filename[:-3], "rU")
    raise IOError("Unknown file: " + filename)


################################################################################
# Hamming distance
# From http://code.activestate.com/recipes/499304-hamming-distance/
################################################################################
def hamm_dist(str1, str2):
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(itertools.imap(ne, str1, str2))


################################################################################
# Select closest barcode with a maximum number of mismatches
# By default it returns barcodes with a maximum of n_max_mismatches=2 mismatches
################################################################################
def select_barcode(seq, barcode_list, n_max_mismatches=1):
    # This compares with all barcodes and selects the one with the smallest hamming distance
    # Before calling this function check if the sequence is already a barcode
    matched_barcodes = list()
    distances = list()
    for barcode in barcode_list:
        h_d = hamm_dist(seq, barcode)
        if h_d <= n_max_mismatches:
            matched_barcodes.append(barcode)
            distances.append(h_d)
    indices = [i for i, x in enumerate(distances) if x == min(distances)]
    return [matched_barcodes[i] for i in indices]


################################################################################
# Mask sequence by quality score
################################################################################
def mask(seq, qual, min_qual=12):

    return "".join((b if (ord(q) - 33) >= min_qual else "N") for b, q in itertools.izip(seq, qual))


################################################################################
# get the reverse-complement DNA sequence
################################################################################
def reverse_complement(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return "".join([seq_dict[base] for base in reversed(seq)])


################################################################################
# Create umi dict
################################################################################
def create_umi_dict(filename):

    umi_file = open_fastq_or_gz(filename)

    umi_dict = dict()

    umi_reads = itertools.izip(umi_file)

    for header_umi in umi_reads:

        seq_umi = umi_reads.next()
        umi_reads.next()
        qual_umi = umi_reads.next()
        umi_dict[header_umi[0].split()[0][1:]] = [seq_umi[0].rstrip(), qual_umi[0].rstrip()]

    return umi_dict


################################################################################
# create list of output files
################################################################################
def create_filename(dir_sample, N7, N5, filetype):
    main_folder = os.path.join(dir_sample, N7 + '_' + N5)
    if filetype == 'mainfolder':
        return main_folder
    elif filetype == 'amplicons':
        return os.path.join(main_folder, 'amplicons')
    elif filetype == 'R1fastq':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.fastq')
    elif filetype == 'R1fastqgz':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.fastq.gz')
    elif filetype == 'R2fastq':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.fastq')
    elif filetype == 'R2fastqgz':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.fastq.gz')
    elif filetype == 'umifastq':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_umi.fastq')
    elif filetype == 'umifastqgz':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_umi.fastq.gz')
    elif filetype == 'R1trimmed':
        return os.path.join(main_folder, 'cutadapt_files', N7 + '_' + N5 + '_R1.trimmed.fastq.gz')
    elif filetype == 'R2trimmed':
        return os.path.join(main_folder, 'cutadapt_files', N7 + '_' + N5 + '_R2.trimmed.fastq.gz')
    elif filetype == 'trimmed_report':
        return os.path.join(main_folder, 'cutadapt_files', N7 + '_' + N5 + '.trimmed.report.txt')
    elif filetype == 'sam_genome_local':
        return os.path.join(main_folder, 'sam_genome_local_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_genome_local':
        return os.path.join(main_folder, 'sam_genome_local_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_genome_local':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_genome_local':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_genome_local':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'sam_plasmid_local':
        return os.path.join(main_folder, 'sam_plasmid_local_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_plasmid_local':
        return os.path.join(main_folder, 'sam_plasmid_local_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'unmapped_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '_unmapped.bam')
    elif filetype == 'qsorted_unmapped_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '_qsorted_unmapped.bam')
    elif filetype == 'unmapped_plasmid_R1fastq':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R1.fastq')
    elif filetype == 'unmapped_plasmid_R2fastq':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R2.fastq')
    elif filetype == 'unmapped_plasmid_R1fastqgz':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R1.fastq.gz')
    elif filetype == 'unmapped_plasmid_R2fastqgz':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R2.fastq.gz')
    elif filetype == 'sam_amplicons':
        return os.path.join(main_folder, 'sam_amplicon_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_amplicons':
        return os.path.join(main_folder, 'sam_amplicon_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'unmapped_bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '_amplicons_unmapped.bam')
    elif filetype == 'qsorted_unmapped_bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '_qsorted_amplicons_unmapped.bam')
    elif filetype == 'unmapped_amplicons_R1fastq':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files', N7 + '_' + N5 + '_amplicons_unmapped_R1.fastq')
    elif filetype == 'unmapped_amplicons_R2fastq':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files', N7 + '_' + N5 + '_amplicons_unmapped_R2.fastq')
    elif filetype == 'unmapped_amplicons_R1fastqgz':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files',
                            N7 + '_' + N5 + '_amplicons_unmapped_R1.fastq.gz')
    elif filetype == 'unmapped_amplicons_R2fastqgz':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files',
                            N7 + '_' + N5 + '_amplicons_unmapped_R2.fastq.gz')
    elif filetype == 'unmapped_amplicons_report':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files', N7 + '_' + N5 + '.unmapped.report.txt')
    elif filetype == 'sam_genome_global':
        return os.path.join(main_folder, 'sam_genome_global_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_genome_global':
        return os.path.join(main_folder, 'sam_genome_global_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_genome_global':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_genome_global':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_genome_global':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'results_amplicons':
        return os.path.join(main_folder, 'results', N7 + '_' + N5)  # We will append the window size later
    elif filetype == 'results_plasmid':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_plasmid.xlsx')
    elif filetype == 'results_all_amplicons':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_all_amplicons.xlsx')
    elif filetype == 'results_genomewide':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_genomewide.xlsx')
    elif filetype == 'summary_all_alignments':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_summary_all_alignments.xlsx')
    elif filetype == 'read_counts':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_read_counts.xlsx')


############################
#
# Demultiplexer
# Input: folder to demultiplex, with Undetermined fastq files and sample info in sample_info.csv
#
# ##########################
def demultiplex(dir_sample):

    # Read indices
    sample_info_filename = os.path.join(dir_sample, 'sample_info.csv')

    experiments = pd.read_csv(sample_info_filename)

    r1_fastq = os.path.join(dir_sample, 'Undetermined_S0_L001_R1_001.fastq.gz')
    r2_fastq = os.path.join(dir_sample, 'Undetermined_S0_L001_R2_001.fastq.gz')
    i1_fastq = os.path.join(dir_sample, 'Undetermined_S0_L001_I1_001.fastq.gz')
    i2_fastq = os.path.join(dir_sample, 'Undetermined_S0_L001_I2_001.fastq.gz')

    index_i1_list = list(experiments['index_I1'])
    barcode_i1_list = list(experiments['barcode_I1'])
    i1_dict = dict(zip(barcode_i1_list, index_i1_list))
    index_i2_list = list(experiments['index_I2'])
    barcode_i2_list = list(experiments['barcode_I2'])
    i2_dict = dict(zip(barcode_i2_list, index_i2_list))

    index_i1_set = set(index_i1_list)

    good_barcode_pairs = dict()

    for bc in index_i1_set:
        good_barcode_pairs[bc] = list(experiments.loc[bc == experiments['index_I1']]['index_I2'])

    barcode_i2_length = len(barcode_i2_list[0])

    files_out = list()

    # Create all directories if necessary
    N7_N5 = itertools.izip(index_i1_list, index_i2_list)
    for (N7, N5) in N7_N5:
        exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')
        if not os.path.exists(exp_dir):
            os.mkdir(exp_dir)
        if not os.path.exists(os.path.dirname(create_filename(dir_sample, N7, N5, 'R1fastq'))):
            os.mkdir(os.path.dirname(create_filename(dir_sample, N7, N5, 'R1fastq')))
        files_out.append(create_filename(dir_sample, N7, N5, 'R1fastq'))
        files_out.append(create_filename(dir_sample, N7, N5, 'R2fastq'))
        files_out.append(create_filename(dir_sample, N7, N5, 'umifastq'))

    # create dict of output basename files, to map to opened files
    n_file = 0
    files_out_dict = dict()

    for file_selected in files_out:
        files_out_dict[os.path.basename(file_selected)] = n_file
        n_file += 1

    if not os.path.exists(os.path.join(dir_sample, 'mismatched')):
        os.mkdir(os.path.join(dir_sample, 'mismatched'))

    if not os.path.exists(os.path.join(dir_sample, 'reports')):
        os.mkdir(os.path.join(dir_sample, 'reports'))

    file_out_not_in_exp_list_r1 = os.path.join(dir_sample, 'mismatched', 'not_in_exp_list_R1.fastq')
    file_out_not_in_exp_list_r2 = os.path.join(dir_sample, 'mismatched', 'not_in_exp_list_R2.fastq')
    file_out_not_in_exp_list_i1 = os.path.join(dir_sample, 'mismatched', 'not_in_exp_list_I1.fastq')
    file_out_not_in_exp_list_i2 = os.path.join(dir_sample, 'mismatched', 'not_in_exp_list_I2.fastq')

    file_out_mismatched_adapters_r1 = os.path.join(dir_sample, 'mismatched', 'mismatched_adapters_R1.fastq')
    file_out_mismatched_adapters_r2 = os.path.join(dir_sample, 'mismatched', 'mismatched_adapters_R2.fastq')
    file_out_mismatched_adapters_i1 = os.path.join(dir_sample, 'mismatched', 'mismatched_adapters_I1.fastq')
    file_out_mismatched_adapters_i2 = os.path.join(dir_sample, 'mismatched', 'mismatched_adapters_I2.fastq')

    # We open all output files
    ref_files = [open(filename, "w") for filename in files_out]

    ref_file_out_not_in_exp_list_r1 = open(file_out_not_in_exp_list_r1, "w")
    ref_file_out_not_in_exp_list_r2 = open(file_out_not_in_exp_list_r2, "w")
    ref_file_out_not_in_exp_list_i1 = open(file_out_not_in_exp_list_i1, "w")
    ref_file_out_not_in_exp_list_i2 = open(file_out_not_in_exp_list_i2, "w")

    ref_file_out_mismatched_adapters_r1 = open(file_out_mismatched_adapters_r1, "w")
    ref_file_out_mismatched_adapters_r2 = open(file_out_mismatched_adapters_r2, "w")
    ref_file_out_mismatched_adapters_i1 = open(file_out_mismatched_adapters_i1, "w")
    ref_file_out_mismatched_adapters_i2 = open(file_out_mismatched_adapters_i2, "w")

    file_read_counts = [0] * len(files_out)

    # We open r1,r2,i1,i2 files and distribute reads
    with open_fastq_or_gz(r1_fastq) as r1_file, open_fastq_or_gz(r2_fastq) as r2_file, open_fastq_or_gz(i1_fastq) as i1_file, open_fastq_or_gz(i2_fastq) as i2_file:
        # Add counters for all reads

        reads_in_experiment_list_count = 0

        reads_not_in_experiment_list_count = 0

        mismatch_count = 0
        mismatch_count_i1 = 0
        mismatch_count_i2 = 0

        mismatch_dict_i1 = dict()
        mismatch_dict_i2 = dict()

        r1_r2_i1_i2 = itertools.izip(r1_file, r2_file, i1_file, i2_file)

        for header_r1, header_r2, header_i1, header_i2 in r1_r2_i1_i2:
    
            seq_r1, seq_r2, seq_i1, seq_i2_plus_umi = r1_r2_i1_i2.next()

            r1_r2_i1_i2.next()

            qual_r1, qual_r2, qual_i1, qual_i2_plus_umi = r1_r2_i1_i2.next()

            seq_i1, seq_i2_plus_umi = seq_i1.rstrip(), seq_i2_plus_umi.rstrip()

            qual_i1, qual_i2_plus_umi = qual_i1.rstrip(), qual_i2_plus_umi.rstrip()

            #       We mask with N any bases with scores below or equal to , (11, default in mask)
            seq_i1 = mask(seq_i1, qual_i1)

            seq_i2 = mask(seq_i2_plus_umi[:barcode_i2_length], qual_i2_plus_umi[:barcode_i2_length])

            umi_qual = qual_i2_plus_umi[barcode_i2_length:]

            umi = mask(seq_i2_plus_umi[barcode_i2_length:], umi_qual)

            # change to 1 for reads with perfect indices or match after correction
            is_good_index = 0

            if (seq_i1 in barcode_i1_list) and (seq_i2 in barcode_i2_list):
                # perfect match case
                is_good_index = 1
            else:
                # We look for barcodes with up to two mismatches, default in select_barcode
                seq_i1_match = select_barcode(seq_i1, barcode_i1_list)
                seq_i2_match = select_barcode(seq_i2, barcode_i2_list)
                if len(seq_i2_match) > 0 and len(seq_i1_match) > 0:
                    # match after selecting adapter with up to 2 mismatches (default in select_barcode)
                    is_good_index = 1
                    seq_i1 = seq_i1_match[0]
                    seq_i2 = seq_i2_match[0]

            if is_good_index:
                # We test whether the read has on of the combination of indices from our experiment list
                # If not save in a separate file
                if i2_dict[seq_i2] in good_barcode_pairs[i1_dict[seq_i1]]:

                    r1f = create_filename(dir_sample, i1_dict[seq_i1], i2_dict[seq_i2], 'R1fastq')
                    r2f = create_filename(dir_sample, i1_dict[seq_i1], i2_dict[seq_i2], 'R2fastq')
                    umif = create_filename(dir_sample, i1_dict[seq_i1], i2_dict[seq_i2], 'umifastq')

                    print("\n".join([header_r1.rstrip(), seq_r1.rstrip(), "+", qual_r1.rstrip()]),
                          file=ref_files[files_out_dict[os.path.basename(r1f)]])
                    file_read_counts[files_out_dict[os.path.basename(r1f)]] += 1

                    print("\n".join([header_r2.rstrip(), seq_r2.rstrip(), "+", qual_r2.rstrip()]),
                          file=ref_files[files_out_dict[os.path.basename(r2f)]])
                    file_read_counts[files_out_dict[os.path.basename(r2f)]] += 1

                    print("\n".join([header_i2.rstrip(), umi, "+", umi_qual]),
                          file=ref_files[files_out_dict[os.path.basename(umif)]])
                    file_read_counts[files_out_dict[os.path.basename(umif)]] += 1

                    reads_in_experiment_list_count += 1

                else:
                    # We print reads with mismatched labels to our experiments
                    print("\n".join([header_r1.rstrip(), seq_r1.rstrip(), "+", qual_r1.rstrip()]),
                          file=ref_file_out_not_in_exp_list_r1)
                    print("\n".join([header_r2.rstrip(), seq_r2.rstrip(), "+", qual_r2.rstrip()]),
                          file=ref_file_out_not_in_exp_list_r2)
                    print("\n".join([header_i1.rstrip(), seq_i1.rstrip(), "+", qual_i1.rstrip()]),
                          file=ref_file_out_not_in_exp_list_i1)
                    print("\n".join([header_i2.rstrip(), seq_i2_plus_umi.rstrip(), "+",
                                     qual_i2_plus_umi.rstrip()]), file=ref_file_out_not_in_exp_list_i2)

                    reads_not_in_experiment_list_count += 1
            else:
                # We print reads with mismatched adapters
                print("\n".join([header_r1.rstrip(), seq_r1.rstrip(), "+", qual_r1.rstrip()]),
                      file=ref_file_out_mismatched_adapters_r1)
                print("\n".join([header_r2.rstrip(), seq_r2.rstrip(), "+", qual_r2.rstrip()]),
                      file=ref_file_out_mismatched_adapters_r2)
                print("\n".join([header_i1.rstrip(), seq_i1.rstrip(), "+", qual_i1.rstrip()]),
                      file=ref_file_out_mismatched_adapters_i1)
                print("\n".join([header_i2.rstrip(), seq_i2_plus_umi.rstrip(), "+",
                                 qual_i2_plus_umi.rstrip()]), file=ref_file_out_mismatched_adapters_i2)

                if seq_i2 not in i2_dict.keys():
                    if seq_i2 in mismatch_dict_i2.keys():
                        mismatch_dict_i2[seq_i2] += 1
                    else:
                        mismatch_dict_i2[seq_i2] = 1
                    mismatch_count_i2 += 1

                if seq_i1 not in i1_dict.keys():
                    if seq_i1 in mismatch_dict_i1.keys():
                        mismatch_dict_i1[seq_i1] += 1
                    else:
                        mismatch_dict_i1[seq_i1] = 1
                    mismatch_count_i1 += 1

                mismatch_count += 1

    # close all files
    for rf in ref_files:
        rf.close()

    ref_file_out_not_in_exp_list_r1.close()
    ref_file_out_not_in_exp_list_r2.close()
    ref_file_out_not_in_exp_list_i1.close()
    ref_file_out_not_in_exp_list_i2.close()

    ref_file_out_mismatched_adapters_r1.close()
    ref_file_out_mismatched_adapters_r2.close()
    ref_file_out_mismatched_adapters_i1.close()
    ref_file_out_mismatched_adapters_i2.close()

    # print report of counts for individual files
    report_file = os.path.join(dir_sample, 'reports', 'report_individual_files.xls')

    x = np.array(file_read_counts)
    fh = open(report_file, "w")
    print("\t".join(['filename', 'reads_count']), file=fh)

    for i in np.nonzero(x)[0]:
        print("\t".join([os.path.basename(files_out[i]), str(file_read_counts[i])]), file=fh)
    fh.close()

    # print report overall counts
    report_file = os.path.join(dir_sample, 'reports', 'report_overall.xls')

    fh = open(report_file, "w")

    print('Total number of reads:\t' + str(reads_in_experiment_list_count + reads_not_in_experiment_list_count +
                                           mismatch_count) + '\n', file=fh)
    print('Reads without I1 or I2 adapters:\t' + str(mismatch_count) + '\n', file=fh)
    print('Reads without I1 adapters:\t' + str(mismatch_count_i1) + '\n', file=fh)
    print('Reads without I2 adapters:\t' + str(mismatch_count_i2) + '\n', file=fh)
    print('Reads with I1 or I2 adapters in the wrong combination:\t' + str(reads_not_in_experiment_list_count) + '\n',
          file=fh)
    print('Reads with I1 or I2 adapters matching our experiments:\t' + str(reads_in_experiment_list_count) + '\n',
          file=fh)
    fh.close()

    # print report mismatched i1_rc adapters
    report_file = os.path.join(dir_sample, 'reports', 'report_mismatched_adapters_i1.xls')

    fh = open(report_file, "w")

    print('I1_RC\tCount', file=fh)

    for dict_element in mismatch_dict_i1:
        print(dict_element + '\t' + str(mismatch_dict_i1[dict_element]), file=fh)

    fh.close()

    # print report mismatched i2 adapters
    report_file = os.path.join(dir_sample, 'reports', 'report_mismatched_adapters_i2.xls')

    fh = open(report_file, "w")

    print('I2\tCount', file=fh)

    for dict_element in mismatch_dict_i2:
        print(dict_element + '\t' + str(mismatch_dict_i2[dict_element]), file=fh)

    fh.close()

    # gzip fastq files
    for fo in files_out:
        with open(fo) as f_in, gzip.open(fo + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)


################################################
# Function to write reference plasmid sequence
################################################
def create_plasmid_reference(dir_sample, amplicon_info):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')
    amplicon_folder = os.path.join(exp_dir, 'amplicons')
    if not os.path.exists(amplicon_folder):
        os.mkdir(amplicon_folder)

    filename = os.path.join(exp_dir, amplicon_folder, 'plasmid.fa')
    file_handle = open(filename, "w")

    # If several plasmids were given, separate them using ';' in sample_info.csv. Here we remove the ';' so that
    # we just concatenate the sequences
    pl_seq = amplicon_info['plasmid_sequence'].replace(';', '')
    seq1 = Seq(pl_seq, IUPAC.unambiguous_dna)
    record1 = SeqRecord(seq1, 'plasmid', description='')
    SeqIO.write(record1, file_handle, 'fasta')

    file_handle.close()
    # Create index file
    initial_dir = os.getcwd()
    os.chdir(amplicon_folder)
    index_err_file = os.path.join(amplicon_folder, 'index_plasmid.err')
    index_out_file = os.path.join(amplicon_folder, 'index_plasmid.out')

    index_err_fh = open(index_err_file, 'wb')
    index_out_fh = open(index_out_file, 'wb')
    subprocess.call(['bowtie2-build',
                     filename, 'plasmid'], stderr=index_err_fh, stdout=index_out_fh)
    os.chdir(initial_dir)
    index_err_fh.close()
    index_out_fh.close()


#################################################################################
# Function to write reference amplicons with various structural rearrangements
#################################################################################
def write_amplicon(dir_sample, amplicon_info, amplicon_list):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')
    amplicon_folder = os.path.join(exp_dir, 'amplicons')
    if not os.path.exists(amplicon_folder):
        os.mkdir(amplicon_folder)

    filename = os.path.join(exp_dir, amplicon_folder, 'amplicons.fa')
    file_handle = open(filename, "w")

    for amps in amplicon_list:
        seq1 = Seq(amps[1], IUPAC.unambiguous_dna)
        record1 = SeqRecord(seq1, amps[0], description='')
        SeqIO.write(record1, file_handle, 'fasta')

    file_handle.close()
    # Create index file
    initial_dir = os.getcwd()
    os.chdir(amplicon_folder)
    index_err_file = os.path.join(amplicon_folder, 'index.err')
    index_out_file = os.path.join(amplicon_folder, 'index.out')

    index_err_fh = open(index_err_file, 'wb')
    index_out_fh = open(index_out_file, 'wb')
    subprocess.call(['bowtie2-build',
                     filename, 'amplicons'], stderr=index_err_fh, stdout=index_out_fh)
    os.chdir(initial_dir)
    index_err_fh.close()
    index_out_fh.close()


############################
#
#  Function to determine the kind of reaction from the number of cuts and their locations
#
# The function classify the cases:
#   - No cut (just UDiTaS primer, used for controls)
#   - Single cut
#   - Dual cut on same chromosome, generates amplicons with large deletions, etc
#   - Dual cut on different chromosomes. Generates 10 amplicons including translocations
#   - Triple cuts on different chromosomes. Generates 21 amplicons including translocations. NOTE that if two of the
#     cuts are in the same chromosome and close together (less than amplicon_window_around_cut) the results may
#     be incorrect since some reads may be mapped to multiple amplicons, but only counted around the cut in one amplicon
#
############################
def get_reaction_type(amplicon_info):
    has_guide1 = type(amplicon_info['chr_guide_1']) is str or type(amplicon_info['chr_guide_1']) is unicode
    has_guide2 = type(amplicon_info['chr_guide_2']) is str or type(amplicon_info['chr_guide_2']) is unicode
    has_guide3 = type(amplicon_info['chr_guide_3']) is str or type(amplicon_info['chr_guide_3']) is unicode

    if not has_guide1 and not has_guide2 and not has_guide3:
        reaction_type = 'control'
    elif has_guide1 and not has_guide2 and not has_guide3:
        reaction_type = 'single_cut'
    elif has_guide1 and has_guide2 and amplicon_info['chr_guide_1'] == amplicon_info['chr_guide_2'] and not has_guide3:
        reaction_type = 'double_cut_same_chromosome'
    elif has_guide1 and has_guide2 and amplicon_info['chr_guide_1'] != amplicon_info['chr_guide_2'] and not has_guide3:
        reaction_type = 'double_cut_different_chromosomes'
    elif has_guide1 and has_guide2 and has_guide3:
        if (amplicon_info['chr_guide_1'] == amplicon_info['chr_guide_2'] or
                amplicon_info['chr_guide_1'] == amplicon_info['chr_guide_3'] or
                amplicon_info['chr_guide_2'] == amplicon_info['chr_guide_3']):
            raise ReactionTypeError('The reaction with three cuts with at least two in the same chromosome is' +
                                    ' not yet supported by current version of UDiTaS')
        reaction_type = 'triple_cut_different_chromosomes'
    else:
        raise ReactionTypeError('Reaction type not yet supported by current version of UDiTaS')

    return reaction_type


############################
#
# Create amplicon. Creates fasta file with the custom reference amplicons including deletions, inversions, etc...
# Input: dir_sample, directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        file_genome_2bit, 2bit file with the reference genome being used
#        amplicon_window_around_cut, value used to grab sequences around cut sites
#
#   This functions creates a set of reference amplicons built from the expected fragments after cutting
#   The amplicons are built differently depending on the case classified by get_reaction_type()
#
# ##########################
def create_amplicon(dir_sample, amplicon_info, file_genome_2bit, amplicon_window_around_cut=1000):
    # We first check the reaction type
    reaction_type = get_reaction_type(amplicon_info)

    genome = twobitreader.TwoBitFile(file_genome_2bit)  # Load genome. Used for getting the sequences

    amplicon_list = []

    # For all reaction types, we check that we don't go out of boundaries when building the amplicons
    # This is unlikely for hg38 or mm10, but could easily happen in the UDiTaS primer is in a plasmid
    if reaction_type == 'control':
        # Case no guides
        if amplicon_info['strand'] == '+':  # This is the UDiTaS oligo strand
            end_coordinate = int(amplicon_info['start']) + amplicon_window_around_cut
            if end_coordinate > len(genome[amplicon_info['chr']]):
                end_coordinate = len(genome[amplicon_info['chr']])
            amplicon_list.append(['wt', genome[amplicon_info['chr']][int(amplicon_info['start']):end_coordinate]])
        elif amplicon_info['strand'] == '-':
            start_coordinate = int(amplicon_info['end']) - amplicon_window_around_cut
            if start_coordinate < 0:
                start_coordinate = 0
            amplicon_list.append(['wt', genome[amplicon_info['chr']][start_coordinate:int(amplicon_info['end'])]])
        else:
            raise StrandError('strand can only have as values + or -')
    elif reaction_type == 'single_cut':
        # Case one guide
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        start_coordinate = int(cut1 - amplicon_window_around_cut)
        if start_coordinate < 0:
            start_coordinate = 0
        end_coordinate = int(cut1 + amplicon_window_around_cut)
        if end_coordinate > len(genome[amplicon_info['chr_guide_1']]):
            end_coordinate = len(genome[amplicon_info['chr_guide_1']])

        seq_upstream = genome[amplicon_info['chr_guide_1']][start_coordinate:int(cut1)]
        seq_downstream = genome[amplicon_info['chr_guide_1']][int(cut1):end_coordinate]

        amplicon_list.append(['wt', seq_upstream + seq_downstream])
        amplicon_list.append(['1a_1a', seq_upstream + reverse_complement(seq_upstream)])
        amplicon_list.append(['1b_1b', reverse_complement(seq_downstream) + seq_downstream])
    elif reaction_type == 'double_cut_same_chromosome':
        # Case two guides on the same chromosome
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if amplicon_info['strand_guide_2'] == '+':
            cut2 = amplicon_info['end_guide_2'] - 3
        elif amplicon_info['strand_guide_2'] == '-':
            cut2 = amplicon_info['start_guide_2'] + 3
        else:
            raise StrandError('strand_guide_2 can only have as values + or -')

        # We switch the coordinates of cut1 and cut2 if the guides are provided so that cut2 < cut1
        # cut1 it will always be smaller than cut2 and in the results cut1 will be the cut site with
        # smaller genomic coordinate
        # cut1 and cut2 also flipped in get_cut_in_reference_amplicon_df

        if cut2 < cut1:
            (cut1, cut2) = (cut2, cut1)

        start_coordinate = int(cut1 - amplicon_window_around_cut)
        if start_coordinate < 0:
            start_coordinate = 0
        end_coordinate = int(cut2 + amplicon_window_around_cut)
        if end_coordinate > len(genome[amplicon_info['chr_guide_1']]):
            end_coordinate = len(genome[amplicon_info['chr_guide_1']])

        seq_upstream = genome[amplicon_info['chr_guide_1']][start_coordinate:int(cut1)]
        seq_cut1_cut2 = genome[amplicon_info['chr_guide_1']][int(cut1):int(cut2)]
        seq_downstream = genome[amplicon_info['chr_guide_1']][int(cut2):end_coordinate]

        amplicon_list.append(['wt', seq_upstream + seq_cut1_cut2 + seq_downstream])
        amplicon_list.append(['large_deletion', seq_upstream + seq_downstream])
        amplicon_list.append(['large_inversion', seq_upstream + reverse_complement(seq_cut1_cut2) + seq_downstream])
        amplicon_list.append(['1a_1a', seq_upstream + reverse_complement(seq_upstream)])
        amplicon_list.append(['2b_2b', reverse_complement(seq_downstream) + seq_downstream])
    elif reaction_type == 'double_cut_different_chromosomes':
        # Case two guides on different chromosomes
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if amplicon_info['strand_guide_2'] == '+':
            cut2 = amplicon_info['end_guide_2'] - 3
        elif amplicon_info['strand_guide_2'] == '-':
            cut2 = amplicon_info['start_guide_2'] + 3
        else:
            raise StrandError('strand_guide_2 can only have as values + or -')

        start_coordinate1 = int(cut1 - amplicon_window_around_cut)
        if start_coordinate1 < 0:
            start_coordinate1 = 0
        end_coordinate1 = int(cut1 + amplicon_window_around_cut)
        if end_coordinate1 > len(genome[amplicon_info['chr_guide_1']]):
            end_coordinate1 = len(genome[amplicon_info['chr_guide_1']])

        start_coordinate2 = int(cut2 - amplicon_window_around_cut)
        if start_coordinate2 < 0:
            start_coordinate2 = 0
        end_coordinate2 = int(cut2 + amplicon_window_around_cut)
        if end_coordinate2 > len(genome[amplicon_info['chr_guide_2']]):
            end_coordinate2 = len(genome[amplicon_info['chr_guide_2']])

        seq_1a = genome[amplicon_info['chr_guide_1']][start_coordinate1:int(cut1)]
        seq_1b = genome[amplicon_info['chr_guide_1']][int(cut1):end_coordinate1]
        seq_2a = genome[amplicon_info['chr_guide_2']][start_coordinate2:int(cut2)]
        seq_2b = genome[amplicon_info['chr_guide_2']][int(cut2):end_coordinate2]

        amplicon_list.append(['1a_1a', seq_1a + reverse_complement(seq_1a)])
        amplicon_list.append(['1a_1b', seq_1a + seq_1b])
        amplicon_list.append(['1a_2a', seq_1a + reverse_complement(seq_2a)])
        amplicon_list.append(['1a_2b', seq_1a + seq_2b])

        amplicon_list.append(['1b_1b', reverse_complement(seq_1b) + seq_1b])
        amplicon_list.append(['2a_1b', seq_2a + seq_1b])
        amplicon_list.append(['2b_1b', reverse_complement(seq_2b) + seq_1b])

        amplicon_list.append(['2a_2a', seq_2a + reverse_complement(seq_2a)])
        amplicon_list.append(['2a_2b', seq_2a + seq_2b])

        amplicon_list.append(['2b_2b', reverse_complement(seq_2b) + seq_2b])
    elif reaction_type == 'triple_cut_different_chromosomes':
        # Case three guides on different chromosomes

        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if amplicon_info['strand_guide_2'] == '+':
            cut2 = amplicon_info['end_guide_2'] - 3
        elif amplicon_info['strand_guide_2'] == '-':
            cut2 = amplicon_info['start_guide_2'] + 3
        else:
            raise StrandError('strand_guide_2 can only have as values + or -')

        if amplicon_info['strand_guide_3'] == '+':
            cut3 = amplicon_info['end_guide_3'] - 3
        elif amplicon_info['strand_guide_3'] == '-':
            cut3 = amplicon_info['start_guide_3'] + 3
        else:
            raise StrandError('strand_guide_3 can only have as values + or -')

        start_coordinate1 = int(cut1 - amplicon_window_around_cut)
        if start_coordinate1 < 0:
            start_coordinate1 = 0
        end_coordinate1 = int(cut1 + amplicon_window_around_cut)
        if end_coordinate1 > len(genome[amplicon_info['chr_guide_1']]):
            end_coordinate1 = len(genome[amplicon_info['chr_guide_1']])

        start_coordinate2 = int(cut2 - amplicon_window_around_cut)
        if start_coordinate2 < 0:
            start_coordinate2 = 0
        end_coordinate2 = int(cut2 + amplicon_window_around_cut)
        if end_coordinate2 > len(genome[amplicon_info['chr_guide_2']]):
            end_coordinate2 = len(genome[amplicon_info['chr_guide_2']])

        start_coordinate3 = int(cut3 - amplicon_window_around_cut)
        if start_coordinate3 < 0:
            start_coordinate3 = 0
        end_coordinate3 = int(cut3 + amplicon_window_around_cut)
        if end_coordinate3 > len(genome[amplicon_info['chr_guide_3']]):
            end_coordinate3 = len(genome[amplicon_info['chr_guide_3']])

        seq_1a = genome[amplicon_info['chr_guide_1']][start_coordinate1:int(cut1)]
        seq_1b = genome[amplicon_info['chr_guide_1']][int(cut1):end_coordinate1]
        seq_2a = genome[amplicon_info['chr_guide_2']][start_coordinate2:int(cut2)]
        seq_2b = genome[amplicon_info['chr_guide_2']][int(cut2):end_coordinate2]
        seq_3a = genome[amplicon_info['chr_guide_3']][start_coordinate3:int(cut3)]
        seq_3b = genome[amplicon_info['chr_guide_3']][int(cut3):end_coordinate3]

        amplicon_list.append(['1a_1a', seq_1a + reverse_complement(seq_1a)])
        amplicon_list.append(['1a_1b', seq_1a + seq_1b])
        amplicon_list.append(['1a_2a', seq_1a + reverse_complement(seq_2a)])
        amplicon_list.append(['1a_2b', seq_1a + seq_2b])
        amplicon_list.append(['1a_3a', seq_1a + reverse_complement(seq_3a)])
        amplicon_list.append(['1a_3b', seq_1a + seq_3b])

        amplicon_list.append(['1b_1b', reverse_complement(seq_1b) + seq_1b])
        amplicon_list.append(['2a_1b', seq_2a + seq_1b])
        amplicon_list.append(['2b_1b', reverse_complement(seq_2b) + seq_1b])
        amplicon_list.append(['3a_1b', seq_3a + seq_1b])
        amplicon_list.append(['3b_1b', reverse_complement(seq_3b) + seq_1b])

        amplicon_list.append(['2a_2a', seq_2a + reverse_complement(seq_2a)])
        amplicon_list.append(['2a_2b', seq_2a + seq_2b])
        amplicon_list.append(['2a_3a', seq_2a + reverse_complement(seq_3a)])
        amplicon_list.append(['2a_3b', seq_2a + seq_3b])

        amplicon_list.append(['2b_2b', reverse_complement(seq_2b) + seq_2b])
        amplicon_list.append(['3a_2b', seq_3a + seq_2b])
        amplicon_list.append(['3b_2b', reverse_complement(seq_3b) + seq_2b])

        amplicon_list.append(['3a_3a', seq_3a + reverse_complement(seq_3a)])
        amplicon_list.append(['3a_3b', seq_3a + seq_3b])

        amplicon_list.append(['3b_3b', reverse_complement(seq_3b) + seq_3b])

    write_amplicon(dir_sample, amplicon_info, amplicon_list)


############################
#
# Remove adapters in fastq files
# Input: directory to be analyzed
#        dir_sample, name of the directory for the whole run, typically with the name of a miseq run
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        process_AMP_seq_run, set to 1 to trim in read2 the same adapter as in GUIDE-Seq
#
# ##########################
def trim_fastq(dir_sample, amplicon_info, process_AMP_seq_run):

    # UDiTaS adapters
    Nv2F = 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
    SBS12 = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'

    if process_AMP_seq_run == 1:
        i2_adapter = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    else:
        i2_adapter = Nv2F

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    file_R1 = create_filename(dir_sample, N7, N5, 'R1fastqgz')
    file_R2 = create_filename(dir_sample, N7, N5, 'R2fastqgz')

    file_cutadapt_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')
    file_cutadapt_R2 = create_filename(dir_sample, N7, N5, 'R2trimmed')
    file_cutadapt_report = create_filename(dir_sample, N7, N5, 'trimmed_report')

    if not os.path.exists(os.path.dirname(file_cutadapt_R1)):
        os.mkdir(os.path.dirname(file_cutadapt_R1))

    # remove adapters with cutadapt
    cutadapt_command = ['cutadapt',
                        '-m', '10',
                        '-e', '0.33',
                        '-a', reverse_complement(SBS12),
                        '-A', reverse_complement(i2_adapter),
                        '-o', file_cutadapt_R1, '-p', file_cutadapt_R2,
                        file_R1, file_R2]

    handle_cutadapt_report = open(file_cutadapt_report, 'wb')
    subprocess.call(cutadapt_command, stdout=handle_cutadapt_report)
    handle_cutadapt_report.close()


############################
#
# Aligns reads to the plasmid sequence using bowtie2 and local alignment.
#
# Input: directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#
# ##########################
def align_plasmid_local(dir_sample, amplicon_info, ncpu=4):

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    # exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')

    file_cutadapt_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')
    file_cutadapt_R2 = create_filename(dir_sample, N7, N5, 'R2trimmed')

    file_sam_plasmid_local = create_filename(dir_sample, N7, N5, 'sam_plasmid_local')
    file_sam_report_plasmid_local = create_filename(dir_sample, N7, N5, 'sam_report_plasmid_local')

    if not os.path.exists(os.path.dirname(file_sam_plasmid_local)):
        os.mkdir(os.path.dirname(file_sam_plasmid_local))

    file_bam_plasmid_local = create_filename(dir_sample, N7, N5, 'bam_plasmid_local')
    file_sorted_bam_plasmid_local = create_filename(dir_sample, N7, N5, 'sorted_bam_plasmid_local')
    # file_sorted_bai_genome_local = create_filename(dir_sample, N7, N5, 'sorted_bai_genome_local')

    if not os.path.exists(os.path.dirname(file_bam_plasmid_local)):
        os.mkdir(os.path.dirname(file_bam_plasmid_local))

    # local alignment to the genome with bowtie2
    initial_dir = os.getcwd()

    folder_amplicons = create_filename(dir_sample, N7, N5, 'amplicons')

    os.chdir(folder_amplicons)

    bowtie2_command = ['bowtie2', '--local', '-p', str(ncpu),
                       '-X', '5000', '-k', '2', '-x', 'plasmid',
                             '-1', file_cutadapt_R1, '-2', file_cutadapt_R2,
                             '-S', file_sam_plasmid_local]

    handle_sam_report_genome_local = open(file_sam_report_plasmid_local, 'wb')

    subprocess.call(bowtie2_command, stderr=handle_sam_report_genome_local)

    handle_sam_report_genome_local.close()

    # convert sam to bam
    sam_to_bam_plasmid_local_command = ['samtools', 'view', '-Sb', file_sam_plasmid_local]

    handle_file_bam_plasmid_local = open(file_bam_plasmid_local, 'wb')

    subprocess.call(sam_to_bam_plasmid_local_command, stdout=handle_file_bam_plasmid_local)

    # sort bam files
    sort_bam_plasmid_local_command = ['samtools', 'sort', file_bam_plasmid_local, '-o', file_sorted_bam_plasmid_local]

    subprocess.call(sort_bam_plasmid_local_command)

    # Create bam index files
    create_bam_plasmid_local_index_command = ['samtools', 'index', file_sorted_bam_plasmid_local]
    subprocess.call(create_bam_plasmid_local_index_command)

    # Clean up
    os.remove(file_sam_plasmid_local)
    os.remove(file_bam_plasmid_local)

    os.chdir(initial_dir)


#################################################################################
# Function to extract reads that did not align to the plasmid
#################################################################################
def extract_unmapped_reads_plasmid(dir_sample, amplicon_info):

    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    file_sorted_bam_plasmid_local = create_filename(dir_sample, N7, N5, 'sorted_bam_plasmid_local')

    file_unmapped_bam_plasmid = create_filename(dir_sample, N7, N5, 'unmapped_bam_plasmid_local')

    file_qsorted_unmapped_bam_plasmid = create_filename(dir_sample, N7, N5, 'qsorted_unmapped_bam_plasmid_local')

    file_R1_unmapped = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R1fastq')
    file_R2_unmapped = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R2fastq')

    if not os.path.exists(os.path.dirname(file_R1_unmapped)):
        os.mkdir(os.path.dirname(file_R1_unmapped))

    extract_unmapped_bam_command = ['samtools', 'view', '-b', '-f', '0x4', file_sorted_bam_plasmid_local, '-o',
                                    file_unmapped_bam_plasmid]

    subprocess.call(extract_unmapped_bam_command)

    qsort_unmapped_bam_command = ['samtools', 'sort', '-n', file_unmapped_bam_plasmid, '-o',
                                  file_qsorted_unmapped_bam_plasmid]

    subprocess.call(qsort_unmapped_bam_command)

    bamtofastq_command = ['bedtools', 'bamtofastq', '-i', file_qsorted_unmapped_bam_plasmid,
                          '-fq', file_R1_unmapped, '-fq2', file_R2_unmapped]

    file_err = file_R1_unmapped[:-9] + '_err.txt'
    handle_file_err = open(file_err, 'wb')

    subprocess.call(bamtofastq_command, stderr=handle_file_err)

    for fo in [file_R1_unmapped, file_R2_unmapped]:
        with open(fo) as f_in, gzip.open(fo + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)
        os.remove(fo)


#################################################################################
# Function to analyze reads aligned to the plasmid
#################################################################################
def analyze_alignments_plasmid(dir_sample, amplicon_info, min_MAPQ, file_genome_2bit, do_plasmid):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')

    file_UMI = create_filename(dir_sample, N7, N5, 'umifastqgz')
    UMI_dict = create_barcode_dict(file_UMI)

    results_folder = os.path.join(exp_dir, 'results')
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    results_file = create_filename(dir_sample, N7, N5, 'results_plasmid')

    if do_plasmid:
        file_sorted_bam_plasmid_local = create_filename(dir_sample, N7, N5, 'sorted_bam_plasmid_local')

        bam_in_alignment_file = pysam.AlignmentFile(file_sorted_bam_plasmid_local, 'rb')

        bam_in = bam_in_alignment_file.fetch()

        genome = twobitreader.TwoBitFile(file_genome_2bit)  # Load genome. Used for getting the sequences

        length_to_test = 15  # We check this number of bases after the primer
        uditas_primer_length = amplicon_info['end'] - amplicon_info['start']

        if amplicon_info['strand'] == '+':  # This is the UDiTaS oligo strand
            seq_after_uditas_primer = genome[amplicon_info['chr']][amplicon_info['end']:(amplicon_info['end'] +
                                                                                        length_to_test)]
        elif amplicon_info['strand'] == '-':
            seq_after_uditas_primer = reverse_complement(genome[amplicon_info['chr']][(amplicon_info['start'] -
                                                                                       length_to_test):amplicon_info['start']])
        n_max_mismatches = 2  # We allow this number of mismatches between the read and the sequence after the primer

        names_list_plasmid_genome = []
        UMI_list_plasmid_genome = []
        names_list_plasmid_only = []
        UMI_list_plasmid_only = []

        for read in bam_in:
            if read.mapping_quality >= min_MAPQ and not read.is_unmapped and not read.is_secondary:
                if read.is_read2:  # R2 is the UDiTaS primer
                    if read.is_reverse:
                        seq_test = reverse_complement(read.query_sequence)[uditas_primer_length:(uditas_primer_length +
                                                                                            length_to_test)]
                    else:
                        seq_test = read.query_sequence[uditas_primer_length:(uditas_primer_length + length_to_test)]

                    # Sometimes, after cutadapt we have a read shorter than uditas_primer_length + length_to_test
                    # We skip those directly without calculating hamm_dist, which doesn't make sense
                    if (len(seq_test) == len(seq_after_uditas_primer.upper()) and
                        hamm_dist(seq_test, seq_after_uditas_primer.upper()) <= n_max_mismatches):
                        # Reads for which the R2 has genomic sequence after the UDiTaS primer
                        UMI_list_plasmid_genome.append(UMI_dict[read.query_name][0])
                        names_list_plasmid_genome.append(read.query_name)
                    else: # We put those short reads into the plasmid only bucket
                        UMI_list_plasmid_only.append(UMI_dict[read.query_name][0])
                        names_list_plasmid_only.append(read.query_name)

        total_reads_plasmid_genome = len(set(names_list_plasmid_genome))
        total_reads_collapsed_plasmid_genome = len(set(UMI_list_plasmid_genome))
        total_reads_plasmid_only = len(set(names_list_plasmid_only))
        total_reads_collapsed_plasmid_only = len(set(UMI_list_plasmid_only))

        results_df = pd.DataFrame({'target_plus_plasmid_total_reads': [total_reads_plasmid_genome],
                                   'target_plus_plasmid_total_reads_collapsed': [total_reads_collapsed_plasmid_genome],
                                   'plasmid_only_total_reads': [total_reads_plasmid_only],
                                   'plasmid_only_total_reads_collapsed': [total_reads_collapsed_plasmid_only]
                                   },
                                  columns=['target_plus_plasmid_total_reads',
                                           'target_plus_plasmid_total_reads_collapsed',
                                           'plasmid_only_total_reads',
                                           'plasmid_only_total_reads_collapsed'])
    else:
        results_df = pd.DataFrame(index=np.arange(1),
                                  columns=['target_plus_plasmid_total_reads',
                                           'target_plus_plasmid_total_reads_collapsed',
                                           'plasmid_only_total_reads',
                                           'plasmid_only_total_reads_collapsed'])

    results_df.to_excel(results_file)

    return results_df


############################
#
# Aligns reads to the whole genome using bowtie2 and local alignment. This is needed to find translocations using
#   split reads and a program like socrates
#
# Input: directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        assembly, name of the assembly to be used by bowtie2. It is convenient to place it in a folder especified by
#        the environmental variable BOWTIE2_INDEXES
#
# ##########################
def align_genome_local(dir_sample, amplicon_info, assembly, check_plasmid_insertions, ncpu=4):

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    has_plasmid = type(amplicon_info['plasmid_sequence']) is str or type(amplicon_info['plasmid_sequence']) is unicode

    if check_plasmid_insertions == 1 and has_plasmid:
        file_R1 = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R1fastqgz')
        file_R2 = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R2fastqgz')
    else:
        file_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')
        file_R2 = create_filename(dir_sample, N7, N5, 'R2trimmed')

    file_sam_genome_local = create_filename(dir_sample, N7, N5, 'sam_genome_local')
    file_sam_report_genome_local = create_filename(dir_sample, N7, N5, 'sam_report_genome_local')

    if not os.path.exists(os.path.dirname(file_sam_genome_local)):
        os.mkdir(os.path.dirname(file_sam_genome_local))

    file_bam_genome_local = create_filename(dir_sample, N7, N5, 'bam_genome_local')
    file_sorted_bam_genome_local = create_filename(dir_sample, N7, N5, 'sorted_bam_genome_local')
    # file_sorted_bai_genome_local = create_filename(dir_sample, N7, N5, 'sorted_bai_genome_local')

    if not os.path.exists(os.path.dirname(file_bam_genome_local)):
        os.mkdir(os.path.dirname(file_bam_genome_local))

    # local alignment to the genome with bowtie2
    initial_dir = os.getcwd()

    bowtie2_command = ['bowtie2', '--local', '-p', str(ncpu),
                       '-X', '5000', '-k', '2', '-x', assembly,
                             '-1', file_R1, '-2', file_R2,
                             '-S', file_sam_genome_local]

    handle_sam_report_genome_local = open(file_sam_report_genome_local, 'wb')

    subprocess.call(bowtie2_command, stderr=handle_sam_report_genome_local)

    handle_sam_report_genome_local.close()

    # convert sam to bam
    sam_to_bam_genome_local_command = ['samtools', 'view', '-Sb', file_sam_genome_local]

    handle_file_bam_genome_local = open(file_bam_genome_local, 'wb')

    subprocess.call(sam_to_bam_genome_local_command, stdout=handle_file_bam_genome_local)

    # sort bam files
    sort_bam_genome_local_command = ['samtools', 'sort', file_bam_genome_local, '-o', file_sorted_bam_genome_local]

    subprocess.call(sort_bam_genome_local_command)

    # Clean up
    os.remove(file_sam_genome_local)
    os.remove(file_bam_genome_local)

    # Create bam index files
    create_bam_genome_local_index_command = ['samtools', 'index', file_sorted_bam_genome_local]
    subprocess.call(create_bam_genome_local_index_command)

    os.chdir(initial_dir)


############################
#
# Aligns reads globally to amplicon.
# Input: directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        file_genome_2bit, 2bit file with the reference genome being used
#
# ##########################
def align_amplicon(dir_sample, amplicon_info, check_plasmid_insertions, ncpu=4):

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    has_plasmid = type(amplicon_info['plasmid_sequence']) is str or type(amplicon_info['plasmid_sequence']) is unicode

    if check_plasmid_insertions == 1 and has_plasmid:
        file_R1 = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R1fastqgz')
        file_R2 = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R2fastqgz')
    else:
        file_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')
        file_R2 = create_filename(dir_sample, N7, N5, 'R2trimmed')

    if not os.path.exists(os.path.dirname(file_R1)):
        os.mkdir(os.path.dirname(file_R1))

    file_sam_amplicons = create_filename(dir_sample, N7, N5, 'sam_amplicons')
    file_sam_report_amplicons = create_filename(dir_sample, N7, N5, 'sam_report_amplicons')

    if not os.path.exists(os.path.dirname(file_sam_amplicons)):
        os.mkdir(os.path.dirname(file_sam_amplicons))

    file_bam_amplicons = create_filename(dir_sample, N7, N5, 'bam_amplicons')
    file_sorted_bam_amplicons = create_filename(dir_sample, N7, N5, 'sorted_bam_amplicons')

    if not os.path.exists(os.path.dirname(file_bam_amplicons)):
        os.mkdir(os.path.dirname(file_bam_amplicons))

    # global alignment to the amplicons with bowtie2
    initial_dir = os.getcwd()
    folder_amplicons = create_filename(dir_sample, N7, N5, 'amplicons')

    os.chdir(folder_amplicons)
    bowtie2_command = ['bowtie2', '-p', str(ncpu), '--very-sensitive',
                       '-X', '5000', '-k', '2', '-x', 'amplicons',
                       '-1', file_R1, '-2', file_R2,
                       '-S', file_sam_amplicons]

    handle_sam_report_amplicons = open(file_sam_report_amplicons, 'wb')

    subprocess.call(bowtie2_command, stderr=handle_sam_report_amplicons)

    handle_sam_report_amplicons.close()

    # convert sam to bam
    sam_to_bam_amplicons_command = ['samtools', 'view', '-Sb', file_sam_amplicons]

    handle_file_bam_amplicons = open(file_bam_amplicons, 'wb')

    subprocess.call(sam_to_bam_amplicons_command, stdout=handle_file_bam_amplicons)

    # sort bam files
    sort_bam_amplicons_command = ['samtools', 'sort', file_bam_amplicons, '-o', file_sorted_bam_amplicons]

    subprocess.call(sort_bam_amplicons_command)

    # Clean up
    os.remove(file_sam_amplicons)
    os.remove(file_bam_amplicons)

    # Create bam index files
    create_bam_amplicons_index_command = ['samtools', 'index', file_sorted_bam_amplicons]
    subprocess.call(create_bam_amplicons_index_command)

    os.chdir(initial_dir)


#################################################################################
# Function to create barcode dict
#################################################################################
def create_barcode_dict(filename):
    barcode_file = open_fastq_or_gz(filename)

    barcode_dict = dict()

    barcode_reads = itertools.izip(barcode_file)

    for header_barcode in barcode_reads:
        seq_barcode = barcode_reads.next()
        barcode_reads.next()
        qual_barcode = barcode_reads.next()
        barcode_dict[header_barcode[0].split()[0][1:]] = [seq_barcode[0].rstrip(), qual_barcode[0].rstrip()]

    return barcode_dict


################################################################################
# parse_indels
#
# Author: David Kelly
#
# Parse the CIGAR tuples for the positions and lengths of insertions and
# deletions implied by the aligned read.
#
# Input
#  aligned_read:  pysam AlignedSegment object
#
# Output
#  indels:        list of (position, indel_length) tuples. indel_lengths are
#                  positive for insertions, negative for deletions.
################################################################################
def parse_indels(aligned_read):
    indels = []

    if not aligned_read.is_unmapped and not aligned_read.is_secondary:  # We only look at primary alignments
        ref_i = aligned_read.reference_start
        for operation, length in aligned_read.cigartuples:
            if operation == 0:
                ref_i += length
            elif operation == 1:
                indels.append((ref_i, length))
            elif operation == 2:
                indels.append((ref_i, -length))
                ref_i += length
            else:
                print >> sys.stderr, 'Unrecognized CIGAR operation for %s' % aligned_read.query_name

    return indels


################################################################################
# Function to get the number of intersecting bases between two intervals
################################################################################
def get_intersection(region1_begin, region1_end, region2_begin, region2_end):
    list1 = range(int(region1_begin) + 1, int(region1_end) + 1)
    list2 = range(int(region2_begin) + 1, int(region2_end) + 1)
    return len(set(list1).intersection(list2))


################################################################################
# find_indels
#
# Input
#  bam_file:               alignments file to process
#
################################################################################
def find_indels(bam_file, strand, region_chr, region_start, region_end, UMI_dict, min_MAPQ, min_AS):
    bam_in_alignment_file = pysam.AlignmentFile(bam_file, 'rb')

    # We get the reads that overlap the window in which we make the counts
    # fetch will get reads with any overlap with the window
    # For UDiTaS, care must be take to ensure that the read covers the whole window, some short reads may cover just
    # one side of the window, depending on the direction of the UDiTaS primer
    bam_in = bam_in_alignment_file.fetch(region_chr, region_start, region_end)

    names_list = []
    position_list = []
    indel_list = []
    UMI_list = []
    for read in bam_in:
        # We add here a check to make sure the read came from the primer, crossed the cut and covered the whole window
        if read.has_tag('AS'):
            read_AS = read.get_tag('AS')
        # We test first if the read is unmapped, otherwise read_AS would be undefined
        if not read.is_unmapped and (((read.reference_start < region_start) and
                 (read.reference_end > region_end)) and
                    read.mapping_quality >= min_MAPQ and
                    read_AS >= min_AS and  not read.is_secondary):
            read_indels = parse_indels(read)
            # if no indels found, write 0
            if len(read_indels) == 0:
                read_indels.append(('-', 0))

            # print indels to table
            for pos, indel in read_indels:
                names_list.append(read.query_name)
                if pos == '-':
                    position_list.append(-1)
                else:
                    position_list.append(int(pos))

                indel_list.append(indel)
                UMI_list.append(UMI_dict[read.query_name][0])

    df = pd.DataFrame({'read_name': names_list,
                       'position': position_list,
                       'indel': indel_list,
                       'UMI': UMI_list})

    df['position_end'] = df.position + np.abs(df.indel)

    overlap = [get_intersection(df.loc[index]['position'], df.loc[index]['position_end'], region_start, region_end)
               for index in range(df.shape[0])]

    position_filter = np.array(overlap) > 0

    deletion_filter = position_filter & np.array(df.indel < 0)

    insertion_filter = position_filter & np.array(df.indel > 0)

    total_reads_in_region = len(set(df['read_name']))
    total_collapsed_reads_in_region = len(set(df['UMI']))

    total_indels = len(set(df.loc[position_filter]['read_name']))
    total_collapsed_indels = len(set(df.loc[position_filter]['UMI']))

    total_deletions = len(set(df.loc[deletion_filter]['read_name']))
    total_collapsed_deletions = len(set(df.loc[deletion_filter]['UMI']))

    total_insertions = len(set(df.loc[insertion_filter]['read_name']))
    total_collapsed_insertions = len(set(df.loc[insertion_filter]['UMI']))

    return [total_reads_in_region, total_indels, total_deletions, total_insertions,
            total_collapsed_reads_in_region,
            total_collapsed_indels,
            total_collapsed_deletions,
            total_collapsed_insertions]


################################################################################
# helper function to create list of fragment coordinates, useful to get size statistics
################################################################################
def create_segments(iter1, bam_in, min_MAPQ):
    segments = list()
    for read in iter1:
        if read.mapping_quality >= min_MAPQ and read.is_read2 and read.is_paired:
            if read.is_reverse:
                segment_start = read.reference_end + read.tlen  # Note: read.tlen is < 0
                segment_end = read.reference_end  # pysam is 0 based
            else:
                segment_start = read.reference_start  # pysam is 0 based
                segment_end = read.reference_start + read.tlen

            if segment_start < 0:  # We don't want to go below 0
                segment_start = 0

            if segment_end > segment_start:
                segments.append((bam_in.getrname(read.reference_id), segment_start,
                                 segment_end))

    return segments


################################################################################
# Function to analyse fragment sizes
################################################################################
def analyze_fragment_sizes(dir_sample, amplicon_info, min_MAPQ):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    file_sorted_bam_amplicons = create_filename(dir_sample, N7, N5, 'sorted_bam_amplicons')

    # Section to plot
    file_figure_classes = (create_filename(dir_sample, N7, N5, 'results_amplicons') + '_fragment_sizes_classes.pdf')
    file_figure = (create_filename(dir_sample, N7, N5, 'results_amplicons') + '_fragment_sizes.pdf')

    bam_in = pysam.AlignmentFile(file_sorted_bam_amplicons, "rb")

    iter_bam_in = bam_in.fetch()

    segments_bam_in = create_segments(iter_bam_in, bam_in, min_MAPQ)

    df = pd.DataFrame(segments_bam_in, columns=['type', 'begin', 'end'])
    df['length'] = df['end'] - df['begin']

    if df.shape[0] > 0:
        median_size = np.median(df['length'])

        file_fragments = (create_filename(dir_sample, N7, N5, 'results_amplicons') + '_fragment_sizes.xlsx')
        df.to_excel(file_fragments, index=False)

        df2 = df.pivot(columns='type', values='length')

        up_limit = 700
        fs = 20
        plt.rcParams["figure.figsize"] = [20, 10]

        # if df2.shape
        # plot with categories separated by colors
        df2.plot.hist(stacked=True, bins=np.arange(0, up_limit, 20))
        plt.xlim(0, up_limit)
        plt.xlabel('Fragment Size (bp)', fontsize=fs)
        plt.ylabel('Counts', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        if len(plt.gca().get_legend_handles_labels()[0]) > 0:  # To prevent warning from legend
            plt.legend(fontsize=12)
        pylab.savefig(file_figure_classes, bbox_inches='tight')
        plt.close(plt.gcf())

        # plot with all categories with the same color
        plt.hist(df['length'], np.arange(0, up_limit, 20))
        plt.xlim(0, up_limit)
        plt.xlabel('Fragment Size (bp)', fontsize=fs)
        plt.ylabel('Counts', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        if len(plt.gca().get_legend_handles_labels()[0]) > 0:  # To prevent warning from legend
            plt.legend(fontsize=12)
        pylab.savefig(file_figure, bbox_inches='tight')
        plt.close(plt.gcf())
    else:
        median_size = 0

    bam_in.close()

    return median_size


#######################################
# Function to get cut list depending on amplicon type
# Inputs are reaction_type and fasta record
# This function is needed to catch the cases when there are zero (control) or two cuts in the reference amplicon,
# eg when we have two cuts in the same chromosome
#
# The cut positions in the output depend on how the amplicons were constructed in create_amplicon. A change in
# that function must be matched in this function
#
# output: list with pairs of prefixes and cut positions [['cut1', pos1], 'cut2', pos2]]
#######################################
def get_cut_in_reference_amplicon_df(amplicon_info, reaction_type, record, strand, window_size,
                                     amplicon_window_around_cut):
    cut_in_reference_amplicon_df = pd.DataFrame(columns=['cut_type', 'cut_position'])
    if reaction_type == 'control':
        if strand == '+':
            # For the control sample (just primer, no cuts) we will look at counts and indels as if there was a cut
            # in position window_size
            cut_in_reference_amplicon_df.loc[0] = ['', window_size + 1]
        else:  # For the - strand we count reads at the end of the amplicon
            cut_in_reference_amplicon_df.loc[0] = ['', amplicon_window_around_cut - window_size - 2]
    elif reaction_type in ['single_cut', 'double_cut_different_chromosomes', 'triple_cut_different_chromosomes']:
        # We need to get cut1 in case its coordinates are smaller than amplicon_window_around_cut
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if cut1 < amplicon_window_around_cut:
            cut_site = cut1
        else:
            cut_site = amplicon_window_around_cut

        # all amplicons have a single cut at position amplicon_window_around_cut
        cut_in_reference_amplicon_df.loc[0] = ['cut1', cut_site]
    elif reaction_type == 'double_cut_same_chromosome':
        # Case two guides on the same chromosome
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if amplicon_info['strand_guide_2'] == '+':
            cut2 = amplicon_info['end_guide_2'] - 3
        elif amplicon_info['strand_guide_2'] == '-':
            cut2 = amplicon_info['start_guide_2'] + 3
        else:
            raise StrandError('strand_guide_2 can only have as values + or -')

        # We switch the coordinates of cut1 and cut2 if the guides are provided so that cut2 < cut1
        # cut1 it will always be smaller than cut2 and in the results cut1 will be the cut site with
        # smaller genomic coordinate
        # cut1 and cut2 also flipped in create_amplicon

        if cut2 < cut1:
            (cut1, cut2) = (cut2, cut1)
        if cut1 < amplicon_window_around_cut:
            cut_site = cut1
        else:
            cut_site = amplicon_window_around_cut

        # In this case some amplicons (wt, large_inversion) have two cuts, the rest have one
        if record.name in ['wt', 'large_inversion']:
            cut1_cut2_length = cut2 - cut1
            cut_in_reference_amplicon_df.loc[0] = ['cut1', cut_site]
            cut_in_reference_amplicon_df.loc[1] = ['cut2', cut_site + cut1_cut2_length]
        else:
            cut_in_reference_amplicon_df.loc[0] = ['cut1', cut_site]
    else:
        raise ReactionTypeError('Reaction type not yet supported by current version of UDiTaS')

    return cut_in_reference_amplicon_df


################################################################################
#  Function to analyze indels and structural rearrangements from aligned reads to amplicons
################################################################################
def analyze_alignments(dir_sample, amplicon_info, window_size, amplicon_window_around_cut, min_MAPQ, min_AS):

    reaction_type = get_reaction_type(amplicon_info)

    # UDiTaS primer strand
    strand = amplicon_info['strand']

    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')
    results_folder = os.path.join(exp_dir, 'results')
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    results_file = (create_filename(dir_sample, N7, N5, 'results_amplicons') + '_results_amplicon_window_'
                    + str(window_size) + '.xlsx')

    file_sorted_bam_amplicons = create_filename(dir_sample, N7, N5, 'sorted_bam_amplicons')

    file_UMI = create_filename(dir_sample, N7, N5, 'umifastqgz')

    UMI_dict = create_barcode_dict(file_UMI)

    filename_amplicons_fa = os.path.join(exp_dir, 'amplicons', 'amplicons.fa')

    results_df = pd.DataFrame(index=[0])

    # We get the reference amplicon list
    with open(filename_amplicons_fa, "rU") as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    for record in records:
        cut_df = get_cut_in_reference_amplicon_df(amplicon_info, reaction_type, record, strand, window_size,
                                                  amplicon_window_around_cut)

        # We go over cut1 and cut2 and when using cut1 we check if we are in control sample
        for i in cut_df.index:
            cut = cut_df.loc[i, 'cut_type']
            cut_position = cut_df.loc[i, 'cut_position']
            region_chr = record.name

            region_start = cut_position - window_size
            region_end = cut_position + window_size + 1

            results = find_indels(file_sorted_bam_amplicons, strand, region_chr, region_start, region_end, UMI_dict,
                                  min_MAPQ, min_AS)
            # This is to catch the control case, no cut there
            if len(cut) > 0:
                prefix = region_chr + '_' + cut
            else:
                prefix = region_chr

            results_df[prefix + '_total_reads'] = [results[0]]
            results_df[prefix + '_total_indels'] = [results[1]]
            results_df[prefix + '_total_deletions'] = [results[2]]
            results_df[prefix + '_total_insertions'] = [results[3]]
            results_df[prefix + '_total_reads_collapsed'] = [results[4]]
            results_df[prefix + '_total_indels_collapsed'] = [results[5]]
            results_df[prefix + '_total_deletions_collapsed'] = [results[6]]
            results_df[prefix + '_total_insertions_collapsed'] = [results[7]]

    median_size = analyze_fragment_sizes(dir_sample, amplicon_info, min_MAPQ)
    results_df['median_fragment_size'] = [median_size]

    results_df.to_excel(results_file)
    return results_df


################################################################################
# Function to calculate the number of reads and collapsed reads aligned to the reference amplicons
################################################################################
def analyze_alignments_all_amplicons(dir_sample, amplicon_info, min_MAPQ, min_AS):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')

    file_UMI = create_filename(dir_sample, N7, N5, 'umifastqgz')
    UMI_dict = create_barcode_dict(file_UMI)

    results_folder = os.path.join(exp_dir, 'results')
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    results_file = create_filename(dir_sample, N7, N5, 'results_all_amplicons')

    file_sorted_bam_amplicons = create_filename(dir_sample, N7, N5, 'sorted_bam_amplicons')

    bam_in_alignment_file = pysam.AlignmentFile(file_sorted_bam_amplicons, 'rb')

    bam_in_all = bam_in_alignment_file.fetch()

    names_list_amplicons = []
    UMI_list_amplicons = []

    for read in bam_in_all:
        if read.has_tag('AS'):
            read_AS = read.get_tag('AS')
        # We test first if the read is unmapped, otherwise read_AS would be undefined
        if not read.is_unmapped and read.mapping_quality >= min_MAPQ \
                and read_AS >= min_AS and not read.is_secondary:
            UMI_list_amplicons.append(UMI_dict[read.query_name][0])
            names_list_amplicons.append(read.query_name)

    all_amplicons_total_reads = len(set(names_list_amplicons))
    all_amplicons_total_reads_collapsed = len(set(UMI_list_amplicons))

    results_df = pd.DataFrame({'all_amplicons_total_reads': [all_amplicons_total_reads],
                               'all_amplicons_total_reads_collapsed': [all_amplicons_total_reads_collapsed]
                               },
                              columns=['all_amplicons_total_reads',
                                       'all_amplicons_total_reads_collapsed'])

    results_df.to_excel(results_file)

    return results_df


################################################################################
# Function to extract all unmapped reads to the amplicons
################################################################################
def extract_unmapped_reads_amplicons(dir_sample, amplicon_info):

    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    file_sorted_bam_amplicons = create_filename(dir_sample, N7, N5, 'sorted_bam_amplicons')

    file_unmapped_bam_amplicons = create_filename(dir_sample, N7, N5, 'unmapped_bam_amplicons')

    file_qsorted_unmapped_bam_amplicons = create_filename(dir_sample, N7, N5, 'qsorted_unmapped_bam_amplicons')

    file_R1_unmapped = create_filename(dir_sample, N7, N5, 'unmapped_amplicons_R1fastq')
    file_R2_unmapped = create_filename(dir_sample, N7, N5, 'unmapped_amplicons_R2fastq')
    file_unmapped_report = create_filename(dir_sample, N7, N5, 'unmapped_amplicons_report')

    if not os.path.exists(os.path.dirname(file_R1_unmapped)):
        os.mkdir(os.path.dirname(file_R1_unmapped))

    extract_unmapped_bam_command = ['samtools', 'view', '-b', '-f', '0x4', file_sorted_bam_amplicons, '-o',
                                    file_unmapped_bam_amplicons]

    subprocess.call(extract_unmapped_bam_command)

    qsort_unmapped_bam_command = ['samtools', 'sort', '-n', file_unmapped_bam_amplicons, '-o',
                                  file_qsorted_unmapped_bam_amplicons]

    subprocess.call(qsort_unmapped_bam_command)

    bamtofastq_command = ['bedtools', 'bamtofastq', '-i', file_qsorted_unmapped_bam_amplicons,
                          '-fq', file_R1_unmapped, '-fq2', file_R2_unmapped]

    handle_unmapped_report = open(file_unmapped_report, 'wb')
    subprocess.call(bamtofastq_command, stderr=handle_unmapped_report)

    for fo in [file_R1_unmapped, file_R2_unmapped]:
        with open(fo) as f_in, gzip.open(fo + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)
        os.remove(fo)


############################
#
# Aligns reads to the whole genome using bowtie2 and global alignment. This is used to find mispriming events
#  and unmapped reads
#
# Input: directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        assembly, name of the assembly to be used by bowtie2. It is convenient to place it in a folder especified by
#        the environmental variable BOWTIE2_INDEXES
#
# ##########################
def align_genome_global(dir_sample, amplicon_info, assembly, ncpu=4):

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    file_R1 = create_filename(dir_sample, N7, N5, 'unmapped_amplicons_R1fastqgz')
    file_R2 = create_filename(dir_sample, N7, N5, 'unmapped_amplicons_R2fastqgz')

    file_sam_genome_global = create_filename(dir_sample, N7, N5, 'sam_genome_global')
    file_sam_report_genome_global = create_filename(dir_sample, N7, N5, 'sam_report_genome_global')

    if not os.path.exists(os.path.dirname(file_sam_genome_global)):
        os.mkdir(os.path.dirname(file_sam_genome_global))

    file_bam_genome_global = create_filename(dir_sample, N7, N5, 'bam_genome_global')
    file_sorted_bam_genome_global = create_filename(dir_sample, N7, N5, 'sorted_bam_genome_global')

    if not os.path.exists(os.path.dirname(file_bam_genome_global)):
        os.mkdir(os.path.dirname(file_bam_genome_global))

    # global alignment to the genome with bowtie2
    initial_dir = os.getcwd()

    bowtie2_command = ['bowtie2', '--very-sensitive', '-p', str(ncpu),
                       '-X', '5000', '-k', '2', '-x', assembly,
                       '-1', file_R1, '-2', file_R2,
                       '-S', file_sam_genome_global]

    handle_sam_report_genome_global = open(file_sam_report_genome_global, 'wb')

    subprocess.call(bowtie2_command, stderr=handle_sam_report_genome_global)

    handle_sam_report_genome_global.close()

    # convert sam to bam
    sam_to_bam_genome_global_command = ['samtools', 'view', '-Sb', file_sam_genome_global]

    handle_file_bam_genome_global = open(file_bam_genome_global, 'wb')

    subprocess.call(sam_to_bam_genome_global_command, stdout=handle_file_bam_genome_global)

    # sort bam files
    sort_bam_genome_global_command = ['samtools', 'sort', file_bam_genome_global, '-o', file_sorted_bam_genome_global]

    subprocess.call(sort_bam_genome_global_command)

    # Clean up
    os.remove(file_sam_genome_global)
    os.remove(file_bam_genome_global)

    # Create bam index files
    create_bam_genome_global_index_command = ['samtools', 'index', file_sorted_bam_genome_global]
    subprocess.call(create_bam_genome_global_index_command)

    os.chdir(initial_dir)


################################################################################
# Function to analyze global alignments to the genome
################################################################################
def analyze_alignments_genome_global(dir_sample, amplicon_info, min_MAPQ, min_AS,  file_genome_2bit):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    genome = twobitreader.TwoBitFile(file_genome_2bit)

    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')

    file_UMI = create_filename(dir_sample, N7, N5, 'umifastqgz')
    UMI_dict = create_barcode_dict(file_UMI)

    results_folder = os.path.join(exp_dir, 'results')
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    results_file = create_filename(dir_sample, N7, N5, 'results_genomewide')

    file_sorted_bam_genome_global = create_filename(dir_sample, N7, N5, 'sorted_bam_genome_global')

    bam_in_alignment_file = pysam.AlignmentFile(file_sorted_bam_genome_global, 'rb')

    bam_in_all = bam_in_alignment_file.fetch()

    names_list_genome = []
    UMI_list_genome = []

    for read in bam_in_all:
        if read.has_tag('AS'):
            read_AS = read.get_tag('AS')
        # We test first if the read is unmapped, otherwise read_AS would be undefined
        if not read.is_unmapped and (read.mapping_quality >= min_MAPQ
                                     and read_AS >= min_AS and not read.is_secondary):
            UMI_list_genome.append(UMI_dict[read.query_name][0])
            names_list_genome.append(read.query_name)

    names_list_target_only = []
    UMI_list_target_only = []

    fetch_window = 1000
    fetch_chr = amplicon_info['chr']

    if amplicon_info['strand'] == '+':  # This is the UDiTaS oligo strand
        fetch_start = amplicon_info['start']
        fetch_end = amplicon_info['end'] + fetch_window
        if fetch_end > len(genome[amplicon_info['chr']]):
            fetch_end = len(genome[amplicon_info['chr']])

    elif amplicon_info['strand'] == '-':
        fetch_start = amplicon_info['start'] - fetch_window
        fetch_end = amplicon_info['end']
        if fetch_start < 0:
            fetch_start = 0

    bam_in_target = bam_in_alignment_file.fetch(fetch_chr, fetch_start, fetch_end)

    for read in bam_in_target:
        if read.mapping_quality >= min_MAPQ and not read.is_unmapped and not read.is_secondary:
            UMI_list_target_only.append(UMI_dict[read.query_name][0])
            names_list_target_only.append(read.query_name)

    genomewide_total_reads = len(set(names_list_genome))
    genomewide_total_reads_collapsed = len(set(UMI_list_genome))
    genomewide_target_only_reads = len(set(names_list_target_only))
    genomewide_target_only_reads_collapsed = len(set(UMI_list_target_only))

    results_df = pd.DataFrame({'genomewide_total_reads': [genomewide_total_reads],
                               'genomewide_total_reads_collapsed': [genomewide_total_reads_collapsed],
                               'genomewide_target_only_reads': [genomewide_target_only_reads],
                               'genomewide_target_only_reads_collapsed': [genomewide_target_only_reads_collapsed]
                               },
                              columns=['genomewide_total_reads',
                                       'genomewide_total_reads_collapsed',
                                       'genomewide_target_only_reads',
                                       'genomewide_target_only_reads_collapsed'])

    results_df.to_excel(results_file)

    return results_df


################################################################################
#  Simple call to samtools to count mapped & primary reads (thus 0x104)
################################################################################
def samtools_count_primary_mapped(bam_file):
    samtools_out = subprocess.check_output(['samtools', 'view', '-c', '-F', '0x104', bam_file])
    return int(samtools_out)


################################################################################
# Function to summarize counts into all amplicons
################################################################################
def summarize_results(results):

    total_reads_list = [k for k in results.keys() if str(k).endswith('_total_reads')]
    total_reads_collapsed_list = [k for k in results.keys() if str(k).endswith('_total_reads_collapsed')]

    results['total_aligned_junctions'] = results[total_reads_list].sum(axis=1)

    for k in total_reads_list:
        # We add np.finfo(float).eps to prevent dividing by 0
        results[k + '_percent'] = 100 * results[k] / (results['total_aligned_junctions'] +
                                                      np.finfo(float).eps)

    results['total_aligned_junctions_collapsed'] = results[total_reads_collapsed_list].sum(axis=1)

    for k in total_reads_collapsed_list:
        results[k + '_percent'] = 100 * results[k] / (results['total_aligned_junctions_collapsed'] +
                                                      np.finfo(float).eps)

    return results


#######################################################################################
# We pivot the table using melt for easier visualization with other tools like Tableau
#######################################################################################
def melt_results(results_summary_with_experiments):

    melt_list = [k for k in results_summary_with_experiments.keys() if
                 str(k).endswith('_total_reads_collapsed_percent')]

    frozen_list = list(results_summary_with_experiments)

    for el in melt_list:
        frozen_list.remove(el)

    results_out = pd.melt(results_summary_with_experiments,
                          value_vars=melt_list,
                          id_vars=frozen_list, var_name='Type', value_name='Percent Editing')

    return results_out


################################################################################
#  Function to summarize counts into the genome
################################################################################
def summarize_results_genome(read_count, all_results_genome_global_df, results_summary):

    df = pd.concat([read_count, results_summary, all_results_genome_global_df], axis=1)

    df['total_reads_aligned'] = (df['total_reads_genome'] + df['total_aligned_amplicon_reads'] +
                                 df['total_reads_plasmid_only'])

    df['percent_reads_aligned'] = 100 * df['total_reads_aligned'] / df['read_count']

    df['percent_reads_aligned_plasmid'] = 100 * df['total_reads_plasmid_only'] / df['total_reads_aligned']
    df['percent_reads_aligned_amplicon'] = 100 * df['total_aligned_amplicon_reads'] / df['total_reads_aligned']
    df['percent_reads_mispriming'] = 100 * df['total_reads_genome'] / df['total_reads_aligned']

    results_out = df[['read_count', 'total_reads_aligned', 'percent_reads_aligned',
                      'total_reads_plasmid_only',
                      'total_aligned_amplicon_reads',
                      'total_reads_genome',
                      'percent_reads_aligned_plasmid',
                      'percent_reads_aligned_amplicon',
                      'percent_reads_mispriming',
                      'median_fragment_size']]

    return results_out


#######################################################################################
# We pivot the table using melt for easier visualization with other tools like Tableau
#######################################################################################
def melt_results_genome(results_summary_genome):

    melt_list = ['percent_reads_aligned_plasmid', 'percent_reads_aligned_amplicon', 'percent_reads_mispriming']

    frozen_list = list(results_summary_genome)

    for el in melt_list:
        frozen_list.remove(el)

    results_out = pd.melt(results_summary_genome,
                          value_vars=melt_list,
                          id_vars=frozen_list, var_name='Type', value_name='Percent Alignment')

    return results_out


################################################################################
# Fast way to count reads, but only works on unix
################################################################################
def wc_unix(filename):
    cat_out = subprocess.Popen(('zcat', filename), stdout=subprocess.PIPE)
    return int(int(subprocess.check_output(["wc", "-l"], stdin=cat_out.stdout).split()[0])/4.)


################################################################################
# Function to count reads
################################################################################
def count_reads(dir_sample, amplicon_info):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    read_counts_file = create_filename(dir_sample, N7, N5, 'read_counts')
    file_cutadapt_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')

    rc = wc_unix(file_cutadapt_R1)
    df = pd.DataFrame({'read_count': [rc]})

    df.to_excel(read_counts_file)

    return rc


################################################################################
#  Function to get the percentages for the alignment of all reads mapped to ALL plasmid, amplicons and genomewide
################################################################################
def get_summary_all_alignments(dir_sample, amplicon_info, read_count, result_plasmid_df,
               result_reads_in_all_amplicons_df, results_genome_global_df):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    summary_all_alignments_file = create_filename(dir_sample, N7, N5, 'summary_all_alignments')

    summary_all_alignments = pd.concat([read_count, result_plasmid_df,
               result_reads_in_all_amplicons_df, results_genome_global_df], axis=1)

    total_reads_list = [k for k in summary_all_alignments.keys() if str(k).endswith('_total_reads')]

    summary_all_alignments['total_aligned'] = summary_all_alignments[total_reads_list].sum(axis=1)

    total_reads_collapsed_list = [k for k in summary_all_alignments.keys() if str(k).endswith('_total_reads_collapsed')]

    summary_all_alignments['total_aligned_collapsed'] = summary_all_alignments[total_reads_collapsed_list].sum(axis=1)

    summary_all_alignments['percent_aligned'] = 100 * (summary_all_alignments['total_aligned'] /
                                                       summary_all_alignments['read_count'])
    summary_all_alignments['percent_aligned_all_amplicons'] = 100 * (summary_all_alignments['all_amplicons_total_reads'] /
                                                        summary_all_alignments['total_aligned'])
    summary_all_alignments.to_excel(summary_all_alignments_file, index=False)

    return summary_all_alignments


#######################################################################################
# We pivot the table using melt for easier visualization with other tools like Tableau
#######################################################################################
def melt_big_results(big_results):

    melt_list1 = [k for k in big_results.keys() if
                 str(k).endswith('_total_reads_collapsed_percent')]

    melt_list2 = ['percent_aligned', 'percent_aligned_all_amplicons']

    melt_list = melt_list1 + melt_list2

    frozen_list = list(big_results)

    for el in melt_list:
        frozen_list.remove(el)

    results_out = pd.melt(big_results,
                          value_vars=melt_list,
                          id_vars=frozen_list, var_name='Type', value_name='Percent Alignment')

    return results_out


################################################################################
# Final summary
################################################################################
def summarize_big_results(big_results_in):

    results = big_results_in.copy()

    # We have the list of all columns to operate, if not present fill with nan
    # TODO add columns to handle three cuts in different chromosomes
    cols_list = ['target_plus_plasmid_total_reads_collapsed',
                 'wt_cut1_total_reads_collapsed', 'wt_cut1_total_indels_collapsed',
                 'wt_cut2_total_reads_collapsed', 'wt_cut2_total_indels_collapsed',
                 'large_deletion_cut1_total_reads_collapsed',
                 'large_inversion_cut1_total_reads_collapsed', 'large_inversion_cut2_total_reads_collapsed',
                 '1a_1a_cut1_total_reads_collapsed', '1b_1b_cut1_total_reads_collapsed',
                 '2b_2b_cut1_total_reads_collapsed',
                 '1a_1b_cut1_total_reads_collapsed', '1a_1b_cut1_total_indels_collapsed',
                 '2a_2b_cut1_total_reads_collapsed', '2a_2b_cut1_total_indels_collapsed',
                 '1a_2a_cut1_total_reads_collapsed', '1a_2b_cut1_total_reads_collapsed',
                 '2a_1b_cut1_total_reads_collapsed', '2a_2a_cut1_total_reads_collapsed',
                 '2b_1b_cut1_total_reads_collapsed', '2b_2b_cut1_total_reads_collapsed']

    for c in cols_list:
        if c not in results.keys():
            results[c] = np.nan

    results['AAV/Plasmid Integrations reads_collapsed'] = results['target_plus_plasmid_total_reads_collapsed']

    results['WT reads_collapsed'] = (results[['wt_cut1_total_reads_collapsed',
                                              'wt_cut2_total_reads_collapsed']].sum(axis=1) -
                                     results[['wt_cut1_total_indels_collapsed',
                                              'wt_cut2_total_indels_collapsed']].sum(axis=1))

    results['Small Indels reads_collapsed'] = (results[['wt_cut1_total_indels_collapsed',
                                                        'wt_cut2_total_indels_collapsed']].sum(axis=1))

    results['Large Deletion reads_collapsed'] = (results['large_deletion_cut1_total_reads_collapsed'])

    results['Inversion reads_collapsed'] = (results[['large_inversion_cut1_total_reads_collapsed',
                                                     'large_inversion_cut2_total_reads_collapsed']].sum(axis=1))

    # Handle two cuts on different chromosomes here
    # We split the indels for the wt amplicons
    results['WT 1a_1b reads_collapsed'] = (results['1a_1b_cut1_total_reads_collapsed'] -
                                           results['1a_1b_cut1_total_indels_collapsed'])
    results['1a_1b indels reads_collapsed'] = results['1a_1b_cut1_total_indels_collapsed']

    results['WT 2a_2b reads_collapsed'] = (results['2a_2b_cut1_total_reads_collapsed'] -
                                           results['2a_2b_cut1_total_indels_collapsed'])
    results['2a_2b indels reads_collapsed'] = results['2a_2b_cut1_total_indels_collapsed']

    # This is the list of all the columns to add to get the total number of aligned reads at the junctions
    total_reads_list = (['AAV/Plasmid Integrations reads_collapsed', 'WT reads_collapsed',
                         'Small Indels reads_collapsed', 'Large Deletion reads_collapsed',
                         'Inversion reads_collapsed',
                         'WT 1a_1b reads_collapsed', '1a_1b indels reads_collapsed',
                         'WT 2a_2b reads_collapsed', '2a_2b indels reads_collapsed',
                         '1a_1a_cut1_total_reads_collapsed', '1b_1b_cut1_total_reads_collapsed',
                         '2b_2b_cut1_total_reads_collapsed',
                         '1a_2a_cut1_total_reads_collapsed', '1a_2b_cut1_total_reads_collapsed',
                         '2a_1b_cut1_total_reads_collapsed', '2a_2a_cut1_total_reads_collapsed',
                         '2b_1b_cut1_total_reads_collapsed'])

    results['total_aligned_junctions_collapsed'] = results[total_reads_list].sum(axis=1)

    # We use eps to prevent errors when dividing by zero. It can happen for samples with zero reads
    eps = np.finfo(float).eps
    results['AAV/Plasmid Integrations'] = (100*results['AAV/Plasmid Integrations reads_collapsed']/
                                         (results['total_aligned_junctions_collapsed'] + eps))

    results['WT'] = (100 * results['WT reads_collapsed'] / (results['total_aligned_junctions_collapsed'] + eps))

    results['Small Indels'] = (100 * results['Small Indels reads_collapsed'] /
                               (results['total_aligned_junctions_collapsed'] + eps))

    results['Large Deletion'] = (100 * results['Large Deletion reads_collapsed'] /
                                 (results['total_aligned_junctions_collapsed'] + eps))

    results['Inversion'] = (100 * results['Inversion reads_collapsed'] /
                            (results['total_aligned_junctions_collapsed']) + eps)

    results['WT 1a_1b'] = (100 * results['WT 1a_1b reads_collapsed'] /
                           (results['total_aligned_junctions_collapsed'] + eps))

    results['1a_1b indels'] = (100 * results['1a_1b indels reads_collapsed'] /
                               (results['total_aligned_junctions_collapsed'] + eps))

    results['WT 2a_2b'] = (100 * results['WT 2a_2b reads_collapsed'] /
                           (results['total_aligned_junctions_collapsed'] + eps))

    results['2a_2b indels'] = (100 * results['2a_2b indels reads_collapsed'] /
                               (results['total_aligned_junctions_collapsed'] + eps))

    results['1a_1a'] = (100 * results['1a_1a_cut1_total_reads_collapsed'] /
                        (results['total_aligned_junctions_collapsed'] + eps))

    results['1b_1b'] = (100 * results['1b_1b_cut1_total_reads_collapsed'] /
                        (results['total_aligned_junctions_collapsed'] + eps))

    results['2b_2b'] = (100 * results['2b_2b_cut1_total_reads_collapsed'] /
                        (results['total_aligned_junctions_collapsed'] + eps))

    results['1a_2a'] = (100 * results['1a_2a_cut1_total_reads_collapsed'] /
                        (results['total_aligned_junctions_collapsed'] + eps))

    results['1a_2b'] = (100 * results['1a_2b_cut1_total_reads_collapsed'] /
                        (results['total_aligned_junctions_collapsed'] + eps))

    results['2a_1b'] = (100 * results['2a_1b_cut1_total_reads_collapsed'] /
                        (results['total_aligned_junctions_collapsed'] + eps))

    results['2a_2a'] = (100 * results['2a_2a_cut1_total_reads_collapsed'] /
                        (results['total_aligned_junctions_collapsed'] + eps))

    results['2b_1b'] = (100 * results['2b_1b_cut1_total_reads_collapsed'] /
                        (results['total_aligned_junctions_collapsed'] + eps))

    results['total_aligned_amplicons'] = results['all_amplicons_total_reads']
    results['total_aligned_amplicons_collapsed'] = results['all_amplicons_total_reads_collapsed']

    list1 = ['read_count',
             'total_aligned', 'percent_aligned', 'total_aligned_collapsed',
             'total_aligned_amplicons', 'percent_aligned_all_amplicons', 'total_aligned_amplicons_collapsed',
             'total_aligned_junctions', 'total_aligned_junctions_collapsed',
             'median_fragment_size']

    list2 = ['AAV/Plasmid Integrations', 'WT', 'Small Indels', 'Large Deletion', 'Inversion', 'WT 1a_1b',
             '1a_1b indels', 'WT 2a_2b', '2a_2b indels', '1a_1a', '1b_1b', '2b_2b', '1a_2a', '1a_2b',
             '2a_1b', '2a_2a', '2b_1b']

    final_list = total_reads_list + list1 + list2

    results_out = results[final_list]

    return results_out


#######################################################################################
# We pivot the table using melt for easier visualization with other tools like Tableau
#######################################################################################
def melt_summarize_big_results(summarized_big_results):

    melt_list1 = ['AAV/Plasmid Integrations', 'WT', 'Small Indels',
                  'Large Deletion', 'Inversion',
                  'WT 1a_1b', '1a_1b indels', 'WT 2a_2b', '2a_2b indels',
                 '1a_1a', '1b_1b', '2b_2b', '1a_2a', '1a_2b', '2a_1b', '2a_2a', '2b_1b']

    melt_list2 = ['percent_aligned', 'percent_aligned_all_amplicons']

    melt_list = melt_list1 + melt_list2

    frozen_list = list(summarized_big_results)

    for el in melt_list:
        frozen_list.remove(el)

    results_out = pd.melt(summarized_big_results,
                          value_vars=melt_list,
                          id_vars=frozen_list, var_name='Editing Type', value_name='Editing Percent')

    return results_out
