#!/usr/bin/env python
"""
    Code to process UDiTaS data
"""

# imports here
from __future__ import print_function

import pandas as pd

import uditas_helpers

import argparse

import os


from _version import __version__

# Metadata
__author__ = "Eugenio Marco"
__credits__ = ["David Kelly", "Hari Jayaram", "Luis Barrera"]
__status__ = "Production"
__copyright__ = "Copyright 2018, Editas Medicine, Inc"


def main():
    # Parse command line options
    parser = argparse.ArgumentParser(description='Process UDiTaS data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("dir_sample", help='Directory with the sample to be processed')
    parser.add_argument("-folder_genome_2bit", help=('Folder containing the 2bit file(s) with the ' +
                                                    'reference genome being used'),
                        default=os.environ.get('GENOMES_2BIT', os.path.join('/', 'reference', 'genomes_2bit')))
    parser.add_argument("-skip_demultiplexing", help='Skip demultiplexing? Options: 0, 1 (skip)', default=0)
    parser.add_argument("-skip_trimming", help='Skip adapter trimming? Options: 0, 1 (skip)', default=0)
    parser.add_argument("-skip_genome_local_alignment", help='Skip genome-wide local alignment? Options: 0 , 1 (skip)',
                        default=1)
    parser.add_argument("-skip_genome_global_alignment", help='Skip genome-wide global alignment? Options: 0 , 1 (skip)',
                        default=0)
    parser.add_argument("-process_amplicon",
                        help='Select row number (0-based) of amplicon to process, set to all to process all amplicons',
                        default='all')
    parser.add_argument("-skip_amplicon_global_alignment", help='Skip amplicon global alignment? Options: 0, 1 (skip)',
                        default=0)
    parser.add_argument("-check_plasmid_insertions",
                        help=('Check for plasmid insertions. Options: 0 (skip), 1 plamid_name ' +
                              'and plasmid_sequence required in sample_info.csv'), default=1)
    parser.add_argument("-skip_plasmid_alignment",
                        help=('Skip plasmid alignment? Note, just alignment. ' +
                              'Counts still evaluated. Options: 0, 1 (skip)'),
                        default=0)
    parser.add_argument("-ncpu", help='Number of CPUs to use', default=4)
    parser.add_argument("-window_size", help='Window size around cut sites used to grab UDiTaS reads', default=15)
    parser.add_argument("-default_amplicon_window_around_cut",
                        help='Window size around cut sites used to create amplicons',
                        default=1000)
    parser.add_argument("-min_MAPQ", help='Minimum mapping quality to include a read', default=5)
    parser.add_argument("-min_AS", help='Minimum alignment score to include a read', default=-180)
    parser.add_argument("-process_AMP_seq_run", help='Set to 1 to process an AMP-seq run using GUIDE-seq adapters', default=0)

    args = parser.parse_args()

    dir_sample = os.path.abspath(os.path.join(os.getcwd(), args.dir_sample))
    folder_genome_2bit = args.folder_genome_2bit
    skip_demultiplexing = int(args.skip_demultiplexing)
    skip_trimming = int(args.skip_trimming)
    skip_genome_local_alignment = int(args.skip_genome_local_alignment)
    skip_genome_global_alignment = int(args.skip_genome_global_alignment)
    process_amplicon = args.process_amplicon
    skip_amplicon_global_alignment = int(args.skip_amplicon_global_alignment)
    check_plasmid_insertions = int(args.check_plasmid_insertions)
    skip_plasmid_alignment = int(args.skip_plasmid_alignment)
    ncpu = int(args.ncpu)
    window_size = int(args.window_size)
    default_amplicon_window_around_cut = int(args.default_amplicon_window_around_cut)
    min_MAPQ = int(args.min_MAPQ)
    min_AS = int(args.min_AS)
    process_AMP_seq_run = int(args.process_AMP_seq_run)

    print('\nUDiTaS version ' + __version__ + '\n\n')

    # Call to demultiplexing step
    if skip_demultiplexing == 0:
        uditas_helpers.demultiplex(dir_sample)

    sample_info_filename = os.path.join(dir_sample, 'sample_info.csv')

    experiments = pd.read_csv(sample_info_filename)
    read_count = pd.DataFrame()
    # Make list of samples to process
    processing_all_amplicons = (process_amplicon == 'all')
    if processing_all_amplicons:
        processing_list = experiments.index
    else:
        processing_list = [int(process_amplicon)]

    print('Processing UDiTaS data in folder ' + dir_sample)

    for i in processing_list:
        amplicon_info = experiments.loc[i]
        assembly = amplicon_info['genome']
        file_genome_2bit = os.path.join(folder_genome_2bit, assembly + '.2bit')
        if skip_trimming == 0:
            uditas_helpers.trim_fastq(dir_sample, amplicon_info, process_AMP_seq_run)
        # We check that we have plasmid_sequence
        has_plasmid = (type(amplicon_info['plasmid_sequence']) is str or
                       type(amplicon_info['plasmid_sequence']) is unicode)
        do_plasmid = check_plasmid_insertions == 1 and has_plasmid
        # Align to plasmid if necessary
        if (do_plasmid and not skip_plasmid_alignment):
            # We align to plasmid reference and extract unmapped reads to map to our amplicons
            uditas_helpers.create_plasmid_reference(dir_sample, amplicon_info)
            uditas_helpers.align_plasmid_local(dir_sample, amplicon_info, ncpu)
            uditas_helpers.extract_unmapped_reads_plasmid(dir_sample, amplicon_info)

        # Get result_plasmid_df, if do_plasmid == False function returns NaN
        result_plasmid_df = uditas_helpers.analyze_alignments_plasmid(dir_sample, amplicon_info, min_MAPQ,
                                                                      file_genome_2bit, do_plasmid)
        result_plasmid_df.index = [i]
        
        # Align locally genomewide
        if skip_genome_local_alignment == 0:
            uditas_helpers.align_genome_local(dir_sample, amplicon_info, assembly, check_plasmid_insertions, ncpu)
        # Map globally to amplicons built out of the expected cuts
        if skip_amplicon_global_alignment == 0:
            uditas_helpers.create_amplicon(dir_sample, amplicon_info, file_genome_2bit,
                                           default_amplicon_window_around_cut)
            uditas_helpers.align_amplicon(dir_sample, amplicon_info, check_plasmid_insertions, ncpu)
            # Now we extract the rest of unmapped reads
            uditas_helpers.extract_unmapped_reads_amplicons(dir_sample, amplicon_info)
        # Get counts AT THE CUT JUNCTIONS in all amplicons
        result_amplicon_df = uditas_helpers.analyze_alignments(dir_sample, amplicon_info, window_size,
                                                               default_amplicon_window_around_cut, min_MAPQ, min_AS)
        result_amplicon_df.index = [i]

        # We get counts of reads aligned to all amplicons, not just at the junctions
        result_reads_in_all_amplicons_df = uditas_helpers.analyze_alignments_all_amplicons(dir_sample, amplicon_info,
                                                                                           min_MAPQ, min_AS)
        result_reads_in_all_amplicons_df.index = [i]

        results_alignments_junction = pd.concat([result_amplicon_df, result_plasmid_df], axis=1)

        if processing_all_amplicons:
            try:
                all_results_alignments_junction = pd.concat([all_results_alignments_junction,
                                                             results_alignments_junction], axis=0)
            except NameError:  # Initialize
                all_results_alignments_junction = results_alignments_junction.copy()

        # We can now align the reads to the genome, with the option of using a genomic_background that is different from
        # assembly, useful in cases of transgenic mice, where the target region and primers are human (assembly), but
        # the genomic background is different (genomic_background). This is used if the field genomic_background is
        # present in sample_info.csv
        has_genomic_background = (type(amplicon_info['genomic_background']) is str or
                                  type(amplicon_info['genomic_background']) is unicode)

        if has_genomic_background:
            genomic_background = amplicon_info['genomic_background']
        else:
            genomic_background = assembly

        if skip_genome_global_alignment == 0:
            uditas_helpers.align_genome_global(dir_sample, amplicon_info, genomic_background, ncpu)
        results_genome_global_df = uditas_helpers.analyze_alignments_genome_global(dir_sample, amplicon_info, min_MAPQ,
                                                                                   min_AS, file_genome_2bit)
        results_genome_global_df.index = [i]

        if processing_all_amplicons:
            try:
                all_results_genome_global_df = pd.concat([all_results_genome_global_df, results_genome_global_df], axis=0)
            except NameError:  # Initialize
                all_results_genome_global_df = results_genome_global_df.copy()

        read_count.loc[i,'read_count'] = uditas_helpers.count_reads(dir_sample, amplicon_info)

        # Calculate and write percentages of all alignments
        summary_all_alignments = uditas_helpers.get_summary_all_alignments(dir_sample, amplicon_info,
                                                                           read_count.loc[[i]], result_plasmid_df,
                                                                           result_reads_in_all_amplicons_df,
                                                                           results_genome_global_df)

        try:
            all_samples_summary_alignments = pd.concat([all_samples_summary_alignments, summary_all_alignments], axis=0)
        except NameError:  # Initialize
            all_samples_summary_alignments = summary_all_alignments.copy()

    if processing_all_amplicons:
        all_results_dir = os.path.join(dir_sample, 'all_results')
        if not os.path.exists(all_results_dir):
            os.mkdir(all_results_dir)
        results_all = pd.concat([experiments, all_results_alignments_junction], axis=1)
        results_all['version'] = __version__
        results_all['args'] = repr(vars(args))
        results_all.to_excel(os.path.join(all_results_dir, 'results_combined_detailed.xlsx'))
        results_summary = uditas_helpers.summarize_results(all_results_alignments_junction)
        results_summary_with_experiments = pd.concat([experiments, results_summary], axis=1)

        results_summary_with_experiments['version'] = __version__
        results_summary_with_experiments['args'] = repr(vars(args))
        results_summary_with_experiments.to_excel(os.path.join(all_results_dir, 'results_summary.xlsx'))

        results_pivot = uditas_helpers.melt_results(results_summary_with_experiments)
        results_pivot.to_excel(os.path.join(all_results_dir, 'results_summary_pivot.xlsx'))

        cols_to_use = all_samples_summary_alignments.columns.difference(results_summary_with_experiments.columns)

        big_results = pd.merge(results_summary_with_experiments, all_samples_summary_alignments[cols_to_use],
                               left_index=True, right_index=True, how='outer')

        summarized_big_results = uditas_helpers.summarize_big_results(big_results)
        experiments['version'] = __version__
        experiments['args'] = repr(vars(args))

        dir_sample_basename = os.path.basename(os.path.normpath(dir_sample))
        final_results = pd.concat([experiments, summarized_big_results], axis=1)
        final_results.to_excel(os.path.join(all_results_dir, dir_sample_basename +
                                            '_' + 'big_results.xlsx'), index=None)

        final_results_pivot = uditas_helpers.melt_summarize_big_results(final_results)
        final_results_pivot.to_excel(os.path.join(all_results_dir, dir_sample_basename +
                                            '_' + 'big_results_pivot.xlsx'), index=None)


if __name__ == "__main__":
    main()