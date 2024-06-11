#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:54:08 2022

@author: MG, MK
"""

import os
import sys
import glob
import shutil
import time
import csv
import re
import random
import pandas as pd
import numpy as np
import warnings
from Bio import SeqIO
from multiprocessing import Pool
from Bio import BiopythonWarning
from tqdm import tqdm

syntax = '''-------------------------------------------------------------------------------
3 June 2024 v1.00
Originally developed for Python 3.10.14

 _____ _____ _____ _____ _____ _____ 
|   | |     |   __|  _  | __  |  _  |
| | | |-   -|__   |     |    -|     |
|_|___|_____|_____|__|__|__|__|__|__| v1.0

Network Analysis of Structure Similarity Alignment V1.0

Syntax:
python NISARA_v1.0.py [base_name] [threshold_cutoff] 
[exponent_threshold] [exponent_degree] [threads]

:param base_name: Base name for the analysis.

:param threshold_cutoff:    Cut-off threshold above which TM-score values are
                            retained and considered for connections (th-cu)

:param exponent_threshold:  The threshold above which values are elevated 
                            according to the exponent_degree (ex-th)

:param exponent_degree:     A selector for the exponent used in 
                            the transformation
                            (ex-de)
                            1 -> exponent = 2
                            2 -> exponent = 12
                            3 -> exponent = 24

:param threads: Number of threads used for parallel TM-align comparisons.

Transformation Logic:
    TM-score values are logarithmized, incremented by 0.01 to avoid zeros, and 
    multiplied by -1 to convert to positive values.
    If the TM-score is higher than the specified threshold, it is additionally
    raised to a power specified by exponent_degree.

    def transform_value(value, threshold, exponent):
        if value > threshold:
            return ((np.log(value) + 0.01) * -1) ** exponent
        else:
            return (np.log(value) + 0.01) * -1

    value:      TM-score
    threshold:  The threshold above which values are amplified according to the 
                exponent_degree
    exponent:   The exponent value used in the transformation

How to Use:
    1. Prepare a folder named 'input' containing pdb files.
    2. Execute the command: python NISARA_v1.0.py name_x 0.5 0.6 2 24

Input:
    *.pdb in prepared folder named 'input'
    Try to keep the names short and uncomplicated.
    
Output:
    -> *_th-cu=*_ex-th=*_ex-de=*_all_results_tmalign_unmaped.tsv
    Full table of TM-align comparisons  of protein pairs. Unmaped names
    -> *_th-cu=*_ex-th=*_ex-de=*_all_results_tmalign_filtered_transformed.tsv'
    Table of TM-align comparisons of protein pairs
    after cutoff and transformation
    -> *_th-cu=*_ex-th=*_ex-de=*_structure-network.clans
    CLANS file
    

Limitations:
    Ensure that structures are in a single piece, not fragmented. For example, 
    e4jw1A3 includes a break in the sequence (A:10-104, A:132-308), which can 
    cause TM-align to only align a fragment of the sequence, potentially leading 
    to less accurate results.

Simple Installation of Requirements with Conda:
    conda install conda-forge::biopython=1.78
    conda install bioconda::tmalign=20170708
    conda install anaconda::pandas=23.1
    conda install anaconda::numpy=2.8.4
    conda install conda-forge::multiprocess=6.3.2
    conda install conda-forge::tqdm=6.3.2
    conda install conda-forge::regex=2022.04.01
-------------------------------------------------------------------------------
'''

if len(sys.argv) != 6:
    print(syntax)
    sys.exit()


base_name = str(sys.argv[1])
treshold_cutoff = sys.argv[2]
exponent_treshold = sys.argv[3]
exponent_degree = sys.argv[4]
threads = int(sys.argv[5])

base_name_with_settings = base_name + '_tr-cu=' + \
    str(treshold_cutoff) + '_ex-th=' + \
    str(exponent_treshold) + '_ex-de=' + exponent_degree


def create_folder_if_not_exists(directory):
    """
    Creates folder if not exists.

    :param directory: name of dictionary.
    """
    try:
        os.makedirs(directory, exist_ok=True)
    except OSError as e:
        print(f"Error creating directory {directory}: {e}")


def remove_item(path):
    """
    Deletes file or folder in path.

    :param path: Path to file oe path.
    """
    if os.path.exists(path):
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)


def format_time_duration(start, end):
    """
    Counts execution time.

    :param start: Start time moment.
    :param end: End time moment.
    :return: Execution time.
    """
    hours, rem = divmod(end - start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:02}:{:02}:{:05.2f}".format(int(hours), int(minutes), seconds)


def create_table_filename(filename):
    """
    Generates a new filename for storing table data based on an existing filename.

    This function takes a given filename, removes its file extension, and appends '_table.tsv' 
    to create a new filename suitable for a tab-separated values (TSV) file. This is typically 
    used to create a distinct filename for a table derived from the original data.

    :param filename: The original filename from which to generate the table filename.
    :return: A string representing the new filename with a '.tsv' extension.
    """
    return f"{filename.split('.')[0]}_table.tsv"


def prepare_mapping():
    """
    Prepares the mapping of Protein Data Bank (PDB) files for further analysis.

    The function performs the following tasks:
    1. Ensures the 'input' folder exists for storing PDB files.
    2. Removes any previous mapping data to prevent conflicts.
    3. Initiates the mapping process for each PDB file in the 'input' directory.
    4. Stores sequences and mapping keys in respective files.
    5. Saves random positions based on the number of sequences processed.
    """

    create_folder_if_not_exists('input')
    remove_item('pdb_mapping')
    remove_item('tmalign')
    create_folder_if_not_exists('pdb_mapping')

    pdb_fasta = []
    mapping_key = ['name\tkey\n']
    counter = 0

    for pdb_file in glob.glob('./input/*.pdb'):
        process_pdb_file(pdb_file, counter, pdb_fasta, mapping_key)
        counter += 1

    save_to_file(base_name_with_settings + '_sequences.txt', pdb_fasta)
    save_to_file(base_name_with_settings + '_mapping_key.txt', mapping_key)
    save_random_positions(len(pdb_fasta))


def process_pdb_file(pdb_file, counter, pdb_fasta, mapping_key):
    """
    Processes a single PDB file to extract relevant data for mapping.

    The function performs the following steps:
    1. Reads the PDB file and processes each record (chain) using the SeqIO parser.
    2. Copies the PDB file to a new directory, assigning a unique filename based on a counter to avoid naming conflicts.
    3. Appends the first chain's processed data to the pdb_fasta list if the file contains any chains, for further analysis.

    :param pdb_file: Path to the PDB file to be processed.
    :param counter: An incrementing number used as a unique identifier for each PDB file processed.
    :param pdb_fasta: A list that accumulates sequences from each processed PDB file.
    :param mapping_key: A list that maintains a mapping of original file names to their corresponding unique identifiers.
    """
    with open(pdb_file, 'r') as file:
        pdb_chains = [process_record(record, pdb_file, counter, mapping_key)
                      for record in SeqIO.parse(file, 'pdb-atom')]

    new_name_pdb_with_path = f'./pdb_mapping/{counter}.pdb'
    shutil.copy(pdb_file, new_name_pdb_with_path)
    if pdb_chains:
        pdb_fasta.append(pdb_chains[0])


def process_record(record, pdb_file, counter, mapping_key):
    """
    Processes a single record from a PDB file to extract the sequence and update mapping keys.

    :param record: A single record (chain) from the PDB file being parsed.
    :param pdb_file: Path to the PDB file, used to extract the base name.
    :param counter: A numeric identifier that provides a unique index for each sequence processed.
    :param mapping_key: A list to which the mapping between the base name and counter is appended.
    :return: A string representing the sequence of the record in FASTA format.
    """
    name = pdb_file.split('/')[-1].split('.pdb')[0]
    sequence = f">{name}[{counter}]\n{record.seq}\n"
    mapping_key.append(f"{name}\t{counter}\n")
    return sequence


def save_to_file(filename, data):
    """
    Saves the provided data to a file.

    This function opens a file with the specified filename and writes the data to it. 
    The data should be a list of strings, where each string will be written to a new line in the file.

    :param filename: The name of the file where the data will be saved.
    :param data: A list of strings to be written to the file.
    """
    with open(filename, 'w') as file:
        file.writelines(data)


def save_random_positions(count):
    """
    Generates a specified number of random positions and saves them to a file.

    This function writes a list of random 3D positions to a file, where each position is 
    represented by three coordinates, each a floating-point number between 0 and 1. 
    The count parameter specifies the number of random positions to generate. Each position is 
    saved in the format "index x-coordinate y-coordinate z-coordinate" on a new line.

    :param count: The number of random positions to generate and save.
    """
    with open(base_name_with_settings + '_random_positions.txt', 'w') as file:
        for i in range(count):
            line = f"{i} {random.uniform(0, 1):.3f} {random.uniform(0, 1):.3f} {random.uniform(0, 1):.3f}\n"
            file.write(line)


def TMalign_run(arg):
    """
    Runs the TM-align program to compare a specified protein structure file against all protein 
    structures in the pdb_mapping directory, and aggregates the alignment results.

    This function takes a path to a PDB file, extracts the filename, and uses it to create specific 
    directories for storing individual and aggregated alignment results. It then performs structural 
    alignment of the specified protein with each protein structure in the pdb_mapping directory using 
    TM-align, saves each result in a dedicated directory, and finally aggregates all results into a 
    single file.

    :param arg: Path to the PDB file which will be used as the reference in the TM-align comparisons.
    """
    full_name_split = arg.rsplit('/')
    file_name = full_name_split[-1]
    file_name_clear = file_name.split('.')
    create_folder_if_not_exists('./tmalign/%s/' % (file_name_clear[0]))
    create_folder_if_not_exists('./tmalign/all/')

    for temp in glob.glob('./pdb_mapping/*.pdb'):
        template_split = temp.rsplit('/')
        template_name = template_split[-1]
        temp_name_clear = template_name.split('.')

        out = str(file_name_clear[0] + '_' + temp_name_clear[0] + '.aln')

        os.system('TMalign %s %s -a > ./tmalign/%s/%s' %
                  (arg, temp, file_name_clear[0], out))
    os.system('cat ./tmalign/%s/*.aln > ./tmalign/all/%s_all_comparisons.aln' %
              (file_name_clear[0], file_name_clear[0]))


def prepare_tmalign_results():
    """
    Processes and aggregates TM-align comparison results from all alignment files into TSV format.

    This function reads through each alignment result file in the 'tmalign/all' directory, extracts
    pertinent comparison data, and formats this into a tab-separated value (TSV) file. Each TSV file
    is named according to the original alignment file and saved in the 'tmalign/all_tab/' directory.
    Data includes details like query and template names, alignment scores, RMSD, and sequence identity
    ratios, among other metrics.

    The function handles parsing of complex output formats from TM-align and ensures all relevant data
    is captured and accurately represented for further analysis.
    """
    create_folder_if_not_exists('./tmalign/all_tab/')
    # create_folder('./tmalign/all_tab_sorted/')

    for file in glob.glob('./tmalign/all/*.aln'):
        """
        Transforms TM-score values above a specified threshold.

        :param value: The logarithmic TM-score value.
        :param threshold: Threshold above which values are transformed.
        :return: The transformed TM-score value.
        """

        file_strip_1 = file.strip('./tmalign/all/')
        file_strip_2 = file_strip_1.strip('.aln')
        new_name_aln_tsv = file_strip_2 + '.tsv'

        with open(file, 'r') as Rs:
            with open('./tmalign/all_tab/' + new_name_aln_tsv, 'w') as table:

                new_table = csv.writer(table, delimiter='\t')

                data = [['query', 'template', 'qlen', 'tlen', 'aligned-len', 'RMSD', 'n_identical/n_aligned', 'TM-score-ch1',
                         'TM-score-ch2', 'TM-score-norm-avg-len', 'avg-len-ch', 'TM-score chosen by shorter chain', 'TM-score-ln']]

                results_read = str(Rs.read())

                # results = results_read.split(' *                        TM-align (Version 20190822)                     *')
                results = re.split(
                    ' \*                        TM-align \(Version \d+\)                     \*', results_read)

                for r in results[1:]:
                    query_strip = file.strip('./tmalign/all/')
                    query_split = query_strip.split('_all_comparisons.aln')
                    query = query_split[0]
                    # print(r)
                    # print(query, '################')
                    if re.search('Length of Chain_\d: +\d+', r):
                        qlen_regex = re.search(
                            'Length of Chain_\d: +\d+', r).group(0)
                        qlen_split = qlen_regex.split('Length of Chain_1:')
                        qlen = qlen_split[-1].strip(' ')
                    else:
                        qlen = 'error'
                        print(r, 'error_temp_qlen')
                        break
                    if re.search('Name of Chain_2: \./pdb_mapping/.+\.pdb', r):
                        temp_regex = re.search(
                            'Name of Chain_2: \./pdb_mapping/.+\.pdb', r).group(0)

                        temp_del_prefix = temp_regex.split('./pdb_mapping/')
                        # print(temp_del_prefix)
                        temp_del_prefix2 = temp_del_prefix[1].split(
                            'Name of Chain_2: ')
                        # print(temp_del_prefix2)
                        temp_split = temp_del_prefix2[0].split('.')
                        temp = temp_split[0]
                        # print(temp)
                    else:
                        temp = 'error'
                        # print(r,'error_temp')

                    if re.search('(TM-score= +\d+\.\d+ \(if normalized by average length of chains = +\d+\.\d+\))|(TM-score= +\d+ \(if normalized by average length of chains = +\d+\))', r):
                        score3_regex = re.search(
                            '(TM-score= +\d+\.\d+ \(if normalized by average length of chains = +\d+\.\d+\))|(TM-score= +\d+ \(if normalized by average length of chains = +\d+\))', r).group(0)
                        score3_split = score3_regex.split(' ')
                        score3 = score3_split[1]
                    else:
                        score3 = 'error'
                    if re.search('Length of Chain_\d: +\d+', r):
                        qlen_regex = re.search(
                            'Length of Chain_\d: +\d+', r).group(0)
                        qlen_split = qlen_regex.split('Length of Chain_1:')
                        qlen = qlen_split[-1].strip(' ')
                    else:
                        qlen = 'error'
                        print(r, 'error_temp_qlen')
                        break
                    if re.search('Length of Chain_2: +\d+', r):
                        tlen_regex = re.search(
                            'Length of Chain_2: +\d+', r).group(0)
                        tlen_split = tlen_regex.split('Length of Chain_2:')
                        tlen = tlen_split[-1].strip(' ')
                    else:
                        tlen = 'error'
                        print(r,  'error_tlen')
                        break
                    if re.search('(TM-score= \d+\.\d+ \(if normalized by length of Chain_1\))|(TM-score= \d+ \(if normalized by length of Chain_1\))', r):
                        score1_regex = re.search(
                            '(TM-score= \d+\.\d+ \(if normalized by length of Chain_1\))|(TM-score= \d+ \(if normalized by length of Chain_1\))', r).group(0)
                        score1_strip1 = score1_regex.strip('TM-score= ')
                        score1_split = score1_strip1.split(' ')
                        score1 = score1_split[0]
                    else:
                        score1 = 'error_TM1'
                        print(r, score1)
                        break
                    if re.search('(TM-score= \d+\.\d+ \(if normalized by length of Chain_2\))|(TM-score= \d+ \(if normalized by length of Chain_2\))', r):
                        score2_regex = re.search(
                            '(TM-score= \d+\.\d+ \(if normalized by length of Chain_2\))|(TM-score= \d+ \(if normalized by length of Chain_2\))', r).group(0)
                        score2_strip1 = score2_regex.strip('TM-score= ')
                        score2_split = score2_strip1.split(' ')
                        score2 = score2_split[0]
                    else:
                        score2 = 'error_TM2'
                        print(r, score2)
                        break
                    if re.search('Aligned length= +\d+', r):
                        aligned_regex = re.search(
                            'Aligned length= +\d+', r).group(0)
                        aligned_strip = aligned_regex.strip('Aligned length=')
                        aligned = aligned_strip.strip(' ')
                    else:
                        aligned = 'error'
                        print(r, 'error_al')
                        break
                    if re.search('(TM-score= +\d+\.\d+ \(if normalized by average length of chains = +\d+\.\d+\))|(TM-score= +\d+ \(if normalized by average length of chains = +\d+\))', r):
                        score3_regex = re.search(
                            '(TM-score= +\d+\.\d+ \(if normalized by average length of chains = +\d+\.\d+\))|(TM-score= +\d+ \(if normalized by average length of chains = +\d+\))', r).group(0)
                        score3_split = score3_regex.split(' ')
                        score3 = score3_split[1]
                        score3avg_strip = score3_regex.split(' ')
                        score3avg = score3avg_strip[-1].strip(')')
                    else:
                        score3 = 'error'
                        print(r, 'error3ORAVG')
                        break
                    if re.search('(RMSD= +\d+\.\d+)|(RMSD=   \d+)', r):
                        RMSD_regex = re.search(
                            '(RMSD= +\d+\.\d+)|(RMSD=   \d+)', r).group(0)
                        RMSD_strip = RMSD_regex.strip('RMSD=')
                        RMSD = RMSD_strip.strip(' ')
                    else:
                        RMSD = 'error_rmsd'
                        print(r, RMSD)
                        break
                    if re.search('(Seq_ID=n_identical/n_aligned= \d+\.\d+)|(Seq_ID=n_identical/n_aligned= \d+)', r):
                        identical_regex = re.search(
                            '(Seq_ID=n_identical/n_aligned= \d+\.\d+)|(Seq_ID=n_identical/n_aligned= \d+)', r).group(0)
                        identical = identical_regex.strip(
                            'Seq_ID=n_identical/n_aligned= ')
                    else:
                        identical = 'error_ident'

                        print(r, 'error_ide')
                        break
                    
                    if float(qlen) == float(tlen):
                        query_all = query.split('_')
                        if query_all[0] == temp:
                            score_chosen = 1.0
                            # print(float(qlen), float(tlen))
                    elif float(qlen) > float(tlen):
                        score_chosen = score2
                    elif float(qlen) < float(tlen):
                        score_chosen = score1
                    
                    data.append([query, temp, qlen, tlen, aligned, RMSD, identical,
                                score1, score2, score3, score3avg, score_chosen])
                new_table.writerows(data)


def stick_all_resaults():
    """
    Aggregates all individual TSV files containing TM-align results into a single file and then filters 
    and saves these results based on specific criteria.

    This function consolidates all the TSV files in the './tmalign/all_tab/' directory into a single TSV
    file. It then reads this aggregated file to filter results based on whether the query and template
    are the same, if the TM-score is below a specified threshold, or if the TM-score is zero. The filtered
    results are then saved into another TSV file, structured to only include the query, template, and the
    TM-score.

    :param base_name_with_settings: Base name to prepend to results files for identification.
    :param treshold_cutoff: Numerical threshold below which TM-scores are considered too low to include.
    """
    with open(base_name_with_settings + '_all_results_tmalign.tsv', 'w') as results:

        new_table = csv.writer(results, delimiter='\t')

        data = [['query', 'template', 'qlen', 'tlen', 'aligned-len', 'RMSD', 'n_identical/n_aligned',
                 'TM-score-ch1', 'TM-score-ch2', 'TM-score-norm-avg-len', 'avg-len-ch', 'TM-score chosen by shorter chain']]
        for file in glob.glob('./tmalign/all_tab/*.tsv'):
            # open file in read mode
            with open(file, 'r') as read_obj:

                # pass the file object to reader() to get the reader object
                csv_reader = csv.reader(read_obj, delimiter='\t')
                header = next(csv_reader)
                # Check file as empty
                if header != None:
                    # Iterate over each row after the header in the csv
                    for row in csv_reader:
                        # row variable is a list that represents a row in csv
                        data.append(row)
        new_table.writerows(data)

    with open(base_name_with_settings + '_all_results_tmalign.tsv', 'r') as results:
        table = results.readlines()
        filtered_lines = ['query\ttemplate\tTM-score\n']

        for line in table[1:]:
            # print(line)
            splited_line = line.split('\t')
            # print(splited_line)
            first_split = splited_line[0].split('_')
            # print(first_split)
            first = first_split[0]
            score = splited_line[-1]
            if first == splited_line[1] or splited_line[-1] == 0.0 or float(score) < float(treshold_cutoff):
                continue
            else:
                # score = math.log(float(splited_line[2]))*(-1)
                new_line_to_clans = str(first) + '\t' + \
                    str(splited_line[1]) + '\t' + str(score)
                # print(new_line_to_clans)
                filtered_lines.append(new_line_to_clans)
    with open(base_name_with_settings + '_all_results_tmalign_filtered.tsv', 'w') as save_pairs:
        for e in filtered_lines:
            save_pairs.write(e)


def transform_high_values(value, threshold, exponent_degree):
    """
    Transforms TM-score values above a specified threshold.

    :param value: The logarithmic TM-score value.
    :param threshold: Threshold above which values are transformed.
    :param exponent_degree: A selector for the exponent used in the transformation:
                            1 -> exponent = 2
                            2 -> exponent = 12
                            3 -> exponent = 24
    :return: The transformed TM-score value.
    """
    # Map exponent_degree to the actual exponent value
    exponent_map = {
        1: 2,
        2: 12,
        3: 24
    }
    # Default to 12 if not specified correctly
    exponent = exponent_map.get(exponent_degree, 12)

    if value > threshold:
        return ((np.log(value) + 0.01) * -1) ** exponent
    else:
        return (np.log(value) + 0.01) * -1


def min_max_normalize_scores():
    """
    Normalizes TM-score values by applying a logarithmic transformation and an optional further transformation 
    for high values. The resulting transformed scores are saved to a new TSV file.

    This function reads TM-score data from a filtered results file, ensures the TM-score column is numeric and 
    contains no NaN values, and applies a logarithmic transformation to the scores. An additional transformation 
    function can be applied to scores above a specified threshold to further modify those values based on the 
    given parameters. The processed data is then saved to a new file.
    """
    input_file_path = base_name_with_settings + '_all_results_tmalign_filtered.tsv'
    output_file_path = base_name_with_settings + '_all_results_tmalign_filtered_transformed.tsv'
    input_data = pd.read_csv(input_file_path, sep='\t')
    input_data['TM-score'] = pd.to_numeric(
        input_data['TM-score'], errors='coerce')
    input_data = input_data.dropna(subset=['TM-score'])
    input_data['Log_Transformed_TM-score'] = (input_data['TM-score'].apply(
        transform_high_values, args=(float(exponent_treshold), float(exponent_degree))))
    input_data.to_csv(output_file_path, index=False, sep='\t')


def preparing_scoring():
    """
    Reads a TSV file containing transformed TM-score values and formats them into a specific string format 
    for further processing or visualization.

    This function reads through a file containing transformed TM-scores, extracts the query, template, and score
    from each line, and formats them into a string of the form 'query template:score'. These formatted strings 
    are then written to a new file intended for mapping or additional processing steps.
    """
    with open(base_name_with_settings + '_all_results_tmalign_filtered_transformed.tsv', 'r') as results:
        table = results.readlines()
        filtered_lines = []

        for line in table[1:]:
            # print(line)
            splited_line = line.split('\t')
            query = splited_line[0].strip()
            template = splited_line[1].strip()
            score = float(splited_line[3].strip())

            new_line_to_clans = str(query) + ' ' + \
                str(template) + ':' + str(score) + '\n'
            # print(new_line_to_clans)
            filtered_lines.append(new_line_to_clans)
    with open(base_name_with_settings + '_mapping.txt', 'w') as save_pairs:
        for e in filtered_lines:
            save_pairs.write(e)


def prepare_clans_file():
    """
    Prepares a file formatted for use with the CLANS software by combining sequence data, positional information,
    and interaction mappings into a single file.

    This function aggregates data from several files: sequences, random positions, and mappings. It reads sequences
    from a fasta file, position data from a text file, and mapping information from another text file. Each section is
    appropriately formatted and compiled into a single file that includes the number of sequences, the sequences themselves,
    their positions, and their interactions in a format compatible with CLANS.
    """
    with open(base_name_with_settings + '_sequences.txt', 'r') as fasta:
        sequences = list(SeqIO.parse(fasta, 'fasta'))
        number_of_seq = len(sequences)
        first_line = 'sequences=' + str(number_of_seq) + '\n'
    with open(base_name_with_settings + '_sequences.txt', 'r') as seq:
        sequences_read = seq.read()
        sequences_section = '<seq>\n' + sequences_read + '</seq>\n'
    with open(base_name_with_settings + '_random_positions.txt', 'r') as pos:
        position_read = pos.read()
        position_section = '<pos>\n' + position_read + '</pos>\n'
    with open(base_name_with_settings + '_mapping.txt', 'r') as mapping:
        mapping_read = mapping.read()
        mapping_section = '<hsp>\n' + mapping_read + '</hsp>\n'
    with open(base_name_with_settings + '_structure-network.clans', 'w') as clans_file_write:
        clans_file = first_line + sequences_section + position_section + mapping_section
        clans_file_write.write(clans_file)


def unmap():
    """
    Removes the '_all_comparisons' suffix from the filenames in the TM-align results file and remaps numerical identifiers
    back to their original names based on a mapping key file.

    The function first reads and modifies the results file to eliminate the '_all_comparisons' suffix from the identifiers.
    Then, it loads these results into a DataFrame. It also reads a key mapping file into another DataFrame and converts this 
    to a dictionary for efficient look-up. Using this dictionary, the function replaces numeric identifiers in the 'query' 
    and 'template' columns of the DataFrame with their corresponding original names. Finally, the modified DataFrame is 
    saved back to a TSV file without the numerical identifiers.

    """

    with open(base_name_with_settings + '_all_results_tmalign.tsv', 'r') as tmresults:
        tmread = tmresults.read()
        tmread_replace = tmread.replace('_all_comparisons', '')
        with open(base_name_with_settings + '_all_results_tmalign.tsv', 'w') as tmresults_write:
            tmresults_write.write(tmread_replace)

    # Create a DataFrame from the data
    df_name = base_name_with_settings + '_all_results_tmalign.tsv'
    df = pd.read_csv(df_name, sep='\t', index_col=False)

    # Mapping of query and template numbers to their respective names
    name_key_mapping = pd.read_csv(
        base_name_with_settings + '_mapping_key.txt', sep='\t', index_col=False)

    # Convert the 'name_key_mapping' DataFrame to a dictionary for easy lookup
    mapping_dict = name_key_mapping.set_index('key')['name'].to_dict()

    # Replace values in the 'query' column using the mapping
    df['query'] = df['query'].map(mapping_dict)
    df['template'] = df['template'].map(mapping_dict)

    df_out = base_name_with_settings + '_all_results_tmalign_unmaped.tsv'
    df.to_csv(df_out, sep="\t", index=False)


def main():
    print('\nArguments:')
    print('name: ', sys.argv[1])
    print('treshold_cutoff: ', sys.argv[2])
    print('exponent_treshold: ', sys.argv[3])
    print('exponent_degree: ', sys.argv[4])
    print('threads: ', sys.argv[5], '\n')
    warnings.filterwarnings('ignore', category=BiopythonWarning)
    start = time.time()
    print('Mapping...')
    prepare_mapping()
    processes = glob.glob('./pdb_mapping/*.pdb')
    print('Running comparisons...')

    pool = Pool(processes=threads)

    for _ in tqdm(pool.imap_unordered(TMalign_run, processes), total=len(processes), desc="Aligning PDB files"):
        pass

    print('Preparing TMalign resaults...')
    prepare_tmalign_results()
    stick_all_resaults()
    print('Preparing scoring...')
    min_max_normalize_scores()
    preparing_scoring()
    print('Preparing clans file...')
    prepare_clans_file()
    print('Unmaping...')
    unmap()
    finish = time.time()
    time_measure = format_time_duration(start, finish)
    print(f'Done! The whole procedure took {time_measure}\n')

    remove_item(base_name_with_settings + '_random_positions.txt')
    remove_item(base_name_with_settings + '_mapping_key.txt')
    remove_item(base_name_with_settings + '_mapping.txt')
    remove_item(base_name_with_settings + '_sequences.txt')
    remove_item(base_name_with_settings + '_all_results_tmalign_filtered.tsv')
if __name__ == '__main__':
    main()
