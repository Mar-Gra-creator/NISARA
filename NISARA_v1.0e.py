#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:54:08 20223

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

###############################################################################
#                          Settings and instructions                          #
###############################################################################
syntax = '''-------------------------------------------------------------------------------
30 July 2024 v1.00
Originally developed in Python 3.10.14
Autors: MG, MK

 _____ _____ _____ _____ _____ _____ 
|   | |     |   __|  _  | __  |  _  |
| | | |-   -|__   |     |    -|     |
|_|___|_____|_____|__|__|__|__|__|__| v1.0

Network Analysis of Structure Similarity Alignment V1.0

Syntax:
python NISARA_v1.0e.py [base_name] [threshold_cutoff] 
[exponent_threshold] [exponent_degree] [threads]

:param base_name: Base name for the analysis.

:param threshold_cutoff:    Cut-off threshold above which TM-score values are
                            retained and considered for connections (th-cu)

:param exponent_threshold:  The threshold above which values are elevated 
                            according to the exponent_degree (ex-th)

:param exponent_degree:     A selector for the exponent used in 
                            the transformation
                            (ex-de)
                            0 -> exponent = raw 
                            1 -> exponent = 2
                            2 -> exponent = 12
                            

:param threads: Number of threads used for parallel TM-align comparisons.

Transformation Logic:
    TM-score values are logarithmized, incremented by 0.01 to avoid zeros, and 
    multiplied by -1 to convert to positive values.
    If the TM-score is higher than the specified threshold, it is additionally
    raised to a power specified by exponent_degree.

    def transform_value(value, threshold, exponent):
        if value > threshold:
            return ((e.g. log(value) + 0.01) * -1) ** exponent
        else:
            return (e.g. log(value) + 0.01) * -1

    value:      TM-score
    threshold:  The threshold above which values are amplified according to the 
                exponent_degree
    exponent:   The exponent value used in the transformation

How to Use:
    1. Prepare a folder named 'input' containing pdb files.
    2. Execute the command: python NISARA_v1.0e.py name_x 0.5 0.6 2 24

Input:
    *.pdb in prepared folder named 'input'
    Try to keep the names short and uncomplicated.
    
Output:
    -> *_th-cu=*_ex-th=*_ex-de=*_all_results_tmalign_unmaped.tsv
       Full table of TM-align comparisons of protein pairs. Unmapped names
    -> *_th-cu=*_ex-th=*_ex-de=*_all_results_tmalign_filtered_transformed.tsv
       Table of TM-align comparisons of protein pairs
       after cutoff and transformation
    -> *_th-cu=*_ex-th=*_ex-de=*_structure-network.clans
       CLANS file
    

Limitations:
    Ensure that structures are in a single piece, not fragmented. 
    (If multiple breaks in the chain, TM-align might align only part of it.)


-------------------------------------------------------------------------------
'''

###############################################################################
#                   Argument Validation and Main Variables                    #
###############################################################################
if len(sys.argv) != 6:
    print(syntax)
    sys.exit()

base_name = str(sys.argv[1])            # e.g.  "test"
treshold_cutoff = sys.argv[2]           # e.g.  "0.4"
exponent_treshold = sys.argv[3]         # e.g.  "0.7"
exponent_degree = sys.argv[4]           # e.g.  "2"
threads = int(sys.argv[5])              # e.g.   8

base_name_with_settings = (
    f"{base_name}_tr-cu={treshold_cutoff}_ex-th={exponent_treshold}_ex-de={exponent_degree}"
)


###############################################################################
#                          Auxiliary functions                                #
###############################################################################
def create_folder_if_not_exists(directory):
    try:
        os.makedirs(directory, exist_ok=True)
    except OSError as e:
        print(f"Error creating directory {directory}: {e}")


def remove_item(path):
    if os.path.exists(path):
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)


def format_time_duration(start, end):
    hours, rem = divmod(end - start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:02}:{:02}:{:05.2f}".format(int(hours), int(minutes), seconds)


def save_to_file(filename, data):
    with open(filename, 'w') as file:
        file.writelines(data)


###############################################################################
#                      Creating mapping and input files                       #
###############################################################################
def prepare_mapping():
    """
	1. Clears/creates folders.
	2. Moves pdb files from ./input to ./pdb_mapping, assigning them IDs (sequential numbers).
	3. Creates a sequence file (FASTA) and a mapping key file.
	4. Creates random 3D positions for files for visualization in CLANS.
    """
    create_folder_if_not_exists('input')
    remove_item('pdb_mapping')
    remove_item('tmalign')
    create_folder_if_not_exists('pdb_mapping')

    pdb_fasta = []
    mapping_key = ['name\tkey\n']
    counter = 0

    for pdb_file in glob.glob('./input/*.pdb'):
        with open(pdb_file, 'r') as file:
            # Odczyt sekwencji z PDB
            pdb_chains = [
                process_record(record, pdb_file, counter, mapping_key)
                for record in SeqIO.parse(file, 'pdb-atom')
            ]
        new_name_pdb_with_path = f'./pdb_mapping/{counter}.pdb'
        shutil.copy(pdb_file, new_name_pdb_with_path)

        # Bierzemy 1-szy łańcuch do FASTA (zwykle wystarczające)
        if pdb_chains:
            pdb_fasta.append(pdb_chains[0])
        counter += 1

    save_to_file(base_name_with_settings + '_sequences.txt', pdb_fasta)
    save_to_file(base_name_with_settings + '_mapping_key.txt', mapping_key)
    save_random_positions(len(pdb_fasta))


def process_record(record, pdb_file, counter, mapping_key):
    """
    Creates a FASTA entry from the PDB string and adds the following line to mapping_key: “file name” -> “number”
    """
    name = pdb_file.split('/')[-1].split('.pdb')[0]
    sequence = f">{name}[{counter}]\n{record.seq}\n"
    mapping_key.append(f"{name}\t{counter}\n")
    return sequence


def save_random_positions(count):
    """
    Generates random 3D coordinates (x, y, z) in the range [0,1], saves to file
    """
    with open(base_name_with_settings + '_random_positions.txt', 'w') as file:
        for i in range(count):
            line = f"{i} {random.uniform(0, 1):.3f} {random.uniform(0, 1):.3f} {random.uniform(0, 1):.3f}\n"
            file.write(line)


###############################################################################
#                          Starting TM-align                                  #
###############################################################################
def TMalign_run(arg):
    """
    Runs TMalign for the arg file (e.g., “0.pdb”) vs. all in ./pdb_mapping/*.pdb.
    The results (single .aln) go to ./tmalign/<id>/.
    Then they are combined in ./tmalign/all/<id>_all_comparisons.aln.
    """
    full_name_split = arg.rsplit('/')
    file_name = full_name_split[-1]
    file_name_clear = file_name.split('.')  # e.g.  ["0", "pdb"]

    create_folder_if_not_exists(f'./tmalign/{file_name_clear[0]}/')
    create_folder_if_not_exists('./tmalign/all/')

    for temp in glob.glob('./pdb_mapping/*.pdb'):
        temp_split = temp.rsplit('/', 1)
        temp_name_clear = temp_split[-1].split('.')  # e.g.  ["1", "pdb"]
        out = f"{file_name_clear[0]}_{temp_name_clear[0]}.aln"

        os.system(f"TMalign {arg} {temp} > ./tmalign/{file_name_clear[0]}/{out}")

    # We combine all matches into one file
    os.system(
        f"cat ./tmalign/{file_name_clear[0]}/*.aln > ./tmalign/all/{file_name_clear[0]}_all_comparisons.aln"
    )


###############################################################################
#                          Parsing TM-align results                        #
###############################################################################
def prepare_tmalign_results():
    """
    Loads .aln files from ./tmalign/all/.
    Each file may contain multiple blocks (because cat merges alignments).
    We search for consecutive blocks starting with “Name of Chain_1: ...”.
    We extract 10 values from each block:
      [query, template, qlen, tlen, aligned-len, RMSD, n_identical/n_aligned, TM-score-ch1, TM-score-ch2, TM-score-chosen]
    We save them to a .tsv file in the ./tmalign/all_tab/ directory.
    """
    create_folder_if_not_exists('./tmalign/all_tab/')

    for file in glob.glob('./tmalign/all/*.aln'):
        # Output file name (tsv)
        file_strip_1 = file.replace('./tmalign/all/', '')
        file_strip_2 = file_strip_1.replace('.aln', '')
        new_tsv = f"{file_strip_2}.tsv"

        data_rows = [
            [
                "query",
                "template",
                "qlen",
                "tlen",
                "aligned-len",
                "RMSD",
                "n_identical/n_aligned",
                "TM-score-ch1",
                "TM-score-ch2",
                "TM-score chosen by shorter chain"
            ]
        ]

        # Load the entire .aln file into memory as a list of lines
        with open(file, 'r') as f:
            lines = [l.rstrip('\n') for l in f]

        # We set the query based on the file name: “0_all_comparisons.aln” -> “0”
        query_split = file_strip_2.split('_all_comparisons')
        query_name = query_split[0] if query_split else "unknown_query"

        # We will be block by block (each block starts with “Name of Chain_1:”)
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith("Name of Chain_1:"):
                # New block
                block_lines = []
                # We add this line to block_lines
                block_lines.append(line)
                i += 1
                # We load into block_lines until the next “Name of Chain_1:” or the end
                while i < len(lines) and not lines[i].startswith("Name of Chain_1:"):
                    block_lines.append(lines[i])
                    i += 1

                # We have a complete block in block_lines from which we will extract parameters.
                drow = parse_single_alignment_block(block_lines, query_name)
                if drow:
                    data_rows.append(drow)
            else:
                i += 1  # we continue searching for “Name of Chain_1:”

        # Write to TSV file
        with open(f"./tmalign/all_tab/{new_tsv}", 'w') as out_tsv:
            writer = csv.writer(out_tsv, delimiter='\t')
            writer.writerows(data_rows)


def parse_single_alignment_block(block_lines, query):
    """
    Accepts a list of lines from a single alignment block (starting with “Name of Chain_1:”)
    Returns a list of 10 values [query, template, qlen, tlen, aligned_len, RMSD, ident, score1, score2, chosen].
    If anything goes wrong -> returns None.
    """
    text = "\n".join(block_lines)

    # template
    # Name of Chain_2: ./pdb_mapping/33.pdb
    match_temp = re.search(r'Name of Chain_2:\s+\./pdb_mapping/(.+)\.pdb', text)
    if match_temp:
        template = match_temp.group(1)
    else:
        template = "error"

    # qlen (Length of Chain_1: 257)
    match_qlen = re.search(r'Length of Chain_1:\s+(\d+)', text)
    if match_qlen:
        qlen = match_qlen.group(1)
    else:
        qlen = "error"

    # tlen (Length of Chain_2: 201)
    match_tlen = re.search(r'Length of Chain_2:\s+(\d+)', text)
    if match_tlen:
        tlen = match_tlen.group(1)
    else:
        tlen = "error"

    # aligned-len  (Aligned length= 143)
    match_aln = re.search(r'Aligned length=\s+(\d+)', text)
    if match_aln:
        aligned_len = match_aln.group(1)
    else:
        aligned_len = "error"

    # RMSD (RMSD=   4.17)
    match_rmsd = re.search(r'RMSD=\s+([\d\.]+)', text)
    if match_rmsd:
        RMSD = match_rmsd.group(1)
    else:
        RMSD = "error_rmsd"

    # ident (Seq_ID=n_identical/n_aligned= 0.091)
    match_ident = re.search(r'Seq_ID=n_identical/n_aligned=\s+([\d\.]+)', text)
    if match_ident:
        ident = match_ident.group(1)
    else:
        ident = "error_ident"

    # TM-score-ch1 (TM-score= 0.41325 (if normalized by length of Chain_1, ...))
    match_score1 = re.search(
        r'TM-score=\s*(\d+(?:\.\d+)?)\s*\(if normalized by length of Chain_1',
        text
    )
    if match_score1:
        score1 = match_score1.group(1)
    else:
        score1 = "error_TM1"

    # TM-score-ch2
    match_score2 = re.search(
        r'TM-score=\s*(\d+(?:\.\d+)?)\s*\(if normalized by length of Chain_2',
        text
    )
    if match_score2:
        score2 = match_score2.group(1)
    else:
        score2 = "error_TM2"

    # We choose the TM score of the shorter chain
    try:
        qf = float(qlen)
        tf = float(tlen)
        s1 = float(score1) if score1 != 'error_TM1' else 0
        s2 = float(score2) if score2 != 'error_TM2' else 0

        if qf == tf:
            chosen_score = s1  # e.g.  we take s1 because the strings are equal
        elif qf > tf:
            chosen_score = s2
        else:
            chosen_score = s1
    except ValueError:
        chosen_score = 0

    row = [
        query,
        template,
        qlen,
        tlen,
        aligned_len,
        RMSD,
        ident,
        score1,
        score2,
        chosen_score
    ]
    return row


###############################################################################
#             Merges and filters the results into a single TSV file           #
###############################################################################
def stick_all_resaults():
    """
    1) Reads all .tsv files from ./tmalign/all_tab/
       -> merges them into one file: *_all_results_tmalign.tsv with 10 columns
    2) Filtering to *_all_results_tmalign_filtered.tsv (columns: query, template, TM-score),
       taking into account threshold_cutoff and removing self-hits.
    """
    out_tsv = base_name_with_settings + '_all_results_tmalign.tsv'

    with open(out_tsv, 'w') as f_out:
        writer = csv.writer(f_out, delimiter='\t')
        # Nagłówek 10 kolumn
        writer.writerow([
            'query',
            'template',
            'qlen',
            'tlen',
            'aligned-len',
            'RMSD',
            'n_identical/n_aligned',
            'TM-score-ch1',
            'TM-score-ch2',
            'TM-score chosen by shorter chain'
        ])

        for tsv_file in glob.glob('./tmalign/all_tab/*.tsv'):
            with open(tsv_file, 'r') as f_in:
                reader = csv.reader(f_in, delimiter='\t')
                next(reader)  # pomijamy nagłówek
                for row in reader:
                    writer.writerow(row)

    # Filtering to the _filtered file
    out_filtered = base_name_with_settings + '_all_results_tmalign_filtered.tsv'
    with open(out_tsv, 'r') as rf, open(out_filtered, 'w') as wf:
        lines = rf.readlines()
        wf.write("query\ttemplate\tTM-score\n") 

        for line in lines[1:]:  # we skip the header
            spl = line.strip().split('\t')
            # 0: query
            # 1: template
            # 9: chosen score
            if len(spl) < 10:
                continue
            query_raw = spl[0]
            template = spl[1]
            chosen_score_str = spl[9]

            # We exclude self-hits (e.g.  “14” == “14”) -> we check if query_raw.split(‘_’)[0] == template
            # or, more simply, query_raw == template.
            # It depends on whether the names are in the style of “14” or “14_something”.
            # Most often, it is enough to check:
            query_short = query_raw.split('_')[0]
            if query_short == template:
                continue

            try:
                chosen_val = float(chosen_score_str)
            except ValueError:
                continue

            if chosen_val <= 0:
                continue
            if chosen_val < float(treshold_cutoff):
                continue

            wf.write(f"{query_short}\t{template}\t{chosen_val}\n")


###############################################################################
#           Transformation of TM-score values and preparation of CLANS        #
###############################################################################
def transform_high_values(value, threshold, exponent_degree):
    if exponent_degree == 0:
        return value
    exponent_map = {1:2, 2:12}
    exponent = exponent_map.get(exponent_degree, 12)
    if value > threshold:
        return ((np.log(value)+0.01)*-1) ** exponent
    else:
        return (np.log(value)+0.01) * -1

def min_max_normalize_scores():
    inp  = base_name_with_settings + '_all_results_tmalign_filtered.tsv'
    outp = base_name_with_settings + '_all_results_tmalign_filtered_transformed.tsv'

    df = pd.read_csv(inp, sep='\t')
    df["TM-score"] = pd.to_numeric(df["TM-score"], errors='coerce')
    df.dropna(subset=["TM-score"], inplace=True)

    # 1) first, a radical transformation
    df["raw_transf"] = df["TM-score"].apply(
        transform_high_values,
        args=(float(exponent_treshold), int(exponent_degree))
    )

    # 2) now min–max normalization to [0,1]
    mn = df["raw_transf"].min()
    mx = df["raw_transf"].max()
    df["Log_Transformed_TM-score"] = (df["raw_transf"] - mn) / (mx - mn)

    # 3) record
    df.drop(columns=["raw_transf"], inplace=True)
    df.to_csv(outp, index=False, sep='\t')



def preparing_scoring():
    """
    At the entrance *_all_results_tmalign_filtered_transformed.tsv (4 columns):
      query, template, TM-score, Log_Transformed_TM-score
    We create a mapping.txt file with the following lines:
      query template:Log_Transformed_TM-score
    e.g.  “14 33:-1.2356”
    """
    inp = base_name_with_settings + '_all_results_tmalign_filtered_transformed.tsv'
    out = base_name_with_settings + '_mapping.txt'

    with open(inp, 'r') as f_in, open(out, 'w') as f_out:
        lines = f_in.readlines()
        for i, line in enumerate(lines):
            if i == 0:  
                continue
            spl = line.strip().split('\t')
            # [query, template, TM-score, Log_Transformed_TM-score]
            if len(spl) < 4:
                continue
            query, template, _, transformed = spl
            f_out.write(f"{query} {template}:{transformed}\n")


def prepare_clans_file():
    """
    Creates a *.clans file:
    sequences=<number_of_sequences>
    <seq>...</seq>
    <pos>...</pos>
    <hsp>...</hsp>
    """
    seq_file = base_name_with_settings + '_sequences.txt'
    pos_file = base_name_with_settings + '_random_positions.txt'
    map_file = base_name_with_settings + '_mapping.txt'

    with open(seq_file, 'r') as ff:
        sequences = list(SeqIO.parse(ff, 'fasta'))
        number_of_seq = len(sequences)

    with open(seq_file, 'r') as fseq:
        seq_data = fseq.read()
    with open(pos_file, 'r') as fpos:
        pos_data = fpos.read()
    with open(map_file, 'r') as fmap:
        map_data = fmap.read()

    first_line = f"sequences={number_of_seq}\n"
    sequences_section = "<seq>\n" + seq_data + "</seq>\n"
    position_section = "<pos>\n" + pos_data + "</pos>\n"
    mapping_section = "<hsp>\n" + map_data + "</hsp>\n"

    with open(base_name_with_settings + '_structure-network.clans', 'w') as out_clans:
        out_clans.write(first_line + sequences_section + position_section + mapping_section)


def unmap():
    """
    Replaces numbers (query, template) in the *_all_results_tmalign.tsv file with original names 
    (from the mapping_key file). Removes “_all_comparisons” if it appears in this file.
    """
    all_results = base_name_with_settings + '_all_results_tmalign.tsv'
    with open(all_results, 'r') as f:
        text = f.read()
    text = text.replace('_all_comparisons', '')
    with open(all_results, 'w') as f:
        f.write(text)

    df = pd.read_csv(all_results, sep='\t', index_col=False)

    mk_file = base_name_with_settings + '_mapping_key.txt'
    mk_df = pd.read_csv(mk_file, sep='\t', index_col=False)
    map_dict = mk_df.set_index('key')['name'].to_dict()

    df['query'] = df['query'].map(map_dict)
    df['template'] = df['template'].map(map_dict)

    df_out = base_name_with_settings + '_all_results_tmalign_unmaped.tsv'
    df.to_csv(df_out, sep="\t", index=False)


###############################################################################
#                                Main function                                #
###############################################################################
def main():
    print('\nArguments:')
    print('name: ', base_name)
    print('treshold_cutoff: ', treshold_cutoff)
    print('exponent_treshold: ', exponent_treshold)
    print('exponent_degree: ', exponent_degree)
    print('threads: ', threads, '\n')

    warnings.filterwarnings('ignore', category=BiopythonWarning)

    start = time.time()
    print('Mapping...')
    prepare_mapping()

    processes = glob.glob('./pdb_mapping/*.pdb')
    print('Running comparisons...')
    pool = Pool(processes=threads)
    for _ in tqdm(pool.imap_unordered(TMalign_run, processes), 
                  total=len(processes), desc="Aligning PDB files"):
        pass

    print('Preparing TMalign results...')
    prepare_tmalign_results()
    stick_all_resaults()

    print('Preparing scoring...')
    min_max_normalize_scores()
    preparing_scoring()

    print('Preparing clans file...')
    prepare_clans_file()

    print('Unmapping...')
    unmap()

    finish = time.time()
    time_measure = format_time_duration(start, finish)
    print(f"Done! The whole procedure took {time_measure}\n")

    # We delete temporary files
    remove_item(base_name_with_settings + '_random_positions.txt')
    remove_item(base_name_with_settings + '_mapping_key.txt')
    remove_item(base_name_with_settings + '_mapping.txt')
    remove_item(base_name_with_settings + '_sequences.txt')
    remove_item(base_name_with_settings + '_all_results_tmalign_filtered.tsv')


if __name__ == '__main__':
    main()
