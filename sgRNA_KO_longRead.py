"""
Used to generate count tables for sgRNA in FASTQ file, using oxford nanopore sequencing results
Date: 12/15/2024
"""
import pandas as pd
import os
from collections import Counter
from Bio import SeqIO
import re


def extract_target(sequence, search_string):
    # Find the position of search_string in the sequence
    pos = sequence.find(search_string)

    if pos != -1:
        # Extract 25 bases before the search_string and 20 bases after
        start = max(pos - 15, 0)  # Ensure the start is not negative
        end = pos + len(search_string)

        # Return the substring
        return sequence[start:end]
    else:
        # If search_string is not found, return NaN or an empty string
        return None



# input data
dataDir = r'C:\Users\neel\Box\cbe-neel\cbe-neel-shared\WeeklyMeetings\Arun\Rosewell collaborations\SLC35A1 Knockouts\S' \
          r'equencing Results\79L5RC_results_FG COLO 357 and ASPC-1\79L5RC_fastq'
ProductSeq ='TCGTCGGCAGCGTCAGTAATGTCTTTGTTGCACGTATTTTCCAGACAATGTCACTTTATTATTCAAGTTATACTGCTTGGCAGTGATGACCCTGATGGCT' \
           'GCAGTCTATACCATAGCTTTAAGATACACAAGGACATCAGACAAAGAACTCTACTTTTCAACCACAGCCGTGTGTATCACAGAAGTTATAAAGTTATTGCTAA' \
           'GTGTGGGAATTTTAGCCCGAGCCCACGAGAC'
sgRNA = 'CAGCCGTGTGTATCACAGAA'
plus_minus = 30
output_file_path = dataDir + '\\' + 'output.xlsx'

position = ProductSeq.find(sgRNA)
if position != -1:
    start_pos1 = max(position - plus_minus, 0)  # Ensure the start position doesn't go below 0
    substring1 = ProductSeq[start_pos1:start_pos1 + 10]

    end_pos2 = position + len(sgRNA) + plus_minus  # plus_minus bases after the end of sgRNA
    substring2 = ProductSeq[end_pos2-10:end_pos2]

fastq_files = [file for file in os.listdir(dataDir) if file.endswith('.fastq')]
for file in fastq_files:
    full_fname = dataDir + '\\' + file
    fq_dict = SeqIO.index(full_fname, "fastq")  # used to index entries based on key-values
    keys1 = list(fq_dict.keys())
    allSeq=[]
    for x in keys1:
        str_temp = str(fq_dict[x].seq)
        if substring1 in str_temp or substring2 in str_temp:
            allSeq.append(str_temp)  # Add to list if match found
    frequency_table = Counter(allSeq)  # make table and then write to pd dataframe
    result_df = pd.DataFrame.from_dict(frequency_table, orient='index').reset_index()
    result_df.columns = ['seq', 'count']
    result_df = result_df.sort_values(by='count', ascending=False)
    result_df = result_df.reset_index(drop=True)
    result_df['target'] = result_df['seq'].apply(lambda x: extract_target(x, sgRNA[-8:]))

    # Check if the file exists
    if os.path.isfile(output_file_path):
        # Open the existing Excel file
        with pd.ExcelWriter(output_file_path, mode='a') as writer:
            # Try to write to the existing sheet, overwrite if it exists
            try:
                result_df.to_excel(writer, sheet_name=file, index=False)
            except ValueError:
                # If the sheet already exists, remove it and then write again
                with pd.ExcelWriter(output_file_path, mode='w') as writer:
                    result_df.to_excel(writer, sheet_name=file, index=False)
    else:
        # If the file does not exist, create it and write to a new sheet
        with pd.ExcelWriter(output_file_path, mode='w') as writer:
            result_df.to_excel(writer, sheet_name=file, index=False)
