"""
Used to generate count tables for sgRNA in FASTQ file
Date: 05/09/2024
"""
import pandas as pd
import os
from collections import Counter
from Bio import SeqIO
import re

# Function to find matching substrings and append corresponding values
def find_matching_sgRNA(row, sgRNA_df):
    matches = sgRNA_df[sgRNA_df['seq'].str.contains(row['seq'], regex=False)]
    if not matches.empty:
        return pd.Series([matches.iloc[0]['sgRNA_name'], matches.iloc[0]['lib']])
    else:
        return pd.Series([None, None])


# main program
# Read expected sgRNA sequences
dir =r'C:\Users\neel\Box\cbe-neel\cbe-neel-shared\JointProjects\SingleCellGlycomics\CRISPR_sgRNA' +'\\'
fname = r'glycoCRISPR_subpoolLibrary_order.xlsx'
sheet = 'glycoCRISPR_subpoolLibrary_orde'
sgRNA_df = pd.read_excel(dir+fname, sheet_name=sheet, usecols=[1, 4, 12], header=None,
                   names=['lib', 'seq', 'sgRNA_name'])
sgRNA_df = sgRNA_df.groupby(['seq', 'sgRNA_name'], as_index=False).agg({'lib': list})

# read Fastq files
dataDir =r'C:\Users\neel\Box\cbe-neel\cbe-neel-shared\DATA FOLDER\NGS data\RPCI_050724-RQ025534\RQ025534-Neelamegham\Arun'
fastq_files = [file for file in os.listdir(dataDir) if 'R1' in file and file.endswith('.fastq')]
for file in fastq_files:
    full_fname = dataDir + '\\' + file
    fq_dict = SeqIO.index(full_fname, "fastq")  # used to index entries based on key-values
    keys1 = list(fq_dict.keys())
    allSeq=[]
#    QC = []
    substring = "AAACACC"   # string just before the sgRNA sequence
    for x in keys1:
        str_temp = str(fq_dict[x].seq)
        match = re.search(substring, str_temp)
        if match:
            start_pos = match.end() + 2
            str_temp = str_temp[start_pos:start_pos + 17]   # read sgRNA sequence and append into allSeq lits
            allSeq.append(str_temp)
    frequency_table = Counter(allSeq)       # make table and then write to pd dataframe
    result_df = pd.DataFrame.from_dict(frequency_table, orient='index').reset_index()
    result_df.columns = ['seq', 'count']
    result_df = result_df.sort_values(by='count', ascending=False)
    result_df = result_df.reset_index(drop=True)
    # add sgRNA name to files
    result_df[['sgRNA_name', 'lib']] = result_df.apply(find_matching_sgRNA, args=(sgRNA_df,), axis=1)
    #        QC.append(str(fq_dict[x].letter_annotations))  # read the sequence and convert into string, then place into list

    # Construct the full file path for the output Excel file
    output_file_path = dataDir + '\\' +'output.xlsx'

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

a=1
