# This program contains three subroutines.
# main() sets up a series of BLAST searches, performing two identical searches for
# each query string, one -ungapped and one -gapped. Other parameters can be set 
# as desired.
# The length of the query strings are set with the {sample_length} variable.
# The {input_fasta_file} contains a set of sequences in FASTA format.
# In this version, a subset of the contigs used in Tomkins 2018 will be used.
# The sequence list is contained in the file 125contigs.csv.
# main() will choose a random starting location on each sequence and create a
# subsequence that is {sample_length} nucleotides long.
# These are sent to blast_sequence(). This subroutine performs the searches and
# saves the results in {blast_results_directory}.
# The combine_blast_results() subroutine takes the reports in {blast_results_directory} 
# and combines all reports whose names include the test string {model} into a
# single file and saves this in the {save_directory}. In that file, The line for 
# each query sequence contains the output data for both -gapped and -ungapped searches.
# If the reports contain more than one line, combine_blast_results() will only 
# retrieve the first hit.
# This is a working program, meaning it is subject to modification and change. I 
# make no promises as to readability or efficiency, but I can attest that it works
# as designed.
# Robert Carter, 13 Sep 2023.

import subprocess
import os
import io
from io import StringIO
import random
from Bio import SeqIO
import csv

input_fasta_file = "E:/chimpcontigs/sequences/contigs.fna"                   # A single file that contains all the sequences you want to examine
blast_db = "E:\panTro6\BlastDB\panTro6_db"                                   # The location and name of the blast database, \ is required on a PC, exact caps and lc also required
blast_results_directory = "e:/panTro6/BlastResults/"                         # Where the blast reports will be saved
log_file = "E:/chimpcontigs/Blast_log_file.csv"                              # Where the log file will be saved
save_directory = "E:/panTro6/"                                               # Where the summary file will be saved
model = "125onPT6"                                                           # Use this to save the results from different models

print ("Blasting 125 on PT6...") # modify this if you like

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def main():
    processed_sequence_count = 1
    contigs_file_path = "e:/chimpcontigs/125contigs.csv"
    allowed_sequence_names = set()
    with open(contigs_file_path, "r") as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            allowed_sequence_names.add(row[0].strip())

    # Load and select file names to process
    for record in SeqIO.parse(input_fasta_file, "fasta"):
        id = record.id.strip()
        if id.lower() in (name.lower() for name in allowed_sequence_names):
           sequence = str(record.seq)
           sequence_length = len(sequence)
           id = f"{id}_{model}_0_{sequence_length}"
           print(processed_sequence_count, id)
   
           output_file_ungapped = os.path.join(blast_results_directory, f"{id}_Ungapped.csv")
           output_file_gapped = os.path.join(blast_results_directory, f"{id}_Gapped.csv")
   
           if os.path.exists(output_file_ungapped):
               print(f"{processed_sequence_count} {output_file_ungapped} already exists. Skipping.")
               processed_sequence_count += 1
           else:
               print("   -ungapped...")               
               blast_sequence(sequence, output_file_ungapped, blast_db, "nodust", "nosoft", "ungapped" )
   
           if os.path.exists(output_file_gapped):
               print(f"{processed_sequence_count} {output_file_gapped} already exists. Skipping.")
           else:
               print("   -gapped...")
               blast_sequence(sequence, output_file_gapped, blast_db, "nodust", "nosoft", "gapped", )
      
        processed_sequence_count += 1

    combine_blast_results()
    bell = '\a'
    print(bell)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def blast_sequence(query_sequence, output_file, blast_db, dust, soft, gapped):

    blastn_executable = "e:/BLAST/blast-2.14.1+/bin/blastn"
    command = [
        blastn_executable,
        "-db", blast_db,
        "-evalue", "0.1",
        "-outfmt", "10 qid qstart qend sseqid sstart send pident nident length mismatch gapopen gaps evalue bitscore",
        "-max_target_seqs", "1",
        "-max_hsps", "1",
#        "-perc_identity", "50",
        "-word_size", "11",
        "-num_threads", "8",
#        "-gapopen", "5",
#        "-gapextend", "3",
#        "-reward", "2"
#        "-penalty", "-3"
#        "-xdrop_ungap", "20",
#        "-xdrop_gap", "30",
#        "-xdrop_gap_final", "100"
    ]

    # Note: these were experimantal but not used in the final analysis. Any of the parameters could be set using such a system.
    if dust == "nodust":
        command.append("-dust")
        command.append("no")
    else:
        command.append("-dust")
        command.append("yes")

    if soft == "nosoft":
        command.append("-soft_masking")
        command.append("false")
    else:
        command.append("-soft_masking")
        command.append("true")

    if gapped == "ungapped":
        command.append("-ungapped")

    try:
        with open(output_file, "w") as f_out:
            # Run BLAST and pipe the output to the output file
            p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=f_out, stderr=subprocess.PIPE, text=True)
            stdout, stderr = p.communicate(query_sequence)
            with open(log_file, "a") as log:
                log.write("Executed BLAST command:\n")
                log.write(" ".join(command) + "\n\n")
            if p.returncode != 0:
                print(f"Error occurred during BLAST search: {stderr}")
    except Exception as e:
        print(f"Error occurred: {str(e)}")

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def combine_blast_results():
    print ("Combining output data...")
    results = {}
    data_dict = {}

    # Iterate through the sequence files and extract data
    file_list = os.listdir(blast_results_directory)
    filtered_files = [filename for filename in file_list if model in filename]
    filenum = 0
    for filename in filtered_files:
        name_data = filename.split('_')
        seq_id = name_data[0]
        if seq_id not in data_dict:
            data_dict[seq_id] = {}
        start = name_data[2]
        if start not in data_dict[seq_id]:
            data_dict[seq_id][start] = {}
        length = name_data[3]
        if length not in data_dict[seq_id][start]:
            data_dict[seq_id][start][length] = {} 

        if name_data[4] == "Ungapped.csv":
            print (seq_id, start, length, "ungapped")
            file_path = os.path.join(blast_results_directory, filename)
            with open(file_path, 'r') as file:
                data = file.readline()
                ungapped_data = data.rstrip("\n")
                data_dict[seq_id][start][length]['ungapped'] = ungapped_data

        if name_data[4] == "Gapped.csv":
            print (seq_id, start, length, "gapped")
            file_path = os.path.join(blast_results_directory, filename)
            with open(file_path, 'r') as file:
                data = file.readline()
                gapped_data = data.rstrip("\n")
                data_dict[seq_id][start][length]['gapped'] = gapped_data

    results_file_path = os.path.join(save_directory, f"{model} combined_blast_results.csv")
    with open(results_file_path, "w") as output_file:
        output_string = ",,,,,Ungapped,,,,,,,,,,,,,,,Gapped,,,,,,,,,,,,,,\nQID,qstart,qlen,sseqid,sstart,send,pident,nident,length,mismatch,gapopen,gaps,evalue,bitscore,sseqid,sstart,send,pident,nident,length,mismatch,gapopen,gaps,evalue,bitscore\n"
        # Note: The things I was saving changed several times, so I cleaned up the headers in the output files after they were written
        # The output was also designed for another program that searches for random substrings
        output_file.write(output_string)
        for seq_id, start_data in data_dict.items():
            for start, len_data in start_data.items():
                for length, data in len_data.items():
                    ungapped_data = data.get("ungapped", ",,,,,,,,,,,,,,")
                    gapped_data = data.get("gapped", ",,,,,,,,,,,,,,")
                    output_string = f"{seq_id},{start},{length},"
                    output_string += f"{ungapped_data},"
                    output_string += f"{gapped_data}\n"
                    output_file.write(output_string)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

if __name__ == "__main__":

    main()
