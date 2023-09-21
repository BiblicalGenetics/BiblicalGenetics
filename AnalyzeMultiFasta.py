import os
import re
import csv
from Bio import SeqIO

folder_path = f"E:/PT3.1.1/Sequences"
block_lengths = []
output_csv_path = "E:/PT3.1.1/N_Histogram.csv" # Change to suit
output_stats_path = "E:/PT3.1.1/SequenceStats.csv" # Change to suit

def find_n_blocks(sequence):
    n_blocks = [len(match.group()) for match in re.finditer(r'N+', sequence)]
    return n_blocks

def main():
    A = 0
    C = 0
    G = 0
    T = 0
    N = 0
    Nblocks = 0
    stats = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".fna"): # FNA = FASTA nucleotide format; change to 'FASTA' if examing a multi-sequence FASTA file
            print (filename)
            file_path = os.path.join(folder_path, filename)
            for record in SeqIO.parse(file_path, "fasta"):
                sequence_name = record.id
                sequence_length = len(record.seq)
                numAs = record.seq.count("A")
                numCs = record.seq.count("C")
                numGs = record.seq.count("G")
                numTs = record.seq.count("T")
                numas = record.seq.count("a")
                numcs = record.seq.count("c")
                numgs = record.seq.count("g")
                numts = record.seq.count("t")
                numNs = record.seq.count("N")
                sequence = str(record.seq)
                n_blocks = find_n_blocks(sequence)
                block_lengths.extend(n_blocks)
                numnblocks = len(n_blocks)
                print (sequence_name, sequence_length, numAs, numas, numCs, numcs, numGs, numgs, numTs, numts, numNs, numnblocks)
                row_data = [sequence_name,sequence_length,numAs,numas,numCs,numcs,numGs,numgs,numTs,numts,numNs,numnblocks]
                stats.append(row_data)
                A += numAs
                A += numas
                C += numCs
                C += numcs
                G += numGs
                G += numgs
                T += numTs
                T += numts
                N += numNs
                Nblocks += numnblocks
                
    print("Totals:")
    print ("   A: ", A)
    print ("   C: ", C)
    print ("   G: ", G)
    print ("   T: ", T)
    print ("   N: ", N)
    print ("   Nblocks: ", Nblocks)    
    
    with open(output_stats_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = ['Sequence Name', 'Sequence Length', 'NumAs', 'Numas', 'NumCs', 'Numcs', 'NumGs', 'Numgs', 'NumTs', 'Numts', 'NumNs', 'Numnblocks']
        writer.writerow(header)
        for stat_row in stats:
            writer.writerow(stat_row)

    with open(output_csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Block Size'])
        for length in block_lengths:
            writer.writerow([length])

    print(f"Data saved to: {output_csv_path}")

if __name__ == "__main__":
    main()
