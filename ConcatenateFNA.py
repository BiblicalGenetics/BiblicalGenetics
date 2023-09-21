import os

input_dir = r'E:\PT3.1.1\Sequences'
output_file = r'E:\PT3.1.1\Sequences\PT3.fna'

with open(output_file, 'w') as output_fna:
    for filename in os.listdir(input_dir):
        if filename.endswith('.fasta'):
            fasta_file_path = os.path.join(input_dir, filename)
            with open(fasta_file_path, 'r') as input_fasta:
                output_fna.write(input_fasta.read())

print(f'All sequences concatenated to "{output_file}".')
