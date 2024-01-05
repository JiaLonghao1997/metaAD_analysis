from Bio import SeqIO
import sys
 
long_sequences = [] # Setup an empty list
sample = sys.argv[1]
infasta = sys.argv[2]
threshold = sys.argv[3]
fasta_threshold = sys.argv[4]
handle = open(infasta, "r")

for record in SeqIO.parse(handle, "fasta") :
    if len(record.seq) >= int(threshold) :
        # Add this record to our list
        record.id = sample + '+' + record.id
        long_sequences.append(record)
handle.close()
 
#print "Found %i long sequences" % len(long_sequences)
 
output_handle = open(fasta_threshold, "w")
SeqIO.write(long_sequences, output_handle, "fasta")
output_handle.close()
