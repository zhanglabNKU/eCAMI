The input is a fasta file to be clustered. The known label of each sequence can be added after the sequence ID and separated with a '|', e.g., the EC numbers or CAZy subfamilies known for this protein. 

Two files are output from CAMI clustering for each cluster: (i) fasta sequence, and (ii) characteristic k-mer peptides. 

For (ii), the output file has its first line being the cluster number/ID (starting from 0), and the second line being the counts of known labels (e.g., the labels will be the EC numbers or CAZy subfamilies provided in the input fasta file; if there are multiple lables, they are separated with tab spaces), and the rest lines are each k-mer peptide followed by their frequency.
