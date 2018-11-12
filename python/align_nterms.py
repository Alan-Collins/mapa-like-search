from Bio.Align.Applications import ClustalOmegaCommandline
in_file = "to_align.txt"
out_file = "aligned.fasta"
clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
print(clustalomega_cline)