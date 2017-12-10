import sys
import argparse
import Bio
import Bio.SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-f",required=True,metavar="FILE",help="The path to the FASTA file")
parser.add_argument("-r",required=True,metavar="FILE",help="The path to the RSEM file")
parser.add_argument("-p",required=True,metavar="FILE",help="The path to PASA pipeline *pasa_assemblies_described.txt")
parser.add_argument("--TPM_threshold",type=float,default=0.01,metavar="0.01",help="If greater than this TPM threshold, pass the filter (propotion of the max transcript in a cluster, rather than an absolute value)")
args = parser.parse_args()

records_to_print = set()
records_seen = set()

#######
## Load the RSEM results into a quickly queryable format
#######

expression_dict = dict()
handle = open(args.r,"rU")
i=1
for line in handle.readlines():
	splitline = line.split("\t")
	if i == 1:
		##Header check
		assert splitline == ['transcript_id', 'gene_id', 'length', 'effective_length', 'expected_count', 'TPM', 'FPKM', 'IsoPct\n']
		i+=1
		continue
	expression_dict[splitline[0]] = float(splitline[5]) ##Store the TPM in the dictionary
	i+=1	

sys.stderr.write("Done loading expression data.\n")
sys.stderr.flush()

#######
##Iterate over the PASA assembly clusters and save the records to print in records_to_print
#######

handle = open(args.p,"rU")
i=1
for line in handle.readlines():
	splitline = line.split("\t")
	if i == 1:
                ##Header check
                assert splitline == ['#scaffold', 'subcluster_id', 'asmbl_acc', 'cdna_accs', 'alignment_description\n']
                i+=1
                continue
	assembly_id = splitline[2]
	included_transcripts = splitline[3].split(",")

	##Find max expression
	max_expression = 0.0
	for t in included_transcripts:
		current_TPM = expression_dict[t]
		if current_TPM > max_expression:
			max_expression = current_TPM
	
	##Filter by expression
	for t in included_transcripts:
		records_seen.add(t)
		current_TPM = expression_dict[t]
		if current_TPM > max_expression * args.TPM_threshold:
			if t in records_to_print:
				#print("Error. Record:"+t+" should not already be recorded.")
				pass
			else:
				records_to_print.add(t)
		else:
			pass
			##Filter out this transcript, by not printing it

sys.stderr.write("Done iterating over PASA clusters.\n")
sys.stderr.flush()

#######
## Iterate over the FASTA file, and only print it if it is in the list.
#######

sys.stderr.write("Done iterating over FASTA file and printing accepted records to stdout.\n")
sys.stderr.flush()

handle = Bio.SeqIO.parse(args.f,"fasta")
i=0
j=0
k=0
for record in handle:
	if record.id in records_to_print:
		sys.stdout.write(record.format("fasta"))
		i+=1
	elif record.id in records_seen:
		j+=1
	elif record.id not in records_seen:
		k+=1
sys.stderr.write("Done filtering FASTA file.\n")
sys.stderr.write(str(i)+" records passed the filter.\n")
sys.stderr.write(str(j+k)+" records were filtered out.\n")
sys.stderr.write("Of the filtered records, "+str(j)+" were filtered based on their low TPM.\n")
sys.stderr.write("Of the filtered records, "+str(k)+" were filtered because they weren't in the PASA file (likely failed alignment validation or didn't align to genome).\n")
sys.stderr.flush()
