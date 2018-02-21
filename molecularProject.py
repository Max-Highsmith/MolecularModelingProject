
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

FASTA_FILE_NAME = "1u19.fasta"

E_VALUE_THRESH = 0.04

PIR_FILENAME = "pir_file.txt"

fasta_string = open(FASTA_FILE_NAME).read();
result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string);
blast_records = NCBIXML.parse(result_handle)
blast_records = list(blast_records);



fh = open(PIR_FILENAME, "w");

count = 0;
for alignment in blast_records[0].alignments:
	if(count<3):
		for hsp in alignment.hsps:
			coverage = hsp.align_length / blast_records[0].query_length;	
			print(coverage);
			if (hsp.expect < E_VALUE_THRESH) and (coverage < .99):
				print("happening2");	
				fh.write(">P1;number"+str(count)+"\n")
				fh.write("sequences\n");
				fh.write(hsp.match+"\n");
				count+=1;

fh.close();
print("done");


			


print("done");



