from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import primer3
import time

start_time = time.time()

def blastRequest():
	NCBIWWW.email = "tom@thomasredden.com"


	sequence = open("mySequence.txt").read()
	result_handle = NCBIWWW.qblast("blastn", "nt", sequence, entrez_query = "Medicago truncatula[organism]", hitlist_size = 4)


	with open("my_blast.xml", "w") as out_handle:
        	out_handle.write(result_handle.read())

	result_handle.close()

	print("--- Reached Checkpoint 1 in %s seconds ---" % (time.time() - start_time))

#	generatePrimers()

def generatePrimers():
#	primer3.calcTm('GTAAAACGACGGCCAGT')

	for record in SeqIO.parse("sample_genes.fasta", "fasta"):

		sequence = str(record._seq)
		id = record.id
		length = len(sequence)

		primers = primer3.bindings.designPrimers(
		{
			'SEQUENCE_ID': id,
			'SEQUENCE_TEMPLATE': sequence
		},
		{
			'PRIMER_OPT_SIZE': 20,
			'PRIMER_PICK_INTERNAL_OLIGO': 1,
			'PRIMER_INTERNAL_MAX_SELF_END': 8,
			'PRIMER_MIN_SIZE': 18,
			'PRIMER_MAX_SIZE': 25,
			'PRIMER_OPT_TM': 60.0,
			'PRIMER_MIN_TM': 57.0,
			'PRIMER_MAX_TM': 63.0,
			'PRIMER_MIN_GC': 20.0,
			'PRIMER_MAX_GC': 80.0,
			'PRIMER_MAX_POLY_X': 100,
			'PRIMER_INTERNAL_MAX_POLY_X': 100,
			'PRIMER_SALT_MONOVALENT': 50.0,
			'PRIMER_DNA_CONC': 50.0,
			'PRIMER_MAX_NS_ACCEPTED': 0,
			'PRIMER_MAX_SELF_ANY': 12,
			'PRIMER_MAX_SELF_END': 8,
			'PRIMER_PAIR_MAX_COMPL_ANY': 12,
			'PRIMER_PAIR_MAX_COMPL_END': 8,
			'PRIMER_PRODUCT_SIZE_RANGE' : [length - 100,length],
			'PRIMER_NUM_RETURN' : 5
			})
		print(primers['PRIMER_PAIR_NUM_RETURNED'])

generatePrimers()

print("--- Reached Checkpoint 2 in %s seconds ---" % (time.time() - start_time))
