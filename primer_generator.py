from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import primer3
import time
import pandas as pd

start_time = time.time()

def blastRequest():
	NCBIWWW.email = "tom@thomasredden.com"


	sequence = open("mySequence.txt").read()
	result_handle = NCBIWWW.qblast("blastn", "nt", sequence, entrez_query = "Medicago truncatula[organism]", hitlist_size = 4)


	with open("my_blast.xml", "w") as out_handle:
        	out_handle.write(result_handle.read())

	result_handle.close()

def generatePrimers():

	df = pd.DataFrame({'gene_id' : [],
			   'primer1' : [],
			   'primer2' : []
			  })
	for record in SeqIO.parse("input.fasta", "fasta"):

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

		complementarySequence = primers['PRIMER_RIGHT_0_SEQUENCE']
		noncompSequence = primers['PRIMER_LEFT_0_SEQUENCE']
		print("primer1: ", complementarySequence, "primer2: ", noncompSequence)
		df.loc[len(df.index)] = [id, complementarySequence, noncompSequence]

	df.to_csv("result.csv", index=False)

generatePrimers()

print("--- Finished  in %s seconds ---" % (time.time() - start_time))
