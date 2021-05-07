import os, subprocess
from fasformat import format_upep, format_smProt, format_query,format_ncEP


ncEP_path = os.path.realpath('..\\Raw_db\\ncEP.txt')
smProt_path = os.path.realpath('..\\Raw_db\\smProt.txt')
upep_path = os.path.realpath('..\\Raw_db\\upep.txt')
query_path = os.path.realpath('..\\query.txt')

# Format text file into FASTA file
format_ncEP(ncEP_path)
format_smProt(smProt_path)
format_upep(upep_path)

# Generate BLAST database from each original databse
make_ncEP_blastdb = ['makeblastdb', '-in', 'ncEP.fasta', '-parse_seqids', '-dbtype', 'prot']
make_smProt_blastdb = ['makeblastdb', '-in', 'smProt.fasta', '-parse_seqids', '-dbtype', 'prot']
make_upep_blastdb = ['makeblastdb', '-in', 'upep.fasta', '-parse_seqids', '-dbtype', 'prot']

subprocess.call(make_ncEP_blastdb)
subprocess.call(make_smProt_blastdb)
subprocess.call(make_upep_blastdb)


print("DBS SET!!!")
