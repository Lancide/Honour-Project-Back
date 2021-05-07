import pandas as p

def format_ncEP(in_file, separator = '\t'):
    ncEP = p.read_csv(in_file, sep=separator)
    modified = ncEP.dropna(subset = ['sequence'])
    modified = modified.drop_duplicates(subset = ['peptide'])
    modified.to_csv('nc_m.csv', index=False)
    with open("ncEP.fasta", 'wt') as outfile:
        for row in modified.itertuples():
            outfile.write(">" + str(row.peptide) + "\n" + str(row.sequence) + '\n')


def format_smProt(in_file, separator = "\t"):
    smProt = p.read_csv(in_file, sep=separator, encoding='unicode_escape', low_memory=False)
    modified = smProt.dropna(subset = ['Sequence'])
    modified.to_csv('sm_m.csv', index=False)
    with open("smProt.fasta", 'wt') as outfile:
        for row in modified.itertuples():
            outfile.write(">" + str(row.SmProtID) + "\n" + str(row.Sequence) + "\n")


def format_upep(in_file, separator = "~"):
    upep = p.read_csv(in_file, sep=separator, encoding='unicode_escape')
    modified = upep.dropna(subset= ['expSeq'])
    modified = modified.drop_duplicates(subset = ['accessionId'])
    modified.to_csv('up_m.csv', index=False)
    with open('upep.fasta', "wt") as outfile:
        for row in modified.itertuples():
            outfile.write(">" + str(row.accessionId) + "\n" + str(row.expSeq) + "\n")


def format_query(query_seq):
    query_fa = open("query_fa.fasta","wt")
    query_fa.write(">" + "Query Sequence" + "\n" + str(query_seq) + "\n")
    query_fa.close()




