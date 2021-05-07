from fasformat import format_ncEP
import subprocess
import os
import pandas as p


make_ncEP_blastdb = ['makeblastdb', '-in', 'ncEP.fasta', '-parse_seqids', '-dbtype', 'prot']
search_in_ncEP = ['blastp', '-query', 'query_fa.fasta', '-out', 'nc_out', '-db', 'ncEP.fasta', '-outfmt', '10']

make_smProt_blastdb = ['makeblastdb', '-in', 'smProt.fasta', '-parse_seqids', '-dbtype', 'prot']
search_in_smProt = ['blastp', '-query', 'query_fa.fasta', '-out', 'sm_out', '-db', 'smProt.fasta', '-outfmt', '10']

make_upep_blastdb = ['makeblastdb', '-in', 'upep.fasta', '-parse_seqids', '-dbtype', 'prot']
search_in_upep = ['blastp', '-query', 'query_fa.fasta', '-out', 'upep_out', '-db', 'upep.fasta', '-outfmt', '10']

set_database = ["py", 'dbs\set_dbs.py']


# remove duplication
# ncEP = p.read_csv('ncEP.txt', sep="\t")
# modified = ncEP.dropna(subset = ['sequence'])
# modified = modified.drop_duplicates(subset = ['peptide'])

ncEP = p.read_csv('modified_db/nc_m.csv', index_col=0)

smProt = p.read_csv('modified_db/sm_m.csv', index_col=0,low_memory=False)

upep = p.read_csv('modified_db/up_m.csv', index_col=0)



nc_cols = {'0': 'chr', '1': 'start', '2': 'end', '3': 'sequence', '4': '', '5': 'length',
           '6': 'ncType', '7': 'species', '8': 'function', '9': 'description', '10': 'title'}


class ResultHandler:
    databases = {}
    result_df = p.DataFrame()
    gen_df = p.DataFrame()
    select_res_cln = []
    selected_gen_cln = []
    format_10_cln = ["query", "id", "identity", "align_len", "mismatches", "gaps",
                     "q_start", "q_end", "s_start", "s_end", "eval", "bitscore"]
    db = ''

    upeps_sp = {'1': 'Human', '2': 'Drosophila melanogaster', '3': 'Mouse', '4': 'Arabidopsis thalliana'}
    upeps_type = {'1': ['uORF', 'A sPEP that encoded from an upstream open reading frame in a transcript.'],
                  '2': ["oORF", 'A sPEP that encoded from a main coding sequence in a transcript.'],
                  '3': ['dORF', 'A sPEP that encoded from a downstream open reading frame in a transcript.'],
                  '4': ['AltORF',
                        'A sPEP that encoded from the non-canonical +2/ or +3 open reading frame in a transcript.'],
                  '5': ['ncRNA', 'A sPEP that encoded from a non-coding RNA transcript.']}
    gen_cln_d = {'0': 'chr', '1': 'starts', '2': 'ends',
                 '3': 'lengths', '4': 'sequence',
                 '5': 'mrna', '6': 'type', '7': 'genedescription', '8': 'species', '9': 'functions',
                 '10': 'ref'}

    db_cln = []


    columns = {}

    final_dic = {}
    seq_ids = []

    id_cln = {}
    res_cln = {}
    gen_cln = {}

    def __init__(self, databases, result_file, uni_cols, db):
        self.databases = databases
        self.db = db
        self.gen_df = databases[db]

        self.result_df = p.read_csv(result_file, sep=',',
                                    names=self.format_10_cln, index_col=1)

        self.seq_ids = [str(i) for i in self.result_df.index]

        self.gen_cln_d = uni_cols


        self.id_cln = {'0': 'db', '1': 'seqid'}

        for id_i, col_name in self.id_cln.items():
            self.final_dic[col_name] = []
            self.columns[str(id_i)] = col_name
        for each_r in self.seq_ids:
            self.final_dic[self.id_cln['0']].append(self.db)
            self.final_dic[self.id_cln['1']].append(each_r)

    def set_result_cln(self, selected):
        """
        :param selected: list of int indices (from 0) of which cln is selected from format_10_cln
        (excluding the seq id column, whose cln index is 1)
        :return:
        """
        initial_num_of_cln = len(self.columns)

        for i, c in enumerate(selected):
            col_name = self.format_10_cln[c]
            uni_index = i + initial_num_of_cln
            self.columns[str(uni_index)] = col_name
            self.res_cln[str(uni_index)] = col_name
            self.final_dic[col_name] = []

            for each_r in self.seq_ids:
                info = self.result_df.at[each_r, col_name]
                self.final_dic[col_name].append(info)

    def set_gen_cln(self, col_map):
        '''
        :param col_map: {str(int): int}
            example -
            gen_cln_ref = {'0': 'chr', '1': 'starts', '2': 'ends',
                            '3': 'lengths', '4': 'sequence',
                            '5': 'mrna', '6': 'type', '7': 'genedescription', '8': 'species', '9': 'functions',
                            '10': 'ref'}
            select_nc_cln = {'0': 0, '1': 1, '2': 2, '3': 6, '4': 5, '5': None,
                            '6': 7, '7': 12, '8': 9, '9': 14, '10': 15}
            index column (the acc. of sequence) does not count!!!!
        :return: None
        '''
        initial_num_of_cln = len(self.final_dic)
        self.db_cln = [str(c) for c in self.gen_df.columns]

        for i, k in enumerate(self.gen_cln_d):
            uni_col_name = self.gen_cln_d[k]
            uni_index = initial_num_of_cln + i

            self.columns[str(uni_index)] = str(uni_col_name)
            self.gen_cln[str(uni_index)] = str(uni_col_name)
            self.final_dic[uni_col_name] = []

            db_col_id = col_map[k]

            for each_r in self.seq_ids:

                if db_col_id is not None:
                    db_col_name = self.db_cln[db_col_id]
                    info = self.gen_df.at[each_r, db_col_name]
                    self.final_dic[uni_col_name].append(info)
                else:
                    info = 'Invalid'
                    self.final_dic[uni_col_name].append(info)

        # need customisation in some cases
        if self.db == 'upep':
            for index, tid in enumerate(self.final_dic['type']):
                type = self.upeps_type[str(tid)][0]
                self.final_dic['type'][index] = type
            for index, spid in enumerate(self.final_dic['species']):
                species = self.upeps_sp[str(spid)]
                self.final_dic['species'][index] = species



    def get_res_dict(self):
        d ={}
        for k, val in self.final_dic.items():
            d[k] = val
        return d
dbs_file = {"ncEP": ncEP, "smProt": smProt, "upep": upep}
default_cols = {'0': 'chr', '1': 'starts', '2': 'ends',
                 '3': 'lengths', '4': 'sequence',
                 '5': 'mrna', '6': 'type', '7': 'genedescription', '8': 'species', '9': 'functions',
                 '10': 'ref'}

select_result = [2, 10, 11]
select_up_cln = {'0': 0, '1': 3, '2': 1, '3': None, '4': 5, '5': 9,
                 '6': 12, '7': 2, '8': 11, '9': None, '10': 7}
select_sm_cln = {'0': None, '1': None, '2': None, '3': 6, '4': 1, '5': 2,
                 '6': 7, '7': 3, '8': 0, '9': 9, '10': 5}
select_nc_cln = {'0': 0, '1': 1, '2': 2, '3': 6, '4': 5, '5': None,
                 '6': 7, '7': 12, '8': 9, '9': 14, '10': 15}
cln_names = ['Data Source','Accession ID',
             'Percentage Identity', 'E Value', 'Bit Score',
             'Chromosome Location', 'Start', 'End',
             'Length', 'Peptide Sequence',
             'mRNA Transcript', 'sORF Type', 'Gene Description', 'Species', 'Functions',
             'Literature Reference']










