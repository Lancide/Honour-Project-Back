from tornado import gen, web
from tornado.ioloop import IOLoop
import pandas as p
import asyncio

from tornado.options import define, options
import subprocess

from fasformat import format_upep, format_ncEP, format_smProt, format_query
from sub_blast import search_in_ncEP, search_in_smProt, search_in_upep, ResultHandler

define("init", default=False, type=bool, help="Initiating the text database collection.")
define("port", default=8000, type=int, help="Server listening port")

# Databases & Constants
local_database = {"ncEP": None, "smProt": None, "upep": None}
local_database["ncEP"] = p.read_csv('modified_db/nc_m.csv', index_col=0)
local_database["smProt"] = p.read_csv('modified_db/sm_m.csv', index_col=0,low_memory=False)
local_database['upep'] = p.read_csv('modified_db/up_m.csv', index_col=0)

default_cols = {'0': 'chr', '1': 'starts', '2': 'ends',
                '3': 'lengths', '4': 'sequence',
                '5': 'mrna', '6': 'type', '7': 'genedescription', '8': 'species', '9': 'functions',
                '10': 'ref'}

select_result = [2, 10, 11]
gen_cln_ref = {'0': 'chr', '1': 'starts', '2': 'ends',
               '3': 'lengths', '4': 'sequence',
               '5': 'mrna', '6': 'type', '7': 'genedescription', '8': 'species', '9': 'functions',
               '10': 'ref'}
cln_names = ['Data Source','Accession ID',
             'Percentage Identity', 'E Value', 'Bit Score',
             'Chromosome Location', 'Start', 'End',
             'Length', 'Peptide Sequence',
             'mRNA Transcript', 'sORF Type', 'Gene Description', 'Species', 'Functions',
             'Literature Reference']
select_nc_cln = {'0': 0, '1': 1, '2': 2, '3': 6, '4': 5, '5': None,
                 '6': 7, '7': 12, '8': 9, '9': 14, '10': 15}  # index col does not count!!!!
select_sm_cln = {'0': None, '1': None, '2': None, '3': 6, '4': 1, '5': 2,
                 '6': 7, '7': 3, '8': 0, '9': 9, '10': 5}
select_up_cln = {'0': 0, '1': 3, '2': 1, '3': None, '4': 5, '5': 9,
                 '6': 12, '7': 2, '8': 11, '9': None, '10': 7}


class BaseHandler(web.RequestHandler):
    def set_default_headers(self):
        self.set_header("Access-Control-Allow-Origin", "*")
        self.set_header("Access-Control-Allow-Headers", "Origin, X-Requested-With, Content-Type, Accept")
        self.set_header('Access-Control-Allow-Methods', 'POST, GET, OPTIONS, PUT')

    def options(self):
        self.set_status(204)
        self.finish()


'''
one search in one db
db: str
'''


class DatabaseHandler(BaseHandler):
    def initialize(self, database):
        self.database = database  # str of add of dbs

    def get(self):
        db = self.get_argument("db")
        peptideSeq = self.get_argument("sequence")

        format_query(str(peptideSeq))  # filename: query_fa.fasta
        d = {}
        if str(db) == 'uPEPeroni':
            subprocess.call(search_in_upep)  # now have upep_out file
            handle_results = ResultHandler(local_database, 'upep_out', default_cols, db='upep')
            handle_results.set_result_cln(select_result)
            handle_results.set_gen_cln(select_up_cln)
            d = handle_results.get_res_dict()

        elif str(db) == 'sProt':
            subprocess.call(search_in_smProt)
            handle_results = ResultHandler(local_database, 'sm_out', default_cols, db='smProt')
            handle_results.set_result_cln(select_result)
            handle_results.set_gen_cln(select_sm_cln)
            d = handle_results.get_res_dict()
        elif str(db) == 'ncEP':
            subprocess.call(search_in_ncEP)
            handle_results = ResultHandler(local_database, 'nc_out', default_cols, db='ncEP')
            handle_results.set_result_cln(select_result)
            handle_results.set_gen_cln(select_nc_cln)
            d = handle_results.get_res_dict()

        df = p.DataFrame(data=d)
        col_js = []

        for index, prop in enumerate(df.columns):
            name = cln_names [index]
            print(index, prop)
            print(index, name)
            col_js.append({"prop": prop, "name": name})

        self.write({"columns": col_js, 'data': df.to_json(orient="records"), "db": db})


if __name__ == "__main__":
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

    application = web.Application([
        (r"/api/", DatabaseHandler, dict(database=local_database))
    ])
    print("Application is listening at port {}. Keep the process open.".format(options.port))
    application.listen(options.port)
    IOLoop.current().start()

# # Press the green button in the gutter to run the script.
# if __name__ == '__main__':
#     print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
