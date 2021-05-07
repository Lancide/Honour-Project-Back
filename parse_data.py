from tornado.httpclient import AsyncHTTPClient
from tornado import gen
from tornado.httputil import url_concat
from tornado.escape import url_escape, json_decode
from bs4 import BeautifulSoup
import urllib.parse
import pandas as pd


class BaseParser:
    def __init__(self, url, agents):
        self.url = url
        self.agents = agents
        self.client = AsyncHTTPClient()
        self.response = None

    def set_agents(self, agents):
        self.agents = agents

    @gen.coroutine
    def get_response(self):
        self.response = yield self.client.fetch(self.url)

    @gen.coroutine
    def parse(self):
        pass


class NCEPParser(BaseParser):
    def __init__(self, url, agents):
        super().__init__(url, agents)

    @gen.coroutine
    def parse(self):
        yield self.get_response()
        body = BeautifulSoup(self.response.body, 'html.parser')

        df = []
        for row in body.find_all('td'):
            url = self.url.replace("browse.html", urllib.parse.quote(row.a["href"], safe="./:?="))
            response = yield self.client.fetch(url)
            bs = BeautifulSoup(response.body, 'html.parser')
            row_count = 0
            data = []
            for tr in bs.find_all('tr'):
                if not row_count:
                    row_count += 1
                else:
                    url_parsed = urllib.parse.urlparse(urllib.parse.unquote(tr.find_all('td')[12].a["href"]))
                    r = urllib.parse.parse_qs(url_parsed.query)
                    for k in r:
                        r[k] = r[k][0]

                    data.append(r)
            df.append(pd.DataFrame(data))

            yield gen.sleep(1)
        df = pd.concat(df, ignore_index=True)
        raise gen.Return(df)

def get_ncEP_async(mail):
    p = NCEPParser("http://www.jianglab.cn/ncEP/browse.html", {"email": mail})
    response = yield p.parse()
    response.to_csv("ncEP.txt", sep="\t", index=False)

###
email = 'lanxin.yi@uqconnect.edu.au'
get_ncEP_async(email)