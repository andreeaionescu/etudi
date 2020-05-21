import logging
import json

from tornado.web import RequestHandler
from Bio import Entrez
from utils import display_basic_details


class EntrezConnection(RequestHandler):
    Entrez.email = 'abcdefg@example.com'

    def __init__(self, application, request, db='pubmed', retmax='250', retmode='json', **kwargs):
        self.db = db
        self.retmax = retmax
        self.retmode = retmode
        super().__init__(application, request, **kwargs)

    def entrez_search(self, query, **kwargs):
        return Entrez.esearch(db=self.db, retmax=self.retmax, retmode=self.retmode, **kwargs, term=query)

    def entrez_handle(self, query):
        try:
            initial_search = self.entrez_search(query)
            if initial_search:
                logging.info('Successfully performed search on Entrez.')
                search = json.loads(initial_search.read())
                ids = ','.join(search['esearchresult']['idlist'])
                return Entrez.efetch(db=self.db, retmode='xml', id=ids)
        except Exception as e:
            logging.error('Something went wrong with the Entrez handle:', e)

    def entrez_read(self, query):
        return Entrez.read(self.entrez_handle(query))

    def query(self, stmt):
        return display_basic_details(self.entrez_read(stmt))
