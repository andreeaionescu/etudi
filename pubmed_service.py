import logging
import json

from tornado.web import RequestHandler
from Bio import Entrez
from utils import display_basic_details
from xml_parser import XmlDictConfig
from xml.etree import cElementTree as ElementTree


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

    def entrez_link_pubmed_pmc(self, pubmed_id):
        return Entrez.read(Entrez.elink(dbfrom='pubmed', db='pmc', LinkName="pubmed_pmc", id=pubmed_id))

    def get_pmc_from_pubmed(self, pubmed_id):
        link_db = self.entrez_link_pubmed_pmc(pubmed_id)[0].get('LinkSetDb')
        return link_db[0].get('Link')[0].get('Id') if len(link_db) != 0 else None

    def entrez_fetch_full_text(self, pubmed_id):
        pmc_id = self.get_pmc_from_pubmed(pubmed_id)
        return Entrez.efetch(db='pmc', id=pmc_id).read() if pmc_id else None

    def query_full_text(self, pubmed_id):
        pmc_id = self.get_pmc_from_pubmed(pubmed_id)
        entrez_fetch_full_text = self.entrez_fetch_full_text(pubmed_id)
        if entrez_fetch_full_text:
            root = ElementTree.XML(entrez_fetch_full_text)
            return {
                'pubmed_id': pubmed_id,
                'pmc_id': pmc_id,
                'full_text': XmlDictConfig(root).get('article', None)
                }
        else:
            return {
                'pubmed_id': pubmed_id,
                'pmc_id': pmc_id,
                'full_text': None
            }

    def query(self, stmt):
        return display_basic_details(self.entrez_read(stmt))
