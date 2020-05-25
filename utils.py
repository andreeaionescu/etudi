import datetime
import json


class DateTimeEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, datetime.datetime):
            return o.isoformat()

        return json.JSONEncoder.default(self, o)


def display_title_abstract(response):
    abstracts = {}
    for article in response['PubmedArticle']:
        pmid = article['MedlineCitation']['PMID']
        title = article['MedlineCitation']['Article']['ArticleTitle']
        if article['MedlineCitation']['Article'].get('Abstract'):
            abstract = ''.join(
                [abstract_text for abstract_text in article['MedlineCitation']['Article']['Abstract']['AbstractText']])
        else:
            abstract = 'No abstract provided.'
        abstracts[pmid] = {'title': title, 'abstract': abstract}
    return abstracts


def get_publication_details(journal_issue):
    publication_details = ''
    pub_date = journal_issue.get('PubDate')
    if pub_date:
        pub_date_day = pub_date.get('Day') + ' ' if pub_date.get('Day') else ''
        pub_date_month = pub_date.get('Month') + ' ' if pub_date.get('Month') else ''
        pub_date_year = pub_date.get('Year', '')
        publication_details += f'{pub_date_day}{pub_date_month}{pub_date_year}'
    return publication_details


def display_basic_details(response):
    details = {}
    for article in response['PubmedArticle']:
        pmid = article['MedlineCitation']['PMID']
        title = article['MedlineCitation']['Article']['ArticleTitle']
        copyright = article['MedlineCitation']['Article']['Abstract'].get('CopyrightInformation', '') if article['MedlineCitation']['Article'].get('Abstract') else None
        keywords = article['MedlineCitation'].get('KeywordList', [])
        authors = [author.get('ForeName', '') + ' ' + author.get('LastName', '') for author in article['MedlineCitation']['Article'].get('AuthorList', [])]
        journal_name = article['MedlineCitation']['Article']['Journal']['Title']
        publication_details = get_publication_details(article['MedlineCitation']['Article']['Journal']['JournalIssue'])  #TODO do we want more details on journal such as JournalIssue, pubdated etc.
        abstract = ''.join([
            abstract_text for abstract_text in article['MedlineCitation']['Article']['Abstract']['AbstractText']]) \
            if article['MedlineCitation']['Article'].get('Abstract') \
            else 'No abstract available.'
        details[pmid] = {'Title': title,
                         'Abstract': abstract,
                         'Authors': authors,
                         'Copyright': copyright,
                         'Keywords': keywords,
                         'Journal': {
                             'PublicationDetails': publication_details,
                             'Title': journal_name
                            }
                         }
    return details
