
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
