import requests
from xml.etree import ElementTree
from pprint import pprint


URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def parse_webenv_and_querykey(xml_root):
    return xml_root.find('WebEnv').text, xml_root.find('QueryKey').text


def retrieve_data_by_query(db, search_txt, output_format):
    try:
        url_base = URL + f'esearch.fcgi?db={db}&term={search_txt}&reldate=1&datetype=edat&retmax=10&usehistory=y'
        response = requests.get(url_base)
        if response.status_code == 200:
            root = ElementTree.fromstring(response.content)
            webenv, querykey = parse_webenv_and_querykey(root)

            # assemble the esummary URL
            esummary_url = url_base + f'esummary.fcgi?db={db}&query_key={querykey}&WebEnv={webenv}'
            # post e-summary url
            docs = requests.get(esummary_url)
            docs.raise_for_status()

            # include for ESearch-EFetch
            # assemble the efetch URL
            efetch_url = url_base + f'efetch.fcgi?db={db}&query_key={querykey}&WebEnv={webenv}&rettype=abstract&retmode={output_format}'
            # post the efetch URL
            data = requests.get(efetch_url)
            data.raise_for_status()
            return data.json() # return body with text abstracts instead of response similar to a POST request

    except requests.exceptions.RequestException as e:
        raise SystemExit(e)


if __name__ == '__main__':
    db = 'pubmed'
    search_txt = 'aspirin'
    output_format = 'JSON'  # XML | CSV | PNG | TXT
    aspirin_record = retrieve_data_by_query(db, search_txt, output_format)
    pprint(aspirin_record)  # returns an xml with UIDs to matching entries
