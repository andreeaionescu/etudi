import requests
from pprint import pprint

# https://pubchem.ncbi.nlm.nih.gov/rest/pug/<input specification>/<operation specification>/[<output specification>][?<operation_options>]

URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
records = ['substances', 'compounds', 'bioassays']


def retrive_data_by_input_specification(input_specs, output_specs, output_format):
    try:
        return requests.get(URL + '/{0}/{1}/{2}'.format(input_specs, output_specs, output_format) )
    except Exception as e:
        print(e)


def main():
    aspirin_input_specs = 'compound/cid/2244'
    aspirin_output_specs = 'property'
    output_format = 'JSON'  # XML | CSV | PNG | TXT
    print('LOL')
    # aspirin_record = retrive_data_by_input_specification(aspirin_input_specs, aspirin_output_specs, output_format)
    #pprint(aspirin_record)
