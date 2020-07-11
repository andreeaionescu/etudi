from tornado.web import Application, HTTPError
from tornado.ioloop import IOLoop
from flask import Flask, request
from pubmed_service import EntrezConnection
from utils import DateTimeEncoder
import json

app = Flask(__name__)
CLIENT_URL = "https://develop.d12x3ux7nu3fim.amplifyapp.com"


class ArticleTemplatesHandler(EntrezConnection):

    def get(self):
        self.render('templates/index.html', items=[])

    def set_default_headers(self, *args, **kwargs):
        self.set_header("Access-Control-Allow-Origin", CLIENT_URL)
        self.set_header('Access-Control-Allow-Headers', 'Origin, Content-Type, X-Auth-Token')
        self.set_header('Access-Control-Allow-Methods', 'POST, GET, OPTIONS')

    def options(self):
        pass

    def post(self):
        message = self.get_body_argument('message')
        response = self.query(message)
        titles = [response[key]['abstract'] for key in response.keys()]
        if not response:
            raise HTTPError(404)
        self.render('templates/index.html', items=titles)


class ArticleHandler(EntrezConnection):

    def set_default_headers(self, *args, **kwargs):
        self.set_header('Access-Control-Allow-Origin', CLIENT_URL)
        self.set_header('Access-Control-Allow-Headers', 'Origin, Content-Type, X-Auth-Token')
        self.set_header('Access-Control-Allow-Methods', 'POST, GET, OPTIONS')
        # self.set_header('Access-Control-Allow-Credentials', 'true')
        self.set_header('Content-Type', 'application/json')

    def options(self):
        pass

    def post(self):
        data = json.loads(self.request.body.decode('utf-8'))
        print('Got JSON data:', data)
        text = data.get('search')
        response = self.query(text)
        print('Response: ', response)
        if not response:
            raise HTTPError(404)
        self.write(json.dumps(response, cls=DateTimeEncoder))


class ArticleByIdHandler(EntrezConnection):

    def set_default_headers(self, *args, **kwargs):
        self.set_header('Access-Control-Allow-Origin', CLIENT_URL)
        self.set_header('Access-Control-Allow-Headers', 'Origin, Content-Type, X-Auth-Token')
        self.set_header('Access-Control-Allow-Methods', 'POST, GET, OPTIONS')
        self.set_header('Content-Type', 'application/json')

    def options(self, pubmed_id):
        pass

    def post(self, pubmed_id):
        data = json.loads(self.request.body.decode('utf-8'))
        print('Got JSON data:', data)
        article_id = data.get('id')
        response = self.query_full_text(article_id)
        print('Response: ', response)
        if not response:
            raise HTTPError(404)
        self.write(json.dumps(response, cls=DateTimeEncoder))


def start_tornado_app():
    server = Application([
        (r'/', ArticleTemplatesHandler),
        (r'/pubmed', ArticleHandler),
        (r'/pubmed/([0-9]+)', ArticleByIdHandler)
    ])
    server.listen(8888)
    IOLoop.instance().start()


app = Flask('etudi')

@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE')
    return response


@app.route('/pubmed', methods=["POST"])
def article_handler():
    entrez_connection = EntrezConnection()
    data = request.get_json()
    print('Got JSON data:', data)
    text = data.get('search')
    response = entrez_connection.query(text)
    print('Response: ', response)
    if not response:
        raise HTTPError(404)
    return json.dumps(response, cls=DateTimeEncoder)


@app.route('/pubmed/<id>', methods=["POST"])
def article_handler_by_id(id):
    entrez_connection = EntrezConnection()
    data = request.get_json()
    print('Got JSON data:', data)
    article_id = data.get('id')
    response = entrez_connection.query_full_text(article_id)
    print('Response: ', response)
    if not response:
        raise HTTPError(404)
    return json.dumps(response, cls=DateTimeEncoder)


if __name__ == '__main__':
    app.run()
