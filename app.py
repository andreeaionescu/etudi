from tornado.web import Application, HTTPError
from tornado.ioloop import IOLoop
from pubmed_service import EntrezConnection
from utils import DateTimeEncoder
import json


class ArticleTemplatesHandler(EntrezConnection):

    def get(self):
        self.render('templates/index.html', items=[])

    def set_default_headers(self, *args, **kwargs):
        self.set_header("Access-Control-Allow-Origin", "http://localhost:8080")
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
        self.set_header('Access-Control-Allow-Origin', 'http://localhost:8080')
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
        self.set_header('Access-Control-Allow-Origin', 'http://localhost:8080')
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


if __name__ == '__main__':
    app = Application([
        (r'/', ArticleTemplatesHandler),
        (r'/pubmed', ArticleHandler),
        (r'/pubmed/([0-9]+)', ArticleByIdHandler)
    ])
    app.listen(8888)
    IOLoop.instance().start()