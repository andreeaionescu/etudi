from tornado.web import Application, HTTPError
from tornado.ioloop import IOLoop
from pubmed_service import EntrezConnection


class ArticleHandler(EntrezConnection):

    def get(self):
        self.render('templates/index.html', items=[])

    def post(self):
        message = self.get_body_argument('message')
        response = self.query(message)
        titles = [response[key]['title'] for key in response.keys()]
        if not response:
            raise HTTPError(404)
        self.render('templates/index.html', items=titles)


if __name__ == '__main__':
    app = Application([
        (r'/', ArticleHandler)
    ])
    app.listen(8888)
    IOLoop.instance().start()