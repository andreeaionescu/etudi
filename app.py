from tornado.web import Application, RequestHandler
from tornado.ioloop import IOLoop
from pubmed_service import EntrezConnection


class ArticleHandler(EntrezConnection):

    def get(self):
        self.render('templates/index.html')

    def post(self):
        self.set_header("Content-Type", "text/plain")
        message = self.get_body_argument("message")
        response = self.query(message)
        self.write(response)


if __name__ == '__main__':
    app = Application([
        (r"/", ArticleHandler)
    ])
    app.listen(8888)
    IOLoop.instance().start()