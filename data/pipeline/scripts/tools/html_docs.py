import pandas as pd

class HtmlReport:
    """An HTML document that reports our source_data."""
    def __init__(self, stylesheets=[]):
        self.html='<html>\n'
        for ss in stylesheets:
            self.html += f'<link rel="stylesheet" href="{ss}">\n'
        self.html += '<body>\n'
    def append(self, code):
        self.html += code
        self.html += '\n'
    def heading(self, string, n=1):
        self.append(f'<h{n}> {string} </h{n}>')
    def img(self, path):
        self.append(f'<img src={path}>')
    def table(self, data):
        try:
            datahtml = data.to_html()
        except:
            datahtml = pd.DataFrame(data).to_html()
        self.append(datahtml)
    def export(self):
        self.append('</body>')
        self.append('</html>')
