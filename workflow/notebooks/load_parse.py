import json

class ParseLoader:
    
    def __init__(self, parse_path):
        self.parse_path = parse_path
        
    def load(self):
        with open(self.parse_path, 'r') as infile:
            return {
                int(node):parse for node, parse 
                in json.load(infile).items()
            }