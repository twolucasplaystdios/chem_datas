import numpy as np

import json


tokens = {}

with open("tokens.json") as f:
    tokens = json.load(f) 


def get_token(sym: str):
    global tokens
    
    return tokens[sym]


def merge2(f1, f2):
    f1["atoms"]
    
    return