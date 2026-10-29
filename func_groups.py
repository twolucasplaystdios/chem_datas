import numpy as np

import json


tokens = {}

with open("tokens.json") as f:
    tokens = json.load(f) 


def get_token(sym: str):
    global tokens
    
    return tokens[sym]


# 羥基
HYDROXYL = {
    "atoms": ["1", "O", "H"],
    "edges": [
        (0, 1, 1),
        (1, 2, 1),
    ],
}

