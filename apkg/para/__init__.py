#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json

settings = {}
js = [os.path.join(os.path.dirname(__file__), i)
      for i in os.listdir(os.path.dirname(__file__))
      if os.path.splitext(i)[-1] == '.json']

for file in js:
    with open(file, "r") as f:
        js_obj = json.load(f)
        settings[os.path.basename(file).replace('.json', '')] = js_obj


if __name__ == "__main__":
    pass
