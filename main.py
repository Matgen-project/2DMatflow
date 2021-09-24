#!/usr/bin/env python
# -*- coding: utf-8 -*-

from apkg.energy import *
from apkg.eex import *

if __name__ == '__main__':
    import sys
    import json

    args = sys.argv
    fp = args[1]
    tp = args[2]
    if tp.lower() == 'bulk':
        data = args[3]
        res = run_bulk(fp, data)
    elif tp.lower() == 'eex':
        lfp = args[3]
        res = run_eex(fp, lfp)
    else:
        res = run_layer(fp)
    with open(r"./{}.log".format(tp), 'w') as f:
        json.dump(res, f, indent=4)


