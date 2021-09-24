#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json


def load(fp):
    with open(fp, 'r') as f:
        data = json.load(f)
    return data


def read_bulk(bfp):
    return load(bfp)


def read_layer(lfp):
    return load(lfp)


def get_eex(Eiso, Ebulk, A, m, n=1):
    # unit meV/A2
    return 1000 * (Eiso - n * Ebulk / m) / A


def run_eex(bfp, lfp):
    bd = read_bulk(bfp)
    ld = read_layer(lfp)
    exres = {}
    for k, v in ld.items():
        mbd = bd.get(k)
        if mbd is None:
            continue
        ebulk = mbd.get('energy')
        elayer = v.get('energy')
        m = mbd.get('ln')
        area = mbd.get('area')
        Eex = get_eex(elayer, ebulk, area, m)
        exres.update(
            {k: Eex}
        )
    return exres
        

if __name__ == "__main__":
    pass

