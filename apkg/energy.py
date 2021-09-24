#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from os.path import *
import re
import json
import numpy as np

from monty.re import reverse_readfile
from pymatgen.analysis import dimensionality
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import SupercellTransformation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from apkg.para import settings


def pull_output(pth):
    for root, _, files in os.walk(pth):
        if 'Scf' in root:
            for f in files:
                if f == 'OUTCAR':
                    yield join(root, f)


def get_energy(ofp):
    for info in reverse_readfile(ofp):
        if 'energy  without entropy' in info:
            energy = info.split('=')[-1]
            return energy

    return None


def get_id(ofp):
    rdn = ofp.split(os.sep)[-3]
    code = re.findall(r"(\d+)", rdn)
    return code[0]


def get_layer(struct):
    cc = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    ca = dimensionality.find_connected_atoms(cc, tolerance=0.45,
                                             ldict=settings['covalent_radii'])
    cluster = dimensionality.find_clusters(cc, ca)
    return cluster


def get_direction_from(fp, code, tp='icsd'):
    with open(r"/HOME/nscc-gz_material_1/matgen_dft/mat2d_work/scripts/analysis/2danalysis/apkg/bulk_info.json", "r") as f:
        bulk_data = json.load(f)

    key = tp+ '_' + str(code)
    sin = bulk_data.get(key) 
    print(key)
    print(sin)
    return sin

def get_direction(bulk_struct):
    cc = SpacegroupAnalyzer(bulk_struct).get_conventional_standard_structure()
    _, _, init_cc = get_layer(cc)
    inc = len(init_cc)
    tm = [[2, 1, 1], [1, 2, 1], [1, 1, 2]]
    cnl = []
    for i in tm:
        sc = SupercellTransformation(scaling_matrix=i).apply_transformation(cc)
        ca = dimensionality.find_connected_atoms(sc, tolerance=0.45,
                                                 ldict=settings['covalent_radii'])
        _, _, cluster = dimensionality.find_clusters(sc, ca)
        cnl.append(len(cluster))
    for i, v in enumerate(cnl):
        if v / inc == 2:
            return i

def get_polygon_area_in_plane(struct, direction):
    ps = [(1, 2), (0, 2), (0, 1)]
    index = ps[direction]
    la = list(map(np.asarray, struct.lattice.as_dict().get('matrix')))
    area = np.linalg.norm(np.cross(la[index[0]], la[index[1]]))
    return area


def run_layer(pth):
    result = {}
    count = 0
    for job in pull_output(pth):
        #print(job)
        egy = get_energy(job)
        code = get_id(job)
        dn = dirname(job)
        bdn = dirname(dn)
        poscar = join(bdn, 'Test_spin')
        sf = join(poscar, 'POSCAR')
        struct = Structure.from_file(sf)
        ln = get_layer(struct)
        if ln == [0, 1, 0]:
            with open(r"./ln_false.log", 'a+') as f:
                f.write(sf + '\n')
            continue
        _, _, nn = ln
        #ln = 1
        direction = get_direction(struct)
        #direction = get_direction_from(sf, code)
        #{code: {'energy': float(egy.strip(' ')), 'direction': direction, 'ln': len(nn)}}
        #{code: {'energy': float(egy.strip(' ')), 'direction': direction['dn'], 'ln': direction['ln']}}
        result.update(
            {code: {'energy': float(egy.strip(' ')), 'direction': direction, 'ln': len(nn)}}
        )
        count += 1
    print('total: ', count)

    return result


def get_bulk_direction(code, dd):
    apply = dd.get(code)
    if apply is not None:
        return dd.get(code).get('direction')
    return


def run_bulk(pth, lsh):
    result = {}
    with open(lsh, 'r') as f:
        dd = json.load(f)
    for job in pull_output(pth):
        #print(job)
        dn = dirname(job)
        sf = join(dn, 'CONTCAR')
        struct = Structure.from_file(sf)
        egy = get_energy(job)
        code = get_id(job)
        ln = get_layer(struct)
        if ln == [0, 1, 0]:
           with open(r"./ln_false.log", 'a+') as f:
                f.write(sf + '\n')
           continue
        else:
            _, _, nn = ln
            n = len(nn)
        direction = get_bulk_direction(code, dd)
        if direction is None:
            continue
        plane_area = get_polygon_area_in_plane(struct, direction)
        result.update(
            {code: {'energy': float(egy.strip(' ')), 'ln': n, 'area': plane_area}}
        )

    return result


if __name__ == "__main__":
    pass

