#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import SupercellTransformation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Cell:
    def __init__(self, struct_file, supercell_scale=None):
        print("reading from cif: {}".format(struct_file))
        self.struct_file = struct_file
        _cell = Structure.from_file(self.struct_file)
        self._cell = SpacegroupAnalyzer(_cell).get_conventional_standard_structure()
        if supercell_scale is None:
            self._scs = np.eye(3) * [2, 2, 2]
        elif np.asarray(supercell_scale).shape != (3, 3):
            self._scs = np.eye(3) * np.asarray(supercell_scale)

    @property
    def _primitive_cell(self):
        return self._cell.get_primitive_structure()

    @property
    def primitive_structure(self):
        return self._primitive_cell

    @property
    def primitive_atoms_coords(self):
        return self._primitive_cell.cart_coords

    @property
    def primitive_lattice(self):
        return np.asarray(self._primitive_cell.lattice.as_dict().get("matrix"))

    @staticmethod
    def get_atomic_number(periodic_site_lst):
        psl = []
        for i in periodic_site_lst:
            psl.append(i.as_dict().get('species')[0].get('element'))

        return psl

    @property
    def primitive_atomic_number(self):
        return self._primitive_cell.atomic_numbers

    @property
    def primitive_atoms_number(self):
        return len(self.primitive_atomic_symbol)

    @property
    def primitive_atomic_symbol(self):
        return self.get_atomic_number(self._primitive_cell.sites)

    def get_distance_between_atoms_in_primitive(self):
        return self._primitive_cell.distance_matrix

    def _get_supercell(self):
        return SupercellTransformation(scaling_matrix=self._scs).apply_transformation(structure=self._primitive_cell)

    def __repr__(self):
        return "{\'Atoms\': {supercell: %r}\n{scale: %r}}" \
               % (self._get_supercell(), self._scs)

    @property
    def supercell_atoms_coords(self):
        return self._get_supercell().cart_coords

    @property
    def supercell_lattice(self):
        return self._get_supercell().lattice.as_dict().get("matrix")

    @property
    def supercell_atomic_number(self):
        return self._get_supercell().atomic_numbers

    @property
    def supercell_atoms_number(self):
        return len(self.supercell_atomic_number)

    @property
    def supercell_atomic_symbol(self):
        return self.get_atomic_number(self._get_supercell().sites)

    def get_distance_between_atoms_in_supercell(self):
        return self._get_supercell().distance_matrix


if __name__ == '__main__':
    pass
