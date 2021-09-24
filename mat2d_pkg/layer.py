#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import math
from os.path import *
from collections.abc import Iterable
from functools import lru_cache
from pymatgen.io.cif import CifWriter

from mat2d_pkg.atoms import Bonded
from mat2d_pkg.mat import *
from mat2d_pkg.utils import FloatRange
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

np.set_printoptions(precision=4, suppress=True)


class Layer:
    def __init__(self, file_path, box_scale=None, bond_type=None, tolerance=None):
        self.ofn = basename(file_path)
        if box_scale is None:
            box_scale = [2, 2, 2]
        self._bx = box_scale
        self._cell = Cell(struct_file=file_path, supercell_scale=box_scale)
        if bond_type is None:
            bond_type = 'vdw_radii'
        if isinstance(tolerance, Iterable):
            try:
                tolerance = list(FloatRange(*tolerance))
            except TypeError:
                raise Exception("If the threshold is a range,"
                                " the format must be (start, end, step)")
        elif tolerance is None:
            tolerance = list(FloatRange(1.1, 1.6, 0.1))
        else:
            tolerance = [tolerance, ]

        self.bt = bond_type
        self.tc = tolerance
        if isinstance(self.tc, Iterable):
            self.long = len(tolerance)
        else:
            self.long = 1

    @property
    def atoms_number(self):
        return {
            "primitive": self._cell.primitive_atoms_number,
            "supercell": self._cell.supercell_atoms_number
        }

    @lru_cache()
    def _cluster_obj(self):
        return {
            "primitive": Bonded(elements_symbol=self._cell.primitive_atomic_symbol,
                                distance_array=self._cell.get_distance_between_atoms_in_primitive(),
                                atomic_position=self._cell.primitive_atoms_coords,
                                bonding_scheme=self.bt,
                                tolerance=self.tc),
            "supercell": Bonded(elements_symbol=self._cell.supercell_atomic_symbol,
                                distance_array=self._cell.get_distance_between_atoms_in_supercell(),
                                atomic_position=self._cell.supercell_atoms_coords,
                                bonding_scheme=self.bt,
                                tolerance=self.tc)
        }

    @property
    @lru_cache()
    def primitive_cluster(self):
        return self._cluster_obj().get("primitive").get_cluster()

    @property
    @lru_cache()
    def supercell_cluster(self):
        return self._cluster_obj().get("supercell").get_cluster()

    @property
    @lru_cache()
    def primitive_cluster_number(self):
        return self._cluster_obj().get("primitive").get_cluster_number()

    @property
    @lru_cache()
    def supercell_cluster_number(self):
        return self._cluster_obj().get("supercell").get_cluster_number()

    @property
    @lru_cache()
    def primitive_cluster_atoms_number(self):
        return self._cluster_obj().get("primitive").get_atoms_number_in_cluster()

    @property
    @lru_cache()
    def supercell_cluster_atoms_number(self):
        return self._cluster_obj().get("supercell").get_atoms_number_in_cluster()

    @staticmethod
    def _mat2d_parser(primitive_cell_max, primitive_cell_min, supercell_max, supercell_min, supercell_atoms_number):

        def mat_type_judge(ratio_number):
            return "2D vdW solid" \
                if ratio_number == 4 else "1D vdW solid" \
                if ratio_number == 2 else "Intercalated 1D/2D"

        if primitive_cell_min == 1:
            return "Exclude! has unbonded atoms!"
        ratio = supercell_max / primitive_cell_max
        if primitive_cell_min == supercell_min:
            return "Exclude! has unbonded molecules!," \
                   "Molecule solid: {}".format(mat_type_judge(ratio))
        else:
            if supercell_max == supercell_atoms_number:
                return "Exclude! 3D solid"
            else:
                return mat_type_judge(ratio)

    def get_classification_results(self):
        supercell_atoms_number = self.atoms_number.get('supercell')
        cr = {}
        for i in range(self.long):
            pc_max, pc_min = list(map(self.primitive_cluster_atoms_number[self.tc[i]].get,
                                      ['max', 'min']))
            sc_max, sc_min = list(map(self.supercell_cluster_atoms_number[self.tc[i]].get,
                                      ['max', 'min']))
            mat_parser = self._mat2d_parser(pc_max, pc_min, sc_max, sc_min, supercell_atoms_number)

            cr.update({
                self.tc[i]: mat_parser
            })

        return cr

    @staticmethod
    def all(result_dict, mat_type):
        # '2D vdW solid'
        for _, values in result_dict.items():
            if values != mat_type:
                return False
        return True

    @staticmethod
    def any(result_dict, mat_type):
        # '2D vdW solid'
        for _, values in result_dict.items():
            if values == mat_type:
                return True
        return False

    @classmethod
    def judge_result(cls, result_dict):
        if cls.all(result_dict, '2D vdW solid'):
            return '2D_vdW_solid'
        if cls.all(result_dict, '1D vdW solid'):
            return '1D_vdW_solid'
        if cls.any(result_dict, '2D vdW solid') & cls.any(result_dict, '1D vdW solid'):
            return 'Weird_vdW_solid'
        if cls.any(result_dict, '2D vdW solid'):
            return '2D_vdW_candidate_solid'
        if cls.any(result_dict, '1D vdW solid'):
            return '1D_vdW_candidate_solid'

        return 'Exclude_mats'


class ConstructLayer(Layer):
    def __init__(self, file_path, box_scale, bond_type=None, tolerance=None):
        self.ofn = basename(file_path)
        if box_scale is None:
            box_scale = [3, 3, 3]
        self._bx = box_scale
        self._cell = Cell(struct_file=file_path, supercell_scale=box_scale)
        super(ConstructLayer, self).__init__(file_path, box_scale, bond_type, tolerance)
        self.d = self.get_direction()
        print(self.d)
        self.li = math.floor(self.get_specific_layer_number(0.45) / 2)

    @lru_cache()
    def get_cluster(self):
        return self.supercell_cluster

    def _run_cluster_loop(self, sort=True):
        all_cluster_obj = {}
        total_cluster = self.get_cluster()
        for index, cl in total_cluster.items():
            t = self._get_cluster_obj(index, cl, sort)
            all_cluster_obj.update(t)
        return all_cluster_obj

    def _get_cluster_obj(self, delta, cluster, sort=True):
        cluster_obj_lst = []
        for ic in cluster:
            cluster_coords_lst, cluster_symbol_lst = [], []
            for index in ic:
                cluster_coords_lst.append(self._cell.supercell_atoms_coords[index])
                cluster_symbol_lst.append(self._cell.supercell_atomic_symbol[index])
            ics = Structure(lattice=self._cell.supercell_lattice,
                            coords=np.asarray(cluster_coords_lst),
                            species=cluster_symbol_lst,
                            coords_are_cartesian=True)
            cluster_obj_lst.append(ics)

        # if sort:
        # return {
        #     delta: self._sort(cluster_obj_lst)
        # }

        return {
            delta: cluster_obj_lst
        }

    @property
    def cluster_dict(self, sort=True):
        return self._run_cluster_loop(sort)

    @lru_cache()
    def get_specific_cluster(self, tolerance):
        # from pprint import pprint
        # pprint(self.cluster_dict)
        return self.cluster_dict.get(tolerance)

    def _sort(self, cluster_obj_list):
        length = list(map(np.linalg.norm, self._cell.supercell_lattice))
        max_length_axis = length.index(max(length))

        return sorted(cluster_obj_list,
                      key=lambda x: max(x.cart_coords[:, max_length_axis]))

    def __getitem__(self, item, tolerance):
        return self.get_specific_cluster(tolerance)[item]

    def get_specific_layer_number(self, tolerance):
        return len(self.get_specific_cluster(tolerance))

    def get_layer_cif(self, save_path, tolerance, layer_index=None, all_layer=False, centralize=False, vacuum=30):
        if not exists(save_path):
            os.makedirs(save_path)
        fn, _ = splitext(self.ofn)
        job_path = join(save_path, fn)
        layer_index = self.li
        if not exists(job_path):
            os.makedirs(job_path)
        try:
            if all_layer:
                if not centralize:
                    for i, struct in enumerate(self.get_specific_cluster(tolerance)):
                        des_fn = join(job_path, str(tolerance) + "-{}-layer-{}.cif".format(fn, i))
                        self.get_cif_file(structure=struct, cif_name=des_fn)
                else:
                    logger.warning("If choice all layer centralized, "
                                   "All cifs may be the same")
                    for i, struct in enumerate(self.get_specific_cluster(tolerance)):
                        des_fn = join(job_path,
                                      str(tolerance) + "-{}-centralized-{}.cif".format(fn, i))
                        self.get_cif_file(structure=self.centralize(self.d, struct, vacuum),
                                          cif_name=des_fn)

            else:
                if not centralize:
                    des_fn = join(job_path,
                                  str(tolerance) + "-{}-layer-{}.cif".format(fn, layer_index))
                    struct = self.get_specific_cluster(tolerance)[layer_index]

                else:
                    des_fn = join(job_path,
                                  str(tolerance) + "-{}-centralized-{}.cif".format(fn, layer_index))
                    struct = self.centralize(self.d, self.get_specific_cluster(tolerance)[layer_index], vacuum)
                self.get_cif_file(structure=struct, cif_name=des_fn)
        except TypeError:
            raise Exception("tolerance value not found !")
        return

    @staticmethod
    def get_cif_file(structure, cif_name):
        return CifWriter(struct=structure).write_file(cif_name)

    @lru_cache()
    def get_layer(self, layer_index, tolerance, all_layer=False):
        layer_index = self.li
        if all_layer:
            return self.get_specific_cluster(tolerance)

        return self.get_specific_cluster(tolerance)[layer_index]

    def get_primitive_layered_mat(self, save_path,
                                  tolerance,
                                  layer_index,
                                  centralize=False, vacuum=30):
        layer_index = self.li
        one_layer = self.get_specific_cluster(tolerance)[layer_index]
        fn, suffix = splitext(basename(self.ofn))
        job_path = join(save_path, fn)
        if not exists(job_path):
            os.makedirs(job_path)
        pol = one_layer.get_primitive_structure()
        cn = '{}-layer-primitive-{}{}'.format(fn, layer_index, suffix)
        cif_name = join(job_path, cn)
        if centralize:
            pol = self.centralize(self.d, pol, vacuum).get_primitive_structure()
            cif_name = join(job_path,
                            str(tolerance) + "-centralized-" + cn)
        return self.get_cif_file(structure=pol, cif_name=cif_name)

    @lru_cache()
    def get_direction(self):
        dft_cl = self.get_specific_cluster(tolerance=0.45)[0]
        tm = [[2, 1, 1], [1, 2, 1], [1, 1, 2]]
        count = 0
        cnl = []
        for i in tm:
            sc = SupercellTransformation(scaling_matrix=i).apply_transformation(structure=dft_cl)
            da = Bonded(Cell.get_atomic_number(sc.sites), sc.distance_matrix,
                        sc.cart_coords, 'covalent_radii', 0.45)
            cn = da.get_cluster_number().get(0.45)
            cnl.append(cn)
        print(cnl)

        if cnl.count(2) > 1 or cnl.count(2) == 0:
            with open(r"./Non_parallel.struct", 'a+') as f:
                f.write(self.ofn + '\n')
            return

        return cnl.index(2)
    @staticmethod
    def centralize(direction, struct, vacuum=30):
        a, b, c = [struct.lattice.a,
                   struct.lattice.b,
                   struct.lattice.c]
        lattice_array = np.asarray(struct.as_dict().get('lattice').get('matrix'))
        # print('vacuum will be build on {} axis'.format(['a', 'b', 'c'][direction]))
        # change lattice
        ml = [a, b, c][direction]
        ml_array = lattice_array[direction]
        positions = struct.cart_coords
        thickness = max(positions[:, direction]) - min(positions[:, direction])
        vector_scale = ((2 * vacuum) + thickness) / ml
        vacuum_vector = ml_array * vector_scale
        lattice_array[direction] = vacuum_vector
        # get atoms position
        # print(positions)
        # group_mid_point = positions.mean(axis=0)
        # print(mp)
        # vertical/ horizontal shift
        # looking for the center of the group
        # print(positions)
        gmp = []
        for val in positions.transpose():
            mid = ((max(val) - min(val)) / 2) + min(val)
            gmp.append(mid)
        group_mid_point = np.asarray(gmp)
        box_mid_point = lattice_array.sum(axis=0) / 2
        diff = group_mid_point - box_mid_point
        new_coords = []
        for index, coords in enumerate(positions):
            new_coord = coords - diff
            new_coords.append(new_coord)
        return Structure(coords=np.asarray(new_coords), species=struct.species, lattice=lattice_array,coords_are_cartesian=True)

if __name__ == '__main__':
    pass
