#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections.abc import Iterable
from functools import lru_cache
# from concurrent.futures.thread import ThreadPoolExecutor
# from concurrent.futures import as_completed

from mat2d_pkg.para import settings
from mat2d_pkg.utils import MergeNestedLst

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class Bonded:
    def __init__(self, elements_symbol, distance_array, atomic_position, bonding_scheme=None, tolerance=None):
        self._es = elements_symbol
        self._da = distance_array
        self._ap = atomic_position
        self._bs = bonding_scheme
        self.tolerance = tolerance

    @property
    def uniq_symbol(self):
        us = []
        for i in self._es:
            if i not in us:
                us.append(i)
        return us

    @property
    def elements_radii(self):
        erd = {}
        for ele in self.uniq_symbol:
            erd.update(
                {
                    ele: settings.get(self._bs).get(ele)
                }
            )
        return erd

    def _scheme(self, tolerance):
        frame_dim, _ = self._da.shape
        bonding_pair = []

        for i in range(frame_dim):
            atom_a_radii = self.elements_radii.get(self._es[i])
            for j in range(frame_dim):
                if j <= i:
                    continue
                atom_b_radii = self.elements_radii.get(self._es[j])
                bonding = self._bonding_standard(ra=atom_a_radii,
                                                 rb=atom_b_radii,
                                                 distance=self._da[i][j],
                                                 delta=tolerance)

                if bonding:
                    bonding_pair.append([i, j])
        return {
            tolerance: bonding_pair
        }

    def _bonding_standard(self, ra, rb, distance, delta):
        try:
            if 'vdw' in self._bs.lower():
                return bool(distance < ra + rb - delta)
            return bool(distance < ra + rb + delta)
        except TypeError:
            raise Exception("This radii table may not contain this element, please check")

    def get_bonding_pair(self):
        if not isinstance(self.tolerance, Iterable):
            return self._scheme(self.tolerance)

        total_bonding_type = {}
        for i in self.tolerance:
            total_bonding_type.update(self._scheme(i))
        # ##########################################
        # the competition between threads causes slower calculations
        # with ThreadPoolExecutor(max_workers=5) as t:
        #    bd_obj_lst = []
        #    for i in self.tolerance:
        #        bd_obj_lst.append(t.submit(self._scheme, i))
        #    for f in as_completed(bd_obj_lst):
        #        data = f.result()
        #        total_bonding_type.update(data)

        return total_bonding_type

    @lru_cache()
    def get_cluster(self):
        bonding_pair = self.get_bonding_pair()
        cluster = {}
        # cluster_thread_lst = []
        # delta = list(bonding_pair.keys())
        # executor = ThreadPoolExecutor(max_workers=5)
        # for _, v in bonding_pair.items():
        #     cluster_thread_lst.append(executor.submit(Cluster(v).get_cluster_list, ))
        # for index, q in enumerate(as_completed(cluster_thread_lst)):
        #    cluster.update({
        #        delta[index]: q.result()
        #     })
        for i, v in bonding_pair.items():
            cluster.update({
                i: Cluster(v).get_cluster_list()
            })

        return cluster

    def get_cluster_number(self):
        return dict(zip(self.get_cluster().keys(),
                        [len(i) for i in self.get_cluster().values()]))

    @lru_cache()
    def get_atoms_number_in_cluster(self):
        delta = list(self.get_cluster().keys())
        num = {}
        for k, i in enumerate(self.get_cluster().values()):
            n = []
            for j in i:
                n.append(len(j))
            if not n:
                continue

            num.update(
                {
                    delta[k]: {"max": max(n),
                               "min": min(n),
                               "number": n}
                }
            )
        return num


class Cluster:
    def __init__(self, bonding_pair_lst):
        """
        :param bonding_pair_lst: A list of bonded atoms-> [[a, b],[a, c]]
        """
        self.bpl = bonding_pair_lst

    def _sort(self):
        """
        :return: make sure all atoms pair in order
        """
        return sorted(self.bpl, key=lambda x: x[0])

    @lru_cache()
    def get_cluster_list(self):
        return MergeNestedLst(self._sort()).merge()

    @property
    def cluster_number(self):
        return len(self.get_cluster_list())


if __name__ == "__main__":
    pass
