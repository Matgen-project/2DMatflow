#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import defaultdict
from functools import lru_cache
from decimal import Decimal
import networkx as nx

import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class FloatRange:
    def __init__(self, start, end, step):
        self.s = Decimal(str(start))
        self.e = Decimal(str(end))
        self.step = Decimal(str(step))

    def __iter__(self):
        if self.s < self.e:
            pop_item = self.s
            while pop_item < self.e:
                yield float(pop_item)
                pop_item += self.step
        elif self.s == self.e:
            yield self.s
        else:
            pop_item = self.s
            while pop_item > self.e:
                yield float(pop_item)
                pop_item -= self.step


class MergeNestedLst:
    def __init__(self, nested_list):
        self._nl = nested_list
        self._sdt = self._start_a_dict()
        self._stl = self._start_a_lst()

    def _start_a_dict(self):
        dt = defaultdict(set)
        for i in self._nl:
            dt[i[0]] = dt[i[0]].union(set(i[1:]))

        return dt

    def _start_a_lst(self):
        stl = []
        for k, v in self._sdt.items():
            if k not in v:
                v.add(k)
            stl.append(list(v))

        return stl

    @lru_cache()
    def merge(self):
        G = nx.Graph()
        G.add_nodes_from(sum(self._stl, []))
        q = [[(s[i], s[i+1]) for i in range(len(s) - 1)] for s in self._stl]
        for i in q:
            G.add_edges_from(i)
        return sorted([list(i) for i in nx.connected_components(G)])


if __name__ == '__main__':
    pass

