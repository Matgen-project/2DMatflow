#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from os.path import *
import shutil
# import fcntl
import json

from .layer import Layer
from .layer import ConstructLayer


def get_struct_files(struct_dir):
    for root, _, files in os.walk(struct_dir):
        for f in files:
            if splitext(f)[-1] == '.cif' or 'POSCAR' in f.upper():
                yield join(root, f)


def split_task(task_dir, top=20):
    count, name = 0, 1
    for item in get_struct_files(task_dir):
        des_dir = join(task_dir, str(name))
        if not exists(join(des_dir)):
            os.makedirs(des_dir)
        shutil.move(item, join(des_dir, basename(item)))
        count += 1
        if count == top:
            name += 1
            top += top


def safe_write(result, fn):
    with open(fn, "a+") as f:
        #fcntl.flock(f, fcntl.LOCK_EX)
        f.write(result)
        #fcntl.flock(f, fcntl.LOCK_UN)


def safe_to_json(file_name, result):
    with open(file_name, "w") as f:
        #fcntl.flock(f, fcntl.LOCK_EX)
        json.dump(result, f, indent=4)
        #fcntl.flock(f, fcntl.LOCK_UN)


if __name__ == "__main__":
    pass
