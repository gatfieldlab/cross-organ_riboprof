#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on  Tue Sep 17 17:23:09 2013
Modified on Sun Jan 11 14:29:30 2015

@author: barpat
"""

import sys, pickle

OK_FLAGS = ['OKEY', 'SNGL', 'MIXD']

class fileError(Exception):
    def __init__(self, m):
        self.message = m
class indexError(Exception):
    def __init__(self, tr_id, density_file):
        self.message = "<{}> could not be found in '{}'\n".format(tr_id, density_file)
class extractionError(Exception):
    def __init__(self, m):
        self.message = m

def get_index_tr(density_file, tr_ids, ok_flags=OK_FLAGS):
    try:
        density_file_handle = open(density_file, 'r')
    except IOError as e:
        msg = "Could not read density file!\nI/O error({0}): {1}\n".format(e.errno, e.strerror)
        raise fileError(msg)
    try:
        index_file = open(density_file+'.idx3', 'rb')
        tr_index = pickle.load(index_file)
    except IOError as e:
        msg = "Could not read index file!\nI/O error({0}): {1}\n".format(e.errno, e.strerror)
        raise fileError(msg)
    except pickle.UnpicklingError as e:
        msg = "Could not unpickle the index file - possibly wrong format!\nUnpickling error: {}\n".format(e.message)
        raise fileError(msg)
    except Exception as e:
        msg = "Could not read/unpickle the index file - unknown error! : {}\n".format(e.__repr__())
        raise fileError(msg)
    for tr_id in tr_ids:
        extract = ''
        try:
            index_pos = tr_index[tr_id]
        except KeyError:
            yield (tr_id, extract, indexError(tr_id, density_file))
            continue
        try:
            density_file_handle.seek(index_pos)
            for line in density_file_handle:
                parsed = line.strip().split('\t',3)
                cur_tr = parsed[0].split('|')[1]
                if cur_tr == tr_id:
                    if parsed[3] in ok_flags:
                        extract += line
                else:
                    break
            yield (tr_id, extract, None)
        except Exception as e:
            yield (tr_id, extract, extractionError("Could not extract '{}' indexed at '{}' - unknown error : {}\n".format(tr_id, index_pos, e)))

def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(argv) < 2:
        sys.stderr.write('usage: get_density DENSITY-FILE TR-IDs\n')
        return 1
    try:
        tr_extracts = get_index_tr(argv[1], argv[2:])
    except fileError as e:
        sys.stderr.write(e.message)
        return 1
    try:
        for tr_id, extract, error in tr_extracts:
            if error:
                if isinstance(error, indexError):
                    sys.stderr.write(error.message)
                    continue
                else:
                    raise error
            else:
                sys.stdout.write(extract)
        return 0
    except extractionError as e:
        sys.stderr.write(e.message)
        return 1

if __name__ == "__main__":
    sys.exit(main())

