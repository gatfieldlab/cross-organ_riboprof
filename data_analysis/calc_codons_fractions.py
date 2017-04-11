#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 17:23:09 2013

@author: barpat
"""

import sys, re, operator
import argparse
import numpy as np
import get_density4 as gd
from collections import defaultdict


np.seterr(invalid='raise')

REGIONS = ('5utr','cds','3utr')
SUMMARY_OPERATIONS = ('mean', 'sum', 'median')
MAX_READ_LEN = 100
ASITE_OFF = {i: 15+int(i>30) for i in range(1,MAX_READ_LEN+1)}
OK_FLAGS = ['OKEY', 'SNGL']
REGION_PATTERN = "\*{,1}(" + "|".join(REGIONS)+ ")\*{,1}|\d+|[\+\-]"
REGION_MATCH = re.compile(REGION_PATTERN)
OPS = {'+': operator.add, '-': operator.sub}
class densityError(Exception):
    def __init__(self, m):
        self.message = m


class Region(object):
    '''
    [ar]{,1}\d+
    '''
    def __init__(self, region):
        self.parse_words = []
        parsed_regions = region.split(":")
        for parsed_region in parsed_regions:
            self.parse_words.append([])
            match_iter = REGION_MATCH.finditer(parsed_region)
            for match in match_iter:
                self.parse_words[-1].append(match.group(0))

    def get_limits(self, tr_info):
        limits = []
        for sentence in self.parse_words:
            limits.append([])
            for word in sentence:
                if word.isdigit():
                    limits[-1].append(int(word))
                elif word in ['+','-']:
                    limits[-1].append(word)
                else:
                    left, right = word[0] == '*', word[-1] == '*'
                    if left:
                        word = word[1:]
                    if right:
                        word = word[:-1]
                    if word in REGIONS:
                        if word == '5utr':
                            word_limits = (0, tr_info['cds_start'])
                        if word == 'cds':
                            word_limits = (tr_info['cds_start'], tr_info['cds_end'])
                        if word == '3utr':
                            word_limits = (tr_info['cds_end'], tr_info['tr_len'])
                    else:
                        word_limits = (0, tr_info['tr_len'])
                    if left:
                        limits[-1].append(word_limits[0])
                    if right:
                        limits[-1].append(word_limits[1])
                    if not (left or right):
                        limits[-1].extend(word_limits)
        result = []
        for limit in limits:
            if not (type(limit[0]) is int and type(limit[-1]) is int):
                print(limits)
                raise Exception('Hoolllalalalal')
            while '-' in limit or '+' in limit:
                o_i, op = next((i,o) for i,o in enumerate(limit) if o == '-' or o == '+')
                try:
                    res = OPS[op](limit[o_i-1], limit[o_i+1])
                except KeyError:
                    print(limit)
                    exit()
                limit = limit[:max(0, o_i-1)] + [res] + limit[min(len(limit), o_i+2):]
            result.append(limit)
        lmin = min([min(l) for l in result])
        lmax = max([max(l) for l in result])
        return (lmin, lmax)

class DensityParser(object):
    '''
    A class
    '''

    def __init__(self, tr_db, tr_ids, density_files, region):
        self.tr_db = tr_db
        self.tr_ids = tr_ids
        self.density_files = []
        for density_file in density_files:
            try:
                self.density_files.append(gd.get_index_tr(density_file, self.tr_ids))
            except gd.fileError as e:
                raise densityError(e.message)

        self.region = Region(region)

    def get_densitybin(self, asite_off=ASITE_OFF):
        while True:
            extracts = []
            is_eof = []
            for densities in self.density_files:
                try:
                    extracts.append(next(densities))
                except gd.fileError as e:
                    raise densityError(e.message)
                except StopIteration:
                    is_eof.append(True)
                else:
                    is_eof.append(False)
            if all(is_eof):
                raise StopIteration
            elif any(is_eof):
                raise densityError("Some density files reached EOF before others.\n")
            else:
                tr_ids, extracts, errors = list(zip(*extracts))
                is_error = [error is not None for error in errors]
                trs_same = len(set(tr_ids)) == 1
                if not trs_same:
                    raise densityError("Some density files return TR-IDs out-of-sync.\n")
                if is_error:
                    pass # for now
                cur_err = []
                result = []
                tr_id = tr_ids[0]
                try:
                    tr_info = self.tr_db[tr_id]
                except KeyError:
                    yield (None, None, None, None, [densityError("<{}> could not be found in DB\n".format(tr_id))])
                else:
                    slice_limits = self.region.get_limits(tr_info)
                    if slice_limits[0] < 0 or slice_limits[1] > tr_info['tr_len']:
                        yield (None, None, None, None, [densityError("Slice {} can't be applied to <{}>\n".format(slice_limits, tr_id))])
                    for i, extract in enumerate(extracts):
                        density = []
                        if isinstance(errors[i], gd.indexError):
                            cur_err.append(densityError(errors[i].message))
                        elif extract == '':
                            cur_err.append(densityError("<{}> returned empty list from density file\n".format(tr_id)))
                        else:
                            ids = extract.split('\t',1)[0]
                            cur_gene, cur_tr = ids.split('|')
                            if not cur_tr == tr_id:
                                cur_err.append(densityError("TR-ID from density extract does not match the returned TR-ID\n"))
                            else:
                                for line in extract.strip().split('\n'):
                                    parsed = line.split('\t')
                                    try:
                                        read_pos = int(parsed[1])
                                        read_len = int(parsed[2])
                                    except IndexError:
                                        cur_err.append(densityError("Can't parse the density data returned by get_density from '{}'!\n".format(tr_id))) #changed from .format(density_file))) (vio)
                                        break
                                    except ValueError:
                                        cur_err.append(densityError("Can't interpret density data returned by get_density!\n"))
                                        break
                                    try:
                                        offset = asite_off[read_len]
                                    except KeyError:
                                        sys.stderr.write("Found a read with a size outside of specified range of 1 to {}\nWill continue setting size to {}!\n".format(MAX_READ_LEN, MAX_READ_LEN))
                                        offset = asite_off[MAX_READ_LEN]
                                    #print('adding:',read_pos, offset) # debug
                                    density.append(read_pos + offset)
                            cur_err.append(None)
                        if density:
                            binned_density=np.bincount(density, minlength=tr_info['tr_len'])
                        else:
                            binned_density=np.zeros(tr_info['tr_len'])
                        result.append(binned_density[slice_limits[0]:slice_limits[1]])
                        #print('outgoing', len(density), tr_len, tr_db[cur_tr], sep='\n') #debug
                    assert len(result) == len(cur_err)
                    yield (cur_gene, cur_tr, slice_limits, result, cur_err)

def summarize_densitybins(densitybins, operation='sum', **kwargs):
    if operation in SUMMARY_OPERATIONS:
        try:
            op = getattr(np, operation)
        except AttributeError:
            sys.stderr.write("{} is not a numpy attribute??".operation)
            sys.exit(1)
        try:
            res = op(np.array(densitybins), axis=0, **kwargs)
        except TypeError as e:
            sys.stderr.write("Numpy error during '{}' operation: {}".format(operation, str(e)))
            sys.exit(1)
        else:
            return res
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', metavar='density_file', type=str)
    parser.add_argument('-t', '--tr_ids', type=str, default='tr_ids.txt')
    parser.add_argument('-r', '--region', type=str)
    parser.add_argument('-c', '--cds', type=str, default='prepared_cds.txt')
    parser.add_argument('-l', '--max_read_len', type=int, default=MAX_READ_LEN, help='specify max length of reads, default={}'.format(MAX_READ_LEN))

    args = parser.parse_args()
    tr_ids = []
    try:
        with open(args.tr_ids) as tr_file:
            for line in tr_file:
                tr_ids.append(line.strip())
    except IOError:
        sys.stderr.write('Could not read the TR-IDs from {}\n'.format(args.tr_ids))
        sys.exit(1)

    tr_db = {}
    try:
        with open(args.cds) as cds_file:
            for line in cds_file:
                parsed = line.strip().split('\t')
                if parsed[2] in tr_ids:
                    (tr_len, cds_start, cds_end) = (int(parsed[i]) for i in (5,6,7))
                    tr_db[parsed[2]] = {'tr_len': tr_len, 'cds_start': cds_start, 'cds_end': cds_end} 
    except IOError:
        sys.stderr.write('Could not read the CDS models from {}\n'.format(args.cds))
        sys.exit(1)

    density_extractor = DensityParser(tr_db, tr_ids, args.files, args.region)
    for gid, tr_id, slice_limits, extracts, errors in density_extractor.get_densitybin():
        if len(errors) == 1 and isinstance(errors[0], densityError):
            sys.stderr.write(errors[0].message)
        else:
            extracts_fractions =[]
#            extracts_sums =[]
            for extract in extracts:
                extracts_fractions.append(extract/max(1,np.sum(extract)))
#                extracts_sums.append(extract)
            summarized = summarize_densitybins(extracts_fractions,'mean') # mean of time points
#            summarized = summarize_densitybins(extracts_sums,'sum') # sum of counts (not fractions) of all time points
            #print(gid, tr_id, summarized)
            for i,f in enumerate(summarized):
                print('{}\t{}\t{}'.format(tr_id,i,f)) #instead of print('{}\t{}'.format(i,f)) to add transcript ID. (vio)

if __name__ == "__main__":
    main()
#         densitybins = get_densitybin(tr_ids, density_file)
#         while True:
#             try:
#                 densitybin = next(densitybins)
#                 if isinstance(densitybin, densityError):
#                     raise densitybin
#             except StopIteration:
#                 break
#             except densityError as e:
#                 sys.stderr.write(e.message)
#             else:
#                 cur_tr = densitybin[0][1]
#                 bp = 0
#                 read_sum = np.sum(densitybin[1])
#                 if read_sum == 0:
#                     read_sum = 1
#                 for count in densitybin[1]:
#                     bp += 1
#                     codon = (bp - 1) // 3 + 1
#                     frame = (bp - 1) % 3 + 1
#                     try:
#                         print(cur_tr, bp, codon, frame, count, count/read_sum, sep='\t')
#                     except FloatingPointError:
#                         sys.stderr.write("{}\n".format('\t'.join(['FlotaingPointError', cur_tr, str(codon), str(frame), str(count), str(read_sum)])))
#                 #print(tr_db[cur_tr]) # debug



