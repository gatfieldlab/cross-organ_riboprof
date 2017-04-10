#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:20:48 2013

@author: barpat
"""

import sys
from collections import defaultdict
from accessories import utils
from argparse import ArgumentParser
import fileinput



parser = ArgumentParser()
parser.add_argument('files', nargs='*', help='specify input SAM files')
parser.add_argument('-o', '--output', help='specify output file. The default is stdout')
parser.add_argument('-c', '--command', help='', type=str, choices=['filter', 'passtru', 'basic', 'rrna'], default='filter')
parser.add_argument('-v', '--verbose', help='a little more verbosity', action='store_true')
parser.add_argument('-s', '--size-file', help='size (bp) file for references in SAM')
parser.add_argument('-k', '--keep-filtered', help='keep the filtered alignments in a separate SAM file', action='store_true')
parser.add_argument('-f', '--file-kept', help='filename for --keep-filtered. if empty default is "filtered.sam"', type=str, default='filtered.sam')
args = parser.parse_args()
if args.output and args.output != '-':
    sys.stdout = open(args.output, 'w')
if args.verbose:
    sys.stderr.write('Opened {} for output\n'.format(sys.stdout.name))
    
if args.keep_filtered:
    kept_file = open(args.file_kept, 'w')
    if args.verbose:
        sys.stderr.write('Opened {} for keeping filtered alignments\n'.format(args.file_kept))
        
#op_t = {'i': int, 'Z': str, 'f': float, 'A': str, 'H': str, 'B': str}

def process_cur_reads(cur_reads, line, command):
    if command == 'passtru':
        sys.stdout.write(line)
        return True
    parsed = line.strip().split('\t')
    bit_op = int(parsed[1])
    if bit_op & 0x4:
        if args.keep_filtered:
            kept_file.write(line)
        return cur_reads
    if command == 'basic':
        sys.stdout.write(line)
        return True
    # command == 'filter' or 'rrna'
    read_reverse = bool(bit_op & 0x10)
    if read_reverse and command == 'filter':
        if args.keep_filtered:
            kept_file.write(line)
        return cur_reads
    if command == 'rrna':
        try:
            ref_name = parsed[2]
            ref_len = rrna_sizes[ref_name]
        except KeyError as e:
            sys.stderr.write("SAM seems to include a references for which the size information is not available.\n{}".format(e.message))
            sys.exit(1)
    else:
        ref_name=None
        ref_len=None
    read_id = parsed[0]
#    ops = {}
    op_as = None
    for f in parsed[11:]:
        op,t,val = f.split(':')
        if op == 'AS':
            op_as = int(val)
            break
#        ops[s[0]] = op_t[s[1]](s[2])
    assert op_as != None, "NO 'AS' tag for this read???"
    if (cur_reads and read_id != cur_reads['read_id']) or not cur_reads:
        if cur_reads:
            if command == 'filter':
                sys.stdout.write(cur_reads['lines'])
            if command == 'rrna':
                hist = cur_reads['hist']
                if len(hist) > 1:
                    size_priority = sorted(range(len(hist)), key=lambda k: hist[k], reverse = True)
                    max_len = hist[size_priority[0]][0]
                    curlines = cur_reads['lines'].split('\n')
                    # Size filtered only -- no directions at this point
                    for refindex in size_priority:
                        if hist[refindex][0] == max_len:
                            sys.stdout.write("{}\n".format(curlines[refindex]))
                        else:
                            if args.keep_filtered:
                                kept_file.write("{}\n".format(curlines[refindex]))
                            else:
                                break
                    # The following is not a complete code, rather it contains
                    # hints for how to approach direction problem, e.g
                    # single read mapping on same ref in both directions with
                    # same match score...
                    '''
                    def_directions = defaultdict(set)
                    for ref in size_priority:
                        if hist[ref][0] == max_len:
                            def_directions[hist[ref][1]].add(hist[ref][2])
                    new_directions = {}
                    show_me = False #len(def_directions) > 1
                    print def_directions
                    for ref in def_directions:
                        if len(def_directions[ref]) == 2:
                            new_directions[ref] = False
                            show_me = True
                        else:
                            new_directions[ref] = list(def_directions[ref])[0]
                    if show_me: #sum(new_directions.values()):
                        print def_directions, new_directions, size_priority
                        print "===="
                    '''
                else:
                    sys.stdout.write(cur_reads['lines'])
        cur_reads = {'read_id': read_id, 'AS': op_as, 'hist': [(ref_len, ref_name, read_reverse)], 
                     'lines': line} # first hit is always kept if not 0x4 (x10)
    else: # cur_reads set and read_id is same -- anything else ??
        if op_as == cur_reads['AS']:
            cur_reads['lines'] += line
            if command == 'rrna':
                cur_reads['hist'].append((ref_len, ref_name, read_reverse))
        else:
            if args.keep_filtered:
                kept_file.write(line)
    return cur_reads

rrna_sizes = {}
if args.command == 'rrna':
    if not args.size_file:
        sys.stderr.write("Can't filter rRNA without the reference size information.")
        sys.exit(1)
    else:
        with open(args.size_file) as sizef:
            for line in sizef:
               parsed = line.strip().split('\t')
               rrna_sizes[parsed[0]] = int(parsed[1])
        if args.verbose:
            sys.stderr.write('Extracted length for {} references from {}\n'.format(len(rrna_sizes), args.size_file))
               
with utils.measureTime('Finished filtering {}'.format(','.join(args.files))):
    samfiles = fileinput.input(args.files)
    cur_reads = None
    while True:
        try:
            line = samfiles.next()
        except StopIteration:
            if args.command == 'filter' and cur_reads:
                sys.stdout.write(cur_reads['lines'])
            break
        if line[0] == '@':
            if not cur_reads:
                sys.stdout.write(line)
            if fileinput.isfirstline():
                sys.stderr.write("Processing input file: '{}'\n".format(fileinput.filename()))
            continue
        else:
            cur_reads = process_cur_reads(cur_reads, line, args.command)
                
