# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:20:48 2013

@author: barpat
"""

import sys
from accessories import utils
from collections import defaultdict
from argparse import ArgumentParser
from itertools import islice
import HTSeq

genes = {}
protein_coding = set()
counter = defaultdict(int)
multigene_counter = defaultdict(int)
gene_count = 0
a_site_readlen = 30
a_site_5offset = (15, 16)
a_site_3offset = 18
map_score = {'cds': 8, '5utr': 4, '3utr': 2, 'db_problem': 1}
sorted_map_keys = sorted(map_score, key=map_score.get, reverse=True)
multiread_count = 0
nonprot_count = 0

track_multi = False
count_multi = False

def set_factory():
    return defaultdict(set)
def int_factory():
    return defaultdict(int)

parser = ArgumentParser()
parser.add_argument('files', nargs='*', help="specify input BAM files")
parser.add_argument('-o', '--output', help="specify output file. The default is stdout")
parser.add_argument('-g', '--gtf', help="specify the GTF file")
parser.add_argument('-s', '--split', help="produce split results for different size ranges", action="store_true")
parser.add_argument('-r', '--ranges', help="size ranges for splitting", type=str, default="26,36")
parser.add_argument('-e', '--expressed', help="specify expressed Gene/Transcript ids")
parser.add_argument('-t', '--type', help="specify the type of reads, def: experiment", default='experiment', choices=['experiment','faux'])
parser.add_argument('-c', '--command', help="specify what you want me to do: 'count', '5prime', 'both'", default='count', choices=['count', '5prime', 'both'])
args = parser.parse_args()
if args.output and args.output != '-':
    sys.stdout = open(args.output, 'w')

with utils.measureTime('CDS gtf models are loaded'):
    with open('all_protein_coding.v38.84.txt') as protlist:
        for ensmusg in protlist:
            protein_coding.add(ensmusg.strip())
    if args.gtf == None:
        sys.stderr.write('Need a GTF file to count the reads against, use -g/--gtf option.\n')
        sys.exit(1)
    with open(args.gtf) as cds:
        for line in cds:
            parsed = line.strip().split('\t')
            if not parsed[0] in genes:
                genes[parsed[0]] = {'type': parsed[1], 'order': gene_count,'tr':{}}
                gene_count += 1
            genes[parsed[0]]['tr'][parsed[2]] = (parsed[3], parsed[4],
                                               int(parsed[5]), int(parsed[6]),
                                               int(parsed[7]))
expressed = {}
if args.expressed:
    with open(args.expressed) as expfile:
        for line in expfile:
            parsed = line.strip().split('\t')
            expressed[parsed[1]] = parsed[0]
if args.split:
    alt_counter = defaultdict(int)
    split_ranges = [int(val) for val in args.ranges.split(',')]
else:
    split_ranges = []

split_size = len(split_ranges) + 1
counter = [defaultdict(int) for i in range(split_size)]
count_by_gene = [defaultdict(int_factory) for i in range(split_size)]
readpos_by_gene = defaultdict(dict)
masked_region = defaultdict(set)

#problem_sets = []
#problem_genes = set([])

if args.command in ('5prime', 'both'):
    handle_5prime = open("{}_5prime_count.txt".format(args.output), 'w')
def process_cur_alns(cur_read, cgenes, split_pos=0):
    global problem_genes
    counter[split_pos]['valid_reads'] += 1
    to_add = None
    if len(cgenes) == 1:
        to_add = cgenes.popitem()
        read_code = 'OKEY'
#        count_by_gene[cgenes.popitem()[0]] += 1
    else:
        expressedg = []
        protcoding = []
        for cgene, trs in cgenes.items():
            etrs = set([tr for tr in trs if tr[0] in expressed])
            if etrs:
                expressedg.append((cgene, etrs))
                if cgene in protein_coding:
                    protcoding.append((cgene, etrs))
        if expressedg:
            if not protcoding:
                counter[split_pos]['nonprot'] += 1
                read_code = 'NPRT'
            elif len(protcoding) > 1:
                if count_multi:
                    multigene_counter[tuple(sorted(protcoding))] += 1
                if track_multi:
                    for cgene in protcoding:
                        masked_region[cgene].add(readpos_by_gene[cgene][cur_read])
                counter[split_pos]['multiread'] += 1
                read_code = 'MULT'
            else:
                #print "Problem: {}, {}".format(protcoding, cgenes)
                if len(expressedg) == 1:
                    counter[split_pos]['single_expressed_prot'] += 1
                    to_add = expressedg.pop()
                    read_code = 'SNGL'
                else:
#                    xgenes = zip(*expressedg)[0]
#                    problem_genes |= set(xgenes)
#                    found = False
#                    for pset in problem_sets:
#                        for cgene in xgenes:
#                            if cgene in pset:
#                                found = True
#                                for ogene in xgenes:
#                                    pset.add(ogene)
#                    if not found:
#                        problem_sets.append(set(xgenes))
                    counter[split_pos]['nonprot_problem'] += 1
                    read_code = 'MIXD'
                    #to_add = protcoding.pop()
        else:
            counter[split_pos]['multi_nonexpressed'] += 1
            read_code = 'MLTN'

    if to_add:
        shall_count = True
        if args.type == 'faux':
            source_id = cur_read[:18]
            if not to_add[0] == source_id:
                #print "{} is not {}".format(cur_read, to_add[0])
                counter[split_pos]['false_faux'] += 1
                shall_count = False
        if shall_count:
            is_cds = set([tr[1] for tr in to_add[1]])
            cds_score = sum([map_score[i] for i in is_cds])
            for what_to_count in sorted_map_keys:
                if cds_score >= map_score[what_to_count]:
                    break
            count_by_gene[split_pos][to_add[0]][what_to_count] += 1
            to_write = dict((to_add,))
    else:
        to_write = cgenes
    if not args.type == 'faux' and args.command in ('5prime', 'both'):
        for alignment in to_write.items():
            for tr in alignment[1]:
                handle_5prime.write("{}|{}\t{}\t{}\t{}\n".format(alignment[0], tr[0], tr[3], tr[2], read_code))


with utils.measureTime('Finished parsing {}'.format(','.join(args.files)), unit='minute'):
    for bamfile in args.files:
        try:
            bam_reader = HTSeq.BAM_Reader(bamfile)
        except:
            continue
        sys.stderr.write("Processing '{}'\n".format(bamfile))
        bam_iter = islice(bam_reader, None)
        cur_read = None
        genes_by_read = defaultdict(set)
        split_pos = 0
        while True:
            try:
                aln = bam_iter.next()
            except StopIteration:
                if cur_read:
                    process_cur_alns(cur_read, genes_by_read)
                break
            read_id = aln.read.name
            if cur_read and not read_id == cur_read:
                process_cur_alns(cur_read, genes_by_read, split_pos)
                genes_by_read = defaultdict(set)
            cur_read = read_id
            read_len = len(aln.read.seq)
            split_pos = 0
            if args.split:
                for split_border in split_ranges:
                    if read_len < split_border:
                        break
                    else:
                        split_pos += 1
            gene_id, tr_id = aln.iv.chrom.split('|',1)
            if args.type == 'faux' and not tr_id in expressed:
                continue
            start = aln.iv.start
            location = '3utr'
            try:
                if start <= max(genes[gene_id]['tr'][tr_id][4]-a_site_3offset,0):
                    location = 'cds'
                if start < max(genes[gene_id]['tr'][tr_id][3]-a_site_5offset[read_len > a_site_readlen],0):
                    location = '5utr'
            except KeyError:
#                        sys.stderr.write('{}|{} is not in the DB\n'.format(gene_id, tr_id)
                location= 'db_problem'
            genes_by_read[gene_id].add((tr_id, location, read_len, start))
                #if track_multi:
                #    cigar = align_utils.parse_cigar(result[5])
                #    tlen = 0
                #    for op in cigar:
                #        if op[-1] in ['M','D','N']:
                #            tlen += int(op[:-1])
                #    readpos_by_gene[gene_id][read_id] = (start, start+tlen)

split_names = ['1_']
for split_border in split_ranges:
    split_names[-1] += str(split_border-1)
    split_names.append(str(split_border)+'_')
split_names[-1] += 'inf'
if args.command in ('count', 'both'):
    for split_pos in range(len(count_by_gene)):
        with open("{}_split_{}.txt".format(args.output, split_names[split_pos]), 'w') as outfile:
            for gene in sorted(genes.keys(), key=lambda x: genes[x]['order']):
                outfile.write("{0}\t{1[5utr]}\t{1[cds]}\t{1[3utr]}\t{1[utr]}\t{1[db_problem]}\n".format(gene, count_by_gene[split_pos][gene]))
if count_multi:
    for genes, count in sorted(multigene_counter.items(), key=lambda (k,v): v):
        print "{}: {}".format(genes, count)
if track_multi:
    for gene in masked_region:
        print gene, masked_region[gene]
if args.command == 'count':
    sys.stdout.write("{}\n".format(counter.__repr__()))
#print "{} problematic multiread genes were found.".format(len(problem_genes))
#for pset in problem_sets:
#    print pset
