# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 11:01:26 2013

@author: barpat
"""

import HTSeq
import sys, argparse
import quicksect
from accessories import utils

parser = argparse.ArgumentParser(description='Prepares transcript models for ribosome profiling analysis')
parser.add_argument('-v','--verbose', help="A little more verbosity helps humanity ;)", action='store_true')
parser.add_argument('-d', '--delimiter', help="delimiter used in output files. Default: '\\t'", type=str,
                    default='\t')
parser.add_argument('-l', '--list', help="file containing a list of GeneID/TranscriptIDs for filtering the GTF transcripts", type=str)
parser.add_argument("gtf", help="the GTF containing all models/exons/CDSs", type=str)

args = parser.parse_args()

try:
    gff_file = HTSeq.GFF_Reader(args.gtf)
except IOError:
    sys.stderr.write("Could not find '{}'\n".format(args.gtf))
    sys.exit(1)

expressed = {}
if args.list:
    try:
        with open(args.list) as expressed_prots:
            next(expressed_prots)
            for line in expressed_prots:
                parsed = line.strip().split('\t')
                if not parsed[0] in expressed:
                    expressed[parsed[0]] = {}
                expressed[parsed[0]][parsed[1]] = parsed[2:]            
    except IOError:
        sys.stderr.write("Could not find '{}'\n".format(args.list))
        sys.exit(1)

with utils.measureTime('GTF file has been parsed'):
    loci = {}
    for feature in gff_file:
        if not (feature.type in ('exon', 'CDS') and feature.source in ('protein_coding', 'nonsense_mediated_decay', 'non_stop_decay')):
            continue
        gene_id = feature.attr['gene_id']
        tr_id = feature.attr['transcript_id']
        if not gene_id in loci:
            loci[gene_id] = {'name': feature.attr['gene_name'], 'biotype': feature.attr['gene_biotype'], 'tr':{}}
        if not tr_id in loci[gene_id]['tr']:
            loci[gene_id]['tr'][tr_id] = {'exon':[], 'CDS': []}
        loci[gene_id]['tr'][tr_id][feature.type].append(feature.iv)


def make_transcript(tr, what='start_stop'): #, tr_id):
    transcript = zip(*[(iv.start, iv.end) for iv in tr['exon']])
    cds = zip(*[(iv.start, iv.end) for iv in tr['CDS']])
    mRNA_length = sum([transcript[1][i] - transcript[0][i] for i in range(len(transcript[0]))])
    cds_length = sum([cds[1][i] - cds[0][i] for i in range(len(cds[0]))])
    try:
        if len(transcript[0]) == 1:
            forward = transcript[1] > transcript[0]
        else:
            forward = transcript[0][-1] > transcript[1][0]
        if forward:
            if len(transcript[0]) == 1:
                first_tr = 0
                last_tr = 0
            elif len(cds[0]) == 1:
                for i in range(len(transcript[0])):
                    if transcript[1][i] >= cds[1][0]:
                        first_tr = i
                        last_tr = i
                        break
            else:
                first_tr = transcript[1].index(cds[1][0])
                last_tr = transcript[0].index(cds[0][-1])
            cds_start = sum([transcript[1][i] - transcript[0][i] for i in range(first_tr)]) + cds[0][0] - transcript[0][first_tr]
    #        cds_stop = mRNA_length - (transcript[1][last_tr] - cds[1][-1])
        else:
            if len(transcript[0]) == 1:
                first_tr = 0
                last_tr = 0
            elif len(cds[0]) == 1:
                for i in range(len(transcript[0])):
                    if transcript[0][i] <= cds[0][0]:
                        first_tr = i
                        last_tr = i
                        break
            else:
                first_tr = transcript[0].index(cds[0][0])
                last_tr = transcript[1].index(cds[1][-1])
            cds_start = sum([transcript[1][i] - transcript[0][i] for i in range(first_tr)]) + transcript[1][first_tr]- cds[1][0]
    #        cds_stop = mRNA_length - (cds[0][-1] - transcript[0][last_tr])
    except Exception, e:
        sys.stderr.write("Problem in '{}' @ {} --- {}\n".format(what, tr, str(e)))
        sys.stderr.write('exons: {}\n'.format(transcript))
        sys.stderr.write('CDSs: {}\n'.format(cds))
        sys.stderr.write('forward: {}\n'.format(forward))
        cds_start=0
        cds_length=0
        mRNA_length=0
        if what=='model':
            return [],[],[],[],[]
    cds_stop = cds_start + cds_length
    #assert (cds_stop - cds_start) % 3 == 0, "{}:\n{} > {}\n{}\n{},{},{}\n".format(tr_id, transcript, cds, first_tr, mRNA_length, cds_start, cds_stop)
    if what == 'start_stop':    
        return mRNA_length, cds_start, cds_stop
    elif what == 'model':
        try:
            if forward:
                utr5 = zip(*[(transcript[0][i], transcript[1][i]) for i in range(first_tr)] + [(transcript[0][first_tr], cds[0][0])])
                utr3 = zip(*[(cds[1][-1], transcript[1][last_tr])] + [(transcript[0][i], transcript[1][i]) for i in range(last_tr+1, len(transcript[0]))])
            else:
                utr5 = zip(*[(transcript[0][i], transcript[1][i]) for i in range(first_tr)] + [(cds[1][0], transcript[1][first_tr])])
                utr3 = zip(*[(transcript[0][last_tr], cds[0][-1])] + [(transcript[0][i], transcript[1][i]) for i in range(last_tr+1, len(transcript[0]))])
            transcript.extend(zip(*[(int(forward)*2 - 1, 'exon')]*len(transcript[0])))
            cds.extend(zip(*[(int(forward)*2 - 1, 'cds')]*len(cds[0]))) 
            utr5.extend(zip(*[(int(forward)*2 - 1, '5utr')]*len(utr5[0])))
            utr3.extend(zip(*[(int(forward)*2 - 1, '3utr')]*len(utr3[0])))
        except Exception, e:
            sys.stderr.write("Problem in '{}' @ {} --- {}\n".format(what, tr, str(e)))
            sys.stderr.write('exons: {}\n'.format(transcript))
            sys.stderr.write('CDSs: {}\n'.format(cds))
            sys.stderr.write('forward: {}\n'.format(forward))
            utr5 = []
            utr3 = []
        return transcript, cds, utr5, utr3
def make_composite(tr_list):
    if not tr_list:
        return 0,0,0
    features = []
    feature_sizes = {'cds':0, '3utr':0, '5utr':0}
    starts, stops, strands, names = range(4)
    all_breaks = []
    for tr in tr_list:
        exon, cds, utr5, utr3 = make_transcript(tr, what='model')
        all_breaks.extend(cds[starts]+utr5[starts]+utr3[starts]+cds[stops]+utr5[stops]+utr3[stops])
        #print cds, utr5, utr3
        #for (start, stop, strand, name) in zip(zip(*cds) + zip(*utr5) + zip(*utr3)):
            #interval = interval.insert(quicksect.Feature(start, stop, strand, name))
        try:
            features += [quicksect.Feature(start, stop-1, strand, name) for (start, stop, strand, name) in zip(*cds) + zip(*utr5) + zip(*utr3)]
        except Exception, e:
            sys.stderr.write('Problem during packing @ {}: {}\n'.format(tr, str(e)))
            sys.stderr.write('exons: {}\n'.format(exon))
            sys.stderr.write('CDSs: {}\n'.format(cds))
            sys.stderr.write('5UTRs: {}\n'.format(utr5))
            sys.stderr.write('3UTRs: {}\n'.format(utr3))
    if not features:
        sys.stderr.write('Problem @ {}\n'.format(tr_list))
        sys.exit(1)
    all_breaks = sorted(list(set(all_breaks)))
    all_intervals = [(all_breaks[i], all_breaks[i+1]-1) for i in range(len(all_breaks)-1)]
    interval_nodes = quicksect.IntervalNode(features[0])
    for feature in features[1:]:
        interval_nodes = interval_nodes.insert(feature)
    for interval in all_intervals:
        overlaps = interval_nodes.intersect(*interval)
#        print "{} : {}".format(interval, overlaps)
        if overlaps:
            label = sorted(set([feature.chr for feature in overlaps]))[-1]
            feature_sizes[label] += interval[1] - interval[0] + 1
            #print "{}:{}:{}".format(interval, overlaps, label)
    return sum(feature_sizes.values()), feature_sizes['5utr'], feature_sizes['5utr']+feature_sizes['cds']
gene_not_expressed = 0
transcript_not_expressed = 0
potential_transcripts = 0
multiprot = 0
superset = 0

with utils.measureTime('GTF has been filtered'):
    for gene_id, locus in loci.items():
        if not gene_id in expressed:
            gene_not_expressed += 1
            continue
        if args.verbose:
            sys.stderr.write('Processing: {}\n'.format(gene_id))
        potential_transcripts += len(locus['tr'])
        cur_trs = {}
        for tr_id in locus['tr']:
            if not tr_id in expressed[gene_id]:
                transcript_not_expressed += 1
                continue
            cur_trs[tr_id] = expressed[gene_id][tr_id]
        gene_type = 'NOT_PROCESSED'
        if len(cur_trs) > 1:
            cur_prots = {'KNOWN': set(), 'PUTATIVE': set(), 'NOVEL': set()}
            for tr_id, tr in cur_trs.items():
                cur_prots[tr[2]].add((tr[0],tuple(locus['tr'][tr_id]['CDS'])))
            if len(cur_prots['KNOWN']) > 1:
                multiprot += 1
                nonoverlapping_set = []
                while cur_prots['KNOWN']:
                    cur_prot = cur_prots['KNOWN'].pop()
                    cur_cds = tuple(sorted([(iv.start, iv.end) for iv in cur_prot[1]]))
                    if not nonoverlapping_set:
                        nonoverlapping_set.append((cur_cds,[cur_prot[0]]))
                        start = (cur_cds,[cur_prot[0]])
                    else:
                        added = False
                        for index, (cds, prot_ids) in enumerate(nonoverlapping_set):
                            a = b = None
                            is_a_set = False
                            if cds[0][0] == cur_cds[0][0] and cds[-1][1] == cur_cds[-1][1]:
                                if len(cds) > len(cur_cds):
                                    a = cds
                                    b = cur_cds
                                else:
                                    a = cur_cds
                                    b = cds
                            elif cds[0][0] <= cur_cds[0][0] and cds[-1][1] >= cur_cds[-1][1]:
                                # check if cur is a subset, if so add prot_id to list
                                a = cds
                                b = cur_cds
                            elif cds[0][0] >= cur_cds[0][0] and cds[-1][1] <= cur_cds[-1][1]:
                                # check if cds is a subset, if so replace tuple at index
                                a = cur_cds
                                b = cds
                            if a:
                                is_a_set = True
                                for exon in b[1:-1]:
                                    if exon in a:
                                        continue
                                    is_a_set = False
                            if is_a_set:
                                prot_ids.append(cur_prot[0])
                                nonoverlapping_set[index] = (a, prot_ids)
                                added = True
                                break
                        if not added:
                            nonoverlapping_set.append((cur_cds,[cur_prot[0]]))                            
                if len(nonoverlapping_set) == 1:
                    gene_type = 'MULTIPLE_INCLUSIVE'
                    superset += 1
                else:
                    gene_type = 'MULTIPLE_NONINCLUSIVE'
                    #print gene_id#, start, nonoverlapping_set
            else:
                gene_type = 'SINGLE_KNOWN'
        else:
            gene_type = 'SINGLE_PROTEIN'
        tr_list = []
        for tr_id, tr in cur_trs.items():
            if tr[2] == 'KNOWN':
                tr_list.append(tr_id)
        com_tr_len, com_cds_start, com_cds_stop = make_composite([locus['tr'][tr_id] for tr_id in tr_list])
        if (com_tr_len + com_cds_start + com_cds_stop) > 0:
            sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_id, gene_type, gene_id, 'composite', 'composite', com_tr_len, com_cds_start, com_cds_stop))
        for tr_id in tr_list:
            tr_len, cds_start, cds_stop = make_transcript(locus['tr'][tr_id]) #, tr_id)
            if (tr_len == 0 and cds_start == 0 and cds_stop == 0):
                sys.stderr.write('PROBLEM: {}:{}\n'.format(gene_id, tr_id))
                continue
            sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_id, gene_type, tr_id, cur_trs[tr_id][2], cur_trs[tr_id][0], tr_len, cds_start, cds_stop))
sys.stderr.write("{} of {} genes are not expressed.\n".format(gene_not_expressed, len(loci)))
sys.stderr.write("{} of {} potential transcripts are not expressed.\n".format(transcript_not_expressed, potential_transcripts))
sys.stderr.write("{} of {} multi-protein coding trascripts has a single superset.\n".format(superset, multiprot))
