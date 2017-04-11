# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 11:01:26 2013

@author: barpat and violeta
v3 looks at diversity (like in v1) or identity (transcripts are identical in specified region and different in at least one of the other regions).
"""

import HTSeq
import sys, argparse
import quicksect
from accessories import utils
from collections import OrderedDict

parser = argparse.ArgumentParser(description='List genes whose transcript diversity comes only from the specified feature (utr5, utr3 or cds)')
parser.add_argument('-v','--verbose', help="A little more verbosity helps humanity ;)", action='store_true')
parser.add_argument('-d', '--delimiter', help="delimiter used in output files. Default: '\\t'", type=str,
                    default='\t')
parser.add_argument('-l', '--list', help="file containing a list of GeneID/TranscriptIDs for filtering the GTF transcripts", type=str)
parser.add_argument("gtf", help="the GTF containing all models/exons/CDSs", type=str)
parser.add_argument('-r','--region', help='Region to look for diversity/identicalness: utr5, utr3 or cds', type=str)
parser.add_argument('-w','--what', help='What to look for: "diversity" (transcripts are different only in the specified region) or "identical" (transcripts are identical only in the specified region)', type=str)
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
        if not (feature.type in ('exon', 'CDS')  and feature.source in ('protein_coding', 'nonsense_mediated_decay', 'non_stop_decay')):
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


gene_not_expressed = 0
transcript_not_expressed = 0
potential_transcripts = 0
multiprot = 0
superset = 0
binned_trs = 0
binned_g = 0
binned_genes = []
utr5_diversity = []
utr3_diversity = []
cds_diversity = []
utr5_identical = []
utr3_identical = []
cds_identical = []
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
        tr_list = []
        for tr_id, tr in cur_trs.items():
            if tr[2] == 'KNOWN':
                tr_list.append(tr_id)
        cds_models = OrderedDict()
        utr3_models = OrderedDict()
        utr5_models = OrderedDict()
        for tr_id in tr_list:
            tr_exon, tr_cds, tr_utr5, tr_utr3 = make_transcript(locus['tr'][tr_id], what='model')
            if (tr_utr5[0],tr_utr5[1]) not in utr5_models:
                utr5_models[tr_utr5[0],tr_utr5[1]] = [tr_id]
            else:
                utr5_models[tr_utr5[0],tr_utr5[1]].append(tr_id)
            if (tr_cds[0],tr_cds[1]) not in cds_models:
                cds_models[tr_cds[0],tr_cds[1]] = [tr_id]
            else:
                cds_models[tr_cds[0],tr_cds[1]].append(tr_id)
            if (tr_utr3[0],tr_utr3[1]) not in utr3_models:
                utr3_models[tr_utr3[0],tr_utr3[1]] = [tr_id]
            else:
                utr3_models[tr_utr3[0],tr_utr3[1]].append(tr_id)
        if len(cds_models) == 1 and len(utr3_models) == 1 and len(utr5_models) != 1:
            utr5_diversity.append(gene_id)
        elif len(cds_models) == 1 and len(utr5_models) == 1 and len(utr3_models) != 1:
            utr3_diversity.append(gene_id)
        elif len(utr5_models) == 1 and len(utr3_models) == 1 and len(cds_models) != 1:
            cds_diversity.append(gene_id)
        if len(cds_models) == 1 and (len(utr3_models) !=1 or len(utr5_models) != 1):
            cds_identical.append(gene_id)
        if len(utr3_models) == 1 and (len(cds_models) !=1 or len(utr5_models) != 1):
            utr3_identical.append(gene_id)
        if len(utr5_models) == 1 and (len(utr3_models) != 1 or len(cds_models) != 1):
            utr5_identical.append(gene_id)

if args.region == 'utr5':
    if args.what == 'diversity':
        for gene in utr5_diversity:
            sys.stdout.write("{}\n".format(gene))
        sys.stderr.write("{} genes differ only in their 5'UTR.\n".format(len(utr5_diversity)))
    elif args.what == 'identical':
        for gene in utr5_identical:
            sys.stdout.write("{}\n".format(gene))
        sys.stderr.write("{} genes are identical in their 5'UTR and different in either or both CDS and 3'UTR.\n".format(len(utr5_identical)))
elif args.region == 'utr3':
    if args.what == 'diversity':
        for gene in utr3_diversity:
            sys.stdout.write("{}\n".format(gene))
        sys.stderr.write("{} genes differ only in their 3'UTR.\n".format(len(utr3_diversity)))
    elif args.what == 'identical':
            for gene in utr3_identical:
                sys.stdout.write("{}\n".format(gene))
            sys.stderr.write("{} genes are identical in their 3'UTR and different in either or both CDS and 5'UTR.\n".format(len(utr3_identical)))
elif args.region == 'cds':
    if args.what == 'diversity':
        for gene in cds_diversity:
            sys.stdout.write("{}\n".format(gene))
        sys.stderr.write("{} genes differ only in their CDS.\n".format(len(cds_diversity)))
    elif args.what == 'identical':
        for gene in cds_identical:
            sys.stdout.write("{}\n".format(gene))
        sys.stderr.write("{} genes genes are identical in their CDS and different in either or both 3'UTR and 5'UTR.\n".format(len(cds_identical)))
