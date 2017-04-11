#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 2016

@author: violeta castelo-szekely
"""

import argparse, pysam, sys, subprocess, re
import numpy as np

transcripts= {}
tr_list = set()
sequences= {}

def get_length(sequence):
    return len(sequence)

def get_gc_content(sequence):
    if len(sequence) < 1:
        return ('NA')
    num_g = sequence.count('G')
    num_c = sequence.count('C')
    num_a = sequence.count('A')
    num_t = sequence.count('T')
    return (num_g + num_c) / float((num_g + num_c + num_a + num_t))

def get_kozak_score(sequence):
    # for incomplete kozak sequences (from transcripts with zero-length or <6 5'UTR), setting the score to -1 to indicate oddball status
    if len(sequence) < 10:
        kozak_score = 'NA'
# kozak consensous is: GccA/GccATGG, where A is the +1 position. Score:
#'G' at -6 (+3 points); 'A'/'G' at -3 (3 points); G at +4 (3 points.)
#'C' at positions -1, -2, -4, -5 (1 point each)
    else:
        kozak_score = 0
        kozak_score += 3 * sum((sequence[0] == "G", sequence[3] == "A" or sequence[3] == "G", sequence[-1] == "G"))
        kozak_score += sum ((sequence[1] == "C", sequence[2] == "C", sequence[4] == "C", sequence[5] == "C"))
    return kozak_score

# MFE (minimum folding energy. Using RNAfold from viennaRNA.)
def execute_external(program):
    process = subprocess.Popen(program,
                         stdout = subprocess.PIPE,
                         shell = True)
    return iter(process.stdout.readline, b'')
def get_energy(sequence):
    if len(sequence) < 1:
        return ('NA')
    lines = execute_external("echo {} | RNAfold --noPS -p".format(sequence))
    lines.next() #skip sequence line
    MFE = float(re.sub('[()]','', lines.next().split()[-1]))
    return (MFE)



parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gtf', help="specify the GTF file")
parser.add_argument('-r', '--region', help="specify the region. Choices: 'utr5', 'cds', 'utr3', 'transcript', 'kozak'.", type=str)
parser.add_argument('-f', '--fasta', help="specify the cDNA fasta.", type=str)
parser.add_argument('-w', '--what', help="feature of interest. Implemented options: 'length', 'gc' (GC content), 'kozak' (Kozak context), 'mfe' (minimum folding energy).", type=str)
parser.add_argument('-l', '--list', help="file containing the list of GeneID and TranscriptIDs of interest", type=str)
args = parser.parse_args()

try:
    with open(args.gtf) as gtf_file:
        for line in gtf_file:
            parsed = line.strip().split('\t')
            if parsed[4] == 'composite' :
                continue
            gid = parsed[0]
            trid = parsed[2]
            if not "|".join([gid,trid]) in transcripts:
                transcripts["|".join([gid,trid])] = {'transcript': (0, int(parsed[5])), 'utr5': (0, int(parsed[6])-1),'cds': (int(parsed[6]), int(parsed[7])), 'utr3': (1+int(parsed[7]), int(parsed[5])), 'kozak': (int(parsed[6])-6, int(parsed[6])+4)}
                #this is to handle transcripts with no 5'utr (cds starting at 0):
                if transcripts["|".join([gid,trid])]['utr5'][1] < 0:
                    transcripts["|".join([gid,trid])]['utr5'] = (0,0)
                    transcripts["|".join([gid,trid])]['kozak'] = (0,4)
except IOError:
    sys.stderr.write("Could not read the GTF file '{}'".format(args.gtf))
    sys.exit(1)

if args.fasta == None:
    sys.stderr.write('Need a cDNA FASTA file, use -f/--fasta option.\n')
    sys.exit(1)
else:
    FA_FILE = args.fasta

if args.region == None:
    sys.stderr.write('Need to specify a region, use -r/--region option. For Kozak score, choose "-r kozak".\n')
    sys.exit(1)
else:
    region = args.region

if args.list == None:
    sys.stderr.write("Could not read the GeneID, TranscriptID file {}".format(args.list))
    sys.exit(1)
with open(args.list) as trs_file:
        for line in trs_file:
            parsed = line.strip().split('\t')
            gid = parsed[0]
            trid = parsed[1]
            seq_obj = pysam.faidx(FA_FILE, "|".join([gid, trid]))
            if not "|".join([gid,trid]) in sequences:
                sequences["|".join([gid,trid])] = "".join([l.strip() for l in seq_obj[1:]])[transcripts["|".join([gid,trid])][region][0]:transcripts["|".join([gid,trid])][region][1]]
for gid in sequences:
    if args.what == 'kozak':
        sys.stdout.write('{}\t{}\t{}\n'.format(gid.split('|')[0], gid.split('|')[1], get_kozak_score(sequences[gid])))
    if args.what == 'mfe':
        mfe = get_energy(sequences[gid])
        sys.stdout.write('{}\t{}\t{}\t{}\t{}\n'.format(gid.split('|')[0], gid.split('|')[1],mfe[0],mfe[1],mfe[2]))
    if args.what == 'length':
        sys.stdout.write('{}\t{}\t{}\n'.format(gid.split('|')[0], gid.split('|')[1], get_length(sequences[gid])))
    if args.what == 'gc':
        sys.stdout.write('{}\t{}\t{}\n'.format(gid.split('|')[0], gid.split('|')[1], get_gc_content(sequences[gid])))
    if args.what == None:
        sys.stderr.write('Need to specify a feature property to analyze, use -w/--what option.')
