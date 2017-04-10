#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#v2 includes orf_composite
import argparse, pysam, sys
import calc_codons as cc
import numpy as np
import bisect
import multiprocessing
# from scipy import stats
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna

FA_DIR = '/home/violeta/projects/database/fasta_sequences/'
FA_FILE = 'Mmusculus.GRCm38.75.cdna.ensembl.fa'
TR_ASITE_OFF = {i: i // 2 for i in range(1,cc.MAX_READ_LEN+1)}
TR_5P_OFF = {i: 0 for i in range(1,cc.MAX_READ_LEN+1)}
START_CODONS = ('ATG')
STOP_CODONS = ('TAG', 'TAA', 'TGA')
DB_FILE, DB_REP, DB_FLAG, DB_SAMPLE = (4, 2, 6, 1)


def get_orfs(sequence, min_pep_len=10, frames=range(3)):
    orfs = []
    main_orf_frame = (len(sequence) + 1) % 3
    for frame in frames:
        starts = []
        stops  = []
        pos = frame
        seq = sequence[frame:]
        while (len(seq) > 2):
            codon = seq[:3]
            seq = seq[3:]
            if codon in START_CODONS:
                starts.append(pos)
            elif codon in STOP_CODONS:
                stops.append(pos)
            pos += 3
        if starts and stops:
            for start_pos in starts:
                can_skip_size_control = False
                smallest_stop = bisect.bisect_right(stops, start_pos)
                if smallest_stop == len(stops):
                    if frame == main_orf_frame:
                        break
                    else:
                        stops = [len(sequence)+frame-2]
                        can_skip_size_control = True
                else:
                    stops = stops[smallest_stop:]
                if can_skip_size_control or (stops[0]-start_pos) >= min_pep_len*3:
                    orfs.append((start_pos, stops[0]))
                else:
                    stops.pop(0)
    return(sorted(orfs))

def make_orf_composite(orfs):
    orf_composite = [(orfs[0])]
    for pos in range(len(orfs)-1):
        if orfs[pos+1][0] in range(orf_composite[-1][0],orf_composite[-1][1]) and orfs[pos+1][1] in range(orf_composite[-1][0],orf_composite[-1][1]): #cur orf contained in previous
            orf_composite = orf_composite
        elif orfs[pos+1][0] in range(orf_composite[-1][0],orf_composite[-1][1]): #partially overlapping (not fully contained)
            orf_composite[-1] = (orf_composite[-1][0], orfs[pos+1][1])
        elif orfs[pos+1][0] > orfs[pos][1]:  #non-overlapping
            orf_composite.append((orfs[pos+1][0], orfs[pos+1][1]))
    return orf_composite


def workerstar(args):
    return worker(*args)

def worker(gid, tr_id, slice_limits, extracts, errors):
    result = ""
    if len(errors) == 1 and isinstance(errors[0], cc.densityError):
        sys.stderr.write(errors[0].message)
    else:
        # apply normalization here
        extracts = extracts / np.array(norm_factors)[:,None]
        summarized = cc.summarize_densitybins(extracts)
        seq_obj = pysam.faidx(FA_DIR+FA_FILE, "|".join([gid, tr_id]))
        seq = "".join([l.strip() for l in seq_obj[1:]])[slice_limits[0]:slice_limits[1]]
        orfs = make_orf_composite(get_orfs(seq, min_pep_len))
        for orf_start, orf_stop in orfs:
            frame_count = [0,0,0]
            coverage = 0
            for pos in range(orf_start, orf_stop):
                pos_count = summarized[pos]
                frame_count[pos % 3] += pos_count
                if pos_count:
                    coverage += 1
            result += "\t".join([gid, tr_id] + [str(i) for i in [orf_start, orf_stop, coverage] +frame_count]) + "\n"
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='*', metavar='density_file', type=str, help='density files to process. if present overwrites --database and --normalization arguments.')
    parser.add_argument('-t', '--tr_ids', type=str, default='robust_single_trs_kidney.txt')
    parser.add_argument('-r', '--region', type=str, default='5utr')
    parser.add_argument('-d', '--database', type=str)
    parser.add_argument('-w', '--what', default='FP', choices=['RP', 'FP', 'TR'])
    parser.add_argument('-c', '--cds', type=str, default='prepared_cds_v38.75_filtered_kidney_vFinal.txt')
    parser.add_argument('-l', '--max_read_len', type=int, default=cc.MAX_READ_LEN, help='specify max length of reads, default={}'.format(cc.MAX_READ_LEN))
    parser.add_argument('-n', '--normalization', type=str) # also enable to supply factors on the fly

    args = parser.parse_args()
    tr_ids = []

    global norm_factors
    global min_pep_len

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

    if args.normalization:
        norm_factor_db = {}
        with open(args.normalization) as normfile:
            for line in normfile:
                parsed = line.strip().split('\t')
                norm_factor_db["{0}ZT{1:0>2}_{2}".format(*parsed)] = float(parsed[3])
    if args.files:
        density_files = args.files
        norm_factors = [1 for i in range(len(density_files))]
    else:
        if not args.database:
            sys.stderr.write('When no density files are specified in command line, a database need to be set with -d/--database option.\n')
            sys.exit(1)
        try:
            density_files = []
            norm_factors = []
            with open(args.database) as db_file:
                for line in db_file:
                    parsed = line.strip().split('\t')
                    if parsed[DB_FLAG] == 'Y' and parsed[DB_FILE][:2] == args.what:
#                        density_files.append("{}_vBAM_vFinal2_5prime_count_sorted.txt".format(parsed[DB_FILE])) #for liver
                        density_files.append("{}_5prime_sorted.txt".format(parsed[DB_FILE])) #for kidney
                        norm_factors.append(norm_factor_db[parsed[DB_SAMPLE] + '_' + parsed[DB_REP]])
        except IOError:
            sys.stderr.write("Could not read the database file '{}'".format(args.database))
            sys.exit(1)
    density_extractor = cc.DensityParser(tr_db, tr_ids, density_files, args.region)
    output_handle = sys.stdout
    min_pep_len = 6
    PROCESSES = 4
    with multiprocessing.Pool(PROCESSES) as pool:
        TASKS = density_extractor.get_densitybin()
        results = pool.imap_unordered(workerstar, TASKS)
        for line in results:
            output_handle.write(line)

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
