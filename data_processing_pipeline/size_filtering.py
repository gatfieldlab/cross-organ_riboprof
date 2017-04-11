# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 11:56:17 2014

@author: violeta
"""
from __future__ import division
import sys
import itertools

try:
    MIN_SEQ_LENGTH = int(sys.argv[3])
except:
    MIN_SEQ_LENGTH = 21 #for RP, 26
try:
    MAX_SEQ_LENGTH = int(sys.argv[4])
except:
    MAX_SEQ_LENGTH = 60 #for RP, 35
    
total_count = 0
filtered_count = 0

# to account for newline character
EFF_MIN_SEQ_LENGTH = MIN_SEQ_LENGTH + 1
EFF_MAX_SEQ_LENGTH = MAX_SEQ_LENGTH + 1

outfile = open(sys.argv[2], 'a')
with open(sys.argv[1], 'r') as trimmed_samples:
    for sample in trimmed_samples:
        sample = sample.strip()
        with open("../trimmed_data/{0}.fastq".format(sample)) as f:
            for identifier,seq,plus,quality in itertools.izip_longest(*[f]*4):
                total_count += 1
                if len(seq) < EFF_MIN_SEQ_LENGTH or len(seq) > EFF_MAX_SEQ_LENGTH:
                    filtered_count += 1
                else:
                    outfile.write(identifier+seq+plus+quality)
                    
outfile.close()
print 'Total number of sequences:', total_count
print 'Number of filtered out sequences:', filtered_count
print 'Number of kept sequences:', total_count - filtered_count
print 'Percentage of kept sequences:', ((total_count - filtered_count)/ total_count)*100
print 'Filter size - min:{}, max:{}'.format(MIN_SEQ_LENGTH, MAX_SEQ_LENGTH)         
        
