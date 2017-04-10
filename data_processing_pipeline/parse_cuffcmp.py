#!/usr/bin/env python

import sys, operator
op = operator
#inc = range(6,34)
inc = range (4,27)
FMI, FPKM, LOW95, HI95, COV = range(5)
cutoff = (2.0, 0.1, 0.05, None, None)
ops = (op.gt, op.gt, op.gt, None, None)

use_for_filter = (FMI, FPKM, LOW95)
count = 0
pass_count = 0
with open(sys.argv[1]) as f:
	for l in f:
                count += 1
		p = l.strip().split('\t')
		#p = l.strip().split()
		passed = []
		for i in inc:
			s = p[i].split('|')
			g = (s[0].split(':')[1], s[1])
			r = [float(v) for v in s[2:7]]
			passed.append(all([ops[feat](r[feat], cutoff[feat]) for feat in use_for_filter]))
		if sum(passed) > 3:
			pass_count += 1
			sys.stdout.write("{}\n".format("\t".join(g)))
sys.stderr.write("{} out of {} passed the filter\n".format(pass_count, count))
			

			
			
