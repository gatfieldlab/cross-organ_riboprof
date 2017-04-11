#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 17:23:09 2013

@author: barpat
"""

import sys, re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import get_density2_py3 as gd
from brewer2mpl import qualitative
from collections import defaultdict
from accessories import align_utils, utils, color_utils


#blue_gradient = color_utils.color_gradient(['navy','skyblue']).gradient
#red_gradient = color_utils.color_gradient(['red','pink']).gradient
blue_gradient = color_utils.color_gradient(['#0066cc','#5aadff']).gradient
red_gradient = color_utils.color_gradient(['#cc6600','#ffad5a']).gradient

REGIONS = ('5utr','cds','3utr')
MIN_WIDTH = 20
MAX_ROW = 6
SUMMARY_OPERATIONS = ('mean', 'sum', 'median')
DB_FILE, DB_REP, DB_FLAG = (4, 2, 6)
MAX_READ_LEN = 100
ASITE_OFF = {i: 15+int(i>30) for i in range(1,MAX_READ_LEN+1)}
RP_TR = "(RP|TR)"
ZT = "ZT(\d+)"
#COLORS = {'primary':{'fg': qualitative.Dark2, 'bg': qualitative.Set2}, 'secondary':{'fg': qualitative.Dark2, 'bg': qualitative.Set2}}
COLORS = {'primary':{'fg': blue_gradient, 'bg': blue_gradient}, 'secondary':{'fg': red_gradient, 'bg': red_gradient}}
LAYERS = {0: 'primary', 1: 'secondary', 'primary': 0, 'secondary': 1}
MIN_PALETTE = 2

def get_densitybin(tr_id, density_file, tr_len, fg_bg, asite_off=ASITE_OFF):
    try:
        extract = gd.get_index_tr(density_file, tr_id)
    except gd.extractionError as e:
        sys.stderr.write(e.message)
        sys.exit(1)
    f_density = []
    b_density = []
    foreground, background = fg_bg.split(',')
    for line in extract.strip().split('\n'):
        parsed = line.split('\t')
        try:
            read_pos = int(parsed[1])
            read_len = int(parsed[2])
        except IndexError:
            sys.stderr.write("Can't parse the density data returned by get_density from '{}'!\n".format(density_file))
            sys.exit(1)
        except ValueError:
            sys.stderr.write("Can't interpret density data returned by get_density!\n")
            sys.exit(1)
        if not foreground == 'none':
            offset = 0
            if foreground == 'Asite':
                try:
                    offset = asite_off[read_len]
                except KeyError:
                    sys.stderr.write("Found a read with a size outside of specified range of 1 to {}\nWill continue setting size to {}!\n".format(MAX_READ_LEN, MAX_READ_LEN))
                    offset = asite_off[MAX_READ_LEN]
            f_density.append(read_pos + offset)
        if not background == 'none': # bg can be only 'none' or otherwise pileup!!
            b_density.extend([i for i in range(read_pos, read_pos + read_len) if i < tr_len])
    return (np.bincount(f_density, minlength=tr_len), np.bincount(b_density, minlength=tr_len))

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

def int_factory():
    return defaultdict(int)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--tr_id', type=str)


parser.add_argument('-f', '--files', nargs='+', type=str)
parser.add_argument('-t', '--timepoints', nargs='+', type=int, default=range(0,24,2))
parser.add_argument('-d', '--database', type=str)
parser.add_argument('-z', '--zt', type=str, default=ZT)
parser.add_argument('-o', '--order', type=str, default=RP_TR)
parser.add_argument('-r', '--region', type=str)
parser.add_argument('-p', '--boxpos', type=str)
parser.add_argument('-e', '--extension', type=int, default=20)

group1 = parser.add_mutually_exclusive_group()
group1.add_argument('-c', '--cds', type=str, default='prepared_cds.txt')
group1.add_argument('-m', '--tr-model', type=str)

parser.add_argument('-n', '--normalization', type=str) # also enable to supply factors on the fly
parser.add_argument('-g', '--gene_name',type=str, default='')
parser.add_argument('-a', '--max_rows', type=int, default=MAX_ROW, help='specify max number of rows in the plot')
parser.add_argument('-l', '--max_read_len', type=int, default=MAX_READ_LEN, help='specify max length of reads, default={}'.format(MAX_READ_LEN))
parser.add_argument('-1', '--primary', type=str, default='Asite,none')
parser.add_argument('-2', '--secondary', type=str, default='none,pileup')
parser.add_argument('-s', '--scale-log', type=str, default='F,F')
parser.add_argument('-b', '--secondary_bg', action='store_true')
parser.add_argument('-C', '--color', type=str)
group2 = parser.add_mutually_exclusive_group()
group2.add_argument('-y', '--local_scale_y', default='', type=str, help='comma delimited string list of primary, secondary, tertiary to set the respective y-axis limits to a local (time point) specific value rather than global')
group2.add_argument('-x', '--fixed_scale_y', default='', type=str, help='example: "primary:25000" or "primary:300,secondary:500"')

#parser.add_argument('-s', '--sum', type=)

args = parser.parse_args()

cds_start = cds_end = tr_len = None
local_scale_layers = args.local_scale_y.split(',')
fixed_scale_layers = {}
for layval in args.fixed_scale_y.split(','):
    if layval == '':
        continue
    try:
        layer, value = layval.split(':')
        value = int(value)
    except ValueError:
        sys.stderr.write("Could not parse the --fixed_scale_y option. Please refer to help\n")
        sys.exit(1)
    fixed_scale_layers[layer] = value

if args.tr_model:
    tr_sizes = re.split(',|;|-|:|/', args.tr_model)
    if len(tr_sizes) < 2:
        sys.stderr.write("Need at least two sizes for 'tr_model'")
        sys.exit(1)
    try:
        cds_start = int(tr_sizes[0])
        cds_end = int(tr_sizes[1])
        tr_len = int(tr_sizes[2])
    except ValueError:
        sys.stderr.write("'tr_model' should contain only delimited numbers")
        sys.exit(1)
    except KeyError:
        sys.stderr.write("'tr_model' does not explicitly define transcript size. Will assume transcript end is same as CDS end!")
        tr_len = cds_end
else:
    try:
        with open(args.cds) as cds_file:
            for line in cds_file:
                parsed = line.strip().split('\t')
                if parsed[2] == args.tr_id:
                    (tr_len, cds_start, cds_end) = (int(parsed[i]) for i in (5,6,7))
                    break
    except IOError:
        sys.stderr.write('Could not read the CDS models from {}\n'.format(args.cds))
        sys.exit(1)

if not tr_len:
    sys.stderr.write('Could not set transcript size and CDS coordinates for {}\n'.format(args.tr_id))
    sys.exit(1)

plot_files = defaultdict(lambda: {'primary':[], 'secondary':[]})
zt_regex = re.compile(args.zt)
prm_scn_regex = re.compile(args.order)
prmscn = args.order[1:-1].split('|')
prm_scn_order = {prmscn[0]: 'primary', prmscn[1]: 'secondary'}

if not args.files:
    if args.timepoints and args.database:
        try:
            with open(args.database) as db_file:
                for line in db_file:
                    parsed = line.strip().split('\t')
                    if parsed[DB_FLAG] == 'Y':
                        zt = zt_regex.search(parsed[1])
                        prm_scn = prm_scn_regex.search(parsed[1])
                        if zt and prm_scn and int(zt.group(1)) in args.timepoints:
                            plot_files[int(zt.group(1))][prm_scn_order[prm_scn.group(1)]].append((int(parsed[DB_REP]), "{}_5prime_sorted.txt".format(parsed[DB_FILE])))
        except IOError:
            sys.stderr.write("Could not read the database file '{}'".format(args.database))
            sys.exit(1)
else: # Needs to be fixed !!!!
    reps_by_zt = defaultdict(int)
    zt_regex = re.compile(args.zt)
    for file in args.files:
        print("Searching: {}".format(file))
        zt = zt_regex.search(file)
        rp_tr = RP_TR_REGEX.search(file)
        if zt and rp_tr:
            reps_by_zt[zt.group(1)] += 1
            plot_files[int(zt.group(1))].append((reps_by_zt[zt.group(1)],file))

norm_factors = defaultdict(lambda: {'primary':{}, 'secondary':{}})
if args.normalization:
    with open(args.normalization) as normfile:
        for line in normfile:
            parsed = line.strip().split('\t')
            norm_factors[int(parsed[1])][prm_scn_order[parsed[0]]][int(parsed[2])] = float(parsed[3])
else:
    for zt in plot_files:
        #norm_factors[zt] = []
        for rep, dump in plot_files[zt]:
            norm_factors[zt][rep] = 1

# parse log-scale argument
scales_log = {'primary': False, 'secondary': False}
if args.scale_log:
    parsed = args.scale_log.split(',',1)
    for layer in scales_log:
        try:
            scales_log[layer] = utils.check_true(parsed[LAYERS[layer]])
        except KeyError:
            pass
print(scales_log)

# Calculating the plot region
region_start = 0
region_end = np.inf
if args.region:
    region_match = re.match('^(\d+)[-:_ ](\d+)$', args.region)
    if region_match:
        region_start = int(region_match.group(1)) - 1
        region_end = int(region_match.group(2)) - 1
    elif args.region in REGIONS:
        if args.region == '5utr':
            region_end = cds_start + args.extension
        elif args.region == 'cds':
            region_start = cds_start - args.extension
            region_end = cds_end + args.extension
        elif args.region == '3utr':
            region_start = cds_end - args.extension
    else:
        sys.stderr.write('Supplied region <{}> could not be interpreted... Whole transcript will be printed.\n')

# canvas calculations
numcols = int((len(plot_files)-1) / MAX_ROW) + 1
numrows = max((len(plot_files) >= MAX_ROW) * MAX_ROW, len(plot_files) % MAX_ROW)
gs = gridspec.GridSpec(numrows, numcols)
print("Canvas created with {} rows and {} columns.".format(numrows, numcols))

# tr and cds calculations
start = max(0, region_start)
end = min(tr_len, region_end)
box_start = max(start, cds_start)
box_end = min(end, cds_end)
box_length = max(0, box_end-box_start)

# parsing and setting manual boxes
manual_boxes = []
if args.boxpos:
    for box_coords in args.boxpos.split(','):
        try:
            mbox_start, mbox_end = (int(i) for i in box_coords.split(':',1))
        except:
            print("Can't parse part of 'boxpos' argument: '{}'".format(box_coords))
            sys.exit(1)
        else:
            manual_boxes.append((mbox_start, mbox_end))

fig = plt.figure(1)
numplot = 0
plot_metadata = {}
maxx_y = { 'primary': 0 , 'secondary': 0 }
for timepoint, density_meta in sorted(plot_files.items(), key=lambda k_v: k_v[0]):
    numplot += 1
    plot_metadata[numplot] = {'primary':{}, 'secondary':{}}
    for layer in ['primary', 'secondary']:
        fg_bg = getattr(args, layer)
        try:
            reps, density_files = zip(*density_meta[layer])
        except TypeError:
            print(layer, density_meta)
            sys.exit(1)
        rep_set = sorted(set(reps))
        plot_metadata[numplot][layer] = {'tpoint': timepoint, 'max_y': 0, 'num_plots': len(rep_set),
                                         'densitybins': {'fg':[], 'bg':[]}, 'colors': {'fg':[], 'bg':[]}}
        cn = max(MIN_PALETTE, min(len(rep_set), 11))
        colorscheme = {ground: COLORS[layer][ground](cn) for ground in ['fg', 'bg']}
        max_y = 0
        for ind, rep in enumerate(rep_set):
            norm_factor = norm_factors[timepoint][layer][rep]
            if norm_factor == 0:
                norm_factor = 1
            indices = [i for i, x in enumerate(reps) if x == rep]
            rep_density_files = [density_files[i] for i in indices]
            densities = list(zip(*[get_densitybin(args.tr_id, density_file, tr_len, fg_bg) for density_file in rep_density_files]))
            densitybin = {}
            curcolor = {}
            densitybin['fg'] = summarize_densitybins(densities[0], 'sum') / norm_factor
            densitybin['bg'] = summarize_densitybins(densities[1], 'sum') / norm_factor
            curcolor['fg'] = colorscheme['fg'][ind]
            curcolor['bg'] = colorscheme['bg'][ind]
            print(layer, timepoint, rep, args.tr_id, norm_factor, rep_density_files, fg_bg, curcolor)
            ground_flags = fg_bg.split(',')
            ground_mode = {'fg': ground_flags[0], 'bg': ground_flags[1]}
            for ground in ['fg', 'bg']:
                plot_metadata[numplot][layer]['densitybins'][ground].append((densitybin[ground], ground_mode[ground]))
                plot_metadata[numplot][layer]['colors'][ground].append(curcolor[ground])
                max_y = max(max_y, max(densitybin[ground][start:end]))
                maxx_y[layer] = max(maxx_y[layer], max_y)
        plot_metadata[numplot][layer]['max_y'] = max_y

tps = sorted(plot_files.keys())
tpdiff = np.array(tps[1:]) - np.array(tps[:-1])
print_title = len(tpdiff) and not np.all(tpdiff == tpdiff[0])
for numplot, metadata in sorted(plot_metadata.items(), key=lambda k_v: k_v[0]):
    curcol = int((numplot-1) / MAX_ROW)
    currow = (numplot - 1) % MAX_ROW
    inner_grid = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[currow, curcol], wspace=0.0, hspace=0.0)
    ax = {}
    cur_y_lims = {}
    for layer in ['primary', 'secondary']:
        ax[layer] = plt.Subplot(fig, inner_grid[LAYERS[layer]])
        ax[layer].tick_params(axis='y', which='major', labelsize=10)
        ax[layer].axvspan(box_start, box_end+1,facecolor='gray', alpha=0.2, linewidth=0)
        for subplot in range(metadata[layer]['num_plots']):
            for ground in ['bg', 'fg']:
                densitybin = metadata[layer]['densitybins'][ground][subplot]
                curcolor = metadata[layer]['colors'][ground][subplot]
                if not densitybin[1] == 'none':
                    ax[layer].bar(range(start+1,end+1), densitybin[0][start:end], width=1.0, # edgecolor=curcolor,
                                  align='center', color=curcolor, linewidth=0, alpha=0.8, log=scales_log[layer])
        ax[layer].set_xlim(start+1, end+1)
        if not scales_log[layer]:
            tight_axis = 'both'
            min_y = 0
        else:
            tight_axis = 'x'
            min_y = 0.5
        ax[layer].locator_params(nbins=3, axis=tight_axis, tight=True)
        if layer in local_scale_layers:
            cur_max_y = metadata[layer]['max_y']
        elif layer in fixed_scale_layers:
            cur_max_y = fixed_scale_layers[layer]
        else:
            cur_max_y = maxx_y[layer]
        cur_y_lims[layer] = cur_max_y
        if layer == 'primary':
            ax[layer].set_ylim(min_y, cur_max_y)
        else:
            ax[layer].set_ylim(cur_max_y, min_y)
    ax['primary'].axes.get_xaxis().set_tick_params(labelbottom=False)
    if numplot % numrows != 0:
        ax['secondary'].axes.get_xaxis().set_tick_params(labelbottom=False)
    else:
    # Manual boxes are plotted via boxpos, e.g uORFs etc
        if manual_boxes:
            box_pos = cur_y_lims['secondary'] * 1.05
            box_size = cur_y_lims['secondary'] * 0.05
            for mbox_start, mbox_end in manual_boxes:
                #manual_box = patches.Rectangle((mbox_start, maxx_y['secondary']+20), mbox_end - mbox_start, 5)
                manual_box = patches.Rectangle((mbox_start, box_pos), mbox_end - mbox_start, box_size)
                manual_box.set_clip_on(False)
                ax['secondary'].add_patch(manual_box)
    cur_tp = "ZT{}".format(metadata['primary']['tpoint'])
    if print_title or numplot == 1 or numplot % numrows == 1:
        ax['primary'].set_title(cur_tp, fontsize=12, position=(0,1), horizontalalignment='left')
    #ax.text(start,int(maxx_y)+3,"ZT{}".format(metadata['tpoint']))
    fig.add_subplot(ax['primary'])
    fig.add_subplot(ax['secondary'])
    print(numplot, cur_tp, curcol, currow, metadata['primary']['tpoint'], cur_y_lims)
plt.suptitle("{} {}".format(args.gene_name,args.tr_id), x=.5, y=.98, horizontalalignment='center', fontsize=14)
#plt.tight_layout()
print('Finished rendering the plot, now saving as pdf ...')
plt.savefig('{}_{}_{}_density_plot.pdf'.format(args.gene_name, args.tr_id, (args.region or 'Whole')), format='PDF')
