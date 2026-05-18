#!/usr/bin/env python
# coding: utf-8

import re
import argparse
import os
import logging
import sys
import gzip
import pandas as pd
import statsmodels.discrete.discrete_model as sm
from statsmodels.api import add_constant
import numpy as np
import copy
import pybedtools as bt
import glob
from scipy.io import mmread

logging.basicConfig(level=logging.INFO, format = '%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('--peaks', nargs='+', help='Paths to peak calls for your cell types. File name is assumed to be {cell_type}.*')
parser.add_argument('--zhang-matrix', default='/lab/work/porchard/2022-muscle-sn/data/sciatac/cre-by-cell-type/matrix.tsv.gz', help='Path to Zhang matrix.tsv.gz (default: /lab/work/porchard/2022-muscle-sn/data/sciatac/cre-by-cell-type/matrix.tsv.gz)')
parser.add_argument('--zhang-cell-types', default='/lab/work/porchard/2022-muscle-sn/data/sciatac/cre-by-cell-type/celltypes.txt.gz', help='Path to Zhang celltypes.txt.gz (default: /lab/work/porchard/2022-muscle-sn/data/sciatac/cre-by-cell-type/celltypes.txt.gz)')
parser.add_argument('--zhang-peaks', default='/lab/work/porchard/2022-muscle-sn/data/sciatac/cre-by-cell-type/cCREs.bed.gz', help='Path to Zhang cCREs.bed.gz (default: /lab/work/porchard/2022-muscle-sn/data/sciatac/cre-by-cell-type/cCREs.bed.gz)')
parser.add_argument('--chrom-sizes', default='/lab/work/arushiv/muscle-sn/data/chrom_sizes/hg38.chrom_sizes', help='Path to chrom size file (default: /lab/work/arushiv/muscle-sn/data/chrom_sizes/hg38.chrom_sizes)')
parser.add_argument('--tss', default='/lab/work/porchard/reference/tss/hg38.gencode.tss.bed.gz', help='Path to TSS bed file (default: /lab/work/porchard/reference/tss/hg38.gencode.tss.bed.gz)')
parser.add_argument('--prefix', default='logregression-ren.', help='Prefix to use for output file(s) (default: "logregression-ren.")')
args = parser.parse_args()

# read in our peaks
#NARROWPEAK_FILES = glob.glob('/lab/work/arushiv/muscle-sn/analyses_hg38/files_for_peter/atac/peaks-by-cluster/atac-*-fusion/atac-*-fusion_peaks.narrowPeak')
NARROWPEAK_FILES = args.peaks
REN_MATRIX = args.zhang_matrix
REN_PEAKS = args.zhang_peaks
REN_CELL_TYPES = args.zhang_cell_types
TSS = args.tss
CHROM_SIZES = args.chrom_sizes
PREFIX = args.prefix


narrowpeaks = pd.concat([pd.read_csv(f, sep='\t', header=None, usecols=range(4), names=['chrom', 'start', 'end', 'name']).assign(cell_type=os.path.basename(f).split('.')[0]) for f in NARROWPEAK_FILES])
narrowpeaks['name'] = narrowpeaks.cell_type
narrowpeaks = narrowpeaks.drop(columns=['cell_type'])

# create master peaks and filter to TSS-distal
master_peaks = bt.BedTool().from_dataframe(narrowpeaks[['chrom', 'start', 'end']]).sort().merge().to_dataframe()
promoters = bt.BedTool(TSS).sort().slop(l=5000, r=0, s=True, g=CHROM_SIZES).sort().merge()
master_peaks = bt.BedTool().from_dataframe(master_peaks).intersect(promoters, v=True).to_dataframe()

# determine accessibility in each cell type
accessibility = bt.BedTool().from_dataframe(master_peaks).sort().intersect(bt.BedTool().from_dataframe(narrowpeaks[['chrom', 'start', 'end', 'name']]).sort(), wa=True, wb=True).to_dataframe()
accessibility['cell_type'] = accessibility.thickStart
accessibility['accessible'] = 1
accessibility = accessibility[['chrom', 'start', 'end', 'cell_type', 'accessible']].drop_duplicates().pivot(index=['chrom', 'start', 'end'], columns='cell_type', values='accessible').fillna(0).astype(int)
lhs = accessibility.reset_index()


ren_matrix = mmread(REN_MATRIX).tocsc()
ren_peaks = pd.read_csv(REN_PEAKS, sep='\t', header=None, names=['chrom', 'start', 'end'])
ren_cell_types = pd.read_csv(REN_CELL_TYPES, sep='\t', header=None, names=['ren_cell_type'])
REN_CELL_TYPE_INDEX_DICT = dict(zip(ren_cell_types.ren_cell_type, ren_cell_types.index))

# create RHS too
# bed tools, intersecting our peaks with theirs
for_rhs = pd.read_csv(REN_MATRIX, sep=' ', skiprows=3, header=None, names=['peak_idx', 'cell_type_idx'])
tmp = ren_peaks.rename_axis(index='peak_idx').reset_index()
tmp['peak_idx'] = tmp['peak_idx'] + 1
for_rhs = for_rhs.merge(tmp)
tmp = ren_cell_types.rename_axis(index='cell_type_idx').reset_index()
tmp['cell_type_idx'] = tmp['cell_type_idx'] + 1
for_rhs = for_rhs.merge(tmp)
for_rhs = for_rhs[['chrom', 'start', 'end', 'ren_cell_type']]

# now build rhs
lhs_bt = bt.BedTool().from_dataframe(lhs[['chrom', 'start', 'end']]).sort()
rhs_bt = bt.BedTool().from_dataframe(for_rhs).sort()
overlap_wide = lhs_bt.intersect(rhs_bt, loj=True).to_dataframe()
overlap_wide = overlap_wide[['chrom', 'start', 'end', 'thickStart']].drop_duplicates()
overlap_wide['val'] = 1
overlap_wide['thickStart'] = overlap_wide['thickStart'].map(lambda x: x if x != '.' else 'DROP')
overlap_wide = overlap_wide.pivot(index=['chrom', 'start', 'end'], columns='thickStart', values='val').fillna(0).astype(int)
rhs = overlap_wide.reset_index().drop(columns=['DROP'])

OUR_CELL_TYPES = [i for i in lhs.columns if i not in ['chrom', 'start', 'end']]
REN_CELL_TYPES = [i for i in rhs.columns if i not in ['chrom', 'start', 'end']]
logging.info('Originally had {} ren cell types'.format(len(ren_cell_types)))
logging.info('Now have {} ren cell types'.format(len(REN_CELL_TYPES)))

# use our cell type specific peaks only
lhs_cell_type_specific_only = lhs[lhs[OUR_CELL_TYPES].sum(axis=1)==1]
rhs_cell_type_specific_only = rhs[lhs[OUR_CELL_TYPES].sum(axis=1)==1]
assert(all(rhs_cell_type_specific_only.start == lhs_cell_type_specific_only.start))

# now model for each pair
p_threshold = 0.05 / len(REN_CELL_TYPES)

RESULTS = []
for OUR_CELL_TYPE in OUR_CELL_TYPES:
    for REN_CELL_TYPE in REN_CELL_TYPES:
        tmp = pd.DataFrame({'lhs': lhs_cell_type_specific_only[OUR_CELL_TYPE], 'rhs': rhs_cell_type_specific_only[REN_CELL_TYPE]})
        if len(tmp.groupby(['lhs', 'rhs']).size()) < 4:
            logging.warning(f'Skipping {OUR_CELL_TYPE} - {REN_CELL_TYPE}\n')
            continue
        mod = sm.Logit(lhs_cell_type_specific_only[OUR_CELL_TYPE], add_constant(rhs_cell_type_specific_only[[REN_CELL_TYPE]], prepend = False)).fit(disp=0)
        p = mod.pvalues[0]
        coef = mod.params[0]
        converged = mod.mle_retvals['converged']
        RESULTS.append([OUR_CELL_TYPE, REN_CELL_TYPE, p, coef, converged])
results_cell_type_specific_only = pd.DataFrame(RESULTS, columns=['our_cell_type', 'ren_cell_type', 'p', 'coef', 'converged'])
results_cell_type_specific_only['bonferroni_significant'] = (results_cell_type_specific_only.p<=p_threshold)
results_cell_type_specific_only.to_csv(f'{PREFIX}results-cell-type-specific-peaks.tsv', sep='\t', index=False)
