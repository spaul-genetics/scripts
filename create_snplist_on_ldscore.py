# Given score.ld file (output from GCTA) it creates bins of variants.
# Binning is based on MAF and quantile of ldscore_SNP
# Author: Subrata Paul
import argparse
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument('--score', type=str, help = 'The score.ld file.')
parser.add_argument('--chrom', type = int, default = 21, help = 'Enter the chromosome number.')
parser.add_argument('--mafs', nargs='+', type=float, default=[0.001, 0.01, 0.05], help = 'Enter the breaks for MAF bins')
parser.add_argument('--lds', type = int, default = 4, help = 'Enter number of LD bins within each MAF bin')
parser.add_argument('--out', type = str, default = '', help='Enter the output file prefix')
args = parser.parse_args()


def create_snplist(dat, CHR, mafs, lds):
    mafs.insert(0,0)
    mafs.insert(len(mafs),0.5)
    dat['MAF'] = (dat.freq >= 0.5) * (1- dat.freq) + (dat.freq < 0.5) * dat.freq
    for maf_bin in range(1, len(mafs)):
        dat_subset = dat[(dat.MAF> mafs[maf_bin-1]) & (dat.MAF <= mafs[maf_bin])]
        ld_quantiles = np.quantile(dat_subset.ldscore_SNP, np.linspace(0,1,lds+1))
        ld_quantiles[0] = ld_quantiles[0] - 0.000001
        for ld_bin in range(1, (lds+1)):
            out = dat_subset[(dat_subset.ldscore_SNP > ld_quantiles[ld_bin - 1]) & (dat_subset.ldscore_SNP <= ld_quantiles[ld_bin])]
            out[['chr', 'SNP']].to_csv(os.path.join(args.out, 'chr' + str(CHR) + '.MAF'+ str(maf_bin) + '.LD' + str(ld_bin) + '.snps'),
                                      index = False, header = False, sep = '\t')
            out[['SNP']].to_csv(os.path.join(args.out, 'chr' + str(CHR) + '.MAF'+ str(maf_bin) + '.LD' + str(ld_bin) + '.snplist'),
                                      index = False, header = False, sep = '\t')


dat = pd.read_csv(args.score, delim_whitespace=True)
create_snplist(dat = dat, CHR = args.chrom, mafs = args.mafs, lds = args.lds)
