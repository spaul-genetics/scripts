import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--chrom', type = int, default = 21, help = 'Enter the chromosome number.')
parser.add_argument('--mfi_loc', type = str,
                    default = '/work/IBG/UKBiobank/Release3/EGAD00010001474/',
                   help = 'Location of the mfi files.')

parser.add_argument('--mfi_prefix', type = str,
                   default = 'ukb_mfi_chr')

parser.add_argument('--mfi_suffix', type = str,
                   default = '_v3.txt')
parser.add_argument('--mac', type  = int, default = 2,
                   help = 'Any variant with less than allele count of less than mac will be removed.')
parser.add_argument('--INFO', type = float, default = 0.3,
                   help = 'Any variant with info score less than INFO will be removed.')
parser.add_argument('--hwe', type = float, default = 1e-6,
                   help = 'Any variant with Hardy-Weienburg p-value < hwe will be removed')
parser.add_argument('--hardy_loc', type = str, default = '/work/IBG/spaul/greml_ukb/bgen/hardy/',
                   help = 'Location of the --hardy output from plink2. The file names whould be of this format chr*.hardy')
parser.add_argument('--out', type = str, default = '',
                   help = 'Directory to write the SNP lists after QC.')
args = parser.parse_args()

mfi = pd.read_csv(args.mfi_loc + args.mfi_prefix + str(args.chrom) + args.mfi_suffix, sep = '\t',
                  names = ['SNP', 'SNP_RS', 'POS', 'A1', 'A2', 'MAF', 'MA', 'INFO'])



hardy = pd.read_csv(args.hardy_loc + 'chr' + str(args.chrom) + '.hardy', sep = '\t')
remove_hwe_snps = hardy.ID[hardy.P<args.hwe]

rem_dup = mfi.shape[0]
mfi = mfi.drop_duplicates('POS', keep = 'first')
hwe_qc = ~mfi.SNP_RS.isin(remove_hwe_snps)
rem_dup = rem_dup - mfi.shape[0]
print(str(rem_dup) + ' SNPS removed for duplicate positions. The first of the duplicated entries are kept.')
mac2 = mfi.MAF > (args.mac -1)/(2*487411)
#mac2 = mfi.MAF >=0.000003
print('After removing the duplicates ' + str(np.sum(~mac2)) + ' singletons are removed.')
snps_only = np.array([len(x)<2 for x in mfi.A1]) & np.array([len(x)<2 for x in mfi.A2])
print('After removing the duplicates ' +str(np.sum(~snps_only)) + ' indels or multiple nucleotides are removed.')
print('After removing the duplicates ' +str(len(remove_hwe_snps)) + ' snps are removed for hwe p < 1e-6.')
info_qc = mfi.INFO >= args.INFO
print('After removing the duplicates ' + str(np.sum(~info_qc)) + ' variants are remove for INFO < ' + str(args.INFO))
after_qc = mfi[mac2 & snps_only & hwe_qc & info_qc].copy()
after_qc['chr'] = args.chrom
print(str(after_qc.shape[0]) + ' SNPs remaining after QC')

after_qc[['chr','SNP']].to_csv(args.out + 'qc.snps.chr' + str(args.chrom) + '.snplist',
               index = False, sep = '\t', header = False)

after_qc[['chr','SNP_RS']].to_csv(args.out + 'qc.snps.chr' + str(args.chrom) + '.rssnplist',
               index = False, sep = '\t', header = False)


with open(args.out + 'qc.snps.chr' + str(args.chrom) + '.qctoolist', 'w') as file:
    for x in after_qc.POS:
        if(args.chrom < 10):
            CHR = '0' + str(args.chrom)
        else:
            CHR = str(args.chrom)
        xx = CHR + ':' + str(x)
        file.write('%s ' % xx)