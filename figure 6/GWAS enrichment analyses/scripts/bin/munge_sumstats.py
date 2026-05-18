#!/usr/bin/env python
# coding: utf-8

# In[6]:


import pandas
import os
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='Format proinsulin summary stats file for fGWAS')
    parser.add_argument('--sumstats', required=True, help ="""Sumstats tsv""")
    # parser.add_argument('--snp', required=True, help ="""snp id column name""")
    parser.add_argument('--chr', default="#snp_chrom", help ="""Chrom column name""")
    parser.add_argument('--pos', required=True, help ="""Pos column name""")
    parser.add_argument('--N', required=True, help ="""sample size column""")
    parser.add_argument('--frq', required=True, help ="""SNP allele frequency column""")
    parser.add_argument('--effect', required=True, help ="""Effect column""")
    parser.add_argument('--stderr', required=True, help ="""Standard error column""")
    parser.add_argument('--prefix', required=True, help ="""outputfile """)

    args = parser.parse_args()
    return args

    
if __name__ == '__main__':
    args = getOpts()

    df = pandas.read_csv(args.sumstats, sep='\t', dtype={args.pos: int})
    columns = df.columns
    ## don't care about F and N because we have Z and SE for all gwas. Still keep these dummy columns.
    if args.frq not in columns and "maf" in columns and "minor_allele" in columns:
        df[args.frq] = df.apply(lambda x: x['maf'] if x['minor_allele'] == x['ALT'] else 1 - x['maf'], axis=1)
    if args.N not in columns and "n_cases" in columns:
        df[args.N] = "n_cases"
    cols = [args.chr, args.pos, args.frq, args.N, args.effect, args.stderr]
    rename = {args.chr : "CHR",
              args.pos : "POS",
              args.frq : "F",
              args.N : "N",
              args.stderr: "SE"}
    outcols = ['CHR', 'POS', 'F', 'Z', 'N', 'SE', 'SNPID']

    df = df[cols].rename(columns = rename)
    df['Z'] = df[args.effect]/df['SE']
    df['SNPID'] = df.apply(lambda x: f"{x['CHR']}:{x['POS']}", axis=1)

    df.sort_values(['CHR','POS'], inplace=True)

    # Duplicates (sometimes multiple entries per chrom:pos because there are multiple alleles) are a problem - remove all
    df.drop_duplicates(subset = ['CHR','POS'], keep=False, inplace=True)
    
    df[outcols].to_csv(f"{args.prefix}.txt", sep=' ', index=False, na_rep="NA")
    df['START'] = df['POS'] - 1
    
    df[['CHR', 'START', 'POS']].to_csv(f"{args.prefix}.bed", sep='\t', index=False, header=False)
