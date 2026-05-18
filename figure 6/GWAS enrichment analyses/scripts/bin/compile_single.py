#!/usr/bin/env python
# coding: utf-8

# In[6]:


import pandas
import numpy
import os
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='Format proinsulin summary stats file for fGWAS')
    parser.add_argument('--pars', nargs='+', required=True, help ="""params files""")
    parser.add_argument('--llk', nargs='+', required=True, help ="""llk files""")
    parser.add_argument('--trait', required=True, help ="""trait name""")
    parser.add_argument('--output', required=True, help ="""outputfile """)

    args = parser.parse_args()
    return args

def annot(fname, suffix, replace):
    annotation = os.path.basename(fname).replace(suffix, "").replace(replace, "")
    return annotation

def fix(filename, trait):
    """Fix .params output dataframe """
    
    annotation = annot(filename, ".params", f"{trait}__annot__")
    d = pandas.read_csv(filename, delim_whitespace=True)
    d['trait'] = trait
    d['annotation'] = annotation
    d = d.tail(1)
    return d


def fix_likelihood(filename, trait):
    """Fix .llk output dataframe """

    annotation = annot(filename, ".llk", f"{trait}__annot__")
    d = pandas.read_csv(filename, delim_whitespace=True, header=None, names=['type','llk'])
    d = d[d['type'] == "ln(lk):"]
    d['trait'] = trait
    d['annotation'] = annotation
    return d


def fix_string(x):
    """Convert fGWAS output confidence intervals to float """
    try:
        out = float(str(x).replace(">", "").replace("<", ""))
    except ValueError:
        out = numpy.nan
    return out

    
if __name__ == '__main__':
    args = getOpts()

    trait = args.trait

    # Compile results
    out = pandas.concat([fix(filename, trait) for filename in args.pars], ignore_index=True)
    llk = pandas.concat([fix_likelihood(filename, trait) for filename in args.llk], ignore_index=True)
    
    m = pandas.merge(out, llk, how="inner", on=['trait','annotation'])
    print(m.head())
    m['CI_lo'] = m['CI_lo'].map(fix_string)
    m['CI_hi'] = m['CI_hi'].map(fix_string)
    m['sig'] = m.apply(lambda x: 1 if (x['estimate']>0 and x['CI_lo']>0) or (x['estimate']<0 and x['CI_hi']<0) else 0, axis=1)

    print(m.head())
    m.to_csv(args.output, index=False, sep='\t', na_rep="NA")
    
