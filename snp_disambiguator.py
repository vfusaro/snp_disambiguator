#!/usr/bin/env python
# encoding: utf-8
"""
snp_disambiguator.py

Created by Vincent Fusaro on 2013-06-13.
"""

import sys
from collections import defaultdict


def main():
    """
    Determines a set of unique haplotypes that can fully identify a given allele.  
    Also reports which alleles cannot be disambiguated based on the full gene halplotype table
    
    A haplotype is considered ambiguous if it is identical to 1 or more alleles and therefore
    cannot be uniquely identified
    
    input: snp_disambiguator <full allele table> <important alleles>
    output: prints to the screen 
        1: the set of ambiguous alleles and then for each allele the 
        2: for each allele a list of tuples containing the column position (zero-based)
        (from the full allele table) and the REF/SNP.   
    """
    
    # argparse would make this more user friendly and scalable if additional inputs are needed
    try:
        allele_table_file = sys.argv[1]
    except:
        raise Exception('you forgot to provide the full allele table')
    try:
        important_alleles_file = sys.argv[2]
    except:
        raise Exception('you forgot to provide the file of important alleles')
    
    haplotypes_by_allele = defaultdict(list)
    allele_ambiguity_count = defaultdict(int)
    alleles_by_haplotype = {}
    is_reference_row = True
    
    with open(allele_table_file, 'r') as atf:
        for line in atf:
            if line.startswith('B*'):
                snps = line.strip().split('\t')
                allele = snps.pop(0)
                if is_reference_row:
                    ref_allele = allele
                    ref_snps = snps
                    is_reference_row = False
                cols = range(len(snps))
                allele_ambiguity_count[allele]
                haplotypes_by_allele[frozenset(zip(cols, snps))].append(allele)
                
                col_snp = zip(cols, snps)

                # remove unwanted characters
                col_snp = [i for i in col_snp if i[1] not in ['_', '*']]

                alleles_by_haplotype[allele] = set(col_snp)
    
    # updates the allele counts. if > 1 then we can't fully resolve the allele    
    for k,v in haplotypes_by_allele.iteritems():
        for item in v:
            allele_ambiguity_count[item] = len(v)
    
    with open(important_alleles_file, 'r') as iaf:
        important_alleles = []
        for line in iaf:
            important_alleles.append(line.strip())
    
    # splits the alleles into 2 categories: ambiguous and not        
    ambiguous_alleles = []
    not_ambiguous_alleles = []
    for allele in important_alleles:
        if allele_ambiguity_count[allele] > 1:
            ambiguous_alleles.append(allele)
        else:
            not_ambiguous_alleles.append(allele)
    
    print 'ambiguous: ', ambiguous_alleles
    # print 'not ambiguous:', not_ambiguous_alleles
    # print 'ref allele name: ', ref_allele
    
    for allele in not_ambiguous_alleles:
        if allele in alleles_by_haplotype:
            set_diff = set.difference(alleles_by_haplotype[allele], 
                                      alleles_by_haplotype[ref_allele])
            format_output(allele, set_diff, ref_snps)
        else:
            print 'This important allele (%s) was not found in the full halplotype table' % allele
                
                
def format_output(allele, col_snp_set, ref_snps):
    """
    Formats the output to be the allele follwed by a list of tuples with column number and REF/SNP
    ex: B*13090101	[(68, 'G/A'), (70, 'A/C')]
    """
    col_snp_lst = []
    for pair in sorted(col_snp_set):
        col_snp_lst.append((pair[0], ref_snps[pair[0]]+'/'+pair[1]))
    print '%s\t%s' % (allele, col_snp_lst)


if __name__ == '__main__':
    main()

