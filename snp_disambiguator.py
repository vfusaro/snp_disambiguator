#!/usr/bin/env python
# encoding: utf-8
"""
snp_disambiguator.py

Created by Vincent Fusaro on 2013-06-13.
"""

import sys
import argparse
from collections import defaultdict


def main():
    """
    Determines a set of unique haplotypes that can fully identify a given allele.  
    Also reports which alleles cannot be disambiguated based on the full gene halplotype table
    
    A haplotype is considered ambiguous if it is identical to 1 or more alleles and therefore
    cannot be uniquely identified
    
    input: snp_disambiguator <haplotype table> <important alleles>
    output: prints to the screen 
        1: the set of ambiguous alleles and then for each allele the 
        2: for each allele a list of tuples containing the column position (zero-based)
        (from the full allele table) and the REF/SNP.   
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('haplotype_table', 
                        help='specify the full allele haplotype table in tab delimited format')
    parser.add_argument('important_alleles', 
                        help='specify a file with important alleles one per line (single column)')
    args = parser.parse_args()
    
    ### main steps ###
    
    # parse the haplotype table into multiple data structures
    allele_ambiguity_count, haplotypes_by_allele, \
    alleles_by_haplotype, ref_allele, ref_snps = parse_haplotype_table(args.haplotype_table)
    
    # parse the important allele list
    important_alleles = parse_important_alleles(args.important_alleles)
    
    # determine ambiguous and not ambiguous alleles by updating the ambiguity count
    ambiguous_alleles, not_ambiguous_alleles = calc_ambiguity(haplotypes_by_allele, allele_ambiguity_count, 
                                                              important_alleles)
    
    # compute the set difference for each allele
    allele_snp_diff = calc_set_difference(not_ambiguous_alleles, alleles_by_haplotype, ref_allele)
    
    # display the results
    format_output(ambiguous_alleles, allele_snp_diff, ref_snps)

    ### end main steps ###
    
                
def parse_haplotype_table(haplotype_file):
    """
    Parses the haplotype table into multiple data structures for computation downstream.  
    """
    haplotypes_by_allele = defaultdict(list)
    allele_ambiguity_count = defaultdict(int)
    alleles_by_haplotype = {}
    is_reference_row = True
    
    with open(haplotype_file, 'r') as f:
        for line in f:
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
                
    return [allele_ambiguity_count, haplotypes_by_allele, alleles_by_haplotype, ref_allele, ref_snps]


def parse_important_alleles(allele_file):
    """
    Reads the important allele file and returns a list
    """
    
    with open(allele_file, 'r') as f:
        important_alleles = []
        for line in f:
            important_alleles.append(line.strip())
    
    return important_alleles
    

def calc_ambiguity(haplotypes_by_allele, allele_ambiguity_count, important_alleles):
    """
    Determines the uniqueness of each haplotype by counting the number of alleles (value) that
    match the haplotype (key).  Then splits the alleles into ambiguous or not based on the counts
    """
    # updates the allele counts. if > 1 then we can't fully resolve the allele    
    for k,v in haplotypes_by_allele.iteritems():
        for item in v:
            allele_ambiguity_count[item] = len(v)
    
    # splits the alleles into 2 categories: ambiguous and not        
    ambiguous_alleles = []
    not_ambiguous_alleles = []
    for allele in important_alleles:
        if allele_ambiguity_count[allele] > 1:
            ambiguous_alleles.append(allele)
        else:
            not_ambiguous_alleles.append(allele)
    
    return [ambiguous_alleles, not_ambiguous_alleles]
                       

def calc_set_difference(not_ambiguous_alleles, alleles_by_haplotype, ref_allele):
    """
    Uses a set difference to compare the important allele to the haplotype table and returns
    the full set of unique positions/SNPs that can disambiguate between reference and allele
    """
    allele_snp_diff = {}
    for allele in not_ambiguous_alleles:
        if allele in alleles_by_haplotype:
            allele_snp_diff[allele] = set.difference(alleles_by_haplotype[allele],
                                                     alleles_by_haplotype[ref_allele])
        else:
            print 'This important allele (%s) was not found in the full halplotype table' % allele
    
    return allele_snp_diff
            
                       
def format_output(ambiguous_alleles, allele_snp_diff, ref_snps):
    """
    Displays the ambiguous alleles and then formats the output to be the allele follwed by a 
    list of tuples with column number and REF/SNP
    ex: B*13090101	[(68, 'G/A'), (70, 'A/C')]
    """
    print 'Ambiguous: ', ambiguous_alleles
    
    for allele, col_snp_set in allele_snp_diff.iteritems():
        col_snp_lst = []
        for pair in sorted(col_snp_set):
            col_snp_lst.append((pair[0], ref_snps[pair[0]]+'/'+pair[1]))
        print '%s\t%s' % (allele, col_snp_lst)


if __name__ == '__main__':
    main()

