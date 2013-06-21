SNP Disambiguator
=====

The purpose of this program is to read in a tab delimited file containing haplotype data for a given gene along with
another file containing a list of "important" alleles.  The goal is find the set of SNPs that uniquely identify the
important alleles and determine which ones cannot be uniquely identified (or disambiguated) from the unimportant set.
Ultimately the goal is to find the minimal haplotype set that requires genotyping the fewest SNPs.

The overall strategy relies on using a set comparison between the important allele and the reference and then a
permutation strategy to test the minimal combinations of SNPs necessary.

##Assumptions
1. An ambiguous haplotype is identical to 1 or more alleles and therefore cannot be uniquely identified
2. The format of the haplotype table has some distinguishing features that identify the alleles and therefore the rows
needed.  In this case everything that starts with "B*".
3. The first row starting with "B*" is the reference sequence
4. A minimal haplotype set is returned if the total number of unique SNPs is less than 15 for computational time.  If
the number of unique SNPs is greater then the full haplotype set is returned.
5. The output format is the allele plus a list of the position/column (zero-based) and the REF/SNP. 

##Usage
Example running with a small set of important alleles
```python
snp_disambiguator.py HLAB.txt important_alleles.txt
```

Example output
```python
Ambiguous:  []
Minimum haplotype (columns)  (139, 185, 248, 283, 299, 314, 339, 407)
B*15020301	[(407, 'G/C')]
B*15020201	[(139, 'C/G')]
B*95590101	[(185, 'A/G'), (187, 'C/G'), (248, 'T/C'), (283, 'T/C'), (285, 'A/C'), (299, 'T/C'), (339, 'T/C'), (468, 'C/T'), (469, 'T/G')]
B*15020401	[(314, 'G/T')]
```

##Test Cases
Detecting multiple sets of ambiguous cases
```python
python snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_1.txt
Ambiguous:  ['B*1A-ambiguous', 'B*1B-ambiguous', 'B*2A-ambiguous', 'B*2B-ambiguous', 'B*2C-ambiguous']
Minimum haplotype (columns)  []
```

Ignoring the "*" (no call) positions and still returning differentiating positions
```python
snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_2.txt
Ambiguous:  []
Minimum haplotype (columns)  [2, 4, 7]
B*3-nocall	[(2, 'T/C'), (4, 'C/T'), (7, 'A/G')]
```

Two important alleles with a single differentiating position
```python
snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_3.txt
Ambiguous:  []
Minimum haplotype (columns)  [2, 10]
B*5-impt-single	[(10, 'G/A')]
B*4-impt-single	[(2, 'T/A')]
```

One important allele with two differentiating positions
```python
snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_4.txt
Ambiguous:  []
Minimum haplotype (columns)  (1,)
B*6-impt-double	[(1, 'C/T'), (10, 'G/A')]
```

An important allele is not found in the haplotype table
```python
snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_5.txt
This important allele (B*7-NOT-FOUND) was not found in the full halplotype table
Ambiguous:  []
Minimum haplotype (columns)  []
```

Evaluating all test cases except B*7-impt-single simultaneously
```python
snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_6.txt
This important allele (B*7-NOT-FOUND) was not found in the full halplotype table
Ambiguous:  ['B*1A-ambiguous', 'B*1B-ambiguous', 'B*2A-ambiguous', 'B*2B-ambiguous', 'B*2C-ambiguous']
Minimum haplotype (columns)  (1, 2, 10)
B*5-impt-single	[(10, 'G/A')]
B*4-impt-single	[(2, 'T/A')]
B*6-impt-double	[(1, 'C/T'), (10, 'G/A')]
B*3-nocall	[(2, 'T/C'), (4, 'C/T'), (7, 'A/G')]
```

Evaluating all test cases simultaneously.  Notice this time we require an extra column #4 with the inclusion of
B*7 as compared to the example above.
```python
snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_7.txt
Ambiguous:  ['B*1A-ambiguous', 'B*1B-ambiguous', 'B*2A-ambiguous', 'B*2B-ambiguous', 'B*2C-ambiguous']
Minimum haplotype (columns)  (1, 2, 4, 10)
B*5-impt-single	[(10, 'G/A')]
B*4-impt-single	[(2, 'T/A')]
B*7-impt-single	[(4, 'C/G')]
B*6-impt-double	[(1, 'C/T'), (10, 'G/A')]
B*3-nocall	[(2, 'T/C'), (4, 'C/T'), (7, 'A/G')]
```

Running a small example within the actual dataset.  Here, the total haplotype set is 12 but gets reduced to 8.
Note: this takes a few seconds to run.  ;)
```python
snp_disambiguator.py HLAB.txt tests/important_alleles_test_8.txt
Ambiguous:  []
Minimum haplotype (columns)  (139, 185, 248, 283, 299, 314, 339, 407)
B*15020301	[(407, 'G/C')]
B*15020201	[(139, 'C/G')]
B*95590101	[(185, 'A/G'), (187, 'C/G'), (248, 'T/C'), (283, 'T/C'), (285, 'A/C'), (299, 'T/C'), (339, 'T/C'), (468, 'C/T'), (469, 'T/G')]
B*15020401	[(314, 'G/T')]
```

##Copyright
Copyright (c) 2013 Vincent Fusaro. Released under the MIT License.







