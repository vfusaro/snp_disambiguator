SNP Disambiguator
=====

The purpose of this program is to read in a tab delimited file containing haplotype data for a given gene along with another file containing a list of "important" alleles.  The goal is find the set of SNPs that uniquely identify the important alleles and determine which ones cannot be uniquely identified (or disambiguated).

The overall strategy relies on using a set comparison between the important allele and the reference.

##Assumptions
1. An ambiguous haplotype is identical to 1 or more alleles and therefore cannot be uniquely identified
2. The format of the haplotype table has some distinguishing features that identify the alleles and therefore the rows needed.  In this case everything that starts with "B*".
3. The first row starting with "B*" is the reference sequence
4. All unique positions/SNPs are returned for a given allele.  This is probably overkill.  Future work could reduce this to the minimal set
5. The output format is the allele plus a list of the position/column (zero-based) and the REF/SNP. 

##Usage
Example running with a small set of important alleles
```python
snp_disambiguator.py HLAB.txt important_alleles.txt
```

Example output
```python
ambiguous:  ['B*55010301', 'B*44190101N']
B*13090101	[(68, 'G/A'), (70, 'A/C'), (119, 'G/A'), (185, 'A/G'), (187, 'C/G'), (283, 'T/C'), (284, 'C/T'), (285, 'A/T'), (286, 'T/G'), (287, 'C/G'), (292, 'G/C'), (309, 'G/C'), (339, 'T/C'), (342, 'G/A'), (349, 'C/T'), (350, 'C/A'), (416, 'G/C'), (436, 'G/T'), (457, 'A/T'), (489, 'C/G'), (490, 'T/A')]
B*07990101	[(29, 'G/T'), (32, 'A/G'), (68, 'G/T'), (130, 'G/A'), (131, 'A/G'), (132, 'T/A'), (135, 'C/A'), (139, 'C/G'), (148, 'A/G'), (198, 'C/A'), (203, 'A/G'), (206, 'A/C'), (208, 'C/G'), (209, 'A/G'), (218, 'T/G'), (283, 'T/C'), (285, 'A/C'), (293, 'G/C'), (299, 'T/C'), (316, 'C/A'), (339, 'T/C'), (349, 'C/A'), (393, 'A/C'), (407, 'G/C'), (469, 'T/G'), (489, 'C/G'), (490, 'T/A'), (533, 'G/C'), (535, 'C/A'), (540, 'C/G'), (548, 'G/T')]
B*18470101	[(23, 'T/C'), (29, 'G/T'), (32, 'A/G'), (68, 'G/T'), (87, 'A/G'), (107, 'G/T'), (132, 'T/C'), (135, 'C/A'), (148, 'A/G'), (160, 'G/A'), (283, 'T/C'), (285, 'A/C'), (299, 'T/C'), (339, 'T/C'), (416, 'G/C'), (457, 'A/T'), (489, 'C/A'), (490, 'T/C'), (513, 'T/C')]
B*35840101	[(23, 'T/C'), (29, 'G/T'), (32, 'A/G'), (68, 'G/A'), (70, 'A/C'), (132, 'T/C'), (135, 'C/A'), (198, 'C/T'), (309, 'G/C'), (317, 'G/C'), (339, 'T/C'), (349, 'C/T'), (416, 'G/C'), (457, 'A/T')]
```

##Test Cases
Detecting multiple sets of ambiguous cases
```python
snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_1.txt
Ambiguous:  ['B*1A-ambiguous', 'B*1B-ambiguous', 'B*2A-ambiguous', 'B*2B-ambiguous', 'B*2C-ambiguous']
```

Ignoring the "*" (no call) positions and still returning differentiating positions
```python
snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_2.txt
Ambiguous:  []
B*3-nocall	[(2, 'T/C'), (4, 'C/T'), (7, 'A/G')]
```

Two important alleles with a single differentiating position
```python
snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_3.txt
Ambiguous:  []
B*5-impt-single	[(10, 'G/A')]
B*4-impt-single	[(2, 'T/A')]
```

One important allele with two differentiating positions
```python
snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_4.txt
Ambiguous:  []
B*6-impt-double	[(1, 'C/T'), (10, 'G/A')]
```

Evaluating all test cases simultaneously 
```python
snp_disambiguator.py tests/HLAB_test_cases.txt tests/important_alleles_test_5.txt
Ambiguous:  ['B*1A-ambiguous', 'B*1B-ambiguous', 'B*2A-ambiguous', 'B*2B-ambiguous', 'B*2C-ambiguous']
B*5-impt-single	[(10, 'G/A')]
B*4-impt-single	[(2, 'T/A')]
B*6-impt-double	[(1, 'C/T'), (10, 'G/A')]
B*3-nocall	[(2, 'T/C'), (4, 'C/T'), (7, 'A/G')]
```

##Copyright
Copyright (c) 2013 Vincent Fusaro. Released under the MIT License.







