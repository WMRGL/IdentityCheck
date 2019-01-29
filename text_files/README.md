## Text files
The text files used in the identity check analysis are as follows:

### Tab delimited text files with chromosome and position of SNPs for each assembly

1. **hg19.nochr.txt** -- e.g.
```
1    67861520
1    179520506
10   100219314
10   104596924
10   104814162
```
2. **hg38.txt** -- e.g.
```
chr1    67395837
chr1    179551371
chr10   98459557
chr10   102837167
chr10   103054405
```

### Identity probabilities and power calculations
The probability of identity (PI) is an estimate of the average probability that 2 independent samples will have the same genotype:

<img src="http://latex.codecogs.com/gif.latex?PI=2\left(\sum_{i}^{n}p_{i}^{2}\right)^{2}-\sum_{i}^{n}p_{i}^{4}" border="0"/>
where p<sub>i</sub> is the allele frequency at the i<sup>th</sup> allele.

The exclusion power is defined as the probability of having a unique assignment and is calculated as follows:

<img src="http://latex.codecogs.com/gif.latex?Power=1-PI" border="0"/>

The text files, one for each assembly, are tab delimited and they contain allele frequencies,identity probabilities and power calculations for each SNP

1. **hg19_total_af_power.txt** -- e.g.
```
ID	GENE	REF	Ref_AF	ALT	Alt_AF	Probability of Identity	Exclusion power 
rs712700	PAX4	T	0.2674	C	0.7326	0.446665412	0.553334588
rs11080572	AFG3L2	C	0.2855	T	0.7145	0.433711909	0.566288091
rs2229546	IL12RB2	C	0.3716	A	0.6284	0.3931174	0.6068826
```
2. **hg38_total_af_power.txt** -- e.g.
```
ID	GENE	REF	Ref_AF	ALT	Alt_AF	Probability of identity	Exclusion power
rs712700	PAX4	T	0.2675	C	0.7325	0.446588719	0.553411281
rs11080572	AFG3L2	C	0.2843	T	0.7157	0.434514776	0.565485224
rs2229546	IL12RB2	C	0.3702	A	0.6298	0.393551179	0.606448821

```
