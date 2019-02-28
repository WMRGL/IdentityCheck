# IdentityCheck

## Comparison Workflow
![workflow](images/agena_flow.png)


### Identity probabilities and power calculations
The probability of identity (PI) is an estimate of the average probability that 2 independent samples will have the same genotype (_Waits, et.al., 2001_):

<img src="http://latex.codecogs.com/gif.latex?PI=2\left(\sum_{i}^{n}p_{i}^{2}\right)^{2}-\sum_{i}^{n}p_{i}^{4}" border="0"/>
where p<sub>i</sub> is the allele frequency at the i<sup>th</sup> allele.

The exclusion power is defined as the probability of having a unique assignment and is calculated as follows:

<img src="http://latex.codecogs.com/gif.latex?Power=1-PI" border="0"/>

<sup><sub>Waits, Lisette & Luikart, Gordon & Taberlet, Pierre. (2001). Estimating the probability of identity among genotypes in natural populations. _Molecular Ecology Notes 10(1)_</sup></sub>
