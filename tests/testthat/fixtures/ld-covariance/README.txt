Numerical fixture for the compiled LD reader (not a mocked LD matrix).

There are three SNPs and two populations, P1 (4 individuals) and P2 (12).
Both populations belong to EUR. Genotype data, one SNP per line:

0022 002200220022
0022 220022002200
0002 112211221122

index.bgz contains:
rs1 1 10 A G 0.5 0
rs2 1 20 A G 0.5 18
rs3 1 30 A G 0.5 36

Both .bgz files contain one BGZF data block followed by an empty EOF block.
Each block uses raw DEFLATE, the standard BC extra field with total block
size minus one, and a CRC32/ISIZE trailer. Because the data fit in the first
block, index virtual offsets equal uncompressed byte offsets (0, 18, 36).
The index's af1ref column is unused by this LD reader.

SNPs 1/2 have opposite LD in the two populations and identical means.
Their corrected correlation is -8/19, using sample covariances weighted
by population sample proportions. The old m-squared inflation gave
-0.9133858268. SNP3 has different means across populations to exercise
the between-population covariance term as well.
