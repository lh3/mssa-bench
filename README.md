This repo evaluates the performance of three multi-string suffix array (SA)
construction methods, [libsais][libsais], [gSACA-K][gsacak] and [ksa].

Let $`\{T_1,T_2,\ldots,T_n\}`$ be a set of strings over $\Sigma$. Their
concatenation is $`T=T_1\$_1T_2\$_2\cdots T_n\$_n`$ where
$`\$_1<\$_2<\cdots\$_n`$ are smaller than all symbols in $\Sigma$.
We would like to compute the suffix array of $T$.

Ksa was adapted from an ancient version of [Yuta Mori][mori]'s sais in 2011. During
symbol comparisons, it implicitly replaces a sentinel $`\$_i`$ with
$`j-|T|`$ where $j$ is the offset of $`\$_i`$ in $T$. Ksa needs at least 9
bytes per element for 64-bit integers.

With libsais, we generate an integer array $`X`$, such that $`X[k]=i`$
if $`T[k]=\$_i`$, or $`X[k]=T[k]+n`$ otherwise. We apply libsais to the integer
array. This approach needs at least 16 bytes per element.

For gSACA-K, see [its paper][gsacak-paper].

Here is the timing for constructing the [CHM13v2 genome][chm13] on both strand (6.1
billion symbols in total) on a Xeon Gold 6130:

|             | ksa|gsaca-k|sais-t1|sais-t4|sais-t8|
|:------------|---:|------:|------:|------:|------:|
|# threads    |   1|      1|      1|      4|      8|
|Elapsed (s)  |1396|crashed|    588|    386|    260|
|CPU time (s) |1395|crashed|    587|   1152|   1439|
|Peak RSS (GB)|52.3|crashed|   92.9|   92.9|   92.9|

gSACA-K worked on small input but crashed for CHM13.

[libsais]: https://github.com/IlyaGrebnov/libsais
[chm13]: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/analysis_set/
[mori]: https://github.com/y-256
[gsacak]: https://github.com/felipelouza/gsa-is
[gsacak-paper]: https://www.sciencedirect.com/science/article/pii/S0304397517302621
[ksa]: https://github.com/lh3/fermi/blob/master/ksa.c
