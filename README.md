This repo evaluates the performance of portable libraries for constructing the
suffix array (SA) of string sets. These libraries only include a few source
files and have no dependencies. They can be building blocks of larger projects.

There are [different ways][ss-review] to define the SA of a string set. We
here focus on the most common definition as follows.  Let
$`\mathcal{T}=\{T_1,T_2,\ldots,T_n\}`$ be a set of strings over $\Sigma$. Their
concatenation is $`T=T_1\$_1T_2\$_2\cdots T_n\$_n`$ where
$`\$_1<\$_2<\cdots<\$_n`$ are smaller than all symbols in $\Sigma$. The SA of
string set $`\mathcal{T}`$ is defined as the SA of string $T$.

Few SA construction libraries directly support string sets. Nonetheless, we can
achieve the goal for any libraries that support integer alphabets, such as
[libsais][libsais], by converting $T$ to an integer array $X$:
```math
X[k]=\left\{\begin{array}{ll}
i & \mbox{if }T[k]=\$_i \\
T[k]+n & \mbox{otherwise}
\end{array}\right.
```
Then the SA of $X$ will be identical to the SA of $T$. A disadvantage of this
method is that we need to convert 8-bit characters to 32-bit or 64-bit
integers. This increases the memory footprint.

To alleviate the issue, I developed [ksa][ksa] in 2011 by adapting an old
version of [Yuta Mori][mori]'s sais. Briefly, during symbol comparisons, ksa
implicitly replaces a sentinel $`\$_i`$ with $`j-|T|`$ where $j$ is the offset
of $`\$_i`$ in $T$. The comparison between symbols takes more time but we do
not need to convert $T$ to integer arrays anymore and can thus save memory.
This repo includes an updated version named [msais][msais].

[Published in 2017][gsacak-paper], [gSACA-K][gsacak] is another library based
on the linear-time SAIS algorithm. Please read its paper for details.

Here is the timing for constructing the [CHM13v2 genome][chm13] on both strand (6.2
billion symbols in total) on a Xeon Gold 6130:

|             |msais|gsaca-k|sais-t1|sais-t4|sais-t8|sais-t8b|sais-t8c|sais16-t8c|
|:------------|---:|------:|------:|------:|------:|-------:|-------:|---------:|
|# threads    |   1|      1|      1|      4|      8|       8|       8|         8|
|Elapsed (s)  |1211|   3356|    588|    386|    260|     374|     473|       296|
|CPU time (s) |1209|   3349|    587|   1152|   1439|    1895|    2602|      1146|
|Peak RSS (GB)|52.3|   53.5|   92.9|   92.9|   92.9|    92.9|    92.9|      58.4|

Some notes and observations:

* libsais is clearly the fastest even on a single thread and we see noticeable
  speedup with multiple threads. A caveat is that the multi-threading
  performance of libsais appears to have large fluctuation. For example,
  the three sais-t8 were run on different nodes with the same configuration but
  the speed was quite different. sais-t8c and sais16-t8c were run on the same
  machine.

* gSACA-K would crash if compiled with `-fopenmp`. I am not sure why.

* msais is faster than gSACA-K and has the same memory footprint. It would be
  good to apply this msais strategy to libsais to reduce its peak memory.

* We omitted [ropebwt2][rb2] and [BEETL][beetl] because they are slow for
  chromosome-long strings and we omitted [grlBWT][grl] because it writes
  temporary files and is not designed as a library. We did not evaluate
  [eGAP][egap] because gSACA-K appears to be faster in multiple third-party
  benchmarks.

[libsais]: https://github.com/IlyaGrebnov/libsais
[chm13]: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/analysis_set/
[mori]: https://github.com/y-256
[gsacak]: https://github.com/felipelouza/gsa-is
[gsacak-paper]: https://www.sciencedirect.com/science/article/pii/S0304397517302621
[ksa]: https://github.com/lh3/fermi/blob/master/ksa.c
[fermi]: https://github.com/lh3/fermi
[fermi-paper]: https://academic.oup.com/bioinformatics/article/28/14/1838/218887
[ss-review]: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btae333/7681884
[rb2]: https://github.com/lh3/ropebwt2
[grl]: https://github.com/ddiazdom/grlBWT
[beetl]: https://github.com/BEETL/BEETL
[egap]: https://github.com/felipelouza/egap
[msais]: https://github.com/lh3/msais-lite
