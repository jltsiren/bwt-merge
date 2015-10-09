# BWT-merge

This is a tool for merging the Burrows-Wheeler transforms (BWT) of large read collections. Querying the merged BWT is often faster than querying the BWT of each dataset separately. If the datasets are similar (e.g. reads from similar genomes), a run-length encoded merged BWT will usually be smaller than the run-length encoded BWTs of individual datasets.

If the dataset is up to a few hundred gigabytes in size, it is probably more practical to use another tool (see "Algorithms and implementations" below) to build the merged BWT directly.

## Usage

BWT-merge is based on the [Succinct Data Structures Library 2.0 (SDSL)](https://github.com/simongog/sdsl-lite). To compile, set `SDSL_DIR` in the Makefile to point to your SDSL directory. The program should compile with a relatively new g++ on both Linux and OS X. It has not been tested with other compilers. Comment out the line `OUTPUT_FLAGS=-DVERBOSE_STATUS_INFO` if you do not want the merging tool to output status information to `stderr`.

There are three tools in the package:

`bwt_convert [options] input output` reads a run-length encoded BWT built by the [String Graph Assembler](https://github.com/jts/sga) from file `input` and writes it to file `output` in the native format of BWT-merge. The converted file is often a bit smaller than the input, even though it includes rank/select indexes. The input/output formats can be changed with options `-i format` and `-o format`.

`bwt_inspect input1 [input2 ...]` tries to identify the BWT formats of the input files. If successful, it will also display some basic information about the files. Only the native format and the SGA format are currently supported.

`bwt_merge [options] input1 input2 [input3 ...] output` reads the input BWT files, merges them, and writes the merged BWT to file `output`. The sequences from each input file are inserted after the sequences from the BWTs that have already been merged. In most cases, the input files should be given from the largest to the smallest. There are several options:

* `-r N` sets the size of **run buffers** to *N* megabytes (default 128). The unsorted run buffers are thread-specific and contain 16-byte values.
* `-b N` sets the size of **thread buffers** to *N* megabytes (default 512). When the run buffer becomes full, its contents are sorted, compressed, and merged with the thread buffer.
* `-m B` sets the number of **merge buffers** to *N* (default 5). The merge buffers are global and numbered from *0* to *N-1*. When a thread buffer becomes full, its contents are merged with one or more merge buffers. Merge buffer *i* contains *2^i* thread buffers. If there is no room in the merge buffers, all *2^N* thread buffers are merged and written to disk.
* `-t N` sets the number of **threads** to *N*. The default is the maximum number of threads OpenMP is allowed to use.
* `-s N` sets the number of **sequence blocks** to *N* (default 4 per thread). Each block consists of roughly the same number of sequences, and the blocks are assigned dynamically to individual threads.
* `-d directory` sets the **temporary directory** (default: working directory).
* `-v patterns` **verifies** the merged BWT by querying it with patterns and comparing the results with those from the inputs. File `patterns` contains one pattern per line.
* `-i formats` speficies **input formats** (default: `native`). Multiple comma-separated formats can be specified.
* `-o format` specifies the **output format** (default: `native`).

## BWT file formats

Many different BWT file formats exist. Some of them are mutually incompatible, while others are essentially different encodings of the same information. There are two main sources of incompatibility:

* **Alphabet:** Many bioinformatics tools use the alphabet `$ACGTN` (in that order). Tools written by computer scientists often assume that the alphabet is sorted, e.g. `$ACGNT`. Because the alphabetic order between characters `N` and `T` is different, the BWTs are incompatible.
* **Endmarkers:** Some BWTs contain explicit endmarkers (usually denoted by `$`) at the end of each sequence, while others handle the end of sequence implicitly. The endmarker is usually the first character in the alphabet, but it can also be the last character.

The native format of BWT-merge uses numeric values 0-5 as its internal alphabet. By default, these **comp values** are interpreted as `$ACGTN`. Other alphabets of size at most 6 can also be supported, as long as the endmarker is the first character in the alphabet. The following formats are currently supported.

|Format         |Alphabet|Details
|:-------------:|:------:|-------
|`native`       |any     |run-length encoded; includes rank/select support
|`plain_default`|default |array of characters
|`plain_sorted` |sorted  |array of characters
|`rfm`          |sorted  |[Relative FM-index](https://github.com/jltsiren/relative-fm): `int_vector_buffer<8>` of comp values
|`ropebwt`      |default |[RopeBWT](https://github.com/lh3/ropebwt): byte array with 3 bits for the character and 5 bits for the length of the run
|`sdsl`         |sorted  |[SDSL](https://github.com/simongog/sdsl-lite): `int_vector_buffer<8>` of characters
|`sga`          |default |[SGA](https://github.com/jts/sga): byte array with 3 bits for the character and 5 bits for the length of the run

## Performance (v0.2.1)

The input consists of several BWT files from the ReadServer project. Each of the files contains unique (error-corrected, trimmed) reads ending with the bases in the file name.

|File         |     AA|     TT|     AT|     TA|
|-------------|:-----:|:-----:|:-----:|:-----:|
|Size         |438 Gbp|436 Gbp|278 Gbp|359 Gbp|
|Sequences    |  4.69G|  4.67G|  2.98G|  3.84G|
|SGA format   |39.9 GB|40.0 GB|26.9 GB|33.5 GB|
|SGA + index  |48.3 GB|48.4 GB|32.6 GB|40.6 GB|
|Native format|38.5 GB|38.7 GB|26.6 GB|32.7 GB|

|Files          |     AA|   AA + TT| AATT + AT|AATTAT + TA|
|---------------|:-----:|:--------:|:--------:|:---------:|
|Size           |438 Gbp|   874 Gbp|  1.15 Tbp|   1.51 Tbp|
|Sequences      |  4.69G|     9.37G|     12.3G|      16.2G|
|Native format  |38.5 GB|   71.3 GB|   91.5 GB|     117 GB|
|Time (step)    |      –|    12.4 h|    10.0 h|     12.5 h|
|Speed (step)   |      –|9.79 Mbp/s|7.74 Mbp/s| 7.97 Mbp/s|
|Memory (step)  |      –|    150 GB|    171 GB|     202 GB|
|Disk (step)    |      –|    267 GB|    220 GB|     289 GB|
|Time (total)   |      –|    12.4 h|    22.3 h|     34.9 h|
|Speed (average)|      –|9.79 Mbp/s|8.88 Mbp/s| 8.55 Mbp/s|
|Memory (peak)  |      –|    150 GB|    171 GB|     202 GB|
|Disk (peak)    |      –|    267 GB|    267 GB|     289 GB|

The experiments were run on a system with two 16-core AMD Opteron 6378 processors and 256 gigabytes of memory. The parameters were 32 **threads**, 128 **sequence blocks**, 128-megabyte **run buffers**, 512-megabyte **thread buffers**, and 5 **merge buffers**. The measured times include disk I/O and index verification.

Some observations:
* Larger `input2` files increase merging speed, while larger `input1` files decrease it. (A larger `input2` file also means that there is more work to do.)
* Memory usage is nondeterministic. The peak seems to be around **2.5 GB/thread** in addition to the input BWTs.
* Disk usage depends primarily on the size of `input2`, but also on the size of `input1`. A good approximation is **5–7 bits/base** in `input2`.

## Background

Building the BWT is a solved problem for sequence collections of up to hundreds of gigabytes in size. For larger collections, there are several issues to consider:

* **Construction time:** As a rough guideline, an algorithm indexing 1 MB/s is good for up to 100 gigabytes of data, and somewhat useful until 1 terabyte. Larger datasets require faster algorithms.
* **Construction space:** When the datasets are much larger than the amount of memory available, we cannot afford using even a single bit of working space per input character.
* **Available hardware:** In a typical computer cluster, a single node has two CPUs (with up to tens of cores), tens to hundreds of gigabytes of local memory, a limited amount of local disk space, and large amounts of shared (but often slow) disk space. Some tools require fast GPUs or large amounts of fast disk space, which are generally not available.
* **Resource usage:** Merging large BWTs is easy by doing a lot of redundant work on multiple nodes. Because computer clusters generally do not have large amounts of unused capacity, good tools should make an efficient use of resources.

## Algorithms and implementations

There are several BWT construction algorithms based on updating an existing BWT. Some algorithms **extend** sequences already in the collection, updating *BWT(T[i+1,j])* to *BWT(T[i,j])*. Others **insert** new sequences to the collection, updating *BWT(S)* to *BWT(S,T)*. In both cases, the algorithm can either do **batch updates** to a **static** BWT representation, or use a **dynamic** representation for the BWT. Algorithms with a static BWT representation require more working space, while algorithms with a dynamic representation have more space overhead for the BWT.

Some common implementations include:

* [BEETL](https://github.com/BEETL/BEETL): (extend, static) on disk
* [NVBIO](http://nvlabs.github.io/nvbio/): (insert, dynamic) using GPU
* [RLCSA](http://jltsiren.kapsi.fi/rlcsa): (insert, static)
* [ropebwt](https://github.com/lh3/ropebwt): (extend, static) or (insert by extending, dynamic)
* [ropebwt2](https://github.com/lh3/ropebwt2): (insert by extending, dynamic)

As this tool merges existing BWTs, it is (insert, static).

There are also other algorithms for building the BWT for large read collections They are based on partitioning the suffixes, sorting each of the partitions separately, and building the BWT directly.

## Version history

### Current version

* Switched from OpenMP to C++11 threads.
* More space-efficient rank/select construction for the BWT.
* Formats: RopeBWT format, faster writing in SGA format.
* `bwt_merge`: Multiple input files, faster RA/BWT merging, multithreaded verification, adjustable input/output formats and temp directory.

### Version 0.2.1

* Minor performance update.
* Faster construction of the rank array and the rank/select structures.

### Version 0.2

* The second pre-release.
* The native BWT format has a header.
* `bwt_inspect` (former `sga_inspect`): Recognizes files in the native format.
* `bwt_convert`: Converts between multiple BWT formats.
* `bwt_merge`: Optimizations, adjustable merge parameters.

### Version 0.1

* The first pre-release.
* `sga_inspect` for inspecting BWT files in the SGA format.
* `bwt_convert` for converting BWT files from the SGA format to the native format.
* `bwt_merge` for merging BWT files in the native format.

## Future work

* Option to load the BWT into a single array to speed up queries.
* New query: extract a sequence based on the lexicographic rank of a suffix.
* RA/BWT merging using more than two threads.
* `bwt_convert`, `bwt_merge`: Input from stdin / output to stdout.
* `bwt_merge`: Different options for where the sequences from `input2` are inserted:
  * after the sequences from `input1` (current behavior)
  * in reverse lexicographic order
  * by position in the reference
* `bwt_merge`: Option to remove duplicate sequences.
* `bwt_convert`: Build rank/select only when necessary.
* Documentation in the wiki.

## References

Wing-Kai Hon, Tak-Wah Lam, Kunihiko Sadakane, Wing-Kin Sung, Siu-Ming Yiu:
**A space and time efficient algorithm for constructing compressed suffix arrays**.
Algorithmica 48(1): 23-36, 2007.
[DOI: 10.1007/s00453-006-1228-8](http://dx.doi.org/10.1007/s00453-006-1228-8)

Jouni Sirén: **Compressed Full-Text Indexes for Highly Repetitive Collections**.
PhD Thesis, University of Helsinki, June 2012.
[http://jltsiren.kapsi.fi/phd](http://jltsiren.kapsi.fi/phd)

Markus J. Bauer, Anthony J. Cox, and Giovanna Rosone:
**Lightweight algorithms for constructing and inverting the BWT of string collections**.
Theoretical Computer Science 483: 134-148, 2013.
[DOI: 10.1016/j.tcs.2012.02.002](http://dx.doi.org/10.1016/j.tcs.2012.02.002)
