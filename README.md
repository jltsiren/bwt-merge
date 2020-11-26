# BWT-merge

This is a tool for merging the Burrows-Wheeler transforms (BWT) of large read collections. Querying the merged BWT is often faster than querying the BWT of each dataset separately. If the datasets are similar (e.g. reads from similar genomes), a run-length encoded merged BWT will usually be smaller than the run-length encoded BWTs of individual datasets.

If the dataset is up to a few hundred gigabytes in size, it is probably more practical to use another tool ([see the wiki](https://github.com/jltsiren/bwt-merge/wiki/BWT-Construction)) to build the merged BWT directly.

Further documentation can be found [in the wiki](https://github.com/jltsiren/bwt-merge/wiki).

## Notes

* This is a proof-of-concept implementation of the BWT-merge algorithm. There are other implementations for merging [GBWTs](https://github.com/jltsiren/gbwt) and [r-indexes](https://github.com/marco-oliva/rimerge).
* There is an error in Figure 1 of the paper. The last two suffixes of R (in lexicographic order) are in wrong order. This affects the final values in a number of arrays.

## Usage

BWT-merge is based on the [Succinct Data Structures Library 2.0 (SDSL)](https://github.com/simongog/sdsl-lite). To compile, set `SDSL_DIR` in the Makefile to point to your SDSL directory. The program should compile with g++ 4.7 or later on both Linux and OS X. It has not been tested with other compilers. Comment out the line `OUTPUT_FLAGS=-DVERBOSE_STATUS_INFO` if you do not want the merging tool to output status information to `stderr`.

There are three tools in the package:

`bwt_convert [options] input output` reads a run-length encoded BWT built by the [String Graph Assembler](https://github.com/jts/sga) from file `input` and writes it to file `output` in the native format of BWT-merge. The converted file is often a bit smaller than the input, even though it includes rank/select indexes. The input/output formats can be changed with options `-i format` and `-o format`.

`bwt_inspect input1 [input2 ...]` tries to identify the BWT formats of the input files. If successful, it will also display some basic information about the files. Only the native format, the RopeBWT format, and the SGA format are currently supported.

`bwt_merge [options] input1 input2 [input3 ...] output` reads the input BWT files, merges them, and writes the merged BWT to file `output`. The sequences from each input file are inserted after the sequences from the BWTs that have already been merged. In most cases, the input files should be given from the largest to the smallest. There are several options:

* `-r N` sets the size of **run buffers** to *N* megabytes (default 128). The unsorted run buffers are thread-specific and contain 16-byte values.
* `-b N` sets the size of **thread buffers** to *N* megabytes (default 256). When the run buffer becomes full, its contents are sorted, compressed, and merged with the thread buffer.
* `-m N` sets the number of **merge buffers** to *N* (default 6). The merge buffers are global and numbered from *0* to *N-1*. When a thread buffer becomes full, its contents are merged with one or more merge buffers. Merge buffer *i* contains *2^i* thread buffers. If there is no room in the merge buffers, all *2^N* thread buffers are merged and written to disk.
* `-t N` sets the number of **threads** to *N*. The default is the number of hardware contexts (~CPU cores) returned by `std::thread::hardware_concurrency`.
* `-s N` sets the number of **sequence blocks** to *N* (default 4 per thread). Each block consists of roughly the same number of sequences, and the blocks are assigned dynamically to individual threads.
* `-d directory` sets the **temporary directory** (default: working directory).
* `-v patterns` **verifies** the merged BWT by querying it with patterns and comparing the results with those from the inputs. File `patterns` contains one pattern per line.
* `-i formats` speficies **input formats** (default: `native`). Multiple comma-separated formats can be specified.
* `-o format` specifies the **output format** (default: `native`).

The list of supported BWT formats includes `native`, `plain_default`, `plain_sorted`, `rfm`, `ropebwt`, `sdsl`, and `sga`. [See the wiki](https://github.com/jltsiren/bwt-merge/wiki/BWT-Formats) for further information.

## Citation

Jouni Sirén: **Burrows-Wheeler transform for terabases**.
Proc. DCC 2016, IEEE, pp. 211-220, Snowbird, Utah, USA, March 29 - April 1, 2016.
[DOI: 10.1109/DCC.2016.17](https://doi.org/10.1109/DCC.2016.17)

## Other references

Wing-Kai Hon, Tak-Wah Lam, Kunihiko Sadakane, Wing-Kin Sung, Siu-Ming Yiu:
**A space and time efficient algorithm for constructing compressed suffix arrays**.
Algorithmica 48(1): 23-36, 2007.
[DOI: 10.1007/s00453-006-1228-8](https://doi.org/10.1007/s00453-006-1228-8)

Jouni Sirén: **Compressed Suffix Arrays for Massive Data**.
Proc. SPIRE 2009, Springer LNCS 5721, pp. 63-74, Saariselkä, Finland, August 25-27, 2009.
[DOI: 10.1007/978-3-642-03784-9_7](https://doi.org/10.1007/978-3-642-03784-9_7)

Jouni Sirén: **Compressed Full-Text Indexes for Highly Repetitive Collections**.
PhD Thesis, University of Helsinki, June 2012.
[http://jltsiren.kapsi.fi/phd](https://jltsiren.kapsi.fi/phd)

Markus J. Bauer, Anthony J. Cox, and Giovanna Rosone:
**Lightweight algorithms for constructing and inverting the BWT of string collections**.
Theoretical Computer Science 483: 134-148, 2013.
[DOI: 10.1016/j.tcs.2012.02.002](https://doi.org/10.1016/j.tcs.2012.02.002)
