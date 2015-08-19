# BWT-merge

Building the Burrows-Wheeler transform (BWT) is a solved problem for sequence collections of up to hundreds of gigabytes in size. For larger collections, there are several issues to consider:

* **Construction time:** As a rough guideline, an algorithm indexing 1 MB/s is good for up to 100 gigabytes of data, and somewhat useful until 1 terabyte. Larger datasets require faster algorithms.
* **Construction space:** When the datasets are much larger than the amount of memory available, we cannot afford using even a single bit of working space per input character.
* **Available hardware:** In a typical computer cluster, a single node has two CPUs (with up to tens of cores), tens to hundreds of gigabytes of local memory, a limited amount of local disk space, and large amounts of shared (but often slow) disk space. Some tools require fast GPUs or large amounts of fast disk space, which are generally not available.
* **Resource usage:** Merging large BWTs is easy by doing a lot of redundant work on multiple nodes. Because computer clusters generally do not have large amounts of unused capacity, good tools should make an efficient use of resources.

The purpose of this tool is to merge run-length encoded BWTs of read collections on a single node. To speed up the computation, the tool will use multiple threads and take advantage of the run-length encoding when building the merging structures. The structures themselves will be run-length encoded, hopefully requiring similar space as the merged BWT.

## Existing algorithms

There are several BWT construction algorithms based on updating an existing BWT. Some algorithms **extend** sequences already in the collection, updating *BWT(T[i+1,j])* to *BWT(T[i,j])*. Others **insert** new sequences to the collection, updating *BWT(S)* to *BWT(S,T)*. In both cases, the algorithm can either do **batch updates** to a **static** BWT representation, or use a **dynamic** representation for the BWT. Algorithms with a static BWT representation require more working space, while algorithms with a dynamic representation have more space overhead for the BWT.

Some common implementations include:

* [BEETL](https://github.com/BEETL/BEETL): (extend, static) on disk
* [NVBIO](http://nvlabs.github.io/nvbio/): (insert, dynamic) using GPU
* [RLCSA](http://jltsiren.kapsi.fi/rlcsa): (insert, static)
* [ropebwt](https://github.com/lh3/ropebwt): (extend, static) or (insert by extending, dynamic)
* [ropebwt2](https://github.com/lh3/ropebwt2): (insert by extending, dynamic)

As this tool will merge existing BWTs, it will be (insert, static).

There are also other algorithms for building the BWT for large read collections They are based on partitioning the suffixes, sorting each of the partitions separately, and building the BWT directly.

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
