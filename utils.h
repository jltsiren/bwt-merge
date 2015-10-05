/*
  Copyright (c) 2015 Genome Research Ltd.

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#ifndef _BWTMERGE_UTILS_H
#define _BWTMERGE_UTILS_H

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>

#include <sdsl/wavelet_trees.hpp>

#include <omp.h>

namespace bwtmerge
{

//------------------------------------------------------------------------------

typedef std::uint64_t size_type;
typedef std::uint8_t  char_type;
typedef std::uint8_t  comp_type;
typedef std::uint8_t  byte_type;

const size_type WORD_BITS = 64;

const size_type BYTE_BITS     = 8;
const size_type KILOBYTE      = 1024;
const size_type MILLION       = 1000000;
const size_type MEGABYTE      = KILOBYTE * KILOBYTE;
const size_type GIGABYTE      = KILOBYTE * MEGABYTE;

const double BYTE_BITS_DOUBLE = 8.0;
const double KILOBYTE_DOUBLE  = 1024.0;
const double MILLION_DOUBLE   = 1000000.0;
const double MEGABYTE_DOUBLE  = KILOBYTE_DOUBLE * KILOBYTE_DOUBLE;
const double GIGABYTE_DOUBLE  = KILOBYTE_DOUBLE * MEGABYTE_DOUBLE;

//------------------------------------------------------------------------------

/*
  range_type stores a closed range [first, second]. Empty ranges are indicated by
  first > second. The emptiness check uses +1 to handle a common special case
  [0, -1].
*/

typedef std::pair<size_type, size_type> range_type;

struct Range
{
  inline static size_type length(range_type range)
  {
    return range.second + 1 - range.first;
  }

  inline static bool empty(range_type range)
  {
    return (range.first + 1 > range.second + 1);
  }
};

template<class A, class B>
std::ostream& operator<<(std::ostream& stream, const std::pair<A, B>& data)
{
  return stream << "(" << data.first << ", " << data.second << ")";
}

//------------------------------------------------------------------------------

/*
  RunBuffer transforms a sequence of values or runs into a sequence of maximal runs.
  Usage:

    RunBuffer buffer;
    while(...)
    {
      if(buffer.add(...)) { doSomething(buffer.run); }
    }
    buffer.flush(); doSomething(buffer.run);
*/

struct RunBuffer
{
  RunBuffer() : value(0), length(0), run(0, 0) { }

  inline bool add(size_type v, size_type n = 1)
  {
    if(v == value) { length += n; return false; }
    else
    {
      this->flush();
      this->value = v; this->length = n;
      return (this->run.second > 0);
    }
  }

  inline bool add(range_type run) { return this->add(run.first, run.second); }

  inline void flush() { this->run.first = this->value; this->run.second = this->length; }

  size_type  value, length;
  range_type run;
};

//------------------------------------------------------------------------------

template<class IntegerType>
inline size_type
bit_length(IntegerType val)
{
  return sdsl::bits::hi(val) + 1;
}

//------------------------------------------------------------------------------

const size_type FNV_OFFSET_BASIS = 0xcbf29ce484222325UL;
const size_type FNV_PRIME        = 0x100000001b3UL;

inline size_type fnv1a_hash(byte_type b, size_type seed)
{
  return (seed ^ b) * FNV_PRIME;
}

inline size_type fnv1a_hash(size_type val, size_type seed)
{
  byte_type* chars = (byte_type*)&val;
  for(size_type i = 0; i < sizeof(val); i++) { seed = fnv1a_hash(chars[i], seed); }
  return seed;
}

template<class ByteArray>
size_type fnv1a_hash(ByteArray& array)
{
  size_type res = FNV_OFFSET_BASIS;
  for(size_type i = 0; i < array.size(); i++) { res = fnv1a_hash((byte_type)(array[i]), res); }
  return res;
}

//------------------------------------------------------------------------------

inline double
inMegabytes(size_type bytes)
{
  return bytes / MEGABYTE_DOUBLE;
}

inline double
inGigabytes(size_type bytes)
{
  return bytes / GIGABYTE_DOUBLE;
}

inline double
inBPC(size_type bytes, size_type size)
{
  return (BYTE_BITS_DOUBLE * bytes) / size;
}

inline double
inMicroseconds(double seconds)
{
  return seconds * MILLION_DOUBLE;
}

const size_type DEFAULT_INDENT = 18;

void printHeader(const std::string& header, size_type indent = DEFAULT_INDENT);
void printSize(const std::string& header, size_type bytes, size_type data_size, size_type indent = DEFAULT_INDENT);
void printTime(const std::string& header, size_type found, size_type matches, size_type bytes, double seconds, size_type indent = DEFAULT_INDENT);
void printTime(const std::string& header, size_type queries, double seconds, size_type indent = DEFAULT_INDENT);

//------------------------------------------------------------------------------

double readTimer();
size_type memoryUsage(); // Peak memory usage in bytes.

//------------------------------------------------------------------------------

// Returns the total length of the rows, excluding line ends.
size_type readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows);

std::string tempFile(const std::string& name_part);

size_type fileSize(std::ifstream& file);
size_type fileSize(std::ofstream& file);

//------------------------------------------------------------------------------

template<class Iterator, class Comparator>
void
sequentialSort(Iterator first, Iterator last, const Comparator& comp)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, comp, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
sequentialSort(Iterator first, Iterator last)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last);
#endif
}

//------------------------------------------------------------------------------

/*
  Split the range approximately evenly between the blocks. The actual number of blocks
  will not be greater than the length of the range or smaller than 1.
*/

std::vector<range_type> getBounds(range_type range, size_type blocks);

//------------------------------------------------------------------------------

/*
  BWT-merge uses a contiguous byte alphabet [0, sigma - 1] internally. Array C is based on the
  number of occurrences of each character in the BWT.
*/

template<class AlphabetType>
inline bool
hasChar(const AlphabetType& alpha, comp_type comp)
{
  return (alpha.C[comp + 1] > alpha.C[comp]);
}

template<class AlphabetType>
inline range_type
charRange(const AlphabetType& alpha, comp_type comp)
{
  return range_type(alpha.C[comp], alpha.C[comp + 1] - 1);
}

template<class AlphabetType>
comp_type
findChar(const AlphabetType& alpha, size_type bwt_pos)
{
  comp_type comp = 0;
  while(alpha.C[comp + 1] <= bwt_pos) { comp++; }
  return comp;
}

// Returns (LF(i), BWT[i]).
template<class BWTType, class AlphabetType>
inline range_type
LF(const BWTType& bwt, const AlphabetType& alpha, size_type i)
{
  auto temp = bwt.inverse_select(i);
  return range_type(temp.first + alpha.C[temp.second], temp.second);
}

template<class BWTType, class AlphabetType>
inline size_type
LF(const BWTType& bwt, const AlphabetType& alpha, size_type i, comp_type comp)
{
  return alpha.C[comp] + bwt.rank(i, comp);
}

template<class BWTType, class AlphabetType>
inline range_type
LF(const BWTType& bwt, const AlphabetType& alpha, range_type range, comp_type comp)
{
  return range_type(LF(bwt, alpha, range.first, comp), LF(bwt, alpha, range.second + 1, comp) - 1);
}

template<class BWTType, class AlphabetType>
inline size_type
Psi(const BWTType& bwt, const AlphabetType& alpha, size_type i)
{
  comp_type comp = findChar(alpha, i);
  return bwt.select(i + 1 - alpha.C[comp], comp);
}

//------------------------------------------------------------------------------

/*
  Some SDSL extensions.
*/

/*
  Utility methods for writing int_vector_buffer<k> compatible files, for k = 8, 16, 32, 64.
*/
template<class Element>
struct IntVectorBuffer
{
  typedef uint64_t code_type;

  static void writeHeader(std::ofstream& out, size_type elements)
  {
    size_type bits = elements * BYTE_BITS * sizeof(Element);
    sdsl::write_member(bits, out);
  }

  static size_type readHeader(std::ifstream& in)
  {
    size_type bits = 0;
    sdsl::read_member(bits, in);
    return bits / (BYTE_BITS * sizeof(Element));
  }

  static void writeData(std::ofstream& out, Element* data, size_type elements)
  {
    size_type bytes = elements * sizeof(Element);
    if(bytes % sizeof(code_type) != 0) { bytes += sizeof(code_type) - bytes % sizeof(code_type); }
    char* ptr = (char*)data;
    for(size_type i = elements * sizeof(Element); i < bytes; i++) { ptr[i] = 0; }
    out.write(ptr, bytes);
  }

  static void readData(std::ifstream& in, Element* data, size_type elements)
  {
    size_type bytes = elements * sizeof(Element);
    if(bytes % sizeof(code_type) != 0) { bytes += sizeof(code_type) - bytes % sizeof(code_type); }
    in.read((char*)data, bytes);
  }
};

/*
  Generic in-memory construction from int_vector_buffer<8> and size. Not very space-efficient, as it
  duplicates the data.
*/
template<class Type>
void
directConstruct(Type& structure, const sdsl::int_vector<8>& data)
{
  std::string ramfile = sdsl::ram_file_name(sdsl::util::to_string(&structure));
  sdsl::store_to_file(data, ramfile);
  {
    sdsl::int_vector_buffer<8> buffer(ramfile); // Must remove the buffer before removing the ramfile.
    Type temp(buffer, data.size());
    structure.swap(temp);
  }
  sdsl::ram_fs::remove(ramfile);
}

/*
  Builds sd_vector directly from a strictly increasing sequence without storing the sequence
  explicitly. No sanity checks. Uses an ugly hack to access sd_vector internals.
*/
struct SDVectorBuilder
{
  SDVectorBuilder();
  SDVectorBuilder(size_type _size, size_type _capacity);

  inline size_type size() const { return this->bits; }
  inline size_type capacity() const { return this->onebits; }

  inline void add(size_type value)
  {
    size_type high = value >> this->low_width;
    this->high_pos += high - prev_high;
    this->prev_high = high;
    this->low[this->tail] = value; this->tail++;
    this->high[this->high_pos] = 1; this->high_pos++;
  }

  size_type bits, onebits;
  size_type low_width;
  size_type tail, prev_high, high_pos;

  sdsl::int_vector<0> low;
  sdsl::bit_vector    high;
};

//------------------------------------------------------------------------------

} // namespace bwtmerge

namespace sdsl
{

template<>
template<>
sd_vector<>::sd_vector<bwtmerge::SDVectorBuilder*>(
  bwtmerge::SDVectorBuilder* builder,
  bwtmerge::SDVectorBuilder*);

} // namespace sdsl

//------------------------------------------------------------------------------

#endif // _BWTMERGE_UTILS_H
