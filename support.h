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

#ifndef _BWTMERGE_SUPPORT_H
#define _BWTMERGE_SUPPORT_H

#include "utils.h"

namespace bwtmerge
{

//------------------------------------------------------------------------------

template<class ByteVector>
void characterCounts(const ByteVector& sequence, const sdsl::int_vector<8>& char2comp, sdsl::int_vector<64>& counts);

/*
  This replaces the SDSL byte_alphabet. The main improvements are:
    - The alphabet can be built from an existing sequence.
    - The comp order does not need to be the same as character order, as long as \0 is the first character.
*/

class Alphabet
{
public:
  typedef bwtmerge::size_type size_type;
  const static size_type MAX_SIGMA = 256;

  const static sdsl::int_vector<8> DEFAULT_CHAR2COMP;
  const static sdsl::int_vector<8> DEFAULT_COMP2CHAR;

  Alphabet();
  Alphabet(const Alphabet& source);
  Alphabet(Alphabet&& source);
  ~Alphabet();

  /*
    ByteVector only has to support operator[] and size(). If there is a clearly faster way for
    sequential access, function characterCounts() should be specialized.
  */
  template<class ByteVector>
  explicit Alphabet(const ByteVector& sequence,
    const sdsl::int_vector<8>& _char2comp = DEFAULT_CHAR2COMP,
    const sdsl::int_vector<8>& _comp2char = DEFAULT_COMP2CHAR) :
    char2comp(_char2comp), comp2char(_comp2char),
    C(sdsl::int_vector<64>(_comp2char.size() + 1, 0)),
    sigma(_comp2char.size())
  {
    if(sequence.size() == 0) { return; }

    characterCounts(sequence, this->char2comp, this->C);
    for(size_type i = 0, sum = 0; i < this->C.size(); i++)
    {
      size_type temp = this->C[i]; this->C[i] = sum; sum += temp;
    }
  }

  /*
    The counts array holds character counts for all comp values.
  */
  explicit Alphabet(const sdsl::int_vector<64>& counts,
    const sdsl::int_vector<8>& _char2comp = DEFAULT_CHAR2COMP,
    const sdsl::int_vector<8>& _comp2char = DEFAULT_COMP2CHAR);

  void swap(Alphabet& source);
  Alphabet& operator=(const Alphabet& source);
  Alphabet& operator=(Alphabet&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  sdsl::int_vector<8>  char2comp, comp2char;
  sdsl::int_vector<64> C;
  size_type            sigma;

private:
  void copy(const Alphabet& a);
};  // class Alphabet

template<class ByteVector>
void
characterCounts(const ByteVector& sequence, const sdsl::int_vector<8>& char2comp, sdsl::int_vector<64>& counts)
{
  for(size_type c = 0; c < counts.size(); c++) { counts[c] = 0; }
  for(size_type i = 0; i < sequence.size(); i++) { counts[char2comp[sequence[i]]]++; }
}

//------------------------------------------------------------------------------

class BlockArray
{
public:
  typedef bwtmerge::size_type size_type;
  const static size_type BLOCK_SIZE = MEGABYTE;

  BlockArray();
  BlockArray(const BlockArray& source);
  BlockArray(BlockArray&& source);
  ~BlockArray();

  void swap(BlockArray& source);
  BlockArray& operator=(const BlockArray& source);
  BlockArray& operator=(BlockArray&& source);

  inline size_type size() const { return this->bytes; }
  inline size_type blocks() const { return this->data.size(); }

  inline static size_type block(size_type i) { return i / BLOCK_SIZE; }
  inline static size_type offset(size_type i) { return i % BLOCK_SIZE; }

  void clear();
  void clear(size_type _block);

  inline byte_type operator[] (size_type i) const
  {
    return this->data[block(i)][offset(i)];
  }

  inline byte_type& operator[] (size_type i)
  {
    return this->data[block(i)][offset(i)];
  }

  inline void push_back(byte_type value)
  {
    if(offset(this->bytes) == 0) { this->data.push_back(new byte_type[BLOCK_SIZE]); }
    (*this)[this->bytes] = value;
    this->bytes++;
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  std::vector<byte_type*> data;
  size_type               bytes;

private:
  void copy(const BlockArray& b);
};  // class BlockArray

//------------------------------------------------------------------------------

struct Run
{
  const static size_type BLOCK_SIZE = 64; // No run can continue past a block boundary.
  const static size_type SIGMA      = 6;
  const static size_type MAX_RUN    = 256 / SIGMA;  // 42; encoded as 6 * 41

  const static size_type DATA_BITS  = 7;
  const static byte_type DATA_MASK  = 0x7F; // Bits 0-6 contain the data.
  const static byte_type NEXT_BYTE  = 0x80; // Bit 7 tells whether the extension continues.

  /*
    Returns (comp value, run length) and updates i to point past the run.
  */
  template<class ByteArray>
  static range_type read(const ByteArray& array, size_type& i)
  {
    range_type run(array[i] % SIGMA, array[i] / SIGMA + 1); i++;
    if(run.second >= MAX_RUN)
    {
      size_type offset = 0;
      run.second += array[i] & DATA_MASK;
      while(array[i] & NEXT_BYTE)
      {
        i++; offset += DATA_BITS;
        run.second += ((size_type)(array[i] & DATA_MASK)) << offset;
      }
      i++;
    }
    return run;
  }

  inline static byte_type basicRun(comp_type comp, size_type length)
  {
    return comp + SIGMA * (length - 1);
  }

  template<class ByteArray>
  static void write(ByteArray& array, comp_type comp, size_type length);

  /*
    Write an extension to the run using 7 bits/byte, least significant byte first.
    The high-order bit tells whether the extension continues to the next byte.
  */
  template<class ByteArray>
  static void writeExtension(ByteArray& array, size_type length);
};

template<class ByteArray>
void
Run::write(ByteArray& array, comp_type comp, size_type length)
{
  while(length > 0)
  {
    size_type bytes_remaining = BLOCK_SIZE - (array.size() % BLOCK_SIZE) - 1;
    size_type basic_length = std::min(length, (bytes_remaining > 0 ? MAX_RUN : MAX_RUN - 1));
    array.push_back(basicRun(comp, basic_length)); length -= basic_length;
    if(length == 0) { return; }

    if(bytes_remaining > 0)
    {
      size_type extension_length = length;
      if(bit_length(length) > DATA_BITS * bytes_remaining)
      {
        extension_length = sdsl::bits::lo_set[DATA_BITS * bytes_remaining];
      }
      writeExtension(array, extension_length); length -= extension_length;
    }
  }
}

template<class ByteArray>
void
Run::writeExtension(ByteArray& array, size_type length)
{
  while(length > DATA_MASK)
  {
    array.push_back((length & DATA_MASK) | NEXT_BYTE);
    length >>= DATA_BITS;
   }
  array.push_back(length);
}

//------------------------------------------------------------------------------

/*
  This class uses an sd_vector to encode the cumulative sum of an array of integers.
  The array contains sum() items in size() elements. The array uses 0-based indexes.
  Each element is encoded as 'items' 0-bits, followed by an 1-bit.
*/
class CumulativeArray
{
public:
  typedef bwtmerge::size_type size_type;

  CumulativeArray();
  CumulativeArray(const CumulativeArray& source);
  CumulativeArray(CumulativeArray&& source);
  ~CumulativeArray();

  /*
    The IntVector has to support operator[] that returns a non-const reference.
    The input is the original array, which is temporarily modified during the construction.
  */
  template<class IntVector>
  explicit CumulativeArray(IntVector& sequence)
  {
    this->m_size = sequence.size();

    for(size_type i = 1; i < this->size(); i++) { sequence[i] += sequence[i - 1] + 1; }
    this->data = sdsl::sd_vector<>(sequence.begin(), sequence.end());
    for(size_type i = this->size() - 1; i > 0; i--) { sequence[i] -= sequence[i - 1] + 1; }

    sdsl::util::init_support(this->rank, &(this->data));
    sdsl::util::init_support(this->select_1, &(this->data));
    sdsl::util::init_support(this->select_0, &(this->data));
  }

  void swap(CumulativeArray& source);
  CumulativeArray& operator=(const CumulativeArray& source);
  CumulativeArray& operator=(CumulativeArray&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  // The number of elements.
  inline size_type size() const { return this->m_size; }

  // The sum of all elements.
  inline size_type sum() const { return this->data.size() - this->size(); }

  // The sum of the first k elements.
  inline size_type sum(size_type k) const
  {
    if(k == 0) { return 0; }
    if(k > this->size()) { k = this->size(); }
    return this->select_1(k) - k + 1;
  }

  inline size_type operator[](size_type i) const { return this->sum(i + 1) - this->sum(i); }

  // The inverse of sum(). Returns the element for item i.
  inline size_type inverse(size_type i) const
  {
    if(i >= this->sum()) { return this->size(); }
    return this->select_0(i + 1) - i;
  }

  // Is item i the last item in its element.
  inline bool isLast(size_type i) const
  {
    if(i >= this->sum()) { return false; }
    return this->data[this->select_0(i + 1) + 1];
  }

  // A combination of the above.
  inline size_type inverse(size_type i, bool& is_last) const
  {
    if(i >= this->sum()) { is_last = false; return this->size(); }
    uint64_t temp = this->select_0(i + 1);
    is_last = this->data[temp + 1];
    return temp - i;
  }

  sdsl::sd_vector<>                data;
  sdsl::sd_vector<>::rank_1_type   rank;
  sdsl::sd_vector<>::select_1_type select_1;
  sdsl::sd_vector<>::select_0_type select_0;
  size_type                  m_size;  // Number of elements.

private:
  void copy(const CumulativeArray& source);
  void setVectors();
};  // class CumulativeArray

//------------------------------------------------------------------------------

} // namespace bwtmerge

#endif // _BWTMERGE_SUPPORT_H
