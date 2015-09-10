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
    The counts array holds character counts for all comp values.
  */
  explicit Alphabet(const sdsl::int_vector<64>& counts,
    const sdsl::int_vector<8>& _char2comp = DEFAULT_CHAR2COMP,
    const sdsl::int_vector<8>& _comp2char = DEFAULT_COMP2CHAR);

  /*
    Creates an alphabet of given size, where char values are also comp values.
  */
  explicit Alphabet(size_type _sigma);

  void swap(Alphabet& source);
  Alphabet& operator=(const Alphabet& source);
  Alphabet& operator=(Alphabet&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  bool operator== (const Alphabet& another) const;
  bool operator!= (const Alphabet& another) const;

  sdsl::int_vector<8>  char2comp, comp2char;
  sdsl::int_vector<64> C;
  size_type            sigma;

private:
  void copy(const Alphabet& a);
};  // class Alphabet

std::ostream& operator<<(std::ostream& stream, const Alphabet& alpha);

//------------------------------------------------------------------------------

class BlockArray
{
public:
  typedef bwtmerge::size_type size_type;
  typedef bwtmerge::byte_type value_type;
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

  void clear(size_type _block)
  {
    delete[] this->data[_block]; this->data[_block] = 0;
  }

  /*
    Removes the block before block(i).
  */
  void clearUntil(size_type i)
  {
    if(block(i) > 0 && this->data[block(i) - 1] != 0) { this->clear(block(i) - 1); }
  }

  inline value_type operator[] (size_type i) const
  {
    return this->data[block(i)][offset(i)];
  }

  inline value_type& operator[] (size_type i)
  {
    return this->data[block(i)][offset(i)];
  }

  inline void push_back(value_type value)
  {
    if(offset(this->bytes) == 0) { this->data.push_back(new value_type[BLOCK_SIZE]); }
    (*this)[this->bytes] = value;
    this->bytes++;
  }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  std::vector<value_type*> data;
  size_type                bytes;

private:
  void copy(const BlockArray& source);
};  // class BlockArray

//------------------------------------------------------------------------------

/*
  Encodes unsigned integers as byte sequences. Each byte contains 7 bits of data
  and one bit telling whether the encoding continues in the next byte. The data is
  stored in LSB order.
*/

struct ByteCode
{
  typedef bwtmerge::size_type    value_type;
  typedef BlockArray::value_type code_type;

  const static size_type DATA_BITS  = 7;
  const static code_type DATA_MASK  = 0x7F;
  const static code_type NEXT_BYTE  = 0x80;

  /*
    Reads the next value and updates i to point to the byte after the value.
  */
  template<class ByteArray>
  static value_type read(const ByteArray& array, size_type& i)
  {
    size_type offset = 0;
    value_type res = array[i] & DATA_MASK;
    while(array[i] & NEXT_BYTE)
    {
      i++; offset += DATA_BITS;
      res += ((value_type)(array[i] & DATA_MASK)) << offset;
    }
    i++;
    return res;
  }

  /*
    Encodes the value and stores it in the array using push_back().
  */
  template<class ByteArray>
  static void write(ByteArray& array, value_type value)
  {
    while(value > DATA_MASK)
    {
      array.push_back((value & DATA_MASK) | NEXT_BYTE);
      value >>= DATA_BITS;
    }
    array.push_back(value);
  }
};

//------------------------------------------------------------------------------

/*
  A run in BWT.
*/

struct Run
{
  typedef bwtmerge::comp_type    comp_type;
  typedef bwtmerge::size_type    length_type;
  typedef BlockArray::value_type code_type;

  const static size_type   BLOCK_SIZE = 64; // No run can continue past a block boundary.
  const static size_type   SIGMA      = 6;
  const static length_type MAX_RUN    = 256 / SIGMA;  // 42; encoded as 6 * 41

  inline static code_type encodeBasic(comp_type comp, length_type length)
  {
    return comp + SIGMA * (length - 1);
  }

  inline static range_type decodeBasic(code_type code)
  {
    return range_type(code % SIGMA, code / SIGMA + 1);
  }

  /*
    Returns (comp value, run length) and updates i to point past the run.
  */
  template<class ByteArray>
  static range_type read(const ByteArray& array, size_type& i)
  {
    range_type run = decodeBasic(array[i]); i++;
    if(run.second >= MAX_RUN) { run.second += ByteCode::read(array, i); }
    return run;
  }

  /*
    Encodes the run and stores it in the array using push_back(). If the encoding
    would continue past a block boundary, the run is split in two.
  */
  template<class ByteArray>
  static void write(ByteArray& array, comp_type comp, length_type length)
  {
    while(length > 0)
    {
      if(length < MAX_RUN)
      {
        array.push_back(encodeBasic(comp, length));
        return;
      }

      size_type bytes_remaining = BLOCK_SIZE - (array.size() % BLOCK_SIZE);
      length_type basic_length = (bytes_remaining > 1 ? MAX_RUN : MAX_RUN - 1);
      array.push_back(encodeBasic(comp, basic_length)); length -= basic_length;
      bytes_remaining--;

      if(bytes_remaining > 0)
      {
        length_type extension_length = length;
        if(bit_length(length) > ByteCode::DATA_BITS * bytes_remaining)
        {
          extension_length = sdsl::bits::lo_set[ByteCode::DATA_BITS * bytes_remaining];
        }
        ByteCode::write(array, extension_length); length -= extension_length;
      }
    }
  }

  template<class ByteArray>
  inline static void write(ByteArray& array, range_type run) { write(array, run.first, run.second); }
};

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
  typedef bwtmerge::size_type value_type;

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
  inline value_type sum() const { return this->data.size() - this->size(); }

  // The sum of the first k elements.
  inline value_type sum(size_type k) const
  {
    if(k == 0) { return 0; }
    if(k > this->size()) { k = this->size(); }
    return this->select_1(k) - k + 1;
  }

  inline value_type operator[](size_type i) const { return this->sum(i + 1) - this->sum(i); }

  // The inverse of sum(). Returns the element for item i.
  inline size_type inverse(value_type i) const
  {
    if(i >= this->sum()) { return this->size(); }
    return this->select_0(i + 1) - i;
  }

  // Is item i the last item in its element.
  inline bool isLast(value_type i) const
  {
    if(i >= this->sum()) { return false; }
    return this->data[this->select_0(i + 1) + 1];
  }

  // A combination of the above.
  inline size_type inverse(value_type i, bool& is_last) const
  {
    if(i >= this->sum()) { is_last = false; return this->size(); }
    size_type temp = this->select_0(i + 1);
    is_last = this->data[temp + 1];
    return temp - i;
  }

  sdsl::sd_vector<>                data;
  sdsl::sd_vector<>::rank_1_type   rank;
  sdsl::sd_vector<>::select_1_type select_1;
  sdsl::sd_vector<>::select_0_type select_0;
  size_type                        m_size;  // Number of elements.

private:
  void copy(const CumulativeArray& source);
  void setVectors();
};  // class CumulativeArray

//------------------------------------------------------------------------------

/*
  A run-length encoded non-decreasing integer array. The non-const iterator is
  destructive: it deletes the blocks it has progressed past.
*/

class RLIterator;
class RLConstIterator;

class RLArray
{
public:
  typedef BlockArray::size_type size_type;
  typedef bwtmerge::size_type   value_type;
  typedef bwtmerge::size_type   length_type;
  typedef std::pair<value_type, length_type> run_type;

  typedef RLIterator      iterator;
  typedef RLConstIterator const_iterator;

  RLArray();
  RLArray(const RLArray& source);
  RLArray(RLArray&& source);
  ~RLArray();

  /*
    Builds an RLArray from the source vector. The vector is sorted during construction.
  */
  explicit RLArray(std::vector<value_type>& source);
  explicit RLArray(std::vector<run_type>& source);

  /*
    Merges the input arrays and deletes their contents.
  */
  RLArray(RLArray& a, RLArray& b);

  void swap(RLArray& source);
  RLArray& operator=(const RLArray& source);
  RLArray& operator=(RLArray&& source);

  inline size_type size() const { return this->run_count; }
  inline size_type values() const { return this->value_count; }
  inline size_type bytes() const { return this->data.size(); }

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  void clear();

  BlockArray data;
  size_type  run_count, value_count;

private:
  void copy(const RLArray& source);

  inline void addRun(value_type value, value_type& prev, length_type length)
  {
    ByteCode::write(this->data, value - prev); prev = value;
    ByteCode::write(this->data, length);
    this->run_count++; this->value_count += length;
  }
};  // class RLArray

template<class Element>
void
add(RLArray& target, std::vector<Element>& values)
{
  if(values.empty()) { return; }

  RLArray temp(values);
  if(target.size() == 0) { target.swap(temp); }
  else { target = RLArray(target, temp); }
}

void add(RLArray& target, RLArray& source);

class RLIterator
{
public:
  typedef RLArray::size_type  size_type;
  typedef RLArray::value_type value_type;
  typedef RLArray::run_type   run_type;

  inline RLIterator(RLArray& _array) :
    array(_array), pos(0), ptr(0), prev(0), run(0, 0)
  {
    this->read();
  }

  inline run_type operator* () const { return this->run; }
  inline run_type* operator-> () { return &(this->run); }
  inline void operator++ () { this->pos++; this->read(); }
  inline bool end() const { return (this->pos >= this->array.size()); }

  RLArray&   array;
  size_type  pos, ptr;
  value_type prev;
  run_type   run;

private:
  inline void read()
  {
    if(this->end()) { return; }
    this->prev = this->run.first;
    this->run.first += ByteCode::read(this->array.data, this->ptr);
    this->run.second = ByteCode::read(this->array.data, this->ptr);
    this->array.data.clearUntil(this->ptr);
  }
};  // class RLIterator

class RLConstIterator
{
public:
  typedef RLArray::size_type  size_type;
  typedef RLArray::value_type value_type;
  typedef RLArray::run_type   run_type;

  inline RLConstIterator(const RLArray& _array) :
    array(_array), pos(0), ptr(0), prev(0), run(0, 0)
  {
    this->read();
  }

  inline run_type operator* () const { return this->run; }
  inline const run_type* operator-> () const { return &(this->run); }
  inline void operator++ () { this->pos++; this->read(); }
  inline bool end() const { return (this->pos >= this->array.size()); }

  const RLArray&  array;
  size_type       pos, ptr;
  value_type      prev;
  run_type        run;

private:
  inline void read()
  {
    if(this->end()) { return; }
    this->prev = this->run.first;
    this->run.first += ByteCode::read(this->array.data, this->ptr);
    this->run.second = ByteCode::read(this->array.data, this->ptr);
  }
};  // class RLConstIterator

//------------------------------------------------------------------------------

} // namespace bwtmerge

#endif // _BWTMERGE_SUPPORT_H
