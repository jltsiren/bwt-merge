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

#include <cstring>

#include "support.h"

namespace bwtmerge
{

//------------------------------------------------------------------------------

/*
  The default alphabet interprets \0 and $ as endmarkers, ACGT and acgt as ACGT,
  and the and the remaining characters as N.
*/

const sdsl::int_vector<8> Alphabet::DEFAULT_CHAR2COMP =
{
  0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

  5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

const sdsl::int_vector<8> Alphabet::DEFAULT_COMP2CHAR = { '$', 'A', 'C', 'G', 'T', 'N' };

//------------------------------------------------------------------------------

Alphabet::Alphabet() :
  char2comp(DEFAULT_CHAR2COMP), comp2char(DEFAULT_COMP2CHAR),
  C(sdsl::int_vector<64>(DEFAULT_COMP2CHAR.size() + 1, 0)),
  sigma(DEFAULT_COMP2CHAR.size())
{
}

Alphabet::Alphabet(const Alphabet& source)
{
  this->copy(source);
}

Alphabet::Alphabet(Alphabet&& source)
{
  *this = std::move(source);
}

Alphabet::Alphabet(const sdsl::int_vector<64>& counts,
  const sdsl::int_vector<8>& _char2comp, const sdsl::int_vector<8>& _comp2char) :
  char2comp(_char2comp), comp2char(_comp2char),
  C(sdsl::int_vector<64>(_comp2char.size() + 1, 0)),
  sigma(_comp2char.size())
{
  for(size_type i = 0; i < counts.size(); i++) { this->C[i + 1] = this->C[i] + counts[i]; }
}

Alphabet::Alphabet(size_type _sigma)
{
  if(_sigma == 0 || _sigma > MAX_SIGMA)
  {
    std::cerr << "Alphabet::Alphabet(): Invalid alphabet size: " << _sigma << std::endl;
    return;
  }

  this->sigma = _sigma;
  this->char2comp = sdsl::int_vector<8>(MAX_SIGMA, 0);
  this->comp2char = sdsl::int_vector<8>(this->sigma, 0);
  this->C = sdsl::int_vector<64>(this->sigma + 1, 0);
  for(size_type c = 0; c < this->sigma; c++)
  {
    this->char2comp[c] = this->comp2char[c] = c;
  }
  for(size_type c = this->sigma; c < MAX_SIGMA; c++)
  {
    this->char2comp[c] = 0;
  }
}

Alphabet::~Alphabet()
{
}

void
Alphabet::copy(const Alphabet& source)
{
  this->char2comp = source.char2comp;
  this->comp2char = source.comp2char;
  this->C = source.C;
  this->sigma = source.sigma;
}

void
Alphabet::swap(Alphabet& source)
{
  if(this != &source)
  {
    this->char2comp.swap(source.char2comp);
    this->comp2char.swap(source.comp2char);
    this->C.swap(source.C);
    std::swap(this->sigma, source.sigma);
  }
}

Alphabet&
Alphabet::operator=(const Alphabet& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

Alphabet&
Alphabet::operator=(Alphabet&& source)
{
  if(this != &source)
  {
    this->char2comp = std::move(source.char2comp);
    this->comp2char = std::move(source.comp2char);
    this->C = std::move(source.C);
    this->sigma = source.sigma;
  }
  return *this;
}

Alphabet::size_type
Alphabet::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += this->char2comp.serialize(out, child, "char2comp");
  written_bytes += this->comp2char.serialize(out, child, "comp2char");
  written_bytes += this->C.serialize(out, child, "C");
  written_bytes += sdsl::write_member(this->sigma, out, child, "sigma");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
Alphabet::load(std::istream& in)
{
  this->char2comp.load(in);
  this->comp2char.load(in);
  this->C.load(in);
  sdsl::read_member(this->sigma, in);
}

bool
Alphabet::sorted() const
{
  for(size_type i = 1; i < this->sigma; i++)
  {
    if(this->comp2char[i - 1] >= this->comp2char[i]) { return false; }
  }
  return true;
}

bool
Alphabet::operator== (const Alphabet& another) const
{
  if(this->sigma != another.sigma) { return false; }
  for(size_type i = 0; i < this->char2comp.size(); i++)
  {
    if(this->char2comp[i] != another.char2comp[i]) { return false; }
  }
  for(size_type i = 0; i < this->comp2char.size(); i++)
  {
    if(this->comp2char[i] != another.comp2char[i]) { return false; }
  }
  return true;
}

bool
Alphabet::operator!= (const Alphabet& another) const
{
  return !(*this == another);
}

std::ostream&
operator<<(std::ostream& stream, const Alphabet& alpha)
{
  stream << "[ ";
  for(size_type i = 0; i < alpha.sigma; i++)
  {
    if(i > 0) { stream << ", "; }
    stream << "(" << (char)(alpha.comp2char[i]) << ", "
           << charRange(alpha, i) << ")";
  }
  stream << " ]";
  return stream;
}

//------------------------------------------------------------------------------

BlockArray::BlockArray()
{
  this->bytes = 0;
}

BlockArray::BlockArray(const BlockArray& source)
{
  this->copy(source);
}

BlockArray::BlockArray(BlockArray&& source)
{
  *this = std::move(source);
}

BlockArray::~BlockArray()
{
  this->clear();
}

void
BlockArray::copy(const BlockArray& source)
{
  this->clear();
  this->data = std::vector<value_type*>(source.data.size());
  this->bytes = source.bytes;

  for(size_type i = 0; i < this->data.size(); i++)
  {
    if(source.data[i] == 0) { this->data[i] = 0; }
    else
    {
      this->data[i] = new value_type[BLOCK_SIZE];
      std::memcpy((void*)(this->data[i]), (void*)(source.data[i]), BLOCK_SIZE);
    }
  }
}

void
BlockArray::clear()
{
  for(size_type i = 0; i < this->data.size(); i++)
  {
    this->clear(i);
  }
  this->data.clear();
  this->bytes = 0;
}

void
BlockArray::swap(BlockArray& source)
{
  if(this != &source)
  {
    this->data.swap(source.data);
    std::swap(this->bytes, source.bytes);
  }
}

BlockArray&
BlockArray::operator=(const BlockArray& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

BlockArray&
BlockArray::operator=(BlockArray&& source)
{
  if(this != &source)
  {
    this->clear();
    std::swap(this->data, source.data); // The source must not delete the data.
    this->bytes = std::move(source.bytes);
  }
  return *this;
}

BlockArray::size_type
BlockArray::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += sdsl::write_member(this->bytes, out, child, "bytes");
  for(size_type i = 0; i < this->data.size(); i++)
  {
    out.write((char*)(this->data[i]), BLOCK_SIZE);
    written_bytes += BLOCK_SIZE;
  }
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
BlockArray::load(std::istream& in)
{
  this->clear();
  sdsl::read_member(this->bytes, in);
  this->data = std::vector<value_type*>((this->bytes + BLOCK_SIZE - 1) / BLOCK_SIZE, 0);
  for(size_type i = 0; i < this->data.size(); i++)
  {
    this->data[i] = new byte_type[BLOCK_SIZE];
    in.read((char*)(this->data[i]), BLOCK_SIZE);
  }
}

//------------------------------------------------------------------------------

CumulativeArray::CumulativeArray()
{
  this->m_size = 0;
}

CumulativeArray::CumulativeArray(const CumulativeArray& source)
{
  this->copy(source);
}

CumulativeArray::CumulativeArray(CumulativeArray&& source)
{
  *this = std::move(source);
}

CumulativeArray::~CumulativeArray()
{
}

CumulativeArray::CumulativeArray(SDVectorBuilder& builder)
{
  this->m_size = builder.capacity();
  this->data = sdsl::sd_vector<>(&builder, &builder);
  this->buildSupport();
}

void
CumulativeArray::copy(const CumulativeArray& source)
{
  this->data = source.data;
  this->rank = source.rank;
  this->select_1 = source.select_1;
  this->select_0 = source.select_0;
  this->setVectors();
  this->m_size = source.m_size;
}

void
CumulativeArray::buildSupport()
{
  sdsl::util::init_support(this->rank, &(this->data));
  sdsl::util::init_support(this->select_1, &(this->data));
  sdsl::util::init_support(this->select_0, &(this->data));
}

void
CumulativeArray::setVectors()
{
  this->rank.set_vector(&(this->data));
  this->select_1.set_vector(&(this->data));
  this->select_0.set_vector(&(this->data));
}

void
CumulativeArray::swap(CumulativeArray& source)
{
  if(this != &source)
  {
    this->data.swap(source.data);
    sdsl::util::swap_support(this->rank, source.rank, &(this->data), &(source.data));
    sdsl::util::swap_support(this->select_1, source.select_1, &(this->data), &(source.data));
    sdsl::util::swap_support(this->select_0, source.select_0, &(this->data), &(source.data));
    std::swap(this->m_size, source.m_size);
  }
}

CumulativeArray&
CumulativeArray::operator=(const CumulativeArray& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

CumulativeArray&
CumulativeArray::operator=(CumulativeArray&& source)
{
  if(this != &source)
  {
    this->data = std::move(source.data);
    this->rank = std::move(source.rank);
    this->select_1 = std::move(source.select_1);
    this->select_0 = std::move(source.select_0);
    this->setVectors();
    this->m_size = source.m_size;
  }
  return *this;
}

CumulativeArray::size_type
CumulativeArray::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += this->data.serialize(out, child, "data");
  written_bytes += this->rank.serialize(out, child, "rank");
  written_bytes += this->select_1.serialize(out, child, "select_1");
  written_bytes += this->select_0.serialize(out, child, "select_0");
  written_bytes += sdsl::write_member(this->m_size, out, child, "size");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
CumulativeArray::load(std::istream& in)
{
  this->data.load(in);
  this->rank.load(in, &(this->data));
  this->select_1.load(in, &(this->data));
  this->select_0.load(in, &(this->data));
  sdsl::read_member(this->m_size, in);
}

//------------------------------------------------------------------------------

void open(RLArray<sdsl::int_vector_buffer<8>>& array, const std::string filename,
  size_type runs, size_type values)
{
  array.data = sdsl::int_vector_buffer<8>(filename);
  array.run_count = runs;
  array.value_count = values;
}

template<>
void
RLArray<BlockArray>::clear()
{
  this->data.clear();
  this->run_count = 0;
  this->value_count = 0;
}

template<>
void
RLArray<sdsl::int_vector_buffer<8>>::clear()
{
  this->data.close();
  this->run_count = 0;
  this->value_count = 0;
}

template<>
void
RLArray<BlockArray>::write(const std::string& filename)
{
  std::ofstream out(filename.c_str(), std::ios_base::binary);
  if(!out)
  {
    std::cerr << "RLArray::write(): Cannot open output file " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }

  IntVectorBuffer<BlockArray::value_type>::writeHeader(out, this->bytes());
  for(size_type block = 0; block < this->data.blocks(); block++)
  {
    size_type bytes = std::min(BlockArray::BLOCK_SIZE, this->bytes() - block * BlockArray::BLOCK_SIZE);
    IntVectorBuffer<BlockArray::value_type>::writeData(out, this->data.data[block], bytes);
    this->data.clear(block);
  }

  out.close();
}

template<>
void
RLIterator<BlockArray>::read()
{
  if(this->end()) { this->run.first = ~(value_type)0; this->run.second = ~(length_type)0; return; }
  this->run.first += ByteCode::read(this->array->data, this->ptr);
  this->run.second = ByteCode::read(this->array->data, this->ptr);
  this->array->data.clearUntil(this->ptr);
}

//------------------------------------------------------------------------------

RankArray::RankArray()
{
}

RankArray::~RankArray()
{
  this->close();
  for(size_type i = 0; i < this->filenames.size(); i++) { remove(this->filenames[i].c_str()); }
}

void
RankArray::open()
{
  this->close();
  this->inputs = std::vector<array_type>(this->size());
  this->iterators = std::vector<iterator>(this->size());

  for(size_type i = 0; i < this->size(); i++)
  {
    bwtmerge::open(this->inputs[i], this->filenames[i], this->run_counts[i], this->value_counts[i]);
    this->iterators[i] = iterator(this->inputs[i]);
  }

  this->heapify();
}

void
RankArray::close()
{
  this->iterators.clear();
  for(size_type i = 0; i < this->inputs.size(); i++) { this->inputs[i].clear(); }
  this->inputs.clear();
}

void
RankArray::heapify()
{
  if(this->size() <= 1) { return; }

  size_type i = parent(this->size() - 1);
  while(true)
  {
    this->down(i);
    if(i == 0) { break; }
    i--;
  }
}

//------------------------------------------------------------------------------

} // namespace bwtmerge
