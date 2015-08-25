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
Alphabet::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
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

Alphabet
rfmAlphabet()
{
  Alphabet alpha;
  std::swap(alpha.comp2char[4], alpha.comp2char[5]);
  std::swap(alpha.char2comp['N'], alpha.char2comp['T']);
  std::swap(alpha.char2comp['n'], alpha.char2comp['t']);
  return alpha;
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
BlockArray::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
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
    this->data[block(i)] = new byte_type[BLOCK_SIZE];
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
CumulativeArray::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
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

RLArray::RLArray()
{
  this->run_count = 0;
  this->value_count = 0;
}

RLArray::RLArray(const RLArray& source)
{
  this->copy(source);
}

RLArray::RLArray(RLArray&& source)
{
  *this = std::move(source);
}

RLArray::~RLArray()
{
}

RLArray::RLArray(std::vector<RLArray::value_type>& source)
{
  this->run_count = 0; this->value_count = 0;
  if(source.empty()) { return; }

  sequentialSort(source.begin(), source.end());
  value_type prev = 0, curr = source[0];
  length_type length = 1;
  for(size_type i = 1; i < source.size(); i++)
  {
    if(source[i] == curr) { length++; }
    else
    {
      this->addRun(curr, prev, length);
      curr = source[i]; length = 1;
    }
  }
  this->addRun(curr, prev, length);
}

RLArray::RLArray(std::vector<RLArray::run_type>& source)
{
  this->run_count = 0; this->value_count = 0;
  if(source.empty()) { return; }

  sequentialSort(source.begin(), source.end());
  value_type prev = 0;
  for(size_type i = 0; i < source.size(); i++) { this->addRun(source[i].first, prev, source[i].second); }
}

RLArray::RLArray(RLArray& a, RLArray& b)
{
  this->run_count = 0; this->value_count = 0;

  iterator a_iter(a), b_iter(b);
  value_type prev = 0;
  while(!(a_iter.end()) && !(b_iter.end()))
  {
    if(a_iter->first < b_iter->first)
    {
      this->addRun(a_iter->first, prev, a_iter->second); ++a_iter;
    }
    else if(b_iter->first < a_iter->first)
    {
      this->addRun(b_iter->first, prev, b_iter->second); ++b_iter;
    }
    else
    {
      this->addRun(a_iter->first, prev, a_iter->second + b_iter->second); ++a_iter; ++b_iter;
    }
  }
  while(!(a_iter.end()))
  {
    this->addRun(a_iter->first, prev, a_iter->second); ++a_iter;
  }
  while(!(b_iter.end()))
  {
    this->addRun(b_iter->first, prev, b_iter->second); ++b_iter;
  }
}

void
RLArray::copy(const RLArray& source)
{
  this->data = source.data;
  this->run_count = source.run_count;
  this->value_count = source.value_count;
}

void
RLArray::swap(RLArray& source)
{
  if(this != &source)
  {
    this->data.swap(source.data);
    std::swap(this->run_count, source.run_count);
    std::swap(this->value_count, source.value_count);
  }
}

RLArray&
RLArray::operator=(const RLArray& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

RLArray&
RLArray::operator=(RLArray&& source)
{
  if(this != &source)
  {
    this->data = std::move(source.data);
    this->run_count = std::move(source.run_count);
    this->value_count = std::move(source.value_count);
  }
  return *this;
}

RLArray::size_type
RLArray::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += this->data.serialize(out, child, "data");
  written_bytes += sdsl::write_member(this->run_count, out, child, "run_count");
  written_bytes += sdsl::write_member(this->value_count, out, child, "value_count");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
RLArray::load(std::istream& in)
{
  this->data.load(in);
  sdsl::read_member(this->run_count, in);
  sdsl::read_member(this->value_count, in);
}

//------------------------------------------------------------------------------

} // namespace bwtmerge
