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

#include <sstream>

#include "sequence.h"

namespace bwtmerge
{

//------------------------------------------------------------------------------

BWT::BWT()
{
}

BWT::BWT(const BWT& b)
{
  this->copy(b);
}

BWT::BWT(BWT&& b)
{
  *this = std::move(b);
}

BWT::~BWT()
{
}

void
BWT::copy(const BWT& b)
{
  this->data = b.data;
  for(size_type c = 0; c < SIGMA; c++) { this->samples[c] = b.samples[c]; }

  this->block_boundaries = b.block_boundaries;
  this->block_rank = b.block_rank;
  this->block_select = b.block_select;
  this->setVectors();
}

void
BWT::setVectors()
{
  this->block_rank.set_vector(&(this->block_boundaries));
  this->block_select.set_vector(&(this->block_boundaries));
}

void
BWT::swap(BWT& b)
{
  if(this != &b)
  {
    this->data.swap(s.data);
    for(size_type c = 0; c < SIGMA; c++) { this->samples[c].swap(b.samples[c]); }
    this->block_boundaries.swap(b.block_boundaries);
    util::swap_support(this->block_rank, b.block_rank, &(this->block_boundaries), &(b.block_boundaries));
    util::swap_support(this->block_select, b.block_select, &(this->block_boundaries), &(b.block_boundaries));
  }
}

BWT&
BWT::operator=(const BWT& b)
{
  if(this != &b) { this->copy(b); }
  return *this;
}

BWT&
BWT::operator=(BWT&& b)
{
  if(this != &b)
  {
    this->data = std::move(b.data);
    for(size_type c = 0; c < SIGMA; c++) { this->samples[c] = std::move(b.samples[c]); }

    this->block_boundaries = std::move(b.block_boundaries);
    this->block_rank = std::move(b.block_rank);
    this->block_select = std::move(b.block_select);
    this->setVectors();
  }
  return *this;
}

BWT::size_type
BWT::serialize(std::ostream& out, structure_tree_node* s, std::string name) const
{
  structure_tree_node* child = structure_tree::add_child(s, name, util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += write_vector(this->data, out, child, "data");
  for(size_type c = 0; c < SIGMA; c++)
  {
    std::stringstream ss; ss << "samples_" << c;
    written_bytes += this->samples[c].serialize(out, child, ss.str());
  }
  written_bytes += this->block_boundaries.serialize(out, child, "block_boundaries");
  written_bytes += this->block_rank.serialize(out, child, "block_rank");
  written_bytes += this->block_select.serialize(out, child, "block_select");
  structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
BWT::load(std::istream& in)
{
  read_vector(this->data, in);
  for(size_type c = 0; c < SIGMA; c++) { this->samples[c].load(in); }

  this->block_boundaries.load(in);
  this->block_rank.load(in, &(this->block_boundaries));
  this->block_select.load(in, &(this->block_boundaries));
}

//------------------------------------------------------------------------------

template<>
void
characterCounts(const BWT& sequence, sdsl::int_vector<64>& counts)
{
  for(uint64_t c = 0; c < counts.size(); c++) { counts[c] = 0; }

  uint64_t rle_pos = 0, seq_pos = 0;
  while(rle_pos < sequence.bytes())
  {
    range_type run = Run::read(sequence.data, rle_pos);
    counts[run.first] += run.second; seq_pos += run.second;
  }
}

//------------------------------------------------------------------------------

} // namespace bwtmerge
