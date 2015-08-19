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

#include "bwt.h"

namespace bwtmerge
{

//------------------------------------------------------------------------------

BWT::BWT()
{
}

BWT::BWT(const BWT& source)
{
  this->copy(source);
}

BWT::BWT(BWT&& source)
{
  *this = std::move(source);
}

BWT::~BWT()
{
}

void
BWT::copy(const BWT& source)
{
  this->data = source.data;
  for(size_type c = 0; c < SIGMA; c++) { this->samples[c] = source.samples[c]; }

  this->block_boundaries = source.block_boundaries;
  this->block_rank = source.block_rank;
  this->block_select = source.block_select;
  this->setVectors();
}

void
BWT::setVectors()
{
  this->block_rank.set_vector(&(this->block_boundaries));
  this->block_select.set_vector(&(this->block_boundaries));
}

void
BWT::swap(BWT& source)
{
  if(this != &source)
  {
    this->data.swap(source.data);
    for(size_type c = 0; c < SIGMA; c++) { this->samples[c].swap(source.samples[c]); }
    this->block_boundaries.swap(source.block_boundaries);
    sdsl::util::swap_support(this->block_rank, source.block_rank, &(this->block_boundaries), &(source.block_boundaries));
    sdsl::util::swap_support(this->block_select, source.block_select, &(this->block_boundaries), &(source.block_boundaries));
  }
}

BWT&
BWT::operator=(const BWT& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

BWT&
BWT::operator=(BWT&& source)
{
  if(this != &source)
  {
    this->data = std::move(source.data);
    for(size_type c = 0; c < SIGMA; c++) { this->samples[c] = std::move(source.samples[c]); }

    this->block_boundaries = std::move(source.block_boundaries);
    this->block_rank = std::move(source.block_rank);
    this->block_select = std::move(source.block_select);
    this->setVectors();
  }
  return *this;
}

BWT::size_type
BWT::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += this->data.serialize(out, child, "data");
  for(size_type c = 0; c < SIGMA; c++)
  {
    std::stringstream ss; ss << "samples_" << c;
    written_bytes += this->samples[c].serialize(out, child, ss.str());
  }
  written_bytes += this->block_boundaries.serialize(out, child, "block_boundaries");
  written_bytes += this->block_rank.serialize(out, child, "block_rank");
  written_bytes += this->block_select.serialize(out, child, "block_select");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
BWT::load(std::istream& in)
{
  this->data.load(in);
  for(size_type c = 0; c < SIGMA; c++) { this->samples[c].load(in); }

  this->block_boundaries.load(in);
  this->block_rank.load(in, &(this->block_boundaries));
  this->block_select.load(in, &(this->block_boundaries));
}

//------------------------------------------------------------------------------

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
