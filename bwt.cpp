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

BWT::BWT(BWT& a, BWT& b, RankArray& ra)
{
  a.destroy(); b.destroy();
  ra.open();

  size_type a_rle_pos = 0, b_rle_pos = 0;
  size_type a_seq_pos = 0;
  range_type a_run = Run::read(a.data, a_rle_pos); a.data.clearUntil(a_rle_pos);
  range_type b_run = Run::read(b.data, b_rle_pos); b.data.clearUntil(b_rle_pos);

  RunBuffer buffer;
  while(!(ra.iterators.empty()))
  {
    RankArray::iterator iter = ra.iterators.top(); ra.iterators.pop();

    while(a_seq_pos < iter->first)
    {
      size_type length = std::min(iter->first - a_seq_pos, a_run.second);
      if(buffer.add(a_run.first, length)) { Run::write(this->data, buffer.run); }
      a_run.second -= length; a_seq_pos += length;
      if(a_run.second == 0 && a_rle_pos < a.data.size())
      {
        a_run = Run::read(a.data, a_rle_pos); a.data.clearUntil(a_rle_pos);
      }
    }
    while(iter->second > 0)
    {
      size_type length = std::min(iter->second, b_run.second);
      if(buffer.add(b_run.first, length)) { Run::write(this->data, buffer.run); }
      b_run.second -= length; iter->second -= length;
      if(b_run.second == 0 && b_rle_pos < b.data.size())
      {
        b_run = Run::read(b.data, b_rle_pos); b.data.clearUntil(b_rle_pos);
      }
    }

    ++iter;
    if(!(iter.end())) { ra.iterators.push(iter); }
  }
  while(a_run.second > 0)
  {
    if(buffer.add(a_run)) { Run::write(this->data, buffer.run); }
    if(a_rle_pos < a.data.size()) { a_run = Run::read(a.data, a_rle_pos); a.data.clearUntil(a_rle_pos); }
    else { a_run.second = 0; }
  }
  buffer.flush(); Run::write(this->data, buffer.run);
  ra.close();

  this->build();
}

//------------------------------------------------------------------------------

void
BWT::build()
{
  std::vector<size_type> block_ends;
  size_type seq_pos = 0, rle_pos = 0;
  while(rle_pos < this->bytes())
  {
    range_type run = Run::read(this->data, rle_pos); seq_pos += run.second;
    if(rle_pos >= this->bytes() || rle_pos % SAMPLE_RATE == 0) { block_ends.push_back(seq_pos - 1); }
  }

  size_type blocks = block_ends.size();
  {
    this->block_boundaries = sdsl::sd_vector<>(block_ends.begin(), block_ends.end());
    sdsl::util::clear(block_ends);
    sdsl::util::init_support(this->block_rank, &(this->block_boundaries));
    sdsl::util::init_support(this->block_select, &(this->block_boundaries));
  }

  sdsl::int_vector<0> counts[SIGMA];
  for(size_type c = 0; c < SIGMA; c++)
  {
    counts[c] = sdsl::int_vector<0>(blocks, 0, bit_length(this->size()));
  }
  for(size_type block = 0; block < blocks; block++)
  {
    size_type limit = std::min(this->bytes(), (block + 1) * SAMPLE_RATE);
    rle_pos = block * SAMPLE_RATE;
    while(rle_pos < limit)
    {
      range_type run = Run::read(this->data, rle_pos);
      counts[run.first][block] += run.second;
    }
  }
  for(size_type c = 0; c < SIGMA; c++)
  {
    this->samples[c] = CumulativeArray(counts[c]);
    sdsl::util::clear(counts[c]);
  }
}

void
BWT::destroy()
{
  for(size_type c = 0; c < SIGMA; c++) { sdsl::util::clear(this->samples[c]); }
  sdsl::util::clear(this->block_boundaries);
  sdsl::util::clear(this->block_rank);
  sdsl::util::clear(this->block_select);
}

//------------------------------------------------------------------------------

void
BWT::characterCounts(sdsl::int_vector<64>& counts)
{
  counts = sdsl::int_vector<64>(SIGMA, 0);

  size_type rle_pos = 0, seq_pos = 0;
  while(rle_pos < this->bytes())
  {
    range_type run = Run::read(this->data, rle_pos);
    counts[run.first] += run.second; seq_pos += run.second;
  }
}

size_type
BWT::hash() const
{
  size_type res = FNV_OFFSET_BASIS;
  size_type rle_pos = 0;
  while(rle_pos < this->bytes())
  {
    range_type run = Run::read(this->data, rle_pos);
    for(size_type i = 0; i < run.second; i++) { res = fnv1a_hash((byte_type)(run.first), res); }
  }
  return res;
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

  for(size_type i = 0; i < this->size(); i++)
  {
    bwtmerge::open(this->inputs[i], this->filenames[i], this->run_counts[i], this->value_counts[i]);
    this->iterators.push(iterator(this->inputs[i]));
  }
}

void
RankArray::close()
{
  this->iterators = std::priority_queue<iterator, std::vector<iterator>, IterComparator>();
  for(size_type i = 0; i < this->inputs.size(); i++) { this->inputs[i].clear(); }
  this->inputs.clear();
}

//------------------------------------------------------------------------------

} // namespace bwtmerge
