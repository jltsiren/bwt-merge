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

#ifndef _BWTMERGE_SEQUENCE_H
#define _BWTMERGE_SEQUENCE_H

#include "support.h"

namespace bwtmerge
{

//------------------------------------------------------------------------------

struct RankArray;

/*
  A basic run-length encoded sequence for alphabet 0-5.
*/
class BWT
{
public:
  typedef BlockArray::size_type size_type;
  typedef bwtmerge::char_type   char_type;
  typedef Run::comp_type        comp_type;
  typedef Run::length_type      length_type;

  const static size_type SAMPLE_RATE = Run::BLOCK_SIZE;
  const static size_type SIGMA       = Run::SIGMA;

  BWT();
  BWT(const BWT& source);
  BWT(BWT&& source);
  ~BWT();

  void swap(BWT& source);
  BWT& operator=(const BWT& source);
  BWT& operator=(BWT&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& i);

//------------------------------------------------------------------------------

  /*
    This constructor interleaves the source BWTs according to the rank array. All the
    input structures will be destroyed in the process.
  */
  BWT(BWT& a, BWT&b, RankArray& ra);

//------------------------------------------------------------------------------

  template<class Format>
  void serialize(const std::string& filename)
  {
    std::ofstream out(filename.c_str(), std::ios_base::binary);
    if(!out)
    {
      std::cerr << "BWT::serialize(): Cannot open output file " << filename << std::endl;
      return;
    }
    Format::write(out, this->data);
    out.close();
  }

  template<class Format>
  void load(const std::string& filename, sdsl::int_vector<64>& counts)
  {
    std::ifstream in(filename.c_str(), std::ios_base::binary);
    if(!in)
    {
      std::cerr << "BWT::load(): Cannot open input file " << filename << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Format::read(in, this->data, counts);
    in.close();
    this->build();
  }

//------------------------------------------------------------------------------

  inline size_type size() const { return this->block_boundaries.size(); }
  inline size_type bytes() const { return this->data.size(); }
  inline size_type count(comp_type c) const { return this->samples[c].sum(); }

  inline size_type rank(size_type i, comp_type c) const
  {
    if(c >= SIGMA) { return 0; }
    if(i > this->size()) { i = this->size(); }

    size_type block = this->block_rank(i);
    size_type res = this->samples[c].sum(block);
    size_type rle_pos = block * SAMPLE_RATE;
    size_type seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);

    while(seq_pos < i)
    {
      range_type run = Run::read(this->data, rle_pos);
      seq_pos += run.second;  // The starting position of the next run.
      if(run.first == c)
      {
        res += run.second;  // Number of c's before the next run.
        if(seq_pos > i) { res -= seq_pos - i; }
      }
    }

    return res;
  }

  inline size_type select(size_type i, comp_type c) const
  {
    if(c >= SIGMA) { return 0; }
    if(i == 0) { return 0; }
    if(i > this->count(c)) { return this->size(); }

    size_type block = this->samples[c].inverse(i - 1);
    size_type count = this->samples[c].sum(block);
    size_type rle_pos = block * SAMPLE_RATE;
    size_type seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
    while(true)
    {
      range_type run = Run::read(this->data, rle_pos);
      seq_pos += run.second - 1;  // The last position in the run.
      if(run.first == c)
      {
        count += run.second;  // Number of c's up to the end of the run.
        if(count >= i) { return seq_pos + i - count; }
      }
      seq_pos++;  // Move to the first position in the next run.
    }
  }

  inline comp_type operator[](size_type i) const
  {
    if(i >= this->size()) { return 0; }

    size_type block = this->block_rank(i);
    size_type rle_pos = block * SAMPLE_RATE;
    size_type seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
    while(true)
    {
      range_type run = Run::read(this->data, rle_pos);
      seq_pos += run.second;  // The start of the next run.
      if(seq_pos > i) { return run.first; }
    }
  }

  // returns (rank(i, seq[i]), seq[i])
  inline range_type inverse_select(size_type i) const
  {
    range_type run(0, 0);
    if(i >= this->size()) { return run; }

    size_type block = this->block_rank(i);
    size_type rle_pos = block * SAMPLE_RATE;
    size_type seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);

    size_type ranks[SIGMA] = {};
    while(seq_pos <= i)
    {
      run = Run::read(this->data, rle_pos);
      seq_pos += run.second;  // The starting position of the next run.
      ranks[run.first] += run.second; // Number of c's before the next run.
    }

    return range_type(this->samples[run.first].sum(block) + ranks[run.first] - (seq_pos - i), run.first);
  }

//------------------------------------------------------------------------------

  template<class ByteVector>
  void extract(range_type range, ByteVector& buffer) const
  {
    if(Range::empty(range) || range.second >= this->size()) { return; }
    buffer.resize(Range::length(range));

    // Find the first character.
    size_type block = this->block_rank(range.first);
    size_type rle_pos = block * SAMPLE_RATE;
    size_type seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
    range_type run(0, 0);

    while(true)
    {
      run = Run::read(this->data, rle_pos);
      seq_pos += run.second - 1;  // The last position in the run.
      if(seq_pos >= range.first) { break; }
      seq_pos++;  // Move to the first position in the next run.
    }

    // Fill the buffer.
    for(size_type i = range.first; i <= range.second; i++)
    {
      if(i > seq_pos)
      {
        run = Run::read(this->data, rle_pos);
        seq_pos += run.second;
      }
      buffer[i - range.first] = run.first;
    }
  }

  void characterCounts(sdsl::int_vector<64>& counts);

  size_type hash() const;

//------------------------------------------------------------------------------

  BlockArray                       data;
  CumulativeArray                  samples[SIGMA];

  sdsl::sd_vector<>                block_boundaries; // Marks the last sequence position in each block.
  sdsl::sd_vector<>::rank_1_type   block_rank;
  sdsl::sd_vector<>::select_1_type block_select;

private:
  void copy(const BWT& source);
  void setVectors();

  // Builds/destroys the rank/select structures.
  void build();
  void destroy();
};  // class BWT

//------------------------------------------------------------------------------

} // namespace bwtmerge

#endif // _BWTMERGE_SEQUENCE_H
