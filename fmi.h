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

#ifndef _BWTMERGE_FMI_H
#define _BWTMERGE_FMI_H

#include <fstream>
#include <iostream>

#include "bwt.h"
#include "formats.h"

namespace bwtmerge
{

//------------------------------------------------------------------------------

class FMI;

void serialize(const FMI& fmi, const std::string& filename, const std::string& format);
void load(FMI& fmi, const std::string& filename, const std::string& format);

//------------------------------------------------------------------------------

struct MergeParameters
{
  const static size_type RUN_BUFFER_SIZE = 8 * MEGABYTE;      // Runs.
  const static size_type THREAD_BUFFER_SIZE = 512 * MEGABYTE; // Bytes.
  const static size_type MERGE_BUFFERS = 5;
  const static size_type BLOCKS_PER_THREAD = 4;

  size_type run_buffer_size, thread_buffer_size;
  size_type merge_buffers;
  size_type threads, sequence_blocks;

  MergeParameters();
};

std::ostream& operator<< (std::ostream& stream, const MergeParameters& parameters);

//------------------------------------------------------------------------------

class FMI
{
public:
  typedef BWT::size_type size_type;

  const static size_type SHORT_RANGE = 3; // Compute LF(i) individually for ranges up to this length.

  FMI();
  FMI(const FMI& source);
  FMI(FMI&& source);
  ~FMI();

  void swap(FMI& source);
  FMI& operator=(const FMI& source);
  FMI& operator=(FMI&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& i);

//------------------------------------------------------------------------------

  /*
    This constructor merges a and b, destroying them in the process.
  */
  FMI(FMI& a, FMI& b, MergeParameters parameters = MergeParameters());

//------------------------------------------------------------------------------

  template<class Format>
  void serialize(const std::string& filename) const
  {
    if(!compatible(this->alpha, Format::order()))
    {
      std::cerr << "FMI::serialize(): Warning: " << Format::name
                << " is not compatible with " << alphabetName(identifyAlphabet(this->alpha))
                << " alphabets!" << std::endl;
    }
    this->bwt.serialize<Format>(filename);
  }

  template<class Format>
  void load(const std::string& filename)
  {
    sdsl::int_vector<64> counts;
    this->bwt.load<Format>(filename, counts);
    Alphabet temp = createAlphabet(Format::order());
    this->alpha = Alphabet(counts, temp.char2comp, temp.comp2char);
  }

//------------------------------------------------------------------------------

  inline size_type size() const { return this->bwt.size(); }
  inline size_type sequences() const { return this->alpha.C[1]; }

  inline range_type charRange(comp_type comp) const
  {
    return bwtmerge::charRange(this->alpha, comp);
  }

  // Returns (LF(i), BWT[i]).
  inline range_type LF(size_type i) const
  {
    return bwtmerge::LF(this->bwt, this->alpha, i);
  }

  inline size_type LF(size_type i, comp_type comp) const
  {
    return bwtmerge::LF(this->bwt, this->alpha, i, comp);
  }

  inline range_type LF(range_type range, comp_type comp) const
  {
    return bwtmerge::LF(this->bwt, this->alpha, range, comp);
  }

  /*
    Computes LF(i) for comp values 1 to sigma - 1.
  */
  inline void LF(size_type i, BWT::ranks_type& results) const
  {
    this->bwt.ranks(i, results);
    for(size_type c = 1; c < this->alpha.sigma; c++) { results[c] += this->alpha.C[c]; }
  }

  /*
    Computes LF(range) for comp values 1 to sigma - 1.
  */
  inline void LF(range_type range, BWT::ranks_type& sp, BWT::ranks_type& ep) const
  {
    this->bwt.ranks(range.first, sp); this->bwt.ranks(range.second + 1, ep);
    for(size_type c = 1; c < this->alpha.sigma; c++)
    {
      sp[c] += this->alpha.C[c]; ep[c] += this->alpha.C[c] - 1;
    }
  }

  template<class Iterator>
  range_type find(Iterator begin, Iterator end) const
  {
    if(begin == end) { return range_type(0, this->size() - 1); }

    --end;
    range_type range = this->charRange(this->alpha.char2comp[*end]);
    while(!Range::empty(range) && end != begin)
    {
      --end;
      range = this->LF(range, this->alpha.char2comp[*end]);
    }

    return range;
  }

  template<class Container>
  range_type find(const Container& pattern) const
  {
    return this->find(pattern.begin(), pattern.end());
  }

  template<class Element>
  range_type find(const Element* pattern, size_type length) const
  {
    return this->find(pattern, pattern + length);
  }

//------------------------------------------------------------------------------

  BWT      bwt;
  Alphabet alpha;

private:
  void copy(const FMI& source);
};

//------------------------------------------------------------------------------

template<>
void
FMI::serialize<NativeFormat>(const std::string& filename) const;

template<>
void
FMI::load<NativeFormat>(const std::string& filename);

//------------------------------------------------------------------------------

} // namespace bwtmerge

#endif // _BWTMERGE_FMI_H
