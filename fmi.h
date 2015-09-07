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

namespace bwtmerge
{

//------------------------------------------------------------------------------

class FMI
{
public:
  typedef BWT::size_type size_type;

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

  template<class Format>
  void load(const std::string& filename)
  {
    this->bwt.load<Format>(filename);
    Alphabet temp = createAlphabet(Format::order());
    sdsl::int_vector<64> counts;
    this->bwt.characterCounts(counts);
    this->alpha = Alphabet(counts, temp.char2comp, temp.comp2char);
  }

//------------------------------------------------------------------------------

  inline size_type size() const { return this->bwt.size(); }
  inline size_type sequences() const { return this->alpha.C[1]; }

  template<class Iterator>
  range_type find(Iterator begin, Iterator end) const
  {
    if(begin == end) { return range_type(0, this->size() - 1); }

    --end;
    range_type range = charRange(this->alpha, this->alpha.char2comp[*end]);
    while(!Range::empty(range) && end != begin)
    {
      --end;
      range = LF(this->bwt, this->alpha, range, this->alpha.char2comp[*end]);
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

} // namespace bwtmerge

#endif // _BWTMERGE_FMI_H
