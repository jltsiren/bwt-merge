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

#include <stack>

#include "fmi.h"

namespace bwtmerge
{

//------------------------------------------------------------------------------

FMI::FMI()
{
}

FMI::FMI(const FMI& source)
{
  this->copy(source);
}

FMI::FMI(FMI&& source)
{
  *this = std::move(source);
}

FMI::~FMI()
{
}

void
FMI::copy(const FMI& source)
{
  this->bwt = source.bwt;
  this->alpha = source.alpha;
}

void
FMI::swap(FMI& source)
{
  if(this != &source)
  {
    this->bwt.swap(source.bwt);
    this->alpha.swap(source.alpha);
  }
}

FMI&
FMI::operator=(const FMI& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

FMI&
FMI::operator=(FMI&& source)
{
  if(this != &source)
  {
    this->bwt = std::move(source.bwt);
    this->alpha = std::move(source.alpha);
  }
  return *this;
}

FMI::size_type
FMI::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += this->bwt.serialize(out, child, "bwt");
  written_bytes += this->alpha.serialize(out, child, "alpha");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
FMI::load(std::istream& in)
{
  this->bwt.load(in);
  this->alpha.load(in);
}

//------------------------------------------------------------------------------

struct MergePosition
{
  size_type  a_pos;
  range_type b_range;

  MergePosition() : a_pos(0), b_range(0, 0) {}
  MergePosition(size_type a, size_type b) : a_pos(a), b_range(b, b) {}
  MergePosition(size_type pos, range_type range) : a_pos(pos), b_range(range) {}
};

const static size_type RUN_BUFFER_SIZE = MEGABYTE;

FMI::FMI(FMI& a, FMI& b)
{
  if(a.alpha != b.alpha)
  {
    std::cerr << "FMI::FMI(): Cannot merge BWTs with different alphabets." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  RLArray ra;
  std::vector<RLArray::run_type> buffer;
  std::stack<MergePosition> positions;

  positions.push(MergePosition(a.sequences(), range_type(0, b.sequences() - 1)));
  while(!(positions.empty()))
  {
    MergePosition curr = positions.top(); positions.pop();
    buffer.push_back(RLArray::run_type(curr.a_pos, Range::length(curr.b_range)));
    if(buffer.size() >= RUN_BUFFER_SIZE) { add(ra, buffer); buffer.clear(); }

    if(Range::length(curr.b_range) < b.alpha.sigma)
    {
      for(size_type i = curr.b_range.first; i <= curr.b_range.second; i++)
      {
        range_type pred = b.LF(i);
        if(pred.second != 0)
        {
          positions.push(MergePosition(a.LF(curr.a_pos, pred.second), pred.first));
        }
      }
    }
    else
    {
      for(comp_type c = 1; c < b.alpha.sigma; c++)
      {
        range_type prev = b.LF(curr.b_range, c);
        if(!(Range::empty(prev))) { positions.push(MergePosition(a.LF(curr.a_pos, c), prev)); }
      }
    }
  }
  if(buffer.size() > 0) { add(ra, buffer); buffer.clear(); }

  std::cout << "Memory usage with RA: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;

  this->bwt = BWT(a.bwt, b.bwt, ra);
  this->alpha = a.alpha;
  for(size_type c = 0; c <= this->alpha.sigma; c++) { this->alpha.C[c] += b.alpha.C[c]; }
}

//------------------------------------------------------------------------------

} // namespace bwtmerge
