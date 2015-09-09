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
#include "formats.h"

using namespace bwtmerge;

//------------------------------------------------------------------------------

const size_type RUN_BUFFER_SIZE = MEGABYTE;

//------------------------------------------------------------------------------

void loadFMI(FMI& fmi, const std::string& filename, const std::string& name,
  const std::vector<std::string>& patterns, size_type chars);

void merge(FMI& result, FMI& index, FMI& increment,
  const std::vector<std::string>& patterns, size_type chars);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 4)
  {
    std::cerr << "Usage: bwt_merge input1 input2 output [patterns]" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }
  bool verify = (argc > 4);

  std::cout << "BWT-merge" << std::endl;
  std::cout << std::endl;

  std::string index_name = argv[1];
  std::string increment_name = argv[2];
  std::string output_name = argv[3];
  std::string pattern_name = (verify ? argv[4] : "");

  std::cout << "Input 1: " << index_name << std::endl;
  std::cout << "Input 2: " << increment_name << std::endl;
  std::cout << "Output: " << output_name << std::endl;
  if(verify) { std::cout << "Patterns: " << pattern_name << std::endl; }
  std::cout << std::endl;

  std::vector<std::string> patterns;
  size_type chars = 0;
  if(verify)
  {
    chars = readRows(pattern_name, patterns, true);
    std::cout << "Read " << patterns.size() << " patterns of total length " << chars << std::endl;
    std::cout << std::endl;
  }

  FMI fmi1, fmi2;
  loadFMI(fmi1, index_name, "BWT 1", patterns, chars);
  loadFMI(fmi2, increment_name, "BWT 2", patterns, chars);

  FMI merged;
  merge(merged, fmi1, fmi2, patterns, chars);

  std::cout << "Memory usage: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
loadFMI(FMI& fmi, const std::string& filename, const std::string& name,
  const std::vector<std::string>& patterns, size_type chars)
{
  sdsl::load_from_file(fmi, filename);
  printSize(name, sdsl::size_in_bytes(fmi), fmi.size());

  if(chars > 0)
  {
    double start = readTimer();
    size_type found = 0, matches = 0;
    for(auto pattern : patterns)
    {
      range_type range = fmi.find(pattern);
      if(!(Range::empty(range))) { found++; matches += Range::length(range); }
    }
    double seconds = readTimer() - start;
    printTime(name, found, matches, chars, seconds);
  }

  std::cout << std::endl;
}

//------------------------------------------------------------------------------

struct MergePosition
{
  size_type  index_pos;
  range_type increment_range;

  MergePosition() : index_pos(0), increment_range(0, 0) {}
  MergePosition(size_type pos, size_type inc_pos) : index_pos(pos), increment_range(inc_pos, inc_pos) {}
  MergePosition(size_type pos, range_type range) : index_pos(pos), increment_range(range) {}
};

void
merge(FMI& result, FMI& index, FMI& increment,
  const std::vector<std::string>& patterns, size_type chars)
{
  double start = readTimer();

  RLArray ra;
  std::vector<RLArray::run_type> buffer;
  std::stack<MergePosition> positions;

  positions.push(MergePosition(index.sequences(), range_type(0, increment.sequences() - 1)));
  while(!(positions.empty()))
  {
    MergePosition curr = positions.top(); positions.pop();
    buffer.push_back(RLArray::run_type(curr.index_pos, Range::length(curr.increment_range)));
    if(buffer.size() >= RUN_BUFFER_SIZE) { add(ra, buffer); buffer.clear(); }

    if(Range::length(curr.increment_range) < increment.alpha.sigma)
    {
      for(size_type i = curr.increment_range.first; i <= curr.increment_range.second; i++)
      {
        range_type pred = increment.LF(i);
        if(pred.second != 0)
        {
          positions.push(MergePosition(index.LF(curr.index_pos, pred.second), pred.first));
        }
      }
    }
    else
    {
      for(comp_type c = 1; c < increment.alpha.sigma; c++)
      {
        range_type prev = increment.LF(curr.increment_range, c);
        if(!(Range::empty(prev))) { positions.push(MergePosition(index.LF(curr.index_pos, c), prev)); }
      }
    }
  }
  if(buffer.size() > 0) { add(ra, buffer); buffer.clear(); }

  double seconds = readTimer() - start;
  std::cout << "Rank array built in " << seconds << " seconds" << std::endl;
  std::cout << std::endl;

  printHeader("RA"); std::cout << ra.values() << " values in " << ra.size() << " runs" << std::endl;
  printSize("RA", sdsl::size_in_bytes(ra), increment.size());
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
