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

  std::cout << "Memory usage before merging: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;

  FMI merged;
  merge(merged, fmi1, fmi2, patterns, chars);
  sdsl::store_to_file(merged, output_name);

  std::cout << "Final memory usage: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
testFMI(FMI& fmi, const std::string& name, const std::vector<std::string>& patterns, size_type chars)
{
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

void
loadFMI(FMI& fmi, const std::string& filename, const std::string& name,
  const std::vector<std::string>& patterns, size_type chars)
{
  sdsl::load_from_file(fmi, filename);
  testFMI(fmi, name, patterns, chars);
}

void
merge(FMI& result, FMI& index, FMI& increment,
  const std::vector<std::string>& patterns, size_type chars)
{
  double increment_mb = inMegabytes(increment.size());

  double start = readTimer();
  result = FMI(index, increment);
  double seconds = readTimer() - start;
  std::cout << "BWTs merged in " << seconds << " seconds ("
            << (increment_mb / seconds) << " MB/s)" << std::endl;
  std::cout << std::endl;

  testFMI(result, "Merged", patterns, chars);
}

//------------------------------------------------------------------------------
