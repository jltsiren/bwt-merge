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

using namespace bwtmerge;

//------------------------------------------------------------------------------

const size_type RUN_BUFFER_SIZE = MEGABYTE;

//------------------------------------------------------------------------------

void verifyFMI(FMI& fmi, const std::string& name,
  const std::vector<std::string>& patterns, std::vector<size_type>& results);

void merge(FMI& result, FMI& index, FMI& increment);

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
  if(verify)
  {
    size_type chars = readRows(pattern_name, patterns, true);
    std::cout << "Read " << patterns.size() << " patterns of total length " << chars << std::endl;
    std::cout << std::endl;
  }

  FMI index;
  std::vector<size_type> index_results;
  index.load<NativeFormat>(index_name);
  verifyFMI(index, "BWT 1", patterns, index_results);

  FMI increment;
  std::vector<size_type> increment_results;
  increment.load<NativeFormat>(increment_name);
  verifyFMI(increment, "BWT 2", patterns, increment_results);

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "bwt_merge: Memory usage before merging: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
#endif

  FMI merged;
  std::vector<size_type> merged_results;
  merge(merged, index, increment);
  merged.serialize<NativeFormat>(output_name);
  verifyFMI(merged, "Merged", patterns, merged_results);

  if(verify)
  {
    size_type errors = 0;
    for(size_type i = 0; i < patterns.size(); i++)
    {
      if(merged_results[i] != index_results[i] + increment_results[i]) { errors++; }
    }
    if(errors > 0)
    {
      std::cout << "Verification failed for " << errors << " patterns." << std::endl;
    }
    else
    {
      std::cout << "Verification successful." << std::endl;
    }
    std::cout << std::endl;
  }

  std::cout << "Peak memory usage: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
verifyFMI(FMI& fmi, const std::string& name,
  const std::vector<std::string>& patterns, std::vector<size_type>& results)
{
  size_type chars = 0;
  for(size_type i = 0; i < patterns.size(); i++) { chars += patterns[i].length(); }
  results.resize(patterns.size());

  printSize(name, sdsl::size_in_bytes(fmi), fmi.size());

  if(chars > 0)
  {
    double start = readTimer();
    size_type found = 0, matches = 0;
    for(size_type i = 0; i < patterns.size(); i++)
    {
      range_type range = fmi.find(patterns[i]);
      results[i] = Range::length(range);
      if(!(Range::empty(range))) { found++; matches += Range::length(range); }
    }
    double seconds = readTimer() - start;
    printTime(name, found, matches, chars, seconds);
  }

  std::cout << std::endl;
}

void
merge(FMI& result, FMI& index, FMI& increment)
{
  double increment_mb = inMegabytes(increment.size());

  double start = readTimer();
  result = FMI(index, increment);
  double seconds = readTimer() - start;
  std::cout << "BWTs merged in " << seconds << " seconds ("
            << (increment_mb / seconds) << " MB/s)" << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
