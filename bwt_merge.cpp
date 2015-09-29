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

#include <string>
#include <unistd.h>

#include "fmi.h"

using namespace bwtmerge;

//------------------------------------------------------------------------------

const size_type RUN_BUFFER_SIZE = MEGABYTE;

//------------------------------------------------------------------------------

void printUsage();

void verifyFMI(FMI& fmi, const std::string& name,
  const std::vector<std::string>& patterns, std::vector<size_type>& results);

void merge(FMI& result, FMI& index, FMI& increment, const MergeParameters& parameters);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 4)
  {
    printUsage();
    std::exit(EXIT_SUCCESS);
  }

  std::cout << "BWT-merge" << std::endl;
  std::cout << std::endl;

  int c = 0;
  bool verify = false;
  MergeParameters parameters;
  std::string index_name, increment_name, output_name, pattern_name;

  while((c = getopt(argc, argv, "b:m:r:s:t:d:v:")) != -1)
  {
    switch(c)
    {
    case 'b':
      parameters.setTB(std::stoul(optarg));
      break;
    case 'm':
      parameters.setMB(std::stoul(optarg));
      break;
    case 'r':
      parameters.setRB(std::stoul(optarg));
      break;
    case 's':
      parameters.setSB(std::stoul(optarg));
      break;
    case 't':
      parameters.setT(std::stoul(optarg));
      break;
    case 'd':
      parameters.setTemp(optarg);
      break;
    case 'v':
      pattern_name = optarg; verify = true;
      break;
    case '?':
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  parameters.sanitize();
  if(optind < argc) { index_name = argv[optind]; }
  else
  {
    std::cerr << "bwt_merge: Input 1 unspecified!" << std::endl;
  }
  if(optind + 1 < argc) { increment_name = argv[optind + 1]; }
  else
  {
    std::cerr << "bwt_merge: Input 2 unspecified!" << std::endl;
  }
  if(optind + 2 < argc) { output_name = argv[optind + 2]; }
  else
  {
    std::cerr << "bwt_merge: Output file unspecified!" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "Input 1:            " << index_name << std::endl;
  std::cout << "Input 2:            " << increment_name << std::endl;
  std::cout << "Output:             " << output_name << std::endl;
  if(verify)
  {
    std::cout << "Patterns:           " << pattern_name << std::endl;
  }
  std::cout << std::endl;
  std::cout << parameters;
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
  merge(merged, index, increment, parameters);
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
printUsage()
{
  std::cerr << "Usage: bwt_merge [options] input1 input2 output" << std::endl;
  std::cerr << std::endl;

  std::cerr << "Options:" << std::endl;
  std::cerr << "  -b N          Set thread buffer size to N megabytes / thread (default: "
            << MergeParameters::defaultTB() << ")" << std::endl;
  std::cerr << "  -m N          Set the number of merge buffers to N (default: "
            << MergeParameters::defaultMB() << ")" << std::endl;
  std::cerr << "  -r N          Set run buffer size to N megabytes / thread (default: "
            << MergeParameters::defaultRB() << ")" << std::endl;
  std::cerr << "  -s N          Set the number of sequence blocks to N (default: "
            << MergeParameters::defaultSB() << " / thread)" << std::endl;
  std::cerr << "  -t N          Use N parallel threads (default: " << MergeParameters::defaultT()
            << " on this system)" << std::endl;
  std::cerr << std::endl;
  std::cerr << "  -d directory  Use the given directory for temporary files (default: .)" << std::endl;
  std::cerr << "  -v filename   Verify by querying with patterns from the given file" << std::endl;
  std::cerr << std::endl;
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
merge(FMI& result, FMI& index, FMI& increment, const MergeParameters& parameters)
{
  double increment_mb = inMegabytes(increment.size());

  double start = readTimer();
  result = FMI(index, increment, parameters);
  double seconds = readTimer() - start;
  std::cout << "BWTs merged in " << seconds << " seconds ("
            << (increment_mb / seconds) << " MB/s)" << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
