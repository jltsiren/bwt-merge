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

#include <sstream>
#include <string>
#include <unistd.h>

#include "fmi.h"

using namespace bwtmerge;

//------------------------------------------------------------------------------

const size_type RUN_BUFFER_SIZE = MEGABYTE;

//------------------------------------------------------------------------------

void printUsage();

/*
  The number of occurrences for each pattern will be added to the results.
*/
void verifyFMI(FMI& fmi, const std::string& name,
  const std::vector<std::string>& patterns, std::vector<size_type>& results);

void merge(FMI& index, FMI& increment, const MergeParameters& parameters);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    printUsage();
    std::exit(EXIT_SUCCESS);
  }

  double start = readTimer();
  std::cout << "BWT-merge" << std::endl;
  std::cout << std::endl;

  int c = 0;
  bool verify = false;
  MergeParameters parameters;
  std::string pattern_name, output_format;
  std::vector<std::string> input_formats;
  while((c = getopt(argc, argv, "b:m:r:s:t:d:v:i:o:")) != -1)
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
    case 'i':
      tokenize(optarg, input_formats, ',');
      break;
    case 'o':
      output_format = optarg;
      break;
    case '?':
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  int inputs = (argc - 1) - optind;
  if(inputs < 2)
  {
    std::cerr << "bwt_merge: Output file not specified" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if(input_formats.empty())
  {
    for(int i = 0; i < inputs; i++) { input_formats.push_back(NativeFormat::tag); }
  }
  if(input_formats.size() == 1)
  {
    for(int i = 1; i < inputs; i++) { input_formats.push_back(input_formats[0]); }
  }
  if(input_formats.size() != (unsigned)inputs)
  {
    std::cerr << "bwt_merge: Specified " << input_formats.size() << " formats for "
              << inputs << " inputs" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if(output_format.length() == 0) { output_format = NativeFormat::tag; }
  parameters.sanitize();
  Parallel::max_threads = parameters.threads;

  for(int i = optind; i < argc - 1; i++)
  {
    std::cout << "Input:            " << argv[i] << " (" << input_formats[i - optind] << ")" << std::endl;
  }
  std::cout << "Output:           " << argv[argc - 1] << " (" << output_format << ")" << std::endl;
  if(verify)
  {
    std::cout << "Patterns:         " << pattern_name << std::endl;
  }
  std::cout << std::endl;
  std::cout << parameters;
  std::cout << std::endl;

  std::vector<std::string> patterns;
  std::vector<size_type> pre_results, post_results;
  if(verify)
  {
    size_type chars = readRows(pattern_name, patterns, true);
    pre_results = std::vector<size_type>(patterns.size(), 0);
    post_results = std::vector<size_type>(patterns.size(), 0);
    std::cout << "Read " << patterns.size() << " patterns of total length " << chars << std::endl;
    std::cout << std::endl;
  }

  FMI index; load(index, argv[optind], input_formats[0]);
  verifyFMI(index, "Input", patterns, pre_results);

  size_type bytes_added = 0;
  for(int input = 1; input < inputs; input++)
  {
    FMI increment; load(increment, argv[optind + input], input_formats[input]);
    bytes_added += increment.size();
    verifyFMI(increment, "Input", patterns, pre_results);
    merge(index, increment, parameters);
  }

  serialize(index, argv[argc - 1], output_format);
  verifyFMI(index, "Output", patterns, post_results);

  if(verify)
  {
    size_type errors = 0;
    for(size_type i = 0; i < patterns.size(); i++)
    {
      if(pre_results[i] != post_results[i]) { errors++; }
    }
    if(errors > 0)
    {
      std::cout << "Verification failed for " << errors << " patterns" << std::endl;
    }
    else
    {
      std::cout << "Verification successful" << std::endl;
    }
    std::cout << std::endl;
  }

  double seconds = readTimer() - start;
  std::cout << "Total time:       " << seconds << " seconds (" << (inMegabytes(bytes_added) / seconds)
            << " MB/s)" << std::endl;
  std::cout << "Peak memory:      " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage()
{
  std::cerr << "Usage: bwt_merge [options] input1 input2 [input3 ...] output" << std::endl;
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

  std::cerr << "  -i formats    Read the inputs in the given formats (default: native)" << std::endl;
  std::cerr << "                Multiple comma-separated formats can be provided." << std::endl;
  std::cerr << "  -o format     Write the output in the given format (default: native)" << std::endl;
  std::cerr << std::endl;

  printFormats(std::cerr);
}

//------------------------------------------------------------------------------

void
queryFMI(ParallelLoop& loop, const FMI& fmi, const std::vector<std::string>& patterns,
  std::vector<size_type>& results,
  std::atomic<size_type>& total_found, std::atomic<size_type>& total_matches)
{
  while(true)
  {
    range_type range = loop.next();
    if(Range::empty(range)) { return; }

    size_type found = 0, matches = 0;
    for(size_type i = range.first; i <= range.second; i++)
    {
      range_type result = fmi.find(patterns[i]);
      results[i] += Range::length(result);
      if(!(Range::empty(range))) { found++; matches += Range::length(result); }
    }

    total_found += found; total_matches += matches;
  }
}

void
verifyFMI(FMI& fmi, const std::string& name,
  const std::vector<std::string>& patterns, std::vector<size_type>& results)
{
  size_type chars = 0;
  for(size_type i = 0; i < patterns.size(); i++) { chars += patterns[i].length(); }

  printSize(name, sdsl::size_in_bytes(fmi), fmi.size());

  if(chars > 0)
  {
    double start = readTimer();
    std::atomic<size_type> found(0), matches(0);
    {
      ParallelLoop loop(0, patterns.size(), Parallel::max_threads, Parallel::max_threads);
      loop.execute(queryFMI, std::ref(fmi), std::ref(patterns),
        std::ref(results), std::ref(found), std::ref(matches));
    }
    double seconds = readTimer() - start;
    printTime(name, found, matches, chars, seconds);
  }

  std::cout << std::endl;
}

void
merge(FMI& index, FMI& increment, const MergeParameters& parameters)
{
  double increment_mb = inMegabytes(increment.size());

  double start = readTimer();
  FMI temp(index, increment, parameters);
  index.swap(temp);
  double seconds = readTimer() - start;
  std::cout << "BWTs merged in " << seconds << " seconds ("
            << (increment_mb / seconds) << " MB/s)" << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
