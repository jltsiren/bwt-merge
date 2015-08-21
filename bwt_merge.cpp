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

#include "fmi.h"

using namespace bwtmerge;

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: bwt_merge input patterns" << std::endl;
    std::cerr << "  This is a temporary verification program." << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::cout << "BWT-merge" << std::endl;
  std::cout << std::endl;

  std::cout << "Input: " << argv[1] << std::endl;
  std::cout << "Patterns: " << argv[2] << std::endl;
  std::cout << std::endl;

  std::vector<std::string> patterns;
  size_type chars = readRows(argv[2], patterns, true);
  std::cout << "Read " << patterns.size() << " patterns of total length " << chars << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  FMI fmi;
  {
    sdsl::int_vector_buffer<8> input(argv[1]);
    Alphabet comp_alpha = Alphabet(BWT::SIGMA);
    fmi.bwt = BWT(input, comp_alpha);
    sdsl::int_vector<64> counts;
    fmi.bwt.characterCounts(counts);
    Alphabet rfm_alpha = rfmAlphabet();
    fmi.alpha = Alphabet(counts, rfm_alpha.char2comp, rfm_alpha.comp2char);
  }

  double start = readTimer();
  size_type found = 0, matches = 0;
  for(auto pattern : patterns)
  {
    range_type range = fmi.find(pattern);
    if(!(Range::empty(range))) { found++; matches += Range::length(range); }
  }
  double seconds = readTimer() - start;

  printSize("FMI", sdsl::size_in_bytes(fmi), fmi.size());
  printTime("FMI", found, matches, chars, seconds);
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "Memory usage: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------
