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

#define TEST_FMI
#define TEST_RLARRAY

void testFMI(std::string input_name, std::string pattern_name);
void testRLArray();

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
  std::cout << std::endl;

#ifdef TEST_FMI
  testFMI(argv[1], argv[2]);
#endif

#ifdef TEST_RLARRAY
  testRLArray();
#endif

  return 0;
}

//------------------------------------------------------------------------------

void
testFMI(std::string input_name, std::string pattern_name)
{
  std::vector<std::string> patterns;
  size_type chars = readRows(pattern_name, patterns, true);
  std::cout << "Read " << patterns.size() << " patterns of total length " << chars << std::endl;
  std::cout << std::endl;

  FMI fmi;
  {
    sdsl::int_vector_buffer<8> input(input_name);
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

  std::cout << "Memory usage: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

const size_type ARRAY_SIZE = 128 * MEGABYTE;
const size_type MAX_VALUE = GIGABYTE;

void
report(const RLArray& array, std::string name)
{
  std::cout << name << ": " << sdsl::size_in_bytes(array) << " bytes, " << array.size() << " runs, "
            << array.values() << " values, " << array.bytes() << " code bytes" << std::endl;
}

void
testRLArray()
{
  std::mt19937_64 rng(0xDEADBEEF);

  std::vector<RLArray::value_type> a_vec(ARRAY_SIZE), b_vec(ARRAY_SIZE);
  for(size_type i = 0; i < ARRAY_SIZE; i++)
  {
    a_vec[i] = rng() % MAX_VALUE; b_vec[i] = rng() % MAX_VALUE;
  }
  std::cout << "Created two arrays of size " << ARRAY_SIZE << std::endl;
  std::cout << std::endl;

  RLArray a(a_vec), b(b_vec);
  report(a, "A"); report(b, "B");

  double start = readTimer();
  RLArray merged(a, b);
  std::cout << "Merged the arrays in " << (readTimer() - start) << " seconds" << std::endl;
  report(merged, "A+B");
  std::cout << std::endl;

  a_vec.reserve(2 * ARRAY_SIZE);
  a_vec.insert(a_vec.end(), b_vec.begin(), b_vec.end());
  sdsl::util::clear(b_vec);
  parallelQuickSort(a_vec.begin(), a_vec.end());

  std::cout << "Starting verification" << std::endl;
  size_type pos = 0;
  for(RLArray::const_iterator iter(merged); !(iter.end()); ++iter)
  {
    for(size_type i = 0; i < iter->second; i++, pos++)
    {
      if(iter->first != a_vec[pos])
      {
        std::cerr << "Array[" << pos << "] = " << a_vec[pos]
                  << ", RLE[" << pos << "] = " << iter->first << std::endl;
        return;
      }
    }
  }
  std::cout << "Verification " << (pos == 2 * ARRAY_SIZE ? "complete" : "failed") << std::endl;
  std::cout << std::endl;

  std::cout << "Memory usage: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
