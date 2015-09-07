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
#include "formats.h"

using namespace bwtmerge;

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: bwt_convert input output" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::cout << "BWT-converter" << std::endl;
  std::cout << std::endl;

  std::string input = argv[1], output = argv[2];
  std::cout << "Input: " << input << std::endl;
  std::cout << "Output: " << output << std::endl;
  std::cout << std::endl;

  double start = readTimer();
  FMI fmi; fmi.load<SGAFormat>(argv[1]);
  double seconds = readTimer() - start;
  std::cout << "BWT converted in " << seconds << " seconds ("
            << (inMegabytes(fmi.size()) / seconds) << " MB/s)" << std::endl;
  std::cout << std::endl;

  printSize("FMI", sdsl::size_in_bytes(fmi), fmi.size());
  std::cout << std::endl;

  sdsl::store_to_file(fmi, argv[2]);
  std::cout << "FMI written to disk" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------
