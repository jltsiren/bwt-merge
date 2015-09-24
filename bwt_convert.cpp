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

#include <unistd.h>

#include "fmi.h"

using namespace bwtmerge;

//------------------------------------------------------------------------------

void printUsage();

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    printUsage();
    std::exit(EXIT_SUCCESS);
  }

  std::cout << "BWT converter" << std::endl;
  std::cout << std::endl;

  int c = 0;
  std::string input_tag = SGAFormat::tag, output_tag = NativeFormat::tag;
  std::string input_name, output_name;

  while((c = getopt(argc, argv, "i:o:")) != -1)
  {
    switch(c)
    {
    case 'i':
      input_tag = optarg;
      break;
    case 'o':
      output_tag = optarg;
      break;
    case '?':
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  if(optind < argc) { input_name = argv[optind]; }
  else
  {
    std::cerr << "bwt_convert: Input file unspecified!" << std::endl;
  }
  if(optind + 1 < argc) { output_name = argv[optind + 1]; }
  else
  {
    std::cerr << "bwt_convert: Output file unspecified!" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "Input:   " << input_name << " (" << input_tag << ")" << std::endl;
  std::cout << "Output:  " << output_name << " (" << output_tag << ")" << std::endl;
  std::cout << std::endl;

  double start = readTimer();
  FMI fmi; load(fmi, input_name, input_tag);
  size_type size = fmi.size();
  printSize("FMI", sdsl::size_in_bytes(fmi), fmi.size());
  std::cout << std::endl;
  serialize(fmi, output_name, output_tag);
  double seconds = readTimer() - start;

  std::cout << "BWT converted in " << seconds << " seconds ("
            << (inMegabytes(size) / seconds) << " MB/s)" << std::endl;
  std::cout << std::endl;

  std::cout << "Memory usage: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

template<class Format>
void
printFormat(std::ostream& stream)
{
  std::string padding;
  if(Format::tag.length() < 15) { padding = std::string(15 - Format::tag.length(), ' '); }
  stream << "  " << Format::tag << padding << Format::name << std::endl;
}

void
printUsage()
{
  std::cerr << "Usage: bwt_convert [options] input output" << std::endl;
  std::cerr << std::endl;

  std::cerr << "Options:" << std::endl;
  std::cerr << "  -i format      Read the input in the given format (default: sga)" << std::endl;
  std::cerr << "  -o format      Write the output in the given format (default: native)" << std::endl;
  std::cerr << std::endl;

  std::cerr << "Formats supporting any alphabetic order:" << std::endl;
  printFormat<NativeFormat>(std::cerr);
  std::cerr << std::endl;

  std::cerr << "Formats using the default alphabet:" << std::endl;
  printFormat<PlainFormatD>(std::cerr);
  printFormat<SGAFormat>(std::cerr);
  std::cerr << std::endl;

  std::cerr << "Formats using sorted alphabet:" << std::endl;
  printFormat<PlainFormatS>(std::cerr);
  printFormat<RFMFormat>(std::cerr);
  printFormat<SDSLFormat>(std::cerr);
  std::cerr << std::endl;
}

//------------------------------------------------------------------------------
