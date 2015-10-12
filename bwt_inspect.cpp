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

#include "formats.h"

using namespace bwtmerge;

//------------------------------------------------------------------------------

template<class HeaderFormat>
bool inspect(std::ifstream& in, size_type& total_sequences, size_type& total_bases);

template<>
bool inspect<RopeHeader>(std::ifstream& in, size_type& total_sequences, size_type& total_bases);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: bwt_inspect input1 [input2 ...]" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::cout << "Inspecting BWT files" << std::endl;
  std::cout << std::endl;

  size_type total_sequences = 0, total_bases = 0;
  for(int arg = 1; arg < argc; arg++)
  {
    std::cout << argv[arg] << ": "; std::cout.flush();

    std::ifstream in(argv[arg], std::ios_base::binary);
    if(!in)
    {
      std::cerr << "bwt_inspect: Cannot open input file " << argv[arg] << std::endl;
      continue;
    }

    if(inspect<NativeHeader>(in, total_sequences, total_bases)) { continue; }
    if(inspect<SGAHeader>(in, total_sequences, total_bases)) { continue; }
    if(inspect<RopeHeader>(in, total_sequences, total_bases)) { continue; }

    in.close();
    std::cout << "Unknown format" << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Total: " << total_sequences << " sequences, " << total_bases << " bases" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

template<class HeaderFormat>
bool
inspect(std::ifstream& in, size_type& total_sequences, size_type& total_bases)
{
  in.seekg(0);
  HeaderFormat header; header.load(in);
  if(!(header.check())) { return false; }

  total_sequences += header.sequences; total_bases += header.bases;
  in.close();
  std::cout << header << std::endl;
  return true;
}

template<>
bool
inspect<RopeHeader>(std::ifstream& in, size_type&, size_type&)
{
  in.seekg(0);
  RopeHeader header; header.load(in);
  if(!(header.check())) { return false; }

  in.close();
  std::cout << header << std::endl;
  return true;
}

//------------------------------------------------------------------------------
