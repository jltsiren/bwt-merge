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

#include <cstdlib>

#include <sys/resource.h>
#include <unistd.h>

#include "utils.h"

namespace bwtmerge
{

//------------------------------------------------------------------------------

void
printHeader(const std::string& header, size_type indent)
{
  std::string padding;
  if(header.length() + 1 < indent) { padding = std::string(indent - 1 - header.length(), ' '); }
  std::cout << header << ":" << padding;
}

void
printSize(const std::string& header, size_type bytes, size_type data_size, size_type indent)
{
  printHeader(header, indent);
  std::cout << inMegabytes(bytes) << " MB (" << inBPC(bytes, data_size) << " bpc)" << std::endl;
}

void
printTime(const std::string& header, size_type found, size_type matches, size_type bytes, double seconds, size_type indent)
{
  printHeader(header, indent);
  std::cout << "Found " << found << " patterns with " << matches << " occ in "
            << seconds << " seconds (" << (inMegabytes(bytes) / seconds) << " MB/s)" << std::endl;
}

void
printTime(const std::string& header, size_type queries, double seconds, size_type indent)
{
  printHeader(header, indent);
  std::cout << queries << " queries in " << seconds << " seconds ("
            << inMicroseconds(seconds / queries) << " µs/query)" << std::endl;
}

//------------------------------------------------------------------------------

double
readTimer()
{
  return omp_get_wtime();
}

size_type
memoryUsage()
{
  rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#ifdef RUSAGE_IN_BYTES
  return usage.ru_maxrss;
#else
  return KILOBYTE * usage.ru_maxrss;
#endif
}

//------------------------------------------------------------------------------

size_type
readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "readRows(): Cannot open input file " << filename << std::endl;
    return 0;
  }

  size_type chars = 0;
  while(input)
  {
    std::string buf;
    std::getline(input, buf);
    if(skip_empty_rows && buf.length() == 0) { continue; }
    rows.push_back(buf);
    chars += buf.length();
  }

  input.close();
  return chars;
}

std::string
tempFile(const std::string& name_part)
{
  char hostname[32];
  gethostname(hostname, 32); hostname[31] = 0;

  return name_part + '_'
    + std::string(hostname) + '_'
    + sdsl::util::to_string(sdsl::util::pid()) + '_'
    + sdsl::util::to_string(sdsl::util::id());
}

size_type
fileSize(std::ifstream& file)
{
  std::streamoff curr = file.tellg();

  file.seekg(0, std::ios::end);
  std::streamoff size = file.tellg();
  file.seekg(0, std::ios::beg);
  size -= file.tellg();

  file.seekg(curr, std::ios::beg);
  return size;
}

size_type
fileSize(std::ofstream& file)
{
  std::streamoff curr = file.tellp();

  file.seekp(0, std::ios::end);
  std::streamoff size = file.tellp();
  file.seekp(0, std::ios::beg);
  size -= file.tellp();

  file.seekp(curr, std::ios::beg);
  return size;
}

//------------------------------------------------------------------------------

std::vector<range_type>
getBounds(range_type range, size_type blocks)
{
  blocks = std::min(blocks, Range::length(range));
  blocks = std::max(blocks, (size_type)1);

  std::vector<range_type> bounds(blocks);
  for(size_type block = 0, start = range.first; block < blocks; block++)
  {
    bounds[block].first = start;
    if(start <= range.second)
    {
      start += std::max((size_type)1, (range.second + 1 - start) / (blocks - block));
    }
    bounds[block].second = start - 1;
  }

  return bounds;
}

//------------------------------------------------------------------------------

SDVectorBuilder::SDVectorBuilder() :
  bits(0), onebits(0),
  low_width(0),
  tail(0), prev_high(0), high_pos(0)
{
}

SDVectorBuilder::SDVectorBuilder(size_type _size, size_type _capacity) :
  bits(_size), onebits(_capacity),
  low_width(0),
  tail(0), prev_high(0), high_pos(0)
{
  if(this->onebits > this->bits) { this->onebits = this->bits; }
  if(this->onebits == 0) { return; }

  size_type logn = bit_length(this->size()), logm = bit_length(this->capacity());
  if(logm == logn) { logm--; }
  this->low_width = logn - logm;
  this->low = sdsl::int_vector<0>(this->capacity(), 0, this->low_width);
  this->high = sdsl::bit_vector(this->capacity() + (((size_type)1) << logm), 0);
}

//------------------------------------------------------------------------------

} // namespace bwtmerge

namespace sdsl
{

template<>
template<>
sd_vector<>::sd_vector<bwtmerge::SDVectorBuilder*>(
  bwtmerge::SDVectorBuilder* builder,
  bwtmerge::SDVectorBuilder*)
{
  this->m_size = builder->size();
  this->m_wl = builder->low_width;
  this->m_low.swap(builder->low);
  this->m_high.swap(builder->high);

  util::init_support(this->m_high_1_select, &(this->m_high));
  util::init_support(this->m_high_0_select, &(this->m_high));
}

} // namespace sdsl

//------------------------------------------------------------------------------
