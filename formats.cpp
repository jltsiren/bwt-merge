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

#include <cstring>

#include "formats.h"

namespace bwtmerge
{

//------------------------------------------------------------------------------
/*
void
RFMFormat::Read(std::istream& in, BlockArray& data)
{
}
*/
//------------------------------------------------------------------------------

struct SGAHeader
{
  uint16_t tag;
  uint64_t reads;
  uint64_t bases;
  uint64_t runs;
  uint32_t flag;

  const static uint16_t DEFAULT_TAG = 0xCACA;
  const static uint32_t DEFAULT_FLAG = 0;

  SGAHeader(std::istream& in)
  {
    sdsl::read_member(this->tag, in);
    sdsl::read_member(this->reads, in);
    sdsl::read_member(this->bases, in);
    sdsl::read_member(this->runs, in);
    sdsl::read_member(this->flag, in);
  }

  SGAHeader() :
    tag(DEFAULT_TAG), reads(0), bases(0), runs(0), flag(DEFAULT_FLAG)
  {
  }

  void write(std::ostream& out)
  {
    sdsl::write_member(this->tag, out);
    sdsl::write_member(this->reads, out);
    sdsl::write_member(this->bases, out);
    sdsl::write_member(this->runs, out);
    sdsl::write_member(this->flag, out);
  }

  bool check()
  {
    return (this->tag == DEFAULT_TAG && this->flag == DEFAULT_FLAG);
  }
};

void
SGAFormat::read(std::istream& in, BlockArray& data)
{
  data.clear();

  SGAHeader header(in);
  if(!(header.check()))
  {
    std::cerr << "SGAFormat::load(): Invalid header" << std::endl;
    return;
  }

  code_type buffer[MEGABYTE];
  range_type run(0, 0);
  for(size_type offset = 0; offset < header.runs; offset++)
  {
    size_type bytes = std::min(MEGABYTE, header.runs - offset);
    in.read((char*)buffer, bytes);
    for(size_type i = 0; i < bytes; i++)
    {
      if(comp(buffer[i]) != run.first)
      {
        if(run.second > 0) { Run::write(data, run.first, run.second); }
        run.first = comp(buffer[i]); run.second = length(buffer[i]);
      }
      else { run.second += length(buffer[i]); }
    }
  }
  if(run.second > 0) { Run::write(data, run.first, run.second); }
}

void
SGAFormat::write(std::ostream& out, const BlockArray& data)
{
  SGAHeader header;
  size_type rle_pos = 0;
  while(rle_pos < data.size())
  {
    range_type run = Run::read(data, rle_pos);
    if(run.first == 0) { header.reads += run.second; }
    header.bases += run.second;
    header.runs += (run.second + MAX_RUN - 1) / MAX_RUN;
  }
  header.write(out);

  rle_pos = 0;
  std::vector<code_type> buffer;
  while(rle_pos < data.size())
  {
    range_type run = Run::read(data, rle_pos);
    while(run.second > MAX_RUN)
    {
      buffer.push_back(encode(run.first, MAX_RUN));
      run.second -= MAX_RUN;
    }
    buffer.push_back(encode(run.first, run.second));
    if(buffer.size() >= MEGABYTE)
    {
      out.write((char*)(buffer.data()), buffer.size());
      buffer.clear();
    }
  }
  out.write((char*)(buffer.data()), buffer.size());
}

//------------------------------------------------------------------------------

} // namespace bwtmerge
