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

Alphabet
createAlphabet(AlphabeticOrder order)
{
  Alphabet alpha;

  switch(order)
  {
  case AO_DEFAULT:
    break;
  case AO_SORTED:
    std::swap(alpha.comp2char[4], alpha.comp2char[5]);
    std::swap(alpha.char2comp['N'], alpha.char2comp['T']);
    std::swap(alpha.char2comp['n'], alpha.char2comp['t']);
    break;
  default:
    break;
  }

  return alpha;
}

AlphabeticOrder
identifyAlphabet(const Alphabet& alpha)
{
  if(alpha.sorted()) { return AO_SORTED; }

  Alphabet default_alpha;
  if(alpha == default_alpha) { return AO_DEFAULT; }

  return AO_UNKNOWN;
}

std::string
alphabetName(AlphabeticOrder order)
{
  switch(order)
  {
  case AO_DEFAULT:
    return "default";
  case AO_SORTED:
    return "sorted";
  case AO_ANY:
    return "any";
  case AO_UNKNOWN:
  default:
    return "unknown";
  }
}

bool
compatible(const Alphabet& alpha, AlphabeticOrder order)
{
  Alphabet default_alpha;
  switch(order)
  {
  case AO_DEFAULT:
    return (alpha == default_alpha);
  case AO_SORTED:
    return alpha.sorted();
  case AO_ANY:
    return true;
  case AO_UNKNOWN:
  default:
    return false;
  }
}

//------------------------------------------------------------------------------

const std::string NativeFormat::name = "Native format";
const std::string NativeFormat::tag = "native";

const std::string PlainFormatD::name = "Plain format (default alphabet)";
const std::string PlainFormatD::tag = "plain_default";

const std::string PlainFormatS::name = "Plain format (sorted alphabet)";
const std::string PlainFormatS::tag = "plain_sorted";

const std::string RFMFormat::name = "RFM format";
const std::string RFMFormat::tag = "rfm";

const std::string SDSLFormat::name = "SDSL format";
const std::string SDSLFormat::tag = "sdsl";

const std::string SGAFormat::name = "SGA format";
const std::string SGAFormat::tag = "sga";

//------------------------------------------------------------------------------

template<class BufferType>
struct PlainData
{
  typedef bwtmerge::char_type char_type;

  const static size_type BUFFER_SIZE = MEGABYTE;

  static void read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts, const Alphabet& alpha)
  {
    data.clear();
    counts = sdsl::int_vector<64>(alpha.sigma, 0);

    RunBuffer run_buffer;
    size_type bytes = BufferType::readHeader(in);
    char_type* buffer = new char_type[BUFFER_SIZE];
    for(size_type offset = 0; offset < bytes; offset += BUFFER_SIZE)
    {
      size_type buffer_size = std::min(BUFFER_SIZE, bytes - offset);
      BufferType::readData(in, buffer, buffer_size);
      for(size_type i = 0; i < buffer_size; i++)
      {
        if(run_buffer.add(buffer[i]))
        {
          run_buffer.run.first = alpha.char2comp[run_buffer.run.first];
          Run::write(data, run_buffer.run);
          counts[run_buffer.run.first] += run_buffer.run.second;
        }
      }
    }
    run_buffer.flush();
    run_buffer.run.first = alpha.char2comp[run_buffer.run.first];
    Run::write(data, run_buffer.run);
    counts[run_buffer.run.first] += run_buffer.run.second;

    delete[] buffer; buffer = 0;
  }

  static void write(std::ofstream& out, const BlockArray& data, const Alphabet& alpha, const NativeHeader& info)
  {
    BufferType::writeHeader(out, info.bases);

    size_type rle_pos = 0, buffer_pos = 0;
    char_type* buffer = new char_type[BUFFER_SIZE];
    while(rle_pos < data.size())
    {
      range_type run = Run::read(data, rle_pos);
      run.first = alpha.comp2char[run.first];
      while(run.second > 0)
      {
        if(buffer_pos >= BUFFER_SIZE)
        {
          BufferType::writeData(out, buffer, buffer_pos);
          buffer_pos = 0;
        }
        size_type length = std::min(BUFFER_SIZE - buffer_pos, run.second); run.second -= length;
        for(size_type i = 0; i < length; i++, buffer_pos++) { buffer[buffer_pos] = run.first; }
      }
    }
    if(buffer_pos > 0) { BufferType::writeData(out, buffer, buffer_pos); }

    delete[] buffer; buffer = 0;
  }
};

/*
  A plain buffer type with interface compatible with the IntVectorBuffer template.
*/
template<class Element>
struct PlainBuffer
{
  typedef uint64_t code_type;

  static void writeHeader(std::ofstream&, size_type) {}

  static size_type readHeader(std::ifstream& in)
  {
    return fileSize(in);
  }

  static void writeData(std::ofstream& out, const Element* data, size_type elements)
  {
    size_type bytes = elements * sizeof(Element);
    out.write((const char*)data, bytes);
  }

  static void readData(std::ifstream& in, Element* data, size_type elements)
  {
    size_type bytes = elements * sizeof(Element);
    in.read((char*)data, bytes);
  }
};

//------------------------------------------------------------------------------

void
PlainFormatD::read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts)
{
  PlainData<PlainBuffer<char_type>>::read(in, data, counts, createAlphabet(order()));
}

void
PlainFormatD::write(std::ofstream& out, const BlockArray& data, const NativeHeader& info)
{
  PlainData<PlainBuffer<char_type>>::write(out, data, createAlphabet(order()), info);
}

//------------------------------------------------------------------------------

void
PlainFormatS::read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts)
{
  PlainData<PlainBuffer<char_type>>::read(in, data, counts, createAlphabet(order()));
}

void
PlainFormatS::write(std::ofstream& out, const BlockArray& data, const NativeHeader& info)
{
  PlainData<PlainBuffer<char_type>>::write(out, data, createAlphabet(order()), info);
}

//------------------------------------------------------------------------------

struct RFMData
{
  const static size_type SIGMA = 6;
};

void
RFMFormat::read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts)
{
  PlainData<IntVectorBuffer<comp_type>>::read(in, data, counts, Alphabet(RFMData::SIGMA));
}

void
RFMFormat::write(std::ofstream& out, const BlockArray& data, const NativeHeader& info)
{
  PlainData<IntVectorBuffer<comp_type>>::write(out, data, Alphabet(RFMData::SIGMA), info);
}

//------------------------------------------------------------------------------

void
SDSLFormat::read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts)
{
  PlainData<IntVectorBuffer<char_type>>::read(in, data, counts, createAlphabet(order()));
}

void
SDSLFormat::write(std::ofstream& out, const BlockArray& data, const NativeHeader& info)
{
  PlainData<IntVectorBuffer<char_type>>::write(out, data, createAlphabet(order()), info);
}

//------------------------------------------------------------------------------

struct SGAData
{
  typedef bwtmerge::comp_type comp_type;
  typedef bwtmerge::size_type length_type;
  typedef uint8_t             code_type;

  inline static code_type encode(comp_type comp, length_type length) { return (comp << RUN_BITS) | length; }
  inline static comp_type comp(code_type code) { return (code >> RUN_BITS); }
  inline static length_type length(code_type code) { return (code & RUN_MASK); }

  const static code_type RUN_MASK = 0x1F;
  const static size_type RUN_BITS = 5;
  const static size_type MAX_RUN = 31;
  const static size_type SIGMA = 6;
};

void
SGAFormat::read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts)
{
  data.clear();
  counts = sdsl::int_vector<64>(SGAData::SIGMA, 0);

  SGAHeader header(in);
  if(!(header.check()))
  {
    std::cerr << "SGAFormat::load(): Invalid header!" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  SGAData::code_type buffer[MEGABYTE];
  RunBuffer run_buffer;
  for(size_type offset = 0; offset < header.runs; offset += MEGABYTE)
  {
    size_type bytes = std::min(MEGABYTE, header.runs - offset);
    in.read((char*)buffer, bytes);
    for(size_type i = 0; i < bytes; i++)
    {
      if(run_buffer.add(SGAData::comp(buffer[i]), SGAData::length(buffer[i])))
      {
        Run::write(data, run_buffer.run);
        counts[run_buffer.run.first] += run_buffer.run.second;
      }
    }
  }
  run_buffer.flush();
  Run::write(data, run_buffer.run);
  counts[run_buffer.run.first] += run_buffer.run.second;
}

void
SGAFormat::write(std::ofstream& out, const BlockArray& data, const NativeHeader& info)
{
  SGAHeader header; header.bases = info.bases; header.sequences = info.sequences;
  size_type rle_pos = 0;
  while(rle_pos < data.size())
  {
    range_type run = Run::read(data, rle_pos);
    header.runs += (run.second + SGAData::MAX_RUN - 1) / SGAData::MAX_RUN;
  }
  header.write(out);

  rle_pos = 0;
  std::vector<SGAData::code_type> buffer;
  while(rle_pos < data.size())
  {
    range_type run = Run::read(data, rle_pos);
    while(run.second > SGAData::MAX_RUN)
    {
      buffer.push_back(SGAData::encode(run.first, SGAData::MAX_RUN));
      run.second -= SGAData::MAX_RUN;
    }
    buffer.push_back(SGAData::encode(run.first, run.second));
    if(buffer.size() >= MEGABYTE)
    {
      out.write((char*)(buffer.data()), buffer.size());
      buffer.clear();
    }
  }
  out.write((char*)(buffer.data()), buffer.size());
}

//------------------------------------------------------------------------------

NativeHeader::NativeHeader() :
  tag(DEFAULT_TAG), flags(0), sequences(0), bases(0)
{
}

NativeHeader::NativeHeader(std::istream& in)
{
  sdsl::read_member(this->tag, in);
  sdsl::read_member(this->flags, in);
  sdsl::read_member(this->sequences, in);
  sdsl::read_member(this->bases, in);
}

void
NativeHeader::write(std::ostream& out) const
{
  sdsl::write_member(this->tag, out);
  sdsl::write_member(this->flags, out);
  sdsl::write_member(this->sequences, out);
  sdsl::write_member(this->bases, out);
}

bool
NativeHeader::check() const
{
  return (this->tag == DEFAULT_TAG);
}

AlphabeticOrder
NativeHeader::order() const
{
  return static_cast<AlphabeticOrder>((this->flags) & ALPHABET_MASK);
}

void
NativeHeader::setOrder(AlphabeticOrder ao)
{
  this->flags &= ~ALPHABET_MASK;
  this->flags |= static_cast<uint32_t>(ao) & ALPHABET_MASK;
}

std::ostream& operator<<(std::ostream& stream, const NativeHeader& header)
{
  return stream << "Native format: " << header.sequences << " sequences, " << header.bases << " bases, "
                << alphabetName(header.order()) << " alphabet";
}

//------------------------------------------------------------------------------

SGAHeader::SGAHeader() :
  tag(DEFAULT_TAG), sequences(0), bases(0), runs(0), flag(DEFAULT_FLAG)
{
}

SGAHeader::SGAHeader(std::istream& in)
{
  sdsl::read_member(this->tag, in);
  sdsl::read_member(this->sequences, in);
  sdsl::read_member(this->bases, in);
  sdsl::read_member(this->runs, in);
  sdsl::read_member(this->flag, in);
}

void
SGAHeader::write(std::ostream& out) const
{
  sdsl::write_member(this->tag, out);
  sdsl::write_member(this->sequences, out);
  sdsl::write_member(this->bases, out);
  sdsl::write_member(this->runs, out);
  sdsl::write_member(this->flag, out);
}

bool
SGAHeader::check() const
{
  return (this->tag == DEFAULT_TAG && this->flag == DEFAULT_FLAG);
}

std::ostream& operator<<(std::ostream& stream, const SGAHeader& header)
{
  return stream << "SGA format: " << header.sequences << " sequences, "
                << header.bases << " bases, " << header.runs << " runs";
}

//------------------------------------------------------------------------------

} // namespace bwtmerge
