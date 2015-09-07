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

#ifndef _BWTMERGE_FORMATS_H
#define _BWTMERGE_FORMATS_H

#include "support.h"

namespace bwtmerge
{

//------------------------------------------------------------------------------

enum AlphabeticOrder { AO_DEFAULT, AO_SORTED };

Alphabet createAlphabet(AlphabeticOrder order);
void readPlain(std::istream& in, BlockArray& data, sdsl::int_vector<64>& counts, const Alphabet& alpha);
void writePlain(std::ostream& out, const BlockArray& data, const Alphabet& alpha);

//------------------------------------------------------------------------------

/*
  BWT file formats. Note that BWT formats are only compatible with those having the
  same alphabetic order.

    load()      reads the BWT from 'in' and stores it in the native format in 'data'
    write()     writes the BWT stored in the native format in 'data' to 'out'
    order()     returns the alphabetic order

    PlainFormat - BWT as a character array; any alphabetic order
    RFMFormat   - BWT as int_vector<8> of comp values; AO_SORTED
    SDSLFormat  - BWT as int_vector<8> of characters; AO_SORTED
    SGAFormat   - SGA assembler; AO_DEFAULT
*/

template<AlphabeticOrder ao>
struct PlainFormat
{
  static void read(std::istream& in, BlockArray& data, sdsl::int_vector<64>& counts)
  {
    readPlain(in, data, counts, createAlphabet(order()));
  }

  static void write(std::ostream& out, const BlockArray& data)
  {
    writePlain(out, data, createAlphabet(order()));
  }

  inline static AlphabeticOrder order() { return ao; }
};

struct RFMFormat
{
  static void read(std::istream& in, BlockArray& data, sdsl::int_vector<64>& counts);
  static void write(std::ostream& out, const BlockArray& data);
  inline static AlphabeticOrder order() { return AO_SORTED; }
};

struct SDSLFormat
{
  static void read(std::istream& in, BlockArray& data, sdsl::int_vector<64>& counts);
  static void write(std::ostream& out, const BlockArray& data);
  inline static AlphabeticOrder order() { return AO_SORTED; }
};

struct SGAFormat
{
  static void read(std::istream& in, BlockArray& data, sdsl::int_vector<64>& counts);
  static void write(std::ostream& out, const BlockArray& data);
  inline static AlphabeticOrder order() { return AO_DEFAULT; }
};

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

  SGAHeader(std::istream& in);
  SGAHeader();

  void write(std::ostream& out);
  bool check();
};

std::ostream& operator<<(std::ostream& stream, const SGAHeader& header);

//------------------------------------------------------------------------------

} // namespace bwtmerge

#endif // _BWTMERGE_FORMATS_H
