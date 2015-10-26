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

enum AlphabeticOrder { AO_DEFAULT = 0, AO_SORTED = 1, AO_ANY = 254, AO_UNKNOWN = 255 };

Alphabet createAlphabet(AlphabeticOrder order);
AlphabeticOrder identifyAlphabet(const Alphabet& alpha);
std::string alphabetName(AlphabeticOrder order);
bool compatible(const Alphabet& alpha, AlphabeticOrder order);

//------------------------------------------------------------------------------

struct NativeHeader
{
  uint32_t tag;
  uint32_t flags;
  uint64_t sequences;
  uint64_t bases;

  const static uint32_t DEFAULT_TAG = 0x54574221;
  const static uint32_t ALPHABET_MASK = 0xFF;

  NativeHeader();

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);
  bool check() const;

  AlphabeticOrder order() const;
  void setOrder(AlphabeticOrder ao);
};

std::ostream& operator<<(std::ostream& stream, const NativeHeader& header);

//------------------------------------------------------------------------------

/*
  BWT file formats. Note that BWT formats are only compatible with those having the
  same alphabetic order.

  load()        reads the BWT from 'in' and stores it in the native format in 'data'
  write()       writes the BWT stored in the native format in 'data' to 'out'
  order()       returns the alphabetic order

  name          meaningful name of the format
  tag           the name used for specifying the format

  NativeFormat  Native BWT format; any alphabetic order
  PlainFormatD  BWT as a character array; AO_DEFAULT
  PlainFormatS  BWT as a character array; AO_SORTED
  RFMFormat     BWT as int_vector<8> of comp values; AO_SORTED
  SDSLFormat    BWT as int_vector<8> of characters; AO_SORTED
  RopeFormat    RopeBWT; AO_DEFAULT
  SGAFormat     SGA assembler; AO_DEFAULT
*/

struct NativeFormat
{
  inline static AlphabeticOrder order() { return AO_ANY; }

  const static std::string name;
  const static std::string tag;
};

struct PlainFormatD
{
  static void read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts);
  static void write(std::ofstream& out, const BlockArray& data, const NativeHeader&);

  inline static AlphabeticOrder order() { return AO_DEFAULT; }

  const static std::string name;
  const static std::string tag;
};

struct PlainFormatS
{
  static void read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts);
  static void write(std::ofstream& out, const BlockArray& data, const NativeHeader&);

  inline static AlphabeticOrder order() { return AO_SORTED; }

  const static std::string name;
  const static std::string tag;
};

struct RFMFormat
{
  static void read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts);
  static void write(std::ofstream& out, const BlockArray& data, const NativeHeader& info);
  inline static AlphabeticOrder order() { return AO_SORTED; }

  const static std::string name;
  const static std::string tag;
};

struct SDSLFormat
{
  static void read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts);
  static void write(std::ofstream& out, const BlockArray& data, const NativeHeader& info);
  inline static AlphabeticOrder order() { return AO_SORTED; }

  const static std::string name;
  const static std::string tag;
};

struct RopeFormat
{
  static void read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts);
  static void write(std::ofstream& out, const BlockArray& data, const NativeHeader& info);
  inline static AlphabeticOrder order() { return AO_DEFAULT; }

  const static std::string name;
  const static std::string tag;
};

struct SGAFormat
{
  static void read(std::ifstream& in, BlockArray& data, sdsl::int_vector<64>& counts);
  static void write(std::ofstream& out, const BlockArray& data, const NativeHeader& info);
  inline static AlphabeticOrder order() { return AO_DEFAULT; }

  const static std::string name;
  const static std::string tag;
};

//------------------------------------------------------------------------------

bool formatExists(const std::string& format);

void printFormats(std::ostream& stream);

template<class Format>
void
printFormat(std::ostream& stream)
{
  std::string padding;
  if(Format::tag.length() < 15) { padding = std::string(15 - Format::tag.length(), ' '); }
  stream << "  " << Format::tag << padding << Format::name << std::endl;
}

//------------------------------------------------------------------------------

struct RopeHeader
{
  uint32_t tag;

  const static uint32_t DEFAULT_TAG = 0x06454C52;
  const static size_type SIZE = 4;

  RopeHeader();

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);
  bool check() const;
};

std::ostream& operator<<(std::ostream& stream, const RopeHeader& header);

struct SGAHeader
{
  uint16_t tag;
  uint64_t sequences;
  uint64_t bases;
  uint64_t bytes;
  uint32_t flags;

  const static uint16_t DEFAULT_TAG = 0xCACA;
  const static uint32_t DEFAULT_FLAGS = 0;

  SGAHeader();

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);
  bool check() const;
};

std::ostream& operator<<(std::ostream& stream, const SGAHeader& header);

//------------------------------------------------------------------------------

} // namespace bwtmerge

#endif // _BWTMERGE_FORMATS_H
