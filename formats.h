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

/*
  BWT file formats. Note that BWT formats are only compatible with those having the
  same alphabetic order.

    load()      reads the BWT from 'in' and stores it in the native format in 'data'
    write()     writes the BWT stored in the native format in 'data' to 'out'
    alphabet()  returns an Alphabet object compatible with the format
    order()     returns the alphabetic order

    RFMFormat - Relative FM-index (uncompressed BWT); incompatible
    SGAFormat - SGA assembler

  FIXME: Add PlainFormat, SDSLFormat, NativeFormat
*/

enum AlphabeticOrder { AO_DEFAULT, AO_SORTED };

struct RFMFormat
{
  static void read(std::istream& in, BlockArray& data);
  static void write(std::ostream& out, const BlockArray& data);
  static Alphabet alphabet();
  inline static AlphabeticOrder order() { return AO_SORTED; }
};

struct SGAFormat
{
  static void read(std::istream& in, BlockArray& data);
  static void write(std::ostream& out, const BlockArray& data);
  static Alphabet alphabet();
  inline static AlphabeticOrder order() { return AO_DEFAULT; }
};

//------------------------------------------------------------------------------

} // namespace bwtmerge

#endif // _BWTMERGE_FORMATS_H
