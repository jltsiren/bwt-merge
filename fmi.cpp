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

#include <stack>

#include "fmi.h"

namespace bwtmerge
{

//------------------------------------------------------------------------------

FMI::FMI()
{
}

FMI::FMI(const FMI& source)
{
  this->copy(source);
}

FMI::FMI(FMI&& source)
{
  *this = std::move(source);
}

FMI::~FMI()
{
}

void
FMI::copy(const FMI& source)
{
  this->bwt = source.bwt;
  this->alpha = source.alpha;
}

void
FMI::swap(FMI& source)
{
  if(this != &source)
  {
    this->bwt.swap(source.bwt);
    this->alpha.swap(source.alpha);
  }
}

FMI&
FMI::operator=(const FMI& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

FMI&
FMI::operator=(FMI&& source)
{
  if(this != &source)
  {
    this->bwt = std::move(source.bwt);
    this->alpha = std::move(source.alpha);
  }
  return *this;
}

FMI::size_type
FMI::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->bwt.serialize(out, child, "bwt");
  written_bytes += this->alpha.serialize(out, child, "alpha");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
FMI::load(std::istream& in)
{
  this->bwt.load(in);
  this->alpha.load(in);
}

//------------------------------------------------------------------------------

template<>
void
FMI::serialize<NativeFormat>(const std::string& filename) const
{
  std::ofstream out(filename.c_str(), std::ios_base::binary);
  if(!out)
  {
    std::cerr << "FMI::serialize(): Cannot open output file " << filename << std::endl;
    return;
  }
  this->serialize(out);
  out.close();
}

template<>
void
FMI::load<NativeFormat>(const std::string& filename)
{
  std::ifstream in(filename.c_str(), std::ios_base::binary);
  if(!in)
  {
    std::cerr << "FMI::load(): Cannot open input file " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->load(in);
  in.close();
}

//------------------------------------------------------------------------------

struct MergeBuffer
{
  typedef RLArray<BlockArray> buffer_type;
  typedef RLArray<BlockArray>::run_type run_type;

  MergeParameters parameters;

  std::mutex               buffer_lock;
  std::vector<buffer_type> merge_buffers;

  std::mutex ra_lock;
  RankArray  ra;
  size_type  ra_values, ra_bytes;

  size_type  size;

  MergeBuffer(size_type _size, const MergeParameters& _parameters) :
    parameters(_parameters),
    merge_buffers(_parameters.merge_buffers),
    ra_values(0), ra_bytes(0), size(_size)
  {
  }

  ~MergeBuffer() {}

  void write(buffer_type& buffer)
  {
    if(buffer.empty()) { return; }

    std::string filename;
    size_type buffer_values = buffer.values(), buffer_bytes = buffer.bytes();
    {
      std::lock_guard<std::mutex> lock(this->ra_lock);
      filename = tempFile(this->parameters.tempPrefix());
      this->ra.filenames.push_back(filename);
      this->ra.run_counts.push_back(buffer.size());
      this->ra.value_counts.push_back(buffer.values());
    }
    buffer.write(filename); buffer.clear();

#ifdef VERBOSE_STATUS_INFO
    double ra_done, ra_gb;
#endif
    {
      std::lock_guard<std::mutex> lock(this->ra_lock);
      this->ra_values += buffer_values;
      this->ra_bytes += buffer_bytes + sizeof(size_type);
#ifdef VERBOSE_STATUS_INFO
      ra_done = (100.0 * this->ra_values) / this->size;
      ra_gb = inGigabytes(this->ra_bytes);
#endif
    }

#ifdef VERBOSE_STATUS_INFO
    {
      std::lock_guard<std::mutex> lock(Parallel::stderr_access);
      std::cerr << "buildRA(): Thread " << std::this_thread::get_id()
                << ": Added the values to the rank array" << std::endl;
      std::cerr << "buildRA(): " << ra_done << "% done; RA size " << ra_gb << " GB" << std::endl;
    }
#endif
  }

  void flush()
  {
    for(size_type i = 1; i < this->merge_buffers.size(); i++)
    {
      this->merge_buffers[i] = buffer_type(this->merge_buffers[i], this->merge_buffers[i - 1]);
    }
#ifdef VERBOSE_STATUS_INFO
    {
      std::lock_guard<std::mutex> lock(Parallel::stderr_access);
      std::cerr << "buildRA(): Flushing "
                << this->merge_buffers[this->merge_buffers.size() - 1].values()
                << " values to disk" << std::endl;
    }
#endif
    this->write(this->merge_buffers[this->merge_buffers.size() - 1]);
  }
};

void
mergeRA(MergeBuffer& mb, MergeBuffer::buffer_type& thread_buffer,
  std::vector<MergeBuffer::run_type>& run_buffer, bool force)
{
  MergeBuffer::buffer_type temp_buffer(run_buffer); run_buffer.clear();
  thread_buffer = MergeBuffer::buffer_type(thread_buffer, temp_buffer);
  if(!force && thread_buffer.bytes() < mb.parameters.thread_buffer_size) { return; }

#ifdef VERBOSE_STATUS_INFO
  {
    std::lock_guard<std::mutex> lock(Parallel::stderr_access);
    std::cerr << "buildRA(): Thread " << std::this_thread::get_id() << ": Adding "
              << thread_buffer.values() << " values to the merge buffer" << std::endl;
  }
#endif

  for(size_type i = 0; i < mb.merge_buffers.size(); i++)
  {
    bool done = false;
    {
      std::lock_guard<std::mutex> lock(mb.buffer_lock);
      if(mb.merge_buffers[i].empty()) { thread_buffer.swap(mb.merge_buffers[i]); done = true; }
      else { temp_buffer.swap(mb.merge_buffers[i]); }
    }
    if(done)
    {
#ifdef VERBOSE_STATUS_INFO
      std::lock_guard<std::mutex> lock(Parallel::stderr_access);
      std::cerr << "buildRA(): Thread " << std::this_thread::get_id()
                << ": Added the values to buffer " << i << std::endl;
#endif
      return;
    }
    thread_buffer = MergeBuffer::buffer_type(thread_buffer, temp_buffer);
  }

  mb.write(thread_buffer);
}

//------------------------------------------------------------------------------

struct MergePosition
{
  size_type  a_pos;
  range_type b_range;

  MergePosition() : a_pos(0), b_range(0, 0) {}
  MergePosition(size_type a, size_type b) : a_pos(a), b_range(b, b) {}
  MergePosition(size_type pos, range_type range) : a_pos(pos), b_range(range) {}
  MergePosition(size_type pos, size_type sp, size_type ep) : a_pos(pos), b_range(sp, ep) {}
};

void
buildRA(ParallelLoop& loop, const FMI& a, const FMI& b, MergeBuffer& mb)
{
  while(true)
  {
    range_type sequence_range = loop.next();
    if(Range::empty(sequence_range)) { return; }

    MergeBuffer::buffer_type thread_buffer;
    std::vector<MergeBuffer::run_type> run_buffer; run_buffer.reserve(mb.parameters.run_buffer_size);
    std::stack<MergePosition> positions;
    BWT::ranks_type a_pos, b_sp, b_ep;
    BWT::rank_ranges_type b_range;

    positions.push(MergePosition(a.sequences(), sequence_range));
    while(!(positions.empty()))
    {
      MergePosition curr = positions.top(); positions.pop();
      run_buffer.push_back(MergeBuffer::run_type(curr.a_pos, Range::length(curr.b_range)));
      if(run_buffer.size() >= mb.parameters.run_buffer_size)
      {
        mergeRA(mb, thread_buffer, run_buffer, false);
      }

      if(Range::length(curr.b_range) == 1)
      {
        range_type pred = b.LF(curr.b_range.first);
        if(pred.second != 0)
        {
          positions.push(MergePosition(a.LF(curr.a_pos, pred.second), pred.first));
        }
      }
      else if(Range::length(curr.b_range) <= FMI::SHORT_RANGE)
      {
        b.LF(curr.b_range, b_range);
        for(size_type c = 1; c < b.alpha.sigma; c++)
        {
          if(!(Range::empty(b_range[c])))
          {
            positions.push(MergePosition(a.LF(curr.a_pos, c), b_range[c]));
          }
        }
      }
      else
      {
        a.LF(curr.a_pos, a_pos); b.LF(curr.b_range, b_sp, b_ep);
        for(size_type c = 1; c < b.alpha.sigma; c++)
        {
          if(b_sp[c] <= b_ep[c]) { positions.push(MergePosition(a_pos[c], b_sp[c], b_ep[c])); }
        }
      }
    }

    mergeRA(mb, thread_buffer, run_buffer, true);
  #ifdef VERBOSE_STATUS_INFO
    {
      std::lock_guard<std::mutex> lock(Parallel::stderr_access);
      std::cerr << "buildRA(): Thread " << std::this_thread::get_id() << ": Finished block "
                << sequence_range << std::endl;
    }
  #endif
  }
}

FMI::FMI(FMI& a, FMI& b, MergeParameters parameters)
{
  if(a.alpha != b.alpha)
  {
    std::cerr << "FMI::FMI(): Cannot merge BWTs with different alphabets" << std::endl;
    std::exit(EXIT_FAILURE);
  }

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "bwt_merge: " << a.sequences() << " sequences of total length " << a.size() << std::endl;
  std::cerr << "bwt_merge: Adding " << b.sequences() << " sequences of total length " << b.size() << std::endl;
  std::cerr << "bwt_merge: Memory usage before merging: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  double start = readTimer();
#endif

  std::vector<range_type> bounds = getBounds(range_type(0, b.sequences() - 1), parameters.sequence_blocks);

  MergeBuffer mb(b.size(), parameters);
  {
    ParallelLoop loop(0, b.sequences(), parameters.sequence_blocks, parameters.threads);
    loop.execute(buildRA, std::ref(a), std::ref(b), std::ref(mb));
  }
  mb.flush();

#ifdef VERBOSE_STATUS_INFO
  double seconds = readTimer() - start;
  std::cerr << "bwt_merge: RA built in " << seconds << " seconds" << std::endl;
  std::cerr << "bwt_merge: Memory usage with RA: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
#endif

  this->bwt = BWT(a.bwt, b.bwt, mb.ra);
  this->alpha = a.alpha;
  for(size_type c = 0; c <= this->alpha.sigma; c++) { this->alpha.C[c] += b.alpha.C[c]; }
}

//------------------------------------------------------------------------------

void
serialize(const FMI& fmi, const std::string& filename, const std::string& format)
{
  if(format == NativeFormat::tag)
  {
    fmi.serialize<NativeFormat>(filename);
  }
  else if(format == PlainFormatD::tag)
  {
    fmi.serialize<PlainFormatD>(filename);
  }
  else if(format == PlainFormatS::tag)
  {
    fmi.serialize<PlainFormatS>(filename);
  }
  else if(format == RFMFormat::tag)
  {
    fmi.serialize<RFMFormat>(filename);
  }
  else if(format == SDSLFormat::tag)
  {
    fmi.serialize<SDSLFormat>(filename);
  }
  else if(format == RopeFormat::tag)
  {
    fmi.serialize<RopeFormat>(filename);
  }
  else if(format == SGAFormat::tag)
  {
    fmi.serialize<SGAFormat>(filename);
  }
  else
  {
    std::cerr << "serialize(): Invalid BWT format: " << format << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void
load(FMI& fmi, const std::string& filename, const std::string& format)
{
  if(format == NativeFormat::tag)
  {
    fmi.load<NativeFormat>(filename);
  }
  else if(format == PlainFormatD::tag)
  {
    fmi.load<PlainFormatD>(filename);
  }
  else if(format == PlainFormatS::tag)
  {
    fmi.load<PlainFormatS>(filename);
  }
  else if(format == RFMFormat::tag)
  {
    fmi.load<RFMFormat>(filename);
  }
  else if(format == SDSLFormat::tag)
  {
    fmi.load<SDSLFormat>(filename);
  }
  else if(format == RopeFormat::tag)
  {
    fmi.load<RopeFormat>(filename);
  }
  else if(format == SGAFormat::tag)
  {
    fmi.load<SGAFormat>(filename);
  }
  else
  {
    std::cerr << "load(): Invalid BWT format: " << format << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//------------------------------------------------------------------------------

const std::string MergeParameters::DEFAULT_TEMP_DIR = ".";
const std::string MergeParameters::TEMP_FILE_PREFIX = ".bwtmerge";

MergeParameters::MergeParameters() :
  run_buffer_size(RUN_BUFFER_SIZE), thread_buffer_size(THREAD_BUFFER_SIZE),
  merge_buffers(MERGE_BUFFERS),
  threads(Parallel::max_threads), sequence_blocks(threads * BLOCKS_PER_THREAD),
  temp_dir(DEFAULT_TEMP_DIR)
{
}

void
MergeParameters::sanitize()
{
  this->threads = Range::bound(this->threads, 1, Parallel::max_threads);
  this->sequence_blocks = std::max(this->sequence_blocks, (size_type)1);
  this->threads = std::min(this->threads, this->sequence_blocks);
}

void
MergeParameters::setTemp(const std::string& directory)
{
  if(directory.length() == 0) { this->temp_dir = DEFAULT_TEMP_DIR; }
  else if(directory[directory.length() - 1] != '/') { this->temp_dir = directory; }
  else { this->temp_dir = directory.substr(0, directory.length() - 1); }
}

std::string
MergeParameters::tempPrefix() const
{
  return this->temp_dir + '/' + TEMP_FILE_PREFIX;
}

std::ostream&
operator<< (std::ostream& stream, const MergeParameters& parameters)
{
  stream << "Run buffers:      "
    << inMegabytes(parameters.run_buffer_size * sizeof(MergeParameters::run_type)) << " MB" << std::endl;
  stream << "Thread buffers:   " << inMegabytes(parameters.thread_buffer_size) << " MB" << std::endl;
  stream << "Merge buffers:    " << parameters.merge_buffers << std::endl;
  stream << "Threads:          " << parameters.threads << std::endl;
  stream << "Sequence blocks:  " << parameters.sequence_blocks << std::endl;
  stream << "Temp directory:   " << parameters.temp_dir << std::endl;
  return stream;
}

//------------------------------------------------------------------------------

} // namespace bwtmerge
