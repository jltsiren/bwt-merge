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
FMI::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
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

struct MergeBuffer
{
  typedef RLArray<BlockArray> buffer_type;
  typedef RLArray<BlockArray>::run_type run_type;

  std::vector<omp_lock_t>  buffer_locks;
  std::vector<buffer_type> merge_buffers;

  const static std::string TEMP_NAME;

  omp_lock_t ra_lock;
  RankArray  ra;
  size_type  ra_values, ra_bytes;

  size_type  size;

  MergeBuffer(size_type _size, size_type buffers = FMI::MERGE_BUFFERS) :
    buffer_locks(buffers), merge_buffers(buffers), ra_values(0), ra_bytes(0), size(_size)
  {
    for(size_type i = 0; i < buffer_locks.size(); i++) { omp_init_lock(&(this->buffer_locks[i])); }
    omp_init_lock(&(this->ra_lock));
  }

  ~MergeBuffer()
  {
    for(size_type i = 0; i < buffer_locks.size(); i++) { omp_destroy_lock(&(this->buffer_locks[i])); }
    omp_destroy_lock(&(this->ra_lock));
  }

  void write(buffer_type& buffer)
  {
    if(buffer.empty()) { return; }

    omp_set_lock(&(this->ra_lock));
    std::string filename = tempFile(TEMP_NAME);
    this->ra.filenames.push_back(filename);
    this->ra.run_counts.push_back(buffer.size());
    this->ra.value_counts.push_back(buffer.values());
    this->ra_values += buffer.values();
    this->ra_bytes += buffer.bytes() + sizeof(size_type);
#ifdef VERBOSE_STATUS_INFO
    double ra_done = (100.0 * this->ra_values) / this->size;
    double ra_gb = inGigabytes(this->ra_bytes);
#endif
    omp_unset_lock(&(this->ra_lock));

    buffer.write(filename); buffer.clear();
#ifdef VERBOSE_STATUS_INFO
    #pragma omp critical
    {
      std::cerr << "buildRA(): Thread " << omp_get_thread_num()
                << ": Added the values to the rank array." << std::endl;
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
    #pragma omp critical
    {
      std::cerr << "buildRA(): Thread " << omp_get_thread_num() << ": Flushing "
                << this->merge_buffers[this->merge_buffers.size() - 1].values()
                << " values to disk." << std::endl;
    }
    this->write(this->merge_buffers[this->merge_buffers.size() - 1]);
  }
};

const std::string MergeBuffer::TEMP_NAME = ".bwtmerge";

void
mergeRA(MergeBuffer& mb, MergeBuffer::buffer_type& thread_buffer,
  std::vector<MergeBuffer::run_type>& buffer, bool force)
{
  MergeBuffer::buffer_type temp_buffer(buffer); buffer.clear();
  thread_buffer = MergeBuffer::buffer_type(thread_buffer, temp_buffer);
  if(!force && thread_buffer.bytes() < FMI::THREAD_BUFFER_SIZE) { return; }

#ifdef VERBOSE_STATUS_INFO
  #pragma omp critical
  {
    std::cerr << "buildRA(): Thread " << omp_get_thread_num() << ": Adding "
              << thread_buffer.values() << " values to the merge buffer." << std::endl;
  }
#endif

  for(size_type i = 0; i < mb.merge_buffers.size(); i++)
  {
    bool done = false;
    omp_set_lock(&(mb.buffer_locks[i]));
    if(mb.merge_buffers[i].empty()) { thread_buffer.swap(mb.merge_buffers[i]); done = true; }
    else { temp_buffer.swap(mb.merge_buffers[i]); }
    omp_unset_lock(&(mb.buffer_locks[i]));
    if(done)
    {
#ifdef VERBOSE_STATUS_INFO
      #pragma omp critical
      {
        std::cerr << "buildRA(): Thread " << omp_get_thread_num()
                  << ": Added the values to buffer " << i << "." << std::endl;
      }
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
};

void
buildRA(const FMI& a, const FMI& b, MergeBuffer& mb, range_type sequence_range)
{
  if(Range::empty(sequence_range)) { return; }

  MergeBuffer::buffer_type thread_buffer;
  std::vector<MergeBuffer::run_type> buffer;
  std::stack<MergePosition> positions;

  positions.push(MergePosition(a.sequences(), sequence_range));
  while(!(positions.empty()))
  {
    MergePosition curr = positions.top(); positions.pop();
    buffer.push_back(MergeBuffer::run_type(curr.a_pos, Range::length(curr.b_range)));
    if(buffer.size() >= FMI::RUN_BUFFER_SIZE) { mergeRA(mb, thread_buffer, buffer, false); }

    if(Range::length(curr.b_range) < b.alpha.sigma)
    {
      for(size_type i = curr.b_range.first; i <= curr.b_range.second; i++)
      {
        range_type pred = b.LF(i);
        if(pred.second != 0)
        {
          positions.push(MergePosition(a.LF(curr.a_pos, pred.second), pred.first));
        }
      }
    }
    else
    {
      for(comp_type c = 1; c < b.alpha.sigma; c++)
      {
        range_type prev = b.LF(curr.b_range, c);
        if(!(Range::empty(prev))) { positions.push(MergePosition(a.LF(curr.a_pos, c), prev)); }
      }
    }
  }

  mergeRA(mb, thread_buffer, buffer, true);
#ifdef VERBOSE_STATUS_INFO
  #pragma omp critical
  {
    std::cerr << "buildRA(): Thread " << omp_get_thread_num() << ": Finished block "
              << sequence_range << "." << std::endl;
  }
#endif
}

FMI::FMI(FMI& a, FMI& b)
{
  if(a.alpha != b.alpha)
  {
    std::cerr << "FMI::FMI(): Cannot merge BWTs with different alphabets." << std::endl;
    std::exit(EXIT_FAILURE);
  }

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "bwt_merge: " << a.sequences() << " sequences of total length " << a.size() << std::endl;
  std::cerr << "bwt_merge: Adding " << b.sequences() << " sequences of total length " << b.size() << std::endl;
  double start = readTimer();
#endif

  size_type threads = omp_get_max_threads();
  std::vector<range_type> bounds = getBounds(range_type(0, b.sequences() - 1), 4 * threads);

  MergeBuffer mb(b.size());
  #pragma omp parallel for schedule(dynamic, 1)
  for(size_type block = 0; block < bounds.size(); block++)
  {
    buildRA(a, b, mb, bounds[block]);
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

} // namespace bwtmerge
