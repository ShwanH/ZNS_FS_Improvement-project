// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
// Copyright (c) 2019-present, Western Digital Corporation
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).

#if !defined(ROCKSDB_LITE) && !defined(OS_WIN) && defined(LIBZBD)

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <libzbd/zbd.h>
#include <linux/blkzoned.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "io_zenfs.h"
#include "rocksdb/env.h"
#include "util/coding.h"
#include "zbd_zenfs.h"

namespace ROCKSDB_NAMESPACE {

Status ZoneExtent::DecodeFrom(Slice* input) {
  if (input->size() != (sizeof(start_) + sizeof(length_)))
    return Status::Corruption("ZoneExtent", "Error: length missmatch");

  GetFixed64(input, &start_);
  GetFixed32(input, &length_);
  return Status::OK();
}

void ZoneExtent::EncodeTo(std::string* output) {
  PutFixed64(output, start_);
  PutFixed32(output, length_);
}

ZoneExtent::~ZoneExtent() {
  if (zone_ != nullptr) {
    zone_->RemoveZoneFile(file_);
  }
}

enum ZoneFileTag : uint32_t {
  kFileID = 1,
  kFileName = 2,
  kFileSize = 3,
  kWriteLifeTimeHint = 4,
  kExtent = 5,
};

void ZoneFile::EncodeTo(std::string* output, uint32_t extent_start) {
  PutFixed32(output, kFileID);
  PutFixed64(output, file_id_);

  PutFixed32(output, kFileName);
  PutLengthPrefixedSlice(output, Slice(filename_));

  PutFixed32(output, kFileSize);
  PutFixed64(output, fileSize);

  for (uint32_t i = extent_start; i < extents_.size(); i++) {
    std::string extent_str;

    PutFixed32(output, kExtent);
    extents_[i]->EncodeTo(&extent_str);
    PutLengthPrefixedSlice(output, Slice(extent_str));
  }

  /* We're not encoding active zone and extent start
   * as files will always be read-only after mount */
}

Status ZoneFile::DecodeFrom(Slice* input) {
  uint32_t tag = 0;

  GetFixed32(input, &tag);
  if (tag != kFileID || !GetFixed64(input, &file_id_))
    return Status::Corruption("ZoneFile", "File ID missing");

  while (true) {
    Slice slice;
    ZoneExtent* extent;
    Status s;

    if (!GetFixed32(input, &tag)) break;

    switch (tag) {
      case kFileName:
        if (!GetLengthPrefixedSlice(input, &slice))
          return Status::Corruption("ZoneFile", "Filename missing");
        filename_ = slice.ToString();
        if (filename_.length() == 0)
          return Status::Corruption("ZoneFile", "Zero length filename");
        break;
      case kFileSize:
        if (!GetFixed64(input, &fileSize))
          return Status::Corruption("ZoneFile", "Missing file size");
        break;
      case kExtent:
        extent = new ZoneExtent(0, 0, nullptr, this);
        GetLengthPrefixedSlice(input, &slice);
        s = extent->DecodeFrom(&slice);
        if (!s.ok()) {
          delete extent;
          return s;
        }
        extent->zone_ = zbd_->GetIOZone(extent->start_);
        if (!extent->zone_)
          return Status::Corruption("ZoneFile", "Invalid zone extent");
        extent->zone_->used_capacity_ += extent->length_;
        extents_.push_back(extent);
        break;
      default:
        return Status::Corruption("ZoneFile", "Unexpected tag");
    }
  }

  MetadataSynced();
  return Status::OK();
}

Status ZoneFile::MergeUpdate(ZoneFile* update) {
  if (file_id_ != update->GetID())
    return Status::Corruption("ZoneFile update", "ID missmatch");

  Rename(update->GetFilename());
  SetFileSize(update->GetFileSize());

  std::vector<ZoneExtent*> update_extents = update->GetExtents();
  for (long unsigned int i = 0; i < update_extents.size(); i++) {
    ZoneExtent* extent = update_extents[i];
    Zone* zone = extent->zone_;
    zone->used_capacity_ += extent->length_;
    extents_.push_back(
        new ZoneExtent(extent->start_, extent->length_, zone, this));
  }

  MetadataSynced();

  return Status::OK();
}

ZoneFile::ZoneFile(ZonedBlockDevice* zbd, std::string filename,
                   uint64_t file_id)
    : zbd_(zbd),
      active_zone_(NULL),
      extent_start_(0),
      extent_filepos_(0),
      lifetime_(Env::WLTH_NOT_SET),
      fileSize(0),
      filename_(filename),
      file_id_(file_id),
      nr_synced_extents_(0) {}

std::string ZoneFile::GetFilename() { return filename_; }

void ZoneFile::Rename(std::string name) { filename_ = name; }

uint64_t ZoneFile::GetFileSize() { return fileSize; }
void ZoneFile::SetFileSize(uint64_t sz) { fileSize = sz; }

ZoneFile::~ZoneFile() {
  for (auto e = std::begin(extents_); e != std::end(extents_); ++e) {
    Zone* zone = (*e)->zone_;

    assert(zone && zone->used_capacity_ >= (*e)->length_);
    zone->used_capacity_ -= (*e)->length_;
    delete *e;
  }
  CloseWR();
}

void ZoneFile::CloseWR() {
  if (active_zone_) {
    active_zone_->CloseWR();
    active_zone_ = NULL;
  }
}

uint64_t ZoneFile::GetOffset(ZoneExtent* target_extent) {
  uint64_t file_offset = 0;
  for (const auto extent : extents_) {
    if (extent->start_ == target_extent->start_) {
      return file_offset;
    }
    file_offset += extent->length_;
  }
  return std::numeric_limits<uint64_t>::max();
}

ZoneExtent* ZoneFile::GetExtent(uint64_t file_offset, uint64_t* dev_offset) {
  for (unsigned int i = 0; i < extents_.size(); i++) {
    if (file_offset < extents_[i]->length_) {
      *dev_offset = extents_[i]->start_ + file_offset;
      return extents_[i];
    } else {
      file_offset -= extents_[i]->length_;
    }
  }
  return NULL;
}

ZoneExtent* ZoneFile::GetExtent(uint64_t extent_start) {
  for (const auto extent : extents_) {
    if (extent->start_ == extent_start) {
      return extent;
    }
  }
  return nullptr;
}

IOStatus ZoneFile::PositionedRead(uint64_t offset, size_t n, Slice* result,
                                  char* scratch, bool direct) {
  int f = zbd_->GetReadFD();
  int f_direct = zbd_->GetReadDirectFD();
  char* ptr;
  uint64_t r_off;
  size_t r_sz;
  ssize_t r = 0;
  size_t read = 0;
  ZoneExtent* extent;
  uint64_t extent_end;
  IOStatus s;

  if (offset >= fileSize) {
    *result = Slice(scratch, 0);
    return IOStatus::OK();
  }

  r_off = 0;
  extent = GetExtent(offset, &r_off);
  if (!extent) {
    /* read start beyond end of (synced) file data*/
    *result = Slice(scratch, 0);
    return s;
  }
  extent_end = extent->start_ + extent->length_;

  /* Limit read size to end of file */
  if ((offset + n) > fileSize)
    r_sz = fileSize - offset;
  else
    r_sz = n;

  ptr = scratch;

  while (read != r_sz) {
    size_t pread_sz = r_sz - read;

    if ((pread_sz + r_off) > extent_end) pread_sz = extent_end - r_off;

    if (direct) {
      assert((uint64_t)ptr % GetBlockSize() == 0);
      assert(pread_sz % GetBlockSize() == 0);
      assert(r_off % GetBlockSize() == 0);
      r = pread(f_direct, ptr, pread_sz, r_off);
    } else {
      r = pread(f, ptr, pread_sz, r_off);
    }

    if (r <= 0) {
      if (r == -1 && errno == EINTR) {
        continue;
      }
      break;
    }

    pread_sz = (size_t)r;

    ptr += pread_sz;
    read += pread_sz;
    r_off += pread_sz;

    if (read != r_sz && r_off == extent_end) {
      extent = GetExtent(offset + read, &r_off);
      if (!extent) {
        /* read beyond end of (synced) file data */
        break;
      }
      r_off = extent->start_;
      extent_end = extent->start_ + extent->length_;
      assert(((size_t)r_off % zbd_->GetBlockSize()) == 0);
    }
  }

  if (r < 0) {
    s = IOStatus::IOError("pread error\n");
    read = 0;
  }

  *result = Slice((char*)scratch, read);
  return s;
}

ZoneExtent* ZoneFile::ReplaceExtent(ZoneExtent* target, ZoneExtent* top) {
  assert(nullptr != target);
  assert(nullptr != top);

  int64_t target_length = (int64_t)target->length_;

  std::vector<ZoneExtent*>::iterator target_pos;
  std::vector<ZoneExtent*> new_extents;
  ZoneExtent* gc_target = nullptr;
  target_pos = find(extents_.begin(), extents_.end(), target);
  assert(extents_.end() != target_pos);

  for (auto rit = extents_.rbegin(); rit != extents_.rend(); rit++) {
    ZoneExtent* extent = *rit;
    assert(nullptr != extent);
    if (extent == top) {
      break;
    }
    new_extents.insert(new_extents.begin(), extent);
  }

  for (auto it = new_extents.begin(); it != new_extents.end(); it++) {
    ZoneExtent* extent = *it;
    assert(nullptr != extent);
    extents_.insert(target_pos + 1, extent);
    extents_.pop_back();
    fileSize -= extent->length_;

    if (target_length >= extent->length_) {
      target_length -= extent->length_;
    } else {
      extent->length_ = target_length;
      target_length = 0;
    }
  }
  if (target_length != 0) {
    fprintf(stderr, "target error: %ld\n", target_length);
  }
  assert(target_length == 0);

  target_pos = find(extents_.begin(), extents_.end(), target);
  assert(target_pos != extents_.end());
  gc_target = *target_pos;
  extents_.erase(target_pos);
  return gc_target;
}

/* You must call this method after the ReplaceExtent() method calls */
void ZoneFile::RestoreExtent(Zone* active_zone, uint64_t /* extent_start */,
                             uint64_t extent_filepos) {
  assert(active_zone == nullptr);
  extent_start_ = extents_.back()->start_;
  extent_filepos_ = extent_filepos;
}

void ZoneFile::PushExtent() {
  uint64_t length;

  assert(fileSize >= extent_filepos_);

  if (!active_zone_) return;

  length = fileSize - extent_filepos_;
  if (length == 0) return;

  assert(length <= (active_zone_->wp_ - extent_start_));
  extents_.push_back(new ZoneExtent(extent_start_, length, active_zone_, this));

  active_zone_->used_capacity_ += length;
  extent_start_ = active_zone_->wp_;
  extent_filepos_ = fileSize;
}

/* Assumes that data and size are block aligned */
IOStatus ZoneFile::Append(void* data, int data_size, int valid_size,
                          bool is_gc = true) {
  uint32_t left = data_size;
  uint32_t wr_size, offset = 0;
  IOStatus s;

  if (active_zone_ == NULL) {
    active_zone_ = zbd_->AllocateZone(lifetime_, this, NULL, is_gc);
    if (!active_zone_) {
      return IOStatus::NoSpace("Zone allocation failure\n");
    }
    extent_start_ = active_zone_->wp_;
    extent_filepos_ = fileSize;
  }

  while (left) {
    if (active_zone_->capacity_ == 0) {
      PushExtent();

      active_zone_->CloseWR();
      active_zone_ = zbd_->AllocateZone(lifetime_, this, active_zone_, is_gc);
      if (!active_zone_) {
        return IOStatus::NoSpace("Zone allocation failure\n");
      }
      extent_start_ = active_zone_->wp_;
      extent_filepos_ = fileSize;
    }

    wr_size = left;
    if (wr_size > active_zone_->capacity_) wr_size = active_zone_->capacity_;

    fprintf(zbd_->GetZoneLogFile(), "%-10ld%-8s%-8lu%-8lu%-45s%-10u%-10lu\n",
            (long int)((double)clock() / CLOCKS_PER_SEC * 1000), "WRITE",
            (unsigned long)0, active_zone_->GetZoneNr(), filename_.c_str(),
            wr_size, fileSize);
    s = active_zone_->Append((char*)data + offset, wr_size);
    if (!s.ok()) {
      return s;
    }

    fileSize += wr_size;
    left -= wr_size;
    offset += wr_size;
  }

  fileSize -= (data_size - valid_size);
  return IOStatus::OK();
}

IOStatus ZoneFile::SetWriteLifeTimeHint(Env::WriteLifeTimeHint lifetime) {
  lifetime_ = lifetime;
  return IOStatus::OK();
}

ZonedWritableFile::ZonedWritableFile(ZonedBlockDevice* zbd, bool _buffered,
                                     ZoneFile* zoneFile,
                                     MetadataWriter* metadata_writer) {
  wp = zoneFile->GetFileSize();
  assert(wp == 0);

  buffered = _buffered;
  block_sz = zbd->GetBlockSize();
  buffer_sz = block_sz * 256;
  buffer_pos = 0;

  zoneFile_ = zoneFile;

  if (buffered) {
    int ret = posix_memalign((void**)&buffer, block_sz, buffer_sz);

    if (ret) buffer = nullptr;

    assert(buffer != nullptr);
  }

  metadata_writer_ = metadata_writer;
}

ZonedWritableFile::~ZonedWritableFile() {
  zoneFile_->CloseWR();
  if (buffered) free(buffer);
};

ZonedWritableFile::MetadataWriter::~MetadataWriter() {}

IOStatus ZonedWritableFile::Truncate(uint64_t size,
                                     const IOOptions& /*options*/,
                                     IODebugContext* /*dbg*/) {
  zoneFile_->SetFileSize(size);
  return IOStatus::OK();
}

IOStatus ZonedWritableFile::Fsync(const IOOptions& /*options*/,
                                  IODebugContext* /*dbg*/) {
  IOStatus s;

  buffer_mtx_.lock();
  s = FlushBuffer();
  buffer_mtx_.unlock();
  if (!s.ok()) {
    return s;
  }
  zoneFile_->PushExtent();

  return metadata_writer_->Persist(zoneFile_);
}

IOStatus ZonedWritableFile::Sync(const IOOptions& options,
                                 IODebugContext* dbg) {
  return Fsync(options, dbg);
}

IOStatus ZonedWritableFile::Flush(const IOOptions& /*options*/,
                                  IODebugContext* /*dbg*/) {
  return IOStatus::OK();
}

IOStatus ZonedWritableFile::RangeSync(uint64_t offset, uint64_t nbytes,
                                      const IOOptions& options,
                                      IODebugContext* dbg) {
  if (wp < offset + nbytes) return Fsync(options, dbg);

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::Close(const IOOptions& options,
                                  IODebugContext* dbg) {
  Fsync(options, dbg);
  zoneFile_->CloseWR();

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::FlushBuffer() {
  uint32_t align, pad_sz = 0, wr_sz;
  IOStatus s;

  if (!buffer_pos) return IOStatus::OK();

  align = buffer_pos % block_sz;
  if (align) pad_sz = block_sz - align;

  if (pad_sz) memset((char*)buffer + buffer_pos, 0x0, pad_sz);

  wr_sz = buffer_pos + pad_sz;
  s = zoneFile_->Append((char*)buffer, wr_sz, buffer_pos);
  if (!s.ok()) {
    return s;
  }

  wp += buffer_pos;
  buffer_pos = 0;

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::BufferedWrite(const Slice& slice) {
  uint32_t buffer_left = buffer_sz - buffer_pos;
  uint32_t data_left = slice.size();
  char* data = (char*)slice.data();
  uint32_t tobuffer;
  int blocks, aligned_sz;
  int ret;
  void* alignbuf;
  IOStatus s;

  if (buffer_pos || data_left <= buffer_left) {
    if (data_left < buffer_left) {
      tobuffer = data_left;
    } else {
      tobuffer = buffer_left;
    }

    memcpy(buffer + buffer_pos, data, tobuffer);
    buffer_pos += tobuffer;
    data_left -= tobuffer;

    if (!data_left) return IOStatus::OK();

    data += tobuffer;
  }

  if (buffer_pos == buffer_sz) {
    s = FlushBuffer();
    if (!s.ok()) return s;
  }

  if (data_left >= buffer_sz) {
    blocks = data_left / block_sz;
    aligned_sz = block_sz * blocks;

    ret = posix_memalign(&alignbuf, block_sz, aligned_sz);
    if (ret) {
      return IOStatus::IOError("failed allocating alignment write buffer\n");
    }

    memcpy(alignbuf, data, aligned_sz);
    s = zoneFile_->Append(alignbuf, aligned_sz, aligned_sz);
    free(alignbuf);

    if (!s.ok()) return s;

    wp += aligned_sz;
    data_left -= aligned_sz;
    data += aligned_sz;
  }

  if (data_left) {
    memcpy(buffer, data, data_left);
    buffer_pos = data_left;
  }

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::Append(const Slice& data,
                                   const IOOptions& /*options*/,
                                   IODebugContext* /*dbg*/) {
  IOStatus s;

  if (buffered) {
    buffer_mtx_.lock();
    s = BufferedWrite(data);
    buffer_mtx_.unlock();
  } else {
    s = zoneFile_->Append((void*)data.data(), data.size(), data.size());
    if (s.ok()) wp += data.size();
  }

  return s;
}

IOStatus ZonedWritableFile::PositionedAppend(const Slice& data, uint64_t offset,
                                             const IOOptions& /*options*/,
                                             IODebugContext* /*dbg*/) {
  IOStatus s;

  if (offset != wp) {
    assert(false);
    return IOStatus::IOError("positioned append not at write pointer");
  }

  if (buffered) {
    buffer_mtx_.lock();
    s = BufferedWrite(data);
    buffer_mtx_.unlock();
  } else {
    s = zoneFile_->Append((void*)data.data(), data.size(), data.size());
    if (s.ok()) wp += data.size();
  }

  return s;
}

void ZonedWritableFile::SetWriteLifeTimeHint(Env::WriteLifeTimeHint hint) {
  zoneFile_->SetWriteLifeTimeHint(hint);
}

IOStatus ZonedSequentialFile::Read(size_t n, const IOOptions& /*options*/,
                                   Slice* result, char* scratch,
                                   IODebugContext* /*dbg*/) {
  IOStatus s;

  s = zoneFile_->PositionedRead(rp, n, result, scratch, direct_);
  if (s.ok()) rp += result->size();

  return s;
}

IOStatus ZonedSequentialFile::Skip(uint64_t n) {
  if (rp + n >= zoneFile_->GetFileSize())
    return IOStatus::InvalidArgument("Skip beyond end of file");
  rp += n;
  return IOStatus::OK();
}

IOStatus ZonedSequentialFile::PositionedRead(uint64_t offset, size_t n,
                                             const IOOptions& /*options*/,
                                             Slice* result, char* scratch,
                                             IODebugContext* /*dbg*/) {
  return zoneFile_->PositionedRead(offset, n, result, scratch, direct_);
}

IOStatus ZonedRandomAccessFile::Read(uint64_t offset, size_t n,
                                     const IOOptions& /*options*/,
                                     Slice* result, char* scratch,
                                     IODebugContext* /*dbg*/) const {
  return zoneFile_->PositionedRead(offset, n, result, scratch, direct_);
}

size_t ZoneFile::GetUniqueId(char* id, size_t max_size) {
  /* Based on the posix fs implementation */
  if (max_size < kMaxVarint64Length * 3) {
    return 0;
  }

  struct stat buf;
  int fd = zbd_->GetReadFD();
  int result = fstat(fd, &buf);
  if (result == -1) {
    return 0;
  }

  char* rid = id;
  rid = EncodeVarint64(rid, buf.st_dev);
  rid = EncodeVarint64(rid, buf.st_ino);
  rid = EncodeVarint64(rid, file_id_);
  assert(rid >= id);
  return static_cast<size_t>(rid - id);

  return 0;
}

size_t ZonedRandomAccessFile::GetUniqueId(char* id, size_t max_size) const {
  return zoneFile_->GetUniqueId(id, max_size);
}

}  // namespace ROCKSDB_NAMESPACE

#endif  // !defined(ROCKSDB_LITE) && !defined(OS_WIN) && defined(LIBZBD)
