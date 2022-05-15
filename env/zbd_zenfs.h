// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
// Copyright (c) 2019-present, Western Digital Corporation
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).

#pragma once

#if !defined(ROCKSDB_LITE) && defined(OS_LINUX) && defined(LIBZBD)

#include <errno.h>
#include <libzbd/zbd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <bitset>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <limits>
#include <mutex>
#include <queue>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#include <map>
#include "rocksdb/env.h"
#include "rocksdb/io_status.h"
#include "cpu_check.h"
#define ZONE_CUSTOM_DEBUG

// #define ORIGINAL
// #define ORIGINAL_GC
// #define HOT_COLD
#define HOT_COLD_GC

//#define DEPENDENT_GC_THREAD
#define INDEPENDENT_GC_THREAD

#define ORIGINAL_IO
//#define EQUAL_IO

#ifdef ORIGINAL
#define ZONE_NO_GC_THREAD
#define ZONE_GC_ENABLE (0)  // is gc enable
#endif

#ifdef ORIGINAL_GC
#define ZONE_USE_RESET_RATE_LIMITER
#define ZONE_USE_RESET_SIGMOID_LIMITER
#define ZONE_GC_ENABLE (1)  // is gc enable
#endif

#ifdef HOT_COLD
#define ZONE_HOT_COLD_SEP
#define ZONE_NO_GC_THREAD
#define ZONE_GC_ENABLE (0)  // is gc disable
#endif

#ifdef HOT_COLD_GC
#define ZONE_HOT_COLD_SEP
#define ZONE_USE_RESET_RATE_LIMITER
#define ZONE_USE_RESET_SIGMOID_LIMITER
#define ZONE_GC_ENABLE (1)  // is gc enable
#endif

#define ZONE_RESET_TRIGGER (100)  // under 100% of free space, RESET started
#define ZONE_GC_WATERMARK (25)    // under 25% of free space, GC started

#define ZONE_GC_THREAD_TICK (std::chrono::milliseconds(1000))

#define ZONE_MAX_NOTIFY_RETRY (10)
#define CPU_CHECK_THREAD_TICK (std::chrono::seconds(5))
//#define ZONE_TIME_CHECK

#if defined(ZONE_CUSTOM_DEBUG)
#pragma message("ZONE CUSTOM DEBUG mode enabled")
#else
#pragma message("ZONE CUSTOM DEBUG mode disabled")
#endif

#if ZONE_GC_ENABLE == 1
#pragma message("ZONE GC mode enabled")
#else
#pragma message("ZONE GC mode disabled")
#endif

#if defined(ZONE_HOT_COLD_SEP)
#pragma message("ZONE_HOT_COLD_SEP mode enabled")
#else
#pragma message("ZONE_HOT_COLD_SEP mode disabled")
#endif

#if defined(ZONE_MIX)
#pragma message("ZONE_MIX mode enabled")
#else
#pragma message("ZONE_MIX mode disabled")
#endif

#if defined(ZONE_USE_RESET_RATE_LIMITER)
#pragma message("fixed reset rate limiter mode enabled")
#elif defined(ZONE_USE_RESET_SIGMOID_LIMITER)
#pragma message("sigmoid reset rate limiter mode enabled")
#else
#pragma message("no reset rate limiter mode enabled")
#endif

namespace ROCKSDB_NAMESPACE {

enum class ZoneGcState { NOT_GC_TARGET, DO_RESET, NORMAL_EXIT };

enum ZoneInvalidLevel {
  ZONE_INVALID_LOW = 0,
  ZONE_INVALID_MIDDLE,
  ZONE_INVALID_HIGH,
  ZONE_INVALID_MAX,
  NR_ZONE_INVALID_LEVEL,
};


class SubZonedBlockDevice; //physical device
class ZonedBlockDevice; //logical device on all physical devices
class ZoneFile;
class ZoneExtent;

#define ZONE_EXTENT_FIND_FAIL (std::numeric_limits<uint64_t>::max())

class ZoneMapEntry {
 public:
  ZoneFile *file_;

  uint64_t extent_start_;
  uint32_t extent_length_;

  std::string filename_;
  std::atomic<bool> is_invalid_;

  Env::WriteLifeTimeHint lifetime_;

  explicit ZoneMapEntry(ZoneFile *file, ZoneExtent *extent);
};

class Zone {
  SubZonedBlockDevice *s_zbd_;

 public:
  explicit Zone(SubZonedBlockDevice *s_zbd, struct zbd_zone *z);

  uint64_t start_;
  uint64_t capacity_; /* remaining capacity */
  uint64_t max_capacity_;
  uint64_t wp_;
  bool open_for_write_;
  Env::WriteLifeTimeHint lifetime_;
  double total_lifetime_;
  std::atomic<long> used_capacity_;
  std::atomic<uint64_t> reset_counter_;
  std::vector<ZoneMapEntry *> file_map_;
#ifdef ZONE_CUSTOM_DEBUG
  std::atomic<uint64_t> file_map_size_;
#endif
  bool has_meta_;

  IOStatus Reset();
  IOStatus Finish();
  IOStatus Close();

  IOStatus Append(char *data, uint32_t size);
  bool IsUsed();
  bool IsFull();
  bool IsEmpty();
  uint64_t GetZoneNr();
  uint64_t GetCapacityLeft();
  double GetInvalidPercentage();
  double GetZoneValidAvgLevel();
  double GetZoneAvgLevel();

  void SetZoneFile(ZoneFile *file, ZoneExtent *extent);
  void RemoveZoneFile(ZoneFile *file, ZoneExtent *extent);
  void PrintZoneFiles(FILE *fp);

  void CloseWR(); /* Done writing */
  SubZonedBlockDevice* GetSubZBD(){return s_zbd_;}
};

class VictimZoneCompare {
 public:
  bool operator()(Zone *z1, Zone *z2) {
    // true means z2 goes to front
    if (z1->GetZoneValidAvgLevel() == z2->GetZoneValidAvgLevel()) {
      return z1->GetInvalidPercentage() < z2->GetInvalidPercentage();
    }
    return z1->GetZoneValidAvgLevel() < z2->GetZoneValidAvgLevel();
  }
};

class EmptyZoneCompare {
 public:
  bool operator()(Zone *z1, Zone *z2) {
    // true means z2 goes to front
    return z1->reset_counter_ > z2->reset_counter_;
  }
};

class SubZonedBlockDevice{
  private:
    std::string bdevname_;
    std::string filename_;
    uint32_t block_sz_;
    uint32_t zone_sz_;
    uint32_t nr_zones_;
    std::vector<Zone *> io_zones;
    std::recursive_mutex io_zones_mtx;
    std::vector<Zone *> meta_zones;
    int read_f_;
    int read_direct_f_;
    int write_f_;
    time_t start_time_;
    std::shared_ptr<Logger> logger_;
    uint32_t finish_threshold_ = 0;
    long int gc_time_sum = 0;

    std::atomic<long> active_io_zones_;
    std::atomic<long> open_io_zones_;
    std::condition_variable zone_resources_;
    std::mutex zone_resources_mtx_; // Protects active/open io zones 

    std::atomic<bool> gc_force_;
    std::atomic<bool> gc_thread_stop_;
    std::atomic<long> active_gc_;
    std::condition_variable gc_thread_cond_;
    std::condition_variable gc_force_cond_;
    std::mutex gc_force_mtx_;
    std::mutex gc_thread_mtx_;
    std::thread *gc_thread_;
    std::priority_queue<Zone *, std::vector<Zone *>, VictimZoneCompare>
      victim_queue_[NR_ZONE_INVALID_LEVEL];
    uint64_t gc_rate_limiter_;
    uint64_t reset_rate_limiter_;

    FILE *zone_log_file_;
    //gc_log, alloc_log
    FILE *gc_log_file_;
    FILE *alloc_log_file_;
    //gc_total count
    uint64_t gc_total = 0;
    //list for gc done zones num
    //std::vector<int> gc_done_Zones;
    //std::string gc_done_Zone;
    uint64_t count_alloc = 0;
    char *gc_buffer_;

    unsigned int max_nr_active_io_zones_;
    unsigned int max_nr_open_io_zones_;
    bool is_gc_running = false;

    Zone *AllocateZoneImpl(Env::WriteLifeTimeHint lifetime, ZoneFile *file,
                         int *is_empty);
  public:
    explicit SubZonedBlockDevice(std::string bdevname,
                              std::shared_ptr<Logger> logger);
    virtual ~SubZonedBlockDevice();

  #if defined(INDEPENDENT_GC_THREAD)
    std::mutex zonefile_mtx_;
  #else
    std::mutex *files_mtx_;
  #endif
    std::mutex *check_thread_mutex_;
    std::condition_variable *cpu_check_cond_;
    // below 2 member only for zone allocation
    std::vector<Zone *> occupied_zones_[Env::WLTH_EXTREME + 1];
    std::priority_queue<Zone *, std::vector<Zone *>, EmptyZoneCompare>
        empty_zones_queue_;

    IOStatus Open(bool readonly = false);

    Zone *GetIOZone(uint64_t offset);

    Zone *AllocateZone(Env::WriteLifeTimeHint lifetime, ZoneFile *zone_file,
                      Zone *before_zone);
    Zone *AllocateMetaZone();

    uint64_t GetFreeSpace();
    std::string GetFilename();
    std::string GetBdevname();
    uint32_t GetBlockSize();
    uint32_t GetEmptyZones();

    void ResetUnusedIOZones();
    void LogZoneStats();
    void LogZoneUsage();

    int GetReadFD() { return read_f_; }
    int GetReadDirectFD() { return read_direct_f_; }
    int GetWriteFD() { return write_f_; }

    uint32_t GetZoneSize() { return zone_sz_; }
    uint32_t GetNrZones() { return nr_zones_; }
    FILE *GetZoneLogFile() { return zone_log_file_; }
    std::vector<Zone *> GetMetaZones() { return meta_zones; }

    void SetFinishTreshold(uint32_t threshold) { finish_threshold_ = threshold; }

    void NotifyIOZoneFull();
    void NotifyIOZoneClosed();
    void NotifyZoneAllocateAvail();
    void NotifyGarbageCollectionRun();
    void NotifyGCOn();
    long GetOpenIOZone() { return open_io_zones_; }
    unsigned int GetMaxNrOpenIOZone() { return max_nr_open_io_zones_; }
    long GetActiveIOZone() { return active_io_zones_; }
    unsigned int GetMaxNrActiveIOZone() { return max_nr_active_io_zones_; }
    bool IsGcRunning(){return is_gc_running;}
    FILE *GetAllocLog(){return alloc_log_file_;}
    void CountAlloc(){count_alloc++;}
  private:
    void GarbageCollectionThread(void);
    uint32_t GarbageCollection(const bool &is_trigger,
                              const Env::WriteLifeTimeHint lifetime,
                              const bool &is_force,
                              const uint64_t &current_empty_space);
    Slice ReadDataFromExtent(const ZoneMapEntry *item, char *scratch,
                            ZoneExtent **target_extent);
    IOStatus CopyDataToFile(const ZoneMapEntry *item, Slice &source,
                            char *scratch);
    void WaitUntilZoneOpenAvail();
    void WaitUntilZoneAllocateAvail();
    template <class Duration>
    bool GarbageCollectionSchedule(const Duration &duration);
    bool ZoneVictimEnableCheck(Zone *z, const bool &is_force);
    void ZoneSelectVictim(const int &invalid_level, const bool &is_force);
    ZoneGcState ValidDataCopy(Env::WriteLifeTimeHint lifetime, Zone *z);
    uint64_t ZoneResetAndFinish(Zone *z, bool reset_condition,
                                bool finish_condition, Zone **callback_victim);

    int AllocateEmptyZone(unsigned int best_diff, Zone *finish_victim,
                          Zone **allocated_zone, Env::WriteLifeTimeHint lifetime);
    int GetAlreadyOpenZone(Zone **allocated_zone, ZoneFile *file,
                          Env::WriteLifeTimeHint lifetime);
    std::string GetZoneFileExt(const std::string filename);
    Zone *AllocateZone(Env::WriteLifeTimeHint lifetime, ZoneFile *file);

};

class ZonedBlockDevice {
 private:
  std::vector<SubZonedBlockDevice*> s_zbds_;
  unsigned int meta_num_; // Turn number of S_ZBD's meta zone allocation.
  std::string filename_;
  std::shared_ptr<Logger> logger_;
  CPUChecker cpu_checker_;
  FILE* cpu_check_file_;
  std::thread *check_thread_;
  std::atomic<bool> check_thread_stop;
  std::mutex check_thread_mutex_;
  std::condition_variable cpu_check_cond_;
 public:
  explicit ZonedBlockDevice(std::string bdevname,
                            std::shared_ptr<Logger> logger);
  virtual ~ZonedBlockDevice();
  
#if defined(INDEPENDENT_GC_THREAD)
  std::map<SubZonedBlockDevice*,std::mutex*> zonefile_mtxs_;
#else
  std::mutex *files_mtx_;
#endif

  void LockMutex();
  void UnlockMutex();
  std::mutex *GetMtxOnFile(ZoneFile *zonefile); 
  SubZonedBlockDevice* AllocateSubZBD(ZoneFile* zonefile);
  IOStatus Open(bool readonly = false);

  Zone *AllocateMetaZone();

  uint64_t GetFreeSpace();
  std::string GetFilename();
  uint32_t GetBlockSize();

  void ResetUnusedIOZones();
  void LogZoneStats(SubZonedBlockDevice* s_zbd = nullptr);
  void LogZoneUsage(SubZonedBlockDevice* s_zbd = nullptr);
  void ShareFileMtx();

  uint32_t GetZoneSize();
  uint32_t GetNrZones();
  FILE *GetZoneLogFile(ZoneFile* zoneFile);
  std::vector<Zone *> GetMetaZones();

  void SetFinishTreshold(uint32_t threshold);
 
 private:
  unsigned int GetMetaNum();
  void LogCPUCheck();
  template <class Duration>
  void WaitUntilGCOn(const Duration &duration);
};
}  // namespace ROCKSDB_NAMESPACE
#endif  // !defined(ROCKSDB_LITE) && defined(OS_LINUX) && defined(LIBZBD)
