#ifndef MEGAHIT_LOCALASM_LOCAL_SUPPORT_H
#define MEGAHIT_LOCALASM_LOCAL_SUPPORT_H

#include <algorithm>
#include <cstddef>
#include <cstdint>

#include "idba/contig_info.h"

enum class LocalSupport { kSupported, kUnsupported, kUnknown };

inline LocalSupport AssessLocalSupport(const ContigInfo &info, size_t length,
                                       uint32_t minimum) {
  if (minimum == 0) return LocalSupport::kSupported;
  const auto k = info.kmer_size();
  const auto &counts = info.counts();
  if (k == 0 || k > length || counts.size() != length - k + 1) {
    return LocalSupport::kUnknown;
  }
  auto supported = std::count_if(counts.begin(), counts.end(),
                                 [minimum](SequenceCountUnitType count) {
                                   return count >= minimum;
                                 });
  return static_cast<size_t>(supported) >= (counts.size() + 1) / 2
             ? LocalSupport::kSupported
             : LocalSupport::kUnsupported;
}

#endif  // MEGAHIT_LOCALASM_LOCAL_SUPPORT_H
