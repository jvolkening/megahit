#include <iostream>

#include "localasm/local_support.h"

int main() {
  ContigInfo info;
  info.set_kmer_size(3);
  info.set_counts({1, 1, 2, 2});
  if (AssessLocalSupport(info, 6, 2) != LocalSupport::kSupported) {
    std::cerr << "Half-supported contig should be retained\n";
    return 1;
  }
  info.set_counts({1, 1, 1, 2});
  if (AssessLocalSupport(info, 6, 2) != LocalSupport::kUnsupported ||
      AssessLocalSupport(info, 6, 0) != LocalSupport::kSupported) {
    std::cerr << "Singleton support or disabled filtering is incorrect\n";
    return 1;
  }
  info.set_counts({1, 1, 2});
  if (AssessLocalSupport(info, 5, 2) != LocalSupport::kUnsupported) {
    std::cerr << "An odd number of positions requires a majority\n";
    return 1;
  }
  info.set_counts({1, 2, 2});
  if (AssessLocalSupport(info, 5, 2) != LocalSupport::kSupported) {
    std::cerr << "Majority support should be retained\n";
    return 1;
  }
  if (AssessLocalSupport(info, 7, 2) != LocalSupport::kUnknown) {
    std::cerr << "Incomplete support metadata should not be filtered\n";
    return 1;
  }
  info.set_kmer_size(0);
  if (AssessLocalSupport(info, 6, 2) != LocalSupport::kUnknown) {
    std::cerr << "Missing k-mer size should not be filtered\n";
    return 1;
  }
  return 0;
}
