// Interface acceptance test: the parameter and matrix types carry the
// documented defaults and empty states. R1 (velocity_x, kubo_lehmann_ed),
// R2 (current_moments), R3 (optical_conductivity), and R4 (drude_weight)
// are implemented and gated in test_r1..test_r4. Every check prints its
// measured value; any failure exits nonzero.

#include "kubo/kubo.hpp"

#include <cstdio>
#include <string>

namespace {

int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

} // namespace

int main() {
  char buf[256];

  const kubo::KuboParams kp;
  std::snprintf(buf, sizeof buf, "n_cheb=%d n_random=%d a=%.2f b=%.1f seed=%llu",
                kp.n_cheb, kp.n_random, kp.a, kp.b,
                static_cast<unsigned long long>(kp.seed));
  check(kp.n_cheb == 512 && kp.n_random == 32 && kp.a == 0.0 && kp.b == 0.0 &&
            kp.seed == 0x9E3779B9ULL,
        "KuboParams defaults (a=0: unset, caller supplies the bound)", buf);

  const kubo::CsrMatrix empty;
  std::snprintf(buf, sizeof buf, "n=%d nnz=%d", empty.n, empty.nnz());
  check(empty.n == 0 && empty.nnz() == 0, "CsrMatrix empty state", buf);

  const kubo::MomentTable mt;
  std::snprintf(buf, sizeof buf, "n_cheb=%d mu.size=%zu", mt.n_cheb,
                mt.mu.size());
  check(mt.n_cheb == 0 && mt.mu.empty(), "MomentTable empty state", buf);

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
