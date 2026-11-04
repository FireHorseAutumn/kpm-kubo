// Interface acceptance test: the stubs exist, carry the
// documented defaults, and throw std::logic_error until implemented.
// R1 (velocity_x, kubo_lehmann_ed) is implemented and gated in test_r1;
// the R2-R4 stubs remain. Every check prints its measured value; any
// failure exits nonzero.

#include "kubo/kubo.hpp"

#include <cstdio>
#include <stdexcept>
#include <string>

namespace {

int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

template <typename F>
bool throws_not_implemented(F&& f) {
  try {
    f();
  } catch (const std::logic_error& e) {
    return std::string(e.what()).find("not implemented") != std::string::npos;
  }
  return false;
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

  check(throws_not_implemented(
            [&] { kubo::current_moments(empty, empty, kp); }),
        "current_moments stub throws", "DESIGN.md R2");
  check(throws_not_implemented([&] { kubo::optical_conductivity(mt, 0.0); }),
        "optical_conductivity stub throws", "DESIGN.md R3");
  check(throws_not_implemented([&] { kubo::drude_weight(empty, empty, 0.0); }),
        "drude_weight stub throws", "DESIGN.md R4");

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
