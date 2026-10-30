#ifndef TANDEM_H
#define TANDEM_H

#include <string>
#include <algorithm>
#include <cctype>
#include <boost/algorithm/string.hpp>

#include "svanno.h"

namespace torali
{

  // Tandem repeat period from REF
  inline int32_t
  _trFlankPeriod(char const* seq, int32_t const N, int32_t const a) {
    if ((a < 0) || (a > N)) return 0;
    int32_t lb = std::max(0, a - 200);
    std::string flank(seq + lb, seq + a);
    int32_t p = 0;
    if (flank.find('N') == std::string::npos) p = detectTandemRepeat(flank, 100, 0.85f).first;
    if (p > 0) return p;
    int32_t rb = std::min(N, a + 200);
    flank = std::string(seq + a, seq + rb);
    if (flank.find('N') != std::string::npos) return 0;
    return detectTandemRepeat(flank, 100, 0.85f).first;
  }

  // Tandem repeat period from REF window
  inline int32_t
  _trWindowPeriod(char const* seq, int32_t const N, int32_t const b, int32_t const e, float const minFrac) {
    if ((b < 0) || (e > N) || (e - b < 10)) return 0;
    std::string w = boost::to_upper_copy(std::string(seq + b, seq + e));
    if (w.find('N') != std::string::npos) return 0;
    return detectTandemRepeat(w, 100, minFrac).first;
  }

  inline int32_t
  _trFlankPeriod(std::string const& chrSeq, int32_t const a) {
    return _trFlankPeriod(chrSeq.data(), (int32_t) chrSeq.size(), a);
  }

  // Reference repeat tract
  inline bool
  _trTract(char const* seq, int32_t const N, int32_t const a, int32_t const p, int32_t& ts, int32_t& te) {
    if ((p <= 0) || (p > 100) || (a < 0) || (a + p > N)) return false;
    double const minFrac = 0.7;
    int32_t const maxTract = 10000;
    te = a;
    while ((te + 2 * p <= N) && (te - a <= maxTract)) {
      int32_t m = 0;
      for(int32_t i = te; i < te + p; ++i) {
	if ((std::toupper(seq[i]) == std::toupper(seq[i + p])) && (std::toupper(seq[i]) != 'N')) ++m;
      }
      if (m < minFrac * p) break;
      te += p;
    }
    if (te + p <= N) te += p;
    ts = a;
    while ((ts - p >= 0) && (a - ts <= maxTract)) {
      int32_t m = 0;
      for(int32_t i = ts - p; i < ts; ++i) {
	if ((std::toupper(seq[i]) == std::toupper(seq[i + p])) && (std::toupper(seq[i]) != 'N')) ++m;
      }
      if (m < minFrac * p) break;
      ts -= p;
    }
    if (te > N) te = N;
    if (te - ts > maxTract) return false;
    return ((te - ts) >= std::max(10, 2 * p));
  }

  inline bool
  _trTract(std::string const& chrSeq, int32_t const a, int32_t const p, int32_t& ts, int32_t& te) {
    return _trTract(chrSeq.data(), (int32_t) chrSeq.size(), a, p, ts, te);
  }

}

#endif
