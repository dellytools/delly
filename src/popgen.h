#ifndef POPGEN_H
#define POPGEN_H

#include <vector>
#include <cmath>
#include <boost/math/distributions/chi_squared.hpp>

namespace torali
{

  // EM estimate of the bi-allelic allele frequency under HWE
  template<typename TConfig, typename TGlVector, typename TValue>
  inline void
  _estBiallelicAF(TConfig const& c, TGlVector const& glVector, TValue (&hweAF)[2]) {
    if (glVector.empty()) return;
    TValue afprior[2];
    afprior[0] = 0.5;
    afprior[1] = 0.5;
    TValue gtprior[3];
    TValue gt[3];
    TValue p;
    TValue err = 1;
    for(std::size_t count = 0; ((err > c.epsilon) && (count < c.maxiter)); ++count) {
      gtprior[0] = afprior[0] * afprior[0];
      gtprior[1] = 2 * afprior[0] * afprior[1];
      gtprior[2] = afprior[1] * afprior[1];
      hweAF[0] = 0;
      hweAF[1] = 0;
      TValue used = 0;
      for(typename TGlVector::const_iterator itG = glVector.begin(); itG != glVector.end(); ++itG) {
	gt[0] = gtprior[0] * (*itG)[0];
	gt[1] = gtprior[1] * (*itG)[1];
	gt[2] = gtprior[2] * (*itG)[2];
	p = gt[0] + gt[1] + gt[2];
	if (p <= 0) continue;
	gt[0] /= p;
	gt[1] /= p;
	gt[2] /= p;
	hweAF[0] += gt[0] + 0.5 * gt[1];
	hweAF[1] += gt[2] + 0.5 * gt[1];
	used += 1;
      }
      if (used > 0) {
	hweAF[0] /= used;
	hweAF[1] /= used;
      }
      err = (afprior[0]-hweAF[0])*(afprior[0]-hweAF[0]) + (afprior[1]-hweAF[1])*(afprior[1]-hweAF[1]);
      afprior[0] = hweAF[0];
      afprior[1] = hweAF[1];
    }
  }

  // EM estimate of the bi-allelic allele frequency under HWE with diploid and haploid samples
  template<typename TConfig, typename TGlVector, typename TGlPairVector, typename TValue>
  inline void
  _estBiallelicAF(TConfig const& c, TGlVector const& glVector, TGlPairVector const& glPairs, TValue (&hweAF)[2]) {
    if (glPairs.empty()) {
      _estBiallelicAF(c, glVector, hweAF);
      return;
    }
    TValue afprior[2];
    afprior[0] = 0.5;
    afprior[1] = 0.5;
    TValue gtprior[3];
    TValue gt[3];
    TValue p;
    TValue err = 1;
    for(std::size_t count = 0; ((err > c.epsilon) && (count < c.maxiter)); ++count) {
      gtprior[0] = afprior[0] * afprior[0];
      gtprior[1] = 2 * afprior[0] * afprior[1];
      gtprior[2] = afprior[1] * afprior[1];
      TValue refAlleles = 0;
      TValue altAlleles = 0;
      for(typename TGlVector::const_iterator itG = glVector.begin(); itG != glVector.end(); ++itG) {
	gt[0] = gtprior[0] * (*itG)[0];
	gt[1] = gtprior[1] * (*itG)[1];
	gt[2] = gtprior[2] * (*itG)[2];
	p = gt[0] + gt[1] + gt[2];
	if (p <= 0) continue;
	refAlleles += (2 * gt[0] + gt[1]) / p;
	altAlleles += (2 * gt[2] + gt[1]) / p;
      }
      for(typename TGlPairVector::const_iterator itH = glPairs.begin(); itH != glPairs.end(); ++itH) {
	gt[0] = afprior[0] * (*itH)[0];
	gt[1] = afprior[1] * (*itH)[1];
	p = gt[0] + gt[1];
	if (p <= 0) continue;
	refAlleles += gt[0] / p;
	altAlleles += gt[1] / p;
      }
      if (refAlleles + altAlleles > 0) {
	hweAF[0] = refAlleles / (refAlleles + altAlleles);
	hweAF[1] = altAlleles / (refAlleles + altAlleles);
      }
      err = (afprior[0]-hweAF[0])*(afprior[0]-hweAF[0]) + (afprior[1]-hweAF[1])*(afprior[1]-hweAF[1]);
      afprior[0] = hweAF[0];
      afprior[1] = hweAF[1];
    }
  }

  // EM estimate of the GT frequencies.
  template<typename TConfig, typename TGlVector, typename TValue>
  inline void
  _estBiallelicGTFreq(TConfig const& c, TGlVector const& glVector, TValue (&mleGTFreq)[3]) {
    if (glVector.empty()) return;
    TValue prior[3];
    prior[0] = 1.0/3.0;
    prior[1] = 1.0/3.0;
    prior[2] = 1.0/3.0;
    TValue gt[3];
    TValue p;
    TValue err = 1;
    for(std::size_t count = 0; ((err > c.epsilon) && (count < c.maxiter)); ++count) {
      mleGTFreq[0] = 0;
      mleGTFreq[1] = 0;
      mleGTFreq[2] = 0;
      TValue used = 0;
      for(typename TGlVector::const_iterator itG = glVector.begin(); itG != glVector.end(); ++itG) {
	gt[0] = prior[0] * (*itG)[0];
	gt[1] = prior[1] * (*itG)[1];
	gt[2] = prior[2] * (*itG)[2];
	p = gt[0] + gt[1] + gt[2];
	if (p <= 0) continue;
	mleGTFreq[0] += gt[0]/p;
	mleGTFreq[1] += gt[1]/p;
	mleGTFreq[2] += gt[2]/p;
	used += 1;
      }
      if (used > 0) {
	mleGTFreq[0] /= used;
	mleGTFreq[1] /= used;
	mleGTFreq[2] /= used;
      }
      err = (prior[0]-mleGTFreq[0])*(prior[0]-mleGTFreq[0]) + (prior[1]-mleGTFreq[1])*(prior[1]-mleGTFreq[1]) + (prior[2]-mleGTFreq[2])*(prior[2]-mleGTFreq[2]);
      prior[0] = mleGTFreq[0];
      prior[1] = mleGTFreq[1];
      prior[2] = mleGTFreq[2];
    }
  }

  // Inbreeding coefficient
  template<typename TGlVector, typename TValue>
  inline void
  _estBiallelicFIC(TGlVector const& glVector, TValue const (&hweAF)[2], TValue& F) {
    if (glVector.empty()) return;
    TValue hweGT[3];
    hweGT[0] = hweAF[0] * hweAF[0];
    hweGT[1] = 2 * hweAF[0] * hweAF[1];
    hweGT[2] = hweAF[1] * hweAF[1];
    TValue sumGLHet = 0;
    TValue denominator = 0;
    for(typename TGlVector::const_iterator itG = glVector.begin(); itG != glVector.end(); ++itG) {
      TValue p = (*itG)[0] * hweGT[0] + (*itG)[1] * hweGT[1] + (*itG)[2] * hweGT[2];
      if (p <= 0) continue;
      sumGLHet += (((*itG)[1] * hweGT[1]) / p);
      denominator += hweGT[1];
    }
    if (denominator > 0) F = 1 - sumGLHet/denominator;
  }

  // Imputation-quality
  template<typename TGlVector, typename TValue>
  inline void
  _estBiallelicRSQ(TGlVector const& glVector, TValue const (&hweAF)[2], TValue& rsq) {
    if (glVector.empty()) return;
    TValue hweGT[3];
    hweGT[0] = hweAF[0] * hweAF[0];
    hweGT[1] = 2 * hweAF[0] * hweAF[1];
    hweGT[2] = hweAF[1] * hweAF[1];
    TValue post[3];
    TValue p = 0;
    TValue sumD = 0;
    TValue sumD2 = 0;
    TValue used = 0;
    for(typename TGlVector::const_iterator itG = glVector.begin(); itG != glVector.end(); ++itG) {
      post[0] = (*itG)[0] * hweGT[0];
      post[1] = (*itG)[1] * hweGT[1];
      post[2] = (*itG)[2] * hweGT[2];
      p = post[0] + post[1] + post[2];
      if (p <= 0) continue;
      post[0] /= p;
      post[1] /= p;
      post[2] /= p;
      sumD += (post[1] + 2 * post[0]);
      sumD2 += (post[1] + 2 * post[0]) * (post[1] + 2 * post[0]);
      used += 1;
    }
    if ((used > 1) && (hweGT[1] > 0)) {
      TValue meanD = sumD/used;
      sumD2 = (sumD2 - used * meanD * meanD);
      if (sumD2 < 0) sumD2 = 0;
      sumD2 /= (used - 1);
      rsq = sumD2 / hweGT[1];
    }
  }

  // Imputation-quality across diploid and haploid samples (pooled within-ploidy dosage variance)
  template<typename TGlVector, typename TValue>
  inline void
  _estBiallelicRSQ(TGlVector const& glVector, TGlVector const& glPairs, TValue const (&hweAF)[2], TValue& rsq) {
    if (glPairs.empty()) {
      _estBiallelicRSQ(glVector, hweAF, rsq);
      return;
    }
    TValue hweGT[3];
    hweGT[0] = hweAF[0] * hweAF[0];
    hweGT[1] = 2 * hweAF[0] * hweAF[1];
    hweGT[2] = hweAF[1] * hweAF[1];
    if (hweGT[1] <= 0) return;
    TValue sumD[2] = {0, 0};
    TValue sumD2[2] = {0, 0};
    TValue used[2] = {0, 0};
    for(typename TGlVector::const_iterator itG = glVector.begin(); itG != glVector.end(); ++itG) {
      TValue post[3];
      post[0] = (*itG)[0] * hweGT[0];
      post[1] = (*itG)[1] * hweGT[1];
      post[2] = (*itG)[2] * hweGT[2];
      TValue p = post[0] + post[1] + post[2];
      if (p <= 0) continue;
      TValue d = (post[1] + 2 * post[0]) / p;
      sumD[0] += d;
      sumD2[0] += d * d;
      used[0] += 1;
    }
    for(typename TGlVector::const_iterator itG = glPairs.begin(); itG != glPairs.end(); ++itG) {
      TValue post[2];
      post[0] = (*itG)[0] * hweAF[0];
      post[1] = (*itG)[1] * hweAF[1];
      TValue p = post[0] + post[1];
      if (p <= 0) continue;
      TValue d = post[0] / p;
      sumD[1] += d;
      sumD2[1] += d * d;
      used[1] += 1;
    }
    TValue num = 0;
    TValue den = 0;
    for(int32_t k = 0; k < 2; ++k) {
      if (used[k] < 2) continue;
      TValue meanD = sumD[k] / used[k];
      TValue ss = sumD2[k] - used[k] * meanD * meanD;
      if (ss < 0) ss = 0;
      num += ss;
      den += (used[k] - 1) * ((k == 0) ? hweGT[1] : (hweAF[0] * hweAF[1]));
    }
    if (den > 0) rsq = num / den;
  }

  template<typename TDosage>
  inline double
  _dosageR2(TDosage const& a, TDosage const& b, int32_t const minShared) {
    double n = 0;
    double sa = 0;
    double sb = 0;
    double saa = 0;
    double sbb = 0;
    double sab = 0;
    std::size_t m = (a.size() < b.size()) ? a.size() : b.size();
    for(std::size_t i = 0; i < m; ++i) {
      if ((a[i] < 0) || (b[i] < 0)) continue;
      double x = a[i];
      double y = b[i];
      n += 1;
      sa += x;
      sb += y;
      saa += x*x;
      sbb += y*y;
      sab += x*y;
    }
    if (n < (double) minShared) return -1;
    double cov = sab - sa * sb / n;
    double va = saa - sa * sa / n;
    double vb = sbb - sb * sb / n;
    if ((va <= 0) || (vb <= 0)) return -1;
    double r = cov / std::sqrt(va * vb);
    return r * r;
  }

  // HWE likelihood-ratio test (HWE-constrained vs unconstrained GT frequencies)
  template<typename TGlVector, typename TValue>
  inline void
  _estBiallelicHWE_LRT(TGlVector const& glVector, TValue const (&hweAF)[2], TValue const (&mleGTFreq)[3], TValue& pvalue) {
    if (glVector.empty()) return;
    TValue hweGT[3];
    hweGT[0] = hweAF[0] * hweAF[0];
    hweGT[1] = 2 * hweAF[0] * hweAF[1];
    hweGT[2] = hweAF[1] * hweAF[1];
    TValue null = 0;
    TValue alt = 0;
    for(typename TGlVector::const_iterator itG = glVector.begin(); itG != glVector.end(); ++itG) {
      TValue pnull = (*itG)[0] * hweGT[0] + (*itG)[1] * hweGT[1] + (*itG)[2] * hweGT[2];
      TValue palt = (*itG)[0] * mleGTFreq[0] + (*itG)[1] * mleGTFreq[1] + (*itG)[2] * mleGTFreq[2];
      if ((pnull <= 0) || (palt <= 0)) continue;
      null += std::log(pnull);
      alt += std::log(palt);
    }
    TValue lrts = -2 * (null - alt);
    if (lrts < 0) lrts = 0;
    boost::math::chi_squared chisqDist(1);
    pvalue = boost::math::cdf(complement(chisqDist, lrts));
  }

}

#endif
