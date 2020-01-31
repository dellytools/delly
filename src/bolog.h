#ifndef BOLOG_H
#define BOLOG_H

#include <boost/math/special_functions/round.hpp>

namespace torali {

#define SMALLEST_GL -1000

template<typename TPrecision>
struct BoLog {
  typedef TPrecision value_type;

  std::vector<TPrecision> phred2prob;

  BoLog() {
    for(int i = 0; i <= boost::math::round(-10 * SMALLEST_GL); ++i) phred2prob.push_back(std::pow(TPrecision(10), -(TPrecision(i)/TPrecision(10))));
  }
};


 template<typename TBoLog, typename TMapqVector>
 inline void
 _computeGLs(TBoLog const& bl, TMapqVector const& mapqRef, TMapqVector const& mapqAlt, float* gls, int32_t* gqval, int32_t* gts, int const file_c) {
   typedef typename TBoLog::value_type FLP;
   FLP gl[3];

   // Compute genotype likelihoods
   for(unsigned int geno=0; geno<=2; ++geno) gl[geno]=0;
   unsigned int peDepth=mapqRef.size() + mapqAlt.size();
   for(typename TMapqVector::const_iterator mapqRefIt = mapqRef.begin();mapqRefIt!=mapqRef.end();++mapqRefIt) {
     gl[0] += std::log10(bl.phred2prob[*mapqRefIt]);
     gl[1] += std::log10(bl.phred2prob[*mapqRefIt] + (FLP(1) - bl.phred2prob[*mapqRefIt]));
     gl[2] += std::log10(FLP(1) - bl.phred2prob[*mapqRefIt]);
   }
   for(typename TMapqVector::const_iterator mapqAltIt = mapqAlt.begin();mapqAltIt!=mapqAlt.end();++mapqAltIt) {
     gl[0] += std::log10(FLP(1) - bl.phred2prob[*mapqAltIt]);
     gl[1] += std::log10((FLP(1) - bl.phred2prob[*mapqAltIt]) + bl.phred2prob[*mapqAltIt]);
     gl[2] += std::log10(bl.phred2prob[*mapqAltIt]);
   }
   gl[1] += -FLP(peDepth) * std::log10(FLP(2));
   unsigned int glBest=0;
   FLP glBestVal=gl[glBest];
   for(unsigned int geno=1; geno<=2; ++geno) {
     if (gl[geno] > glBestVal) {
       glBestVal=gl[geno];
       glBest = geno;
     }
   }
   // Rescale by best genotype
   for(unsigned int geno=0; geno<=2; ++geno) {
     gl[geno] -= glBestVal;
     // Cap at smallest GL
     gl[geno] = (gl[geno] > SMALLEST_GL) ? gl[geno] : SMALLEST_GL;
   }

   // Phred-scaled genotype likelihoods
   uint32_t pl[3];
   pl[0] = (uint32_t) boost::math::round(-10 * gl[0]);
   pl[1] = (uint32_t) boost::math::round(-10 * gl[1]);
   pl[2] = (uint32_t) boost::math::round(-10 * gl[2]);
   if ((peDepth) && (pl[0] + pl[1] + pl[2] > 0)) {
     FLP likelihood = (FLP) std::log10((1-1/(bl.phred2prob[pl[0]]+bl.phred2prob[pl[1]]+bl.phred2prob[pl[2]])));
     likelihood = (likelihood > SMALLEST_GL) ? likelihood : SMALLEST_GL;
     gqval[file_c] = (int32_t) boost::math::round(-10 * likelihood);
     if (glBest==0) {
       gts[file_c * 2] = bcf_gt_unphased(1);
       gts[file_c * 2 + 1] = bcf_gt_unphased(1);
     } else if (glBest==1) {
       gts[file_c * 2] = bcf_gt_unphased(0);
       gts[file_c * 2 + 1] = bcf_gt_unphased(1);
     } else {
       gts[file_c * 2] = bcf_gt_unphased(0);
       gts[file_c * 2 + 1] = bcf_gt_unphased(0);
     }
   } else {
     gts[file_c * 2] = bcf_gt_missing;
     gts[file_c * 2 + 1] = bcf_gt_missing;
     gqval[file_c] = 0;
   }
   gls[file_c * 3 + 2] = (float) gl[0];
   gls[file_c * 3 + 1] = (float) gl[1];
   gls[file_c * 3] = (float) gl[2];
 }
 

}

#endif
