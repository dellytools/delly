/*
============================================================================
DELLY: Structural variant discovery by integrated PE mapping and SR analysis
============================================================================
Copyright (C) 2012 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#ifndef BOLOG_H
#define BOLOG_H

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace torali {


template<typename TPrecision>
struct BoLog {
  std::vector<TPrecision> phred2prob;

  BoLog() {
    for(int i = 0; i<=256; ++i) phred2prob.push_back(boost::multiprecision::pow(TPrecision(10), -(TPrecision(i)/TPrecision(10))));
  }
};


 template<typename TMapqVector, typename FLP, typename TGQ, typename TGenotype>
   inline bool
   _computeGLs(TMapqVector const& mapqRef, TMapqVector const& mapqAlt, FLP* gl, TGQ& gq, TGenotype& gtype) {
   // Compute genotype likelihoods
   BoLog<FLP> bl;
   
   for(unsigned int geno=0; geno<=2; ++geno) gl[geno]=0;
   unsigned int peDepth=mapqRef.size() + mapqAlt.size();
   for(typename TMapqVector::const_iterator mapqRefIt = mapqRef.begin();mapqRefIt!=mapqRef.end();++mapqRefIt) {
     gl[0] += boost::multiprecision::log10(bl.phred2prob[*mapqRefIt]);
     gl[1] += boost::multiprecision::log10(bl.phred2prob[*mapqRefIt] + (FLP(1) - bl.phred2prob[*mapqRefIt]));
     gl[2] += boost::multiprecision::log10(FLP(1) - bl.phred2prob[*mapqRefIt]);
   }
   for(typename TMapqVector::const_iterator mapqAltIt = mapqAlt.begin();mapqAltIt!=mapqAlt.end();++mapqAltIt) {
     gl[0] += boost::multiprecision::log10(FLP(1) - bl.phred2prob[*mapqAltIt]);
     gl[1] += boost::multiprecision::log10((FLP(1) - bl.phred2prob[*mapqAltIt]) + bl.phred2prob[*mapqAltIt]);
     gl[2] += boost::multiprecision::log10(bl.phred2prob[*mapqAltIt]);
   }
   gl[1] += -FLP(peDepth) * boost::multiprecision::log10(FLP(2));
   unsigned int glBest=0;
   FLP glBestVal=gl[glBest];
   for(unsigned int geno=1; geno<=2; ++geno) {
     if (gl[geno] > glBestVal) {
       glBestVal=gl[geno];
       glBest = geno;
     }
   }
   // Rescale by best genotype and get second best genotype for GQ
   FLP glSecondBestVal=-std::numeric_limits<FLP>::infinity();
   for(unsigned int geno=0; geno<=2; ++geno) {
     if ((gl[geno]>glSecondBestVal) && (gl[geno]<=glBestVal) && (geno!=glBest)) glSecondBestVal=gl[geno];
     gl[geno] -= glBestVal;
   }

   // Phred-scaled genotype likelihoods
   uint32_t pl[3];
   pl[0] = boost::multiprecision::iround(-10 * gl[0]);
   pl[1] = boost::multiprecision::iround(-10 * gl[1]);
   pl[2] = boost::multiprecision::iround(-10 * gl[2]);
   for (unsigned int geno=0; geno<=2; ++geno) pl[geno] = (pl[geno] < 255) ? pl[geno] : 255;
   if ((peDepth) && (pl[0] + pl[1] + pl[2] > 0)) {
     gq = boost::multiprecision::iround(-10 * boost::multiprecision::log10((1-1/(bl.phred2prob[pl[0]]+bl.phred2prob[pl[1]]+bl.phred2prob[pl[2]]))));
     if (glBest==0) gtype = "1/1";
     else if (glBest==1) gtype = "0/1";
     else gtype = "0/0";
   } else {
     gtype = "./.";
     gq = 0;
   }
   return (gq > 0);
 }
 

}

#endif
