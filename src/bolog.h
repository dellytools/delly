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



}

#endif
