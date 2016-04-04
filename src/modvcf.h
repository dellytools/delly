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

#ifndef MODVCF_H
#define MODVCF_H

#include <htslib/sam.h>
#include <htslib/vcf.h>


namespace torali
{


void _remove_info_tag(bcf_hdr_t* hdr, bcf1_t* rec, std::string const& tag) {
  bcf_update_info(hdr, rec, tag.c_str(), NULL, 0, BCF_HT_INT);  // Type does not matter for n = 0
}

void _remove_format_tag(bcf_hdr_t* hdr, bcf1_t* rec, std::string const& tag) {
  bcf_update_format(hdr, rec, tag.c_str(), NULL, 0, BCF_HT_INT);  // Type does not matter for n = 0
}

void _remove_info(bcf_hdr_t* hdr, bcf1_t* rec) {
  std::string tmp[] = {"CT", "PRECISE", "IMPRECISE", "SVTYPE", "SVMETHOD", "CIEND", "CIPOS", "CHR2", "END", "PE", "MAPQ", "SR", "SRQ", "CONSENSUS"};
  std::set<std::string> keepInfo(tmp, tmp + sizeof(tmp)/sizeof(tmp[0]));

  if (!(rec->unpacked & BCF_UN_INFO)) bcf_unpack(rec, BCF_UN_INFO);
  
  for (int i = 0; i < rec->n_info; ++i){
    bcf_info_t* inf = &rec->d.info[i];
    const char* key = bcf_hdr_int2id(hdr, BCF_DT_ID, inf->key);
    if (keepInfo.find(std::string(key)) != keepInfo.end()) continue;
    
    if (inf->vptr_free) {
      free(inf->vptr - inf->vptr_off);
      inf->vptr_free = 0;
    }
    rec->d.shared_dirty |= BCF1_DIRTY_INF;
    inf->vptr = NULL;
  }
}

void _remove_format(bcf_hdr_t* hdr, bcf1_t* rec) {
  if (!(rec->unpacked & BCF_UN_FMT)) bcf_unpack(rec, BCF_UN_FMT);
  
  for(int i = 0; i<rec->n_fmt; ++i) {
    bcf_fmt_t* fmt = &rec->d.fmt[i];
    const char* key = bcf_hdr_int2id(hdr, BCF_DT_ID, fmt->id);
    bcf_update_format(hdr, rec, key, NULL, 0, BCF_HT_INT); // the type is irrelevant for n = 0
    // Keep GT
    //if ((key[0]=='G') && key[1]=='T' && (!key[2])) continue;

    if (fmt->p_free) {
      free(fmt->p - fmt->p_off);
      fmt->p_free = 0;
    }
    rec->d.indiv_dirty = 1;
    fmt->p = NULL;
  }
}

inline int
_getInfoType(bcf_hdr_t const* hdr, std::string const& key) {
  return bcf_hdr_id2type(hdr, BCF_HL_INFO, bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str()));
}

inline int
_getFormatType(bcf_hdr_t const* hdr, std::string const& key) {
  return bcf_hdr_id2type(hdr, BCF_HL_FMT, bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str()));
}

inline bool _missing(bool const value) {
  return !value;
}

inline bool _missing(float const value) {
  return bcf_float_is_missing(value);
}

inline bool _missing(int8_t const value) {
  return (value == bcf_int8_missing);
}

inline bool _missing(int16_t const value) {
  return (value == bcf_int16_missing);
}

inline bool _missing(int32_t const value) {
  return (value == bcf_int32_missing);
}

inline bool _missing(std::string const& value) {
  return ((value.empty()) || (value == "."));
}

inline bool
_isKeyPresent(bcf_hdr_t const* hdr, std::string const& key) {
  return (bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str())>=0);
}


}

#endif
