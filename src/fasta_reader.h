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

#ifndef FASTA_READER_H
#define FASTA_READER_H

#include "alphabet.h"
#include "memory_mapped_file.h"


namespace torali
{

  template<typename TName, typename TLen, typename TAlphabet, typename TSeq, typename TAlignDir>
    struct FastaRecord {
      typedef TSeq TSequence;
      TName chrName;
      TLen chrLen;
      TSeq seq;
      TAlignDir alignDir;

      inline void
      reserve() {
	seq.reserve(0);
      }

      inline void
      clear() {
	chrName = "";
	chrLen = 0;
	seq.clear();
      }
    };

  template<typename TName, typename TLen, typename TAlphabet, typename TSeq>
    struct FastaRecord<TName, TLen, TAlphabet, TSeq, void> {
      typedef TSeq TSequence;
      TName chrName;
      TLen chrLen;
      TSeq seq;

      inline void
      reserve() {
	seq.reserve(0);
      }

      inline void
      clear() {
	chrName = "";
	chrLen = 0;
	seq.clear();
      }
    };

  template<typename TName, typename TLen>
    struct FastaRecord<TName, TLen, void, void, void> {
      typedef void TSequence;
      TName chrName;
      TLen chrLen;

      inline void
      reserve() {}

      inline void 
      clear() {
	chrName = "";
	chrLen = 0;
      }
    };


  template<typename TName, typename TLen, typename TAlphabet, typename TSeq, typename TAlignDir, typename TInputName>
    inline
    void addFastaRecordName(FastaRecord<TName, TLen, TAlphabet, TSeq, TAlignDir>& rec, TInputName const& name) 
  {
    typename TInputName::size_type firstHit=name.find_first_of(" \t\f\v\n\r");
    if (firstHit != name.npos) rec.chrName=name.substr(0, firstHit);
    else rec.chrName=name;
  }

  template<typename TName, typename TLen, typename TAlphabet, typename TSeq, typename TAlignDir, typename TInputLen>
    inline
    void addFastaRecordLen(FastaRecord<TName, TLen, TAlphabet, TSeq, TAlignDir>& rec, TInputLen len) 
  {
    rec.chrLen = len;
  }

  template<typename TName, typename TLen, typename TChar>
    inline
    void addFastaRecordSeq(FastaRecord<TName, TLen, void, void, void>&, TChar) 
  {
  }

  template<typename TName, typename TLen, typename TSeq, typename TAlignDir, typename TChar>
    inline
    void addFastaRecordSeq(FastaRecord<TName, TLen, Dna5Alphabet, TSeq, TAlignDir>& rec, TChar c) 
  {
    rec.seq.push_back(dna5_encode[int(c)]);
  }

  template<typename TName, typename TLen, typename TSeq, typename TAlignDir, typename TChar>
    inline
    void addFastaRecordSeq(FastaRecord<TName, TLen, Dna5GapAlphabet, TSeq, TAlignDir>& rec, TChar c) 
  {
    rec.seq.push_back(dna5gap_encode[int(c)]);
  }


  template<typename TName, typename TLen, typename TSeq, typename TAlignDir, typename TChar>
    inline
    void addFastaRecordSeq(FastaRecord<TName, TLen, DnaAlphabet, TSeq, TAlignDir>& rec, TChar c) 
  {
    rec.seq.push_back(dna_encode[int(c)]);
  }

  template<typename TName, typename TLen, typename TSeq, typename TAlignDir, typename TChar>
    inline
    void addFastaRecordSeq(FastaRecord<TName, TLen, DnaGapAlphabet, TSeq, TAlignDir>& rec, TChar c) 
  {
    rec.seq.push_back(dnagap_encode[int(c)]);
  }

  template<typename TName, typename TLen, typename TAlphabet, typename TSeq, typename TBool>
    inline
    void addAlignDir(FastaRecord<TName, TLen, TAlphabet, TSeq, void>&, TBool const) 
  {
  }

  template<typename TName, typename TLen, typename TAlphabet, typename TSeq, typename TAlignDir, typename TBool>
    inline
    void addAlignDir(FastaRecord<TName, TLen, TAlphabet, TSeq, TAlignDir>& rec, TBool const align) 
  {
    rec.alignDir=align;
  }

  


template<typename TFastaRecord>
class Fasta_reader {
public:
  typedef unsigned long int TLen;
  typedef std::vector<TFastaRecord> TFastaVec;

  enum {
    MAX_BUFFER_SIZE = 4096
  };

  TFastaVec faVec;
  TLen totalSize;

  Fasta_reader() : totalSize(0) {}

  inline void read_fasta(std::string const& in_path) {
    Memory_mapped_file map_file(in_path);
    char buffer[Fasta_reader::MAX_BUFFER_SIZE];
    std::string chr;
    bool header = false;
    TLen seqSize = 0;
    TFastaRecord a;    
    while (map_file.left_bytes() > 0) {
      Memory_mapped_file::TLen readLen = map_file.read_bytes(buffer, Fasta_reader::MAX_BUFFER_SIZE);
      char* p = buffer;
      char* pEnd = buffer + readLen;
      for(; p!=pEnd; ++p) {
	if (*p == '>') {
	  if (seqSize > 0) {
	    addFastaRecordName(a, chr);
	    addFastaRecordLen(a, seqSize);
	    faVec.push_back(a);
	    faVec[faVec.size() - 1].reserve();
	    a.clear();
	    totalSize+=seqSize;
	  }
	  header = true;
	  seqSize = 0;
	  chr.clear();
	  chr.reserve(500);
	} else if (header) {
	  if ((*p == '\n') || (*p == '\r')) header = false;
	  else chr.push_back(*p);
	} else {
	  if ((*p != '\n') && (*p != '\r')) {
	    addFastaRecordSeq(a,*p);
	    ++seqSize;
	  }
	}
      }
    }
    if (seqSize > 0) {
      addFastaRecordName(a, chr);
      addFastaRecordLen(a, seqSize);
      faVec.push_back(a);
      faVec[faVec.size() - 1].reserve();
      a.clear();
      totalSize+=seqSize;
    }
    map_file.close();
  }


  template<typename TFile, typename TBuffer, typename TSize>
    inline void _read_single_line_fasta_with_key(TFile& map_file, TBuffer* buffer, TSize& readLen, std::string const& key) 
  {
    totalSize = 0;
    faVec.clear();    
    std::string chr;
    std::string keyId;
    bool validRead = false;
    bool processing = true;
    TFastaRecord a;    
    int field = 0;
    TLen seqSize = 0;
    while ((processing) && (map_file.left_bytes() > 0)) {
      if (readLen == 0) {
	readLen = map_file.read_bytes(buffer, Fasta_reader::MAX_BUFFER_SIZE);
      }
      char* p = buffer;
      char* pEnd = buffer + readLen;
      for(; p!=pEnd; ++p) {
	if (field == 0) {
	  if (*p == '\t') {
	    ++field;
	    chr.clear();
	    chr.reserve(100);
	    a.clear();
	    seqSize = 0;
	    if (keyId.compare(key) == 0) validRead = true;
	    else {
	      validRead = false;
	      field = 0;
	      if (totalSize > 0) processing = false;
	    }
	    keyId.clear();
	  } else {
	    if ((*p == '\n') || (*p == '\r')) keyId.clear();
	    else keyId.push_back(*p);
	  }
	} else if (validRead) {
	  if (field == 1) {
	    if (*p == '\t') ++field;
	    else chr.push_back(*p);
	  } else if (field == 2) {
	    if (*p == '\t') ++field;
	    else {
	      if (*p == '1') addAlignDir(a, true);
	      else addAlignDir(a, false);
	    }
	  } else {
	    if ((*p != '\n') && (*p != '\r')) {
	      addFastaRecordSeq(a,*p);
	      ++seqSize;
	    } else {
	      if (seqSize > 0) {
		addFastaRecordName(a, chr);
		addFastaRecordLen(a, seqSize);
		faVec.push_back(a);
		totalSize+=seqSize;
	      }
	      field = 0;
	    }
	  }
	} else {
	  if ((*p == '\n') || (*p == '\r')) field = 0;
	}
      }
      if (processing) readLen = 0;
    }
  }


};

}

#endif
