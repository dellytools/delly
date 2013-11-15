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


#ifndef TOKENIZER_H
#define TOKENIZER_H

#include "memory_mapped_file.h"

namespace torali
{

class Tokenizer {
public:
  enum {
    MAX_LINE_LENGTH = Memory_mapped_file::MAX_LINE_LENGTH
  };

  typedef unsigned int TPos;
  char* data;
  TPos cPos;
  TPos cSize;
  char split;

  Tokenizer() : data(0), cPos(0), cSize(0), split('\t') {}

  Tokenizer(char* buffer, TPos total_size) : 
    data(buffer),
    cPos(0),
    cSize(total_size),
    split('\t') {}


  inline void _readNextChunk(char* local) {
    TPos pPos = cPos;
    for(;pPos<cSize;++pPos) {
      if ((data[pPos] == split) || (data[pPos] == '\0')) {
	memcpy(local, data+cPos, pPos);
	local[pPos - cPos] = '\0';
	break;
      }
    }
    cPos = pPos + 1;
  }

  inline void resetBuffer(char* buffer, TPos total_size) {
    data = buffer;
    cPos = 0;
    cSize = total_size;
  }

  inline TPos _getLenNextChunk() {
    TPos pPos = cPos;
    for(;pPos<cSize;++pPos)
      if (data[pPos] == split) 
	return pPos - cPos;
    return 0;
  }

  inline void skipNextChunk() {
    TPos pPos = cPos;
    for(;pPos<cSize;++pPos) {
      if (data[pPos] == split) break;
    }
    cPos = pPos + 1;
  }

  inline void getString(std::string& name) {
    char buffer[MAX_LINE_LENGTH];
    _readNextChunk(buffer);
    name = buffer;
  }

  inline void getBuffer(char buffer[]) {
    _readNextChunk(buffer);
  }
  
  inline int getInt() {    
    char buffer[MAX_LINE_LENGTH];
    _readNextChunk(buffer);
    return atoi(buffer);
  }

  inline unsigned int getUInt() {    
    char buffer[MAX_LINE_LENGTH];
    _readNextChunk(buffer);
    return (unsigned int) atoi(buffer);
  }

  inline unsigned short getUShort() {    
    char buffer[MAX_LINE_LENGTH];
    _readNextChunk(buffer);
    return (unsigned short) atoi(buffer);
  }

  inline short getShort() {    
    char buffer[MAX_LINE_LENGTH];
    _readNextChunk(buffer);
    return (short) atoi(buffer);
  }

  inline double getDouble() {    
    char buffer[MAX_LINE_LENGTH];
    _readNextChunk(buffer);
    return (double) atof(buffer);
  }

  inline float getFloat() {    
    char buffer[MAX_LINE_LENGTH];
    _readNextChunk(buffer);
    return (float) atof(buffer);
  }

};


}

#endif

