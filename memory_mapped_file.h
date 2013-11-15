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

#ifndef MEMORY_MAPPED_FILE_H
#define MEMORY_MAPPED_FILE_H

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/filesystem/convenience.hpp>

namespace torali
{


class Memory_mapped_file {
public:
  typedef unsigned long int TLen;
  enum {
    //MAX_LINE_LENGTH = 512,
    MAX_LINE_LENGTH = 8192,
    //MAX_LINE_LENGTH = 32768,
    NUM_PAGES = 10
  };

  std::string path;
  TLen page_start;
  TLen page_offset;
  TLen page_size;
  TLen mod_value;
  TLen div_value;
  TLen mmap_size;
  TLen file_size;
  TLen still_left;
  boost::iostreams::mapped_file_source map_file;

  Memory_mapped_file(const std::string& in_path) : 
    path(in_path),
    page_start(0), 
    page_offset(0), 
    page_size(map_file.alignment()),
    mod_value(page_size - 1),
    div_value(~(page_size - 1)),
    mmap_size(page_size * Memory_mapped_file::NUM_PAGES), 
    file_size(boost::filesystem::file_size(in_path)),
    still_left(file_size),
    map_file(path, (mmap_size > file_size) ? file_size : mmap_size, 0) {}


  Memory_mapped_file(const std::string& in_path, unsigned int num_pages) : 
    path(in_path),
    page_start(0), 
    page_offset(0), 
    page_size(map_file.alignment()),
    mod_value(page_size - 1),
    div_value(~(page_size - 1)), 
    mmap_size(page_size * num_pages), 
    file_size(boost::filesystem::file_size(in_path)),
    still_left(file_size),
    map_file(path, (mmap_size > file_size) ? file_size : mmap_size, 0) {}
	     
  
  inline TLen read_bytes(void* buffer, TLen num_bytes) {
    if (num_bytes > still_left) num_bytes = still_left;
    TLen end_pointer = page_offset + num_bytes;
    if (end_pointer >= mmap_size) {
      page_start += (page_offset & div_value);
      page_offset = (page_offset & mod_value);
      map_file.close();
      map_file.open(path, mmap_size, page_start);
    }
    memcpy(buffer, map_file.data() + page_offset, num_bytes);
    page_offset += num_bytes;
    still_left -= num_bytes;
    return num_bytes;
  }

  inline TLen read_line(char* buffer, TLen num_bytes = Memory_mapped_file::MAX_LINE_LENGTH) {
    if (num_bytes > still_left) num_bytes = still_left;
    TLen end_pointer = page_offset + num_bytes;
    if (end_pointer >= mmap_size) {
      page_start += (page_offset & div_value);
      page_offset = (page_offset & mod_value);
      map_file.close();
      map_file.open(path, mmap_size, page_start);
    }
    memcpy(buffer, map_file.data() + page_offset, num_bytes);
    char* pointBuf = buffer;
    TLen read_bytes = 0; 
    for(; read_bytes<num_bytes; ++read_bytes, ++pointBuf) {
      if (*pointBuf == '\n') {
	*pointBuf = '\0';
	++pointBuf;
	++read_bytes;
	break;
      } 
    }
    page_offset += read_bytes;
    still_left -= read_bytes;
    return read_bytes;
  }

  inline TLen read_line(std::string& line) {
    char buffer[Memory_mapped_file::MAX_LINE_LENGTH];    
    TLen l=read_line(buffer);
    line = buffer; 
    return l; 
  }

  inline TLen left_bytes() const {
    return still_left;
  }

  inline void close() {
    map_file.close();
  }
};

}

#endif
