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

#ifndef RECORD_H
#define RECORD_H

#include <iostream>

#include "api/BamReader.h"
#include "memory_mapped_file.h"
#include "tokenizer.h"

namespace torali
{
  
  
  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    struct Record {
      T0 f0;
      T1 f1;
      T2 f2;
      T3 f3;
      T4 f4;
      T5 f5;
      T6 f6;
      T7 f7;
      T8 f8;
      T9 f9;
      T10 f10;
      T11 f11;
    };


  template<typename T0, typename T1, typename T2, typename T3, typename T6, typename T7, typename T8, typename T9, typename T10>
    struct Record<T0, T1, T2, T3, void, void, T6, T7, T8, T9, T10, void> {
      T0 f0;
      T1 f1;
      T2 f2;
      T3 f3;
      T6 f6;
      T7 f7;
      T8 f8;
      T9 f9;
      T10 f10;
    };

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T6, typename T7, typename T8, typename T9, typename T10>
    struct Record<T0, T1, T2, T3, T4, void, T6, T7, T8, T9, T10, void> {
      T0 f0;
      T1 f1;
      T2 f2;
      T3 f3;
      T4 f4;
      T6 f6;
      T7 f7;
      T8 f8;
      T9 f9;
      T10 f10;
    };


  template<typename T0, typename T1, typename T2, typename T3, typename T6, typename T7, typename T8, typename T9>
    struct Record<T0, T1, T2, T3, void, void, T6, T7, T8, T9, void, void> {
      T0 f0;
      T1 f1;
      T2 f2;
      T3 f3;
      T6 f6;
      T7 f7;
      T8 f8;
      T9 f9;
    };

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
    struct Record<T0, T1, T2, T3, T4, T5, T6, T7, T8, void, void, void> {
      T0 f0;
      T1 f1;
      T2 f2;
      T3 f3;
      T4 f4;
      T5 f5;
      T6 f6;
      T7 f7;
      T8 f8;
    };

  template<typename T1, typename T2, typename T3, typename T6, typename T7, typename T8, typename T9>
    struct Record<void, T1, T2, T3, void, void, T6, T7, T8, T9, void, void> {
      T1 f1;
      T2 f2;
      T3 f3;
      T6 f6;
      T7 f7;
      T8 f8;
      T9 f9;
    };

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
    struct Record<T0, T1, T2, T3, T4, T5, T6, void, void, void, void, void> {
      T0 f0;
      T1 f1;
      T2 f2;
      T3 f3;
      T4 f4;
      T5 f5;
      T6 f6;
    };

  template<typename T1, typename T2, typename T3, typename T4, typename T5>
    struct Record<void, T1, T2, T3, T4, T5, void, void, void, void, void, void> {
      T1 f1;
      T2 f2;
      T3 f3;
      T4 f4;
      T5 f5;
    };

  template<typename T0, typename T1, typename T2, typename T3, typename T4>
    struct Record<T0, T1, T2, T3, T4, void, void, void, void, void, void, void> {
      T0 f0;
      T1 f1;
      T2 f2;
      T3 f3;
      T4 f4;
    };

  template<typename T1, typename T2, typename T3, typename T4>
    struct Record<void, T1, T2, T3, T4, void, void, void, void, void, void, void> {
      T1 f1;
      T2 f2;
      T3 f3;
      T4 f4;
    };


  template<typename T1,typename T3,  typename T6, typename T7, typename T8>
    struct Record<void, T1, void, T3, void, void, T6, T7, T8, void, void, void> {
      T1 f1;
      T3 f3;
      T6 f6;
      T7 f7;
      T8 f8;
    };


  template<typename T0, typename T1,  typename T2, typename T3, typename T6, typename T7, typename T9>
    struct Record<T0, T1, T2, T3, void, void, T6, T7, void, T9, void, void> {
      T0 f0;
      T1 f1;
      T2 f2;
      T3 f3;
      T6 f6;
      T7 f7;
      T9 f9;
    };


  template<typename T0,typename T1,  typename T2, typename T4>
    struct Record<T0, T1, T2, void, T4, void, void, void, void, void, void, void> {
      T0 f0;
      T1 f1;
      T2 f2;
      T4 f4;
    };

  template<typename T0,typename T1,  typename T2, typename T3>
    struct Record<T0, T1, T2, T3, void, void, void, void, void, void, void, void> {
      T0 f0;
      T1 f1;
      T2 f2;
      T3 f3;
    };

  template<typename T0,typename T1,  typename T2>
    struct Record<T0, T1, T2, void, void, void, void, void, void, void, void, void> {
      T0 f0;
      T1 f1;
      T2 f2;
    };


  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF0(Tokenizer& token, Record<std::string, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    token.getString(rec.f0);
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF0(Tokenizer& token, Record<void, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>&) 
  {
    token.skipNextChunk();
  }

  template<typename T0, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF1(Tokenizer& token, Record<T0, unsigned short, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    rec.f1=token.getUShort();
  }

  template<typename T0, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF1(Tokenizer& token, Record<T0, unsigned int, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    rec.f1=token.getUInt();
  }

  template<typename T0, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF1(Tokenizer& token, Record<T0, void, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>&) 
  {
    token.skipNextChunk();
  }


  template<typename T0, typename T1, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF2(Tokenizer& token, Record<T0, T1, std::string, T3, T4, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    token.getString(rec.f2);
  }

  template<typename T0, typename T1, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF2(Tokenizer& token, Record<T0, T1, unsigned int, T3, T4, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    rec.f2=token.getUInt();
  }

  template<typename T0, typename T1, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF2(Tokenizer& token, Record<T0, T1, void, T3, T4, T5, T6, T7, T8, T9, T10, T11>&) 
  {
    token.skipNextChunk();
  }

  template<typename T0, typename T1, typename T2, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF3(Tokenizer& token, Record<T0, T1, T2, unsigned int, T4, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    rec.f3=token.getUInt();
  }

  template<typename T0, typename T1, typename T2, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF3(Tokenizer& token, Record<T0, T1, T2, int, T4, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    rec.f3=token.getInt();
  }

  template<typename T0, typename T1, typename T2, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF3(Tokenizer& token, Record<T0, T1, T2, std::string, T4, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    token.getString(rec.f3);
  }

  template<typename T0, typename T1, typename T2, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF3(Tokenizer& token, Record<T0, T1, T2, void, T4, T5, T6, T7, T8, T9, T10, T11>&) 
  {
    token.skipNextChunk();
  }


  template<typename T0, typename T1, typename T2, typename T3, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF4(Tokenizer& token, Record<T0, T1, T2, T3, unsigned short, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    rec.f4=token.getUShort();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF4(Tokenizer& token, Record<T0, T1, T2, T3, unsigned int, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    rec.f4=token.getUInt();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF4(Tokenizer& token, Record<T0, T1, T2, T3, int, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    rec.f4=token.getInt();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF4(Tokenizer& token, Record<T0, T1, T2, T3, std::string, T5, T6, T7, T8, T9, T10, T11>& rec) 
  {
    token.getString(rec.f4);
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF4(Tokenizer& token, Record<T0, T1, T2, T3, void, T5, T6, T7, T8, T9, T10, T11>&) 
  {
    token.skipNextChunk();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF5(Tokenizer& token, Record<T0, T1, T2, T3, T4, std::string, T6, T7, T8, T9, T10, T11>& rec) 
  {
    token.getString(rec.f5);
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF5(Tokenizer& token, Record<T0, T1, T2, T3, T4, int, T6, T7, T8, T9, T10, T11>& rec) 
  {
    rec.f5=token.getInt();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF5(Tokenizer& token, Record<T0, T1, T2, T3, T4, double, T6, T7, T8, T9, T10, T11>& rec) 
  {
    rec.f5=token.getDouble();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF5(Tokenizer& token, Record<T0, T1, T2, T3, T4, void, T6, T7, T8, T9, T10, T11>&) 
  {
    token.skipNextChunk();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF6(Tokenizer& token, Record<T0, T1, T2, T3, T4, T5, std::string, T7, T8, T9, T10, T11>& rec) 
  {
    token.getString(rec.f6);
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF6(Tokenizer& token, Record<T0, T1, T2, T3, T4, T5, int, T7, T8, T9, T10, T11>& rec) 
  {
    rec.f6=token.getInt();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T7, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF6(Tokenizer& token, Record<T0, T1, T2, T3, T4, T5, void, T7, T8, T9, T10, T11>&) 
  {
    token.skipNextChunk();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF7(Tokenizer& token, Record<T0, T1, T2, T3, T4, T5, T6, unsigned int, T8, T9, T10, T11>& rec) 
  {
    rec.f7=token.getUInt();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T8, typename T9, typename T10, typename T11>
    inline 
    void addF7(Tokenizer& token, Record<T0, T1, T2, T3, T4, T5, T6, void, T8, T9, T10, T11>&) 
  {
    token.skipNextChunk();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T9, typename T10, typename T11>
    inline 
    void addF8(Tokenizer& token, Record<T0, T1, T2, T3, T4, T5, T6, T7, int, T9, T10, T11>& rec) 
  {
    rec.f8=token.getInt();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T9, typename T10, typename T11>
    inline 
    void addF8(Tokenizer& token, Record<T0, T1, T2, T3, T4, T5, T6, T7, unsigned int, T9, T10, T11>& rec) 
  {
    rec.f8=token.getUInt();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T9, typename T10, typename T11>
    inline 
    void addF8(Tokenizer& token, Record<T0, T1, T2, T3, T4, T5, T6, T7, void, T9, T10, T11>&) 
  {
    token.skipNextChunk();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T10, typename T11>
    inline 
    void addF9(Tokenizer& token, Record<T0, T1, T2, T3, T4, T5, T6, T7, T8, unsigned int, T10, T11>& rec) 
  {
    rec.f9=token._getLenNextChunk();
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T10, typename T11>
    inline 
    void addF9(Tokenizer& token, Record<T0, T1, T2, T3, T4, T5, T6, T7, T8, std::string, T10, T11>& rec) 
  {
    token.getString(rec.f9);
  }

  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T10, typename T11>
    inline 
    void addF9(Tokenizer& token, Record<T0, T1, T2, T3, T4, T5, T6, T7, T8, void, T10, T11>&) 
  {
    token.skipNextChunk();
  }

}

#endif
