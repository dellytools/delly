#ifndef VERSION_H
#define VERSION_H

namespace torali
{


  std::string dellyVersionNumber = "1.2.6";

  inline 
    void printTitle(std::string const& title) 
    {
      std::cerr << "**********************************************************************" << std::endl;
      std::cerr << "Program: Delly" << std::endl;
      std::cerr << "This is free software, and you are welcome to redistribute it under" << std::endl;
      std::cerr << "certain conditions (BSD License); for license details use '-l'." << std::endl;
      std::cerr << "This program comes with ABSOLUTELY NO WARRANTY; for details use '-w'." << std::endl;
      std::cerr <<  std::endl;
      std::cerr <<  title << " (Version: " << dellyVersionNumber << ")" << std::endl;
      std::cerr << "Contact: Tobias Rausch (rausch@embl.de)" << std::endl;
      std::cerr << "**********************************************************************" << std::endl;
      std::cerr << std::endl;
    }

  inline
    void displayWarranty()
    {
      std::cerr << "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, ";
      std::cerr << "THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR ";
      std::cerr << "CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, ";
      std::cerr << "PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, ";
      std::cerr << "WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ";
      std::cerr << "ADVISED OF THE POSSIBILITY OF SUCH DAMAGE." << std::endl;
      std::cerr << std::endl;
    }


    inline void
    bsd() {
      std::cerr << "Copyright (c) 2023 European Molecular Biology Laboratory (EMBL)." << std::endl;
      std::cerr << "All rights reserved." << std::endl;
      std::cerr << std::endl;
      std::cerr << "Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:" << std::endl;
      std::cerr << "    1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer." << std::endl;
      std::cerr << "    2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution." << std::endl;
      std::cerr << "    3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission." << std::endl;
      std::cerr << std::endl;
      std::cerr << "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, ";
      std::cerr << "THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR ";
      std::cerr << "CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, ";
      std::cerr << "PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, ";
      std::cerr << "WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ";
      std::cerr << "ADVISED OF THE POSSIBILITY OF SUCH DAMAGE." << std::endl;
      std::cerr << std::endl;
    }


}

#endif
