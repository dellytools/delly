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

#ifndef SIGNALHANDLER_H
#define SIGNALHANDLER_H

#include <iostream>
#include <stdexcept>
#include <signal.h>
#include <errno.h>

using std::runtime_error;


namespace torali {

  std::string finalDeletionCmd;

  class SignalHandler {
  public:
    SignalHandler() {}

    ~SignalHandler() {}
    
    static void exitSignalHandler(int) {
      std::cerr << "Program performs clean-up and exits!" << std::endl;
      std::cerr << finalDeletionCmd.c_str() << std::endl;
      if (std::system(finalDeletionCmd.c_str())) {
	exit(EXIT_FAILURE);
      }
      exit(-1);
    }

    void setupSignalHandlers(std::string& cmd) {
      if ((signal((int) SIGINT, SignalHandler::exitSignalHandler) == SIG_ERR) || (signal((int) SIGTERM, SignalHandler::exitSignalHandler) == SIG_ERR)) {
	std::cerr << "!!! Error setting up signal handlers !!!" << std::endl;
	exit(-1);
      } else {
	finalDeletionCmd = cmd;
      }
    }

  };

}

#endif

