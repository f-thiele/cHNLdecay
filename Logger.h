// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <sstream>
#include <string>

// LogSink is a backend consuming preformatted messages
// there can be several different instances depending on where
// to send the data
enum class Level : int { DEBUG, INFO, WARNING, ERROR };
extern Level gLOGLEVEL;

class Logger {
public:
  Logger(Level l) { _level = l; }; // Level l, LogSink& ls);

  void operator()(std::string const &message, char const *function,
                  char const *file, int line) {
    if (_level >= gLOGLEVEL)
      std::cout << "[" << file << ":" << line << "] (" << function << ") -- "
                << message << std::endl;
  };

private:
  Level _level;
  // LogSink& _sink;
};
#define LOG(Logger_, Message_)                                                 \
  Logger_(static_cast<std::ostringstream &>(std::ostringstream().flush()       \
                                            << Message_)                       \
              .str(),                                                          \
          __FUNCTION__, __FILE__, __LINE__);

Logger &Debug();
Logger &Info();
Logger &Warning();
Logger &Error();
#define LOG_DEBUG(Message_) LOG(Debug(), Message_);
#define LOG_INFO(Message_) LOG(Info(), Message_);
#define LOG_WARNING(Message_) LOG(Warning(), Message_);
#define LOG_ERROR(Message_) LOG(Error(), Message_);

#endif // end of header definition
