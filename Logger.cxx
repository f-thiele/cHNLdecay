#include "Logger.h"

Logger& Debug() {
  static Logger logger(Level::DEBUG);
  return logger;
};
Logger& Info() {
  static Logger logger(Level::INFO);
  return logger;
};
Logger& Warning() {
  static Logger logger(Level::WARNING);
  return logger;
};
Logger& Error() {
  static Logger logger(Level::ERROR);
  return logger;
};
