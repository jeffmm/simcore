#include "simcore/logger.hpp"

/****************************/
/******** SINGLETON *********/
/****************************/

LoggerService::LoggerService() {}
LoggerService::~LoggerService() {}
FILE *LoggerService::log_file_ = nullptr;

void LoggerService::SetOutput(const char *fname) const {
  LoggerService::log_file_ = fopen(fname, "w");
}

void LoggerService::Info(const char *msg, va_list args) const {
  WriteMsg("INFO ", msg, args);
}
void LoggerService::Warning(const char *msg, va_list args) const {
  WriteMsg("WARN ", msg, args);
}

void LoggerService::Debug(const char *msg, va_list args) const {
#if defined(DEBUG) || defined(TRACE)
  WriteMsg("DEBUG", msg, args);
#endif
}
void LoggerService::Trace(const char *msg, va_list args) const {
#ifdef TRACE
  WriteMsg("TRACE", msg, args);
#endif
}

void LoggerService::Error(const char *msg, va_list args) const {
  WriteMsg("ERROR", msg, args);
}

void LoggerService::WriteMsg(const char *severity, const char *msg, va_list args) const {
  std::lock_guard<std::mutex> lk(mtx_);

  // Create time string.
  char timestr[64];
  time_t t = time(NULL);
  struct tm *p = localtime(&t);
  strftime(timestr, 64, "%F::%T", p);

  // Log time, severity, and message to log file
  if (log_file_ != nullptr) {
    // Copy args for log file
    va_list arg_cpy;
    va_copy(arg_cpy, args);

    fprintf(log_file_, "%s", timestr);
    fprintf(log_file_, " - [%s] - ", severity);
    vfprintf(log_file_, msg, arg_cpy);
    fprintf(log_file_, "\n");
    fflush(log_file_);
    va_end(arg_cpy);
  }

  // Then to stderr
  fprintf(stderr, "%s", timestr);
  fprintf(stderr, " - [%s] - ", severity);
  vfprintf(stderr, msg, args);
  fprintf(stderr, "\n");
}

// Get logger instance
const LoggerService &LoggerService::Get() {
  static LoggerService logger;
  return logger;
}


/****************************/
/******** INTERFACE *********/
/****************************/

void Logger::SetOutput(const char *fname) {
  LoggerService::Get().SetOutput(fname);
}

void Logger::Trace(const char *msg, ...) {
  va_list args;
  va_start(args, msg);
  LoggerService::Get().Trace(msg, args);
  va_end(args);
}

void Logger::Debug(const char *msg, ...) {
  va_list args;
  va_start(args, msg);
  LoggerService::Get().Debug(msg, args);
  va_end(args);
}

void Logger::Info(const char *msg, ...) {
  va_list args;
  va_start(args, msg);
  LoggerService::Get().Info(msg, args);
  va_end(args);
}

void Logger::Warning(const char *msg, ...) {
  va_list args;
  va_start(args, msg);
  LoggerService::Get().Warning(msg, args);
  va_end(args);
}

void Logger::Error(const char *msg, ...) {
  va_list args;
  va_start(args, msg);
  LoggerService::Get().Error(msg, args);
  va_end(args);
  // Kill the program on error
  exit(1);
}
