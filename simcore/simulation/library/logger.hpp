#ifndef _SIMCORE_LOGGER_H_
#define _SIMCORE_LOGGER_H_

#include <csignal>
#include <cstdio>
#include <ctime>
#include <mutex>

// Logger singleton
class LoggerService {
private:
  void WriteMsg(const char *severity, const char *msg, va_list args) const;
  // Log file
  static FILE *log_file_;
  // Make logging thread-safe
  mutable std::mutex mtx_;

public:
  LoggerService();
  ~LoggerService();
  // No copies allowed
  LoggerService(LoggerService const &that) = delete;
  LoggerService &operator=(LoggerService const &that) = delete;

  // Initialize log file name
  void SetOutput(const char *fname) const;
  // Logging functions
  void Trace(const char *msg, va_list args) const;
  void Debug(const char *msg, va_list args) const;
  void Info(const char *msg, va_list args) const;
  void Warning(const char *msg, va_list args) const;
  void Error(const char *msg, va_list args) const;

  // Get logging instance
  static const LoggerService &Get();
};


// Logger interface
class Logger {
public:
  // Initialize log file name
  static void SetOutput(const char *fname);
  // Logging functions.
  static void Trace(const char *msg, ...);
  static void Debug(const char *msg, ...);
  static void Info(const char *msg, ...);
  static void Warning(const char *msg, ...);
  static void Error(const char *msg, ...);
};

#endif
