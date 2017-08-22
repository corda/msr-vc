#pragma once

enum verbosity { SILENT = 0, WARNING, INFO, DEBUG };
#define WARN(...)  if (Config::verbose_level >= WARNING) { printf(__VA_ARGS__); }
#define INFO(...)  if (Config::verbose_level >= INFO)    { printf(__VA_ARGS__); }
#define DEBUG(...) if (Config::verbose_level >= DEBUG)   { printf(__VA_ARGS__); }

class Config {
public:
  static int max_mem;
  static bool precomp_free;
  static int verbose_level;
  static bool machine_readable;
};