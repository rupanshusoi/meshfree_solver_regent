#include <array>
#include <deque>
#include <iostream>
#include <fstream>
#include <regex>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mappers/default_mapper.h"
#include "mappers/logging_wrapper.h"

using namespace Legion;
using namespace Legion::Mapping;

static void create_mappers(Machine machine,
                           Runtime* rt,
                           const std::set<Processor>& local_procs) {
  for (Processor proc : local_procs) {
    rt->replace_default_mapper(new LoggingWrapper(new DefaultMapper(rt, machine, proc)));
  }
}

void register_mappers() {
  Runtime::add_registration_callback(create_mappers);
}
