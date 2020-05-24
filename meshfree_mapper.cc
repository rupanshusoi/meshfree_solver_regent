#include "mappers/default_mapper.h"
#include "mappers/logging_wrapper.h"
#include "realm/logging.h"


#include "meshfree_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

static void create_mappers(Machine machine,
                           Runtime* rt,
                           const std::set<Processor>& local_procs) {
  for (Processor proc : local_procs) {
    rt->replace_default_mapper(new LoggingWrapper(new DefaultMapper(rt->get_mapper_runtime(), machine, proc)), proc);
  }
}

void register_mappers() {
  Runtime::add_registration_callback(create_mappers);
}
