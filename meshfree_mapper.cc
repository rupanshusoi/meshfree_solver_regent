#include "mappers/default_mapper.h"
//#include "mappers/logging_wrapper.h"
//#include "realm/logging.h"


#include "meshfree_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

void default_policy_select_target_processors(
                                    MapperContext ctx,
                                    const Task &task,
                                    std::vector<Processor> &target_procs)
{
  target_procs.push_back(task.target_proc);
}

static void create_mappers(Machine machine,
                           Runtime* rt,
                           const std::set<Processor>& local_procs) {
  for (Processor proc : local_procs) {
    rt->replace_default_mapper(new DefaultMapper(rt->get_mapper_runtime(), machine, proc), proc);
  }
}

void register_mappers() {
  Runtime::add_registration_callback(create_mappers);
}
