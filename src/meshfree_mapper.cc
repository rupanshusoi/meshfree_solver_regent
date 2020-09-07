#include "mappers/default_mapper.h"
//#include "mappers/logging_wrapper.h"
//#include "realm/logging.h"


#include "meshfree_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

class MeshfreeMapper : public DefaultMapper
{
public:
  MeshfreeMapper(MapperRuntime *rt, Machine machine, Processor local);
  virtual void default_policy_select_target_processors(MapperContext ctx, const Task &task, std::vector<Processor> &target_procs);
  virtual Memory default_policy_select_target_memory(MapperContext ctx, Processor target_proc, const RegionRequirement &req);
};

MeshfreeMapper::MeshfreeMapper(MapperRuntime *rt, Machine machine, Processor local)
  : DefaultMapper(rt, machine, local)
{
}

Memory MeshfreeMapper::default_policy_select_target_memory(MapperContext ctx,
                                               Processor target_proc,
                                               const RegionRequirement &req)
{
  bool prefer_rdma = ((req.tag & DefaultMapper::PREFER_RDMA_MEMORY) != 0);

  // TODO: deal with the updates in machine model which will
  //       invalidate this cache
  std::map<Processor,Memory>::iterator it;
  if (prefer_rdma)
  {
    it = cached_rdma_target_memory.find(target_proc);
    if (it != cached_rdma_target_memory.end()) return it->second;
  } else {
    it = cached_target_memory.find(target_proc);
    if (it != cached_target_memory.end()) return it->second;
  }

  // Find the visible memories from the processor for the given kind
  Machine::MemoryQuery visible_memories(machine);
  visible_memories.has_affinity_to(target_proc);
  if (visible_memories.count() == 0)
  {
    //log_mapper.error("No visible memories from processor " IDFMT "! "
    //                 "This machine is really messed up!", target_proc.id);
    assert(false);
  }
  // Use zero copy memory for every processor
  unsigned best_bandwidth = 0;
  Memory best_memory = Memory::NO_MEMORY;
  std::vector<Machine::ProcessorMemoryAffinity> affinity(1);
  for (Machine::MemoryQuery::iterator it = visible_memories.begin();
        it != visible_memories.end(); it++)
  {
    affinity.clear();
    machine.get_proc_mem_affinity(affinity, target_proc, *it,
    			      false /*not just local affinities*/);
    assert(affinity.size() == 1);
    if (affinity[0].m.kind() == Memory::Z_COPY_MEM && affinity[0].bandwidth > best_bandwidth) {
      best_memory = *it;
      best_bandwidth = affinity[0].bandwidth;
    }    
  }
  assert(best_memory.exists());

  cached_target_memory[target_proc] = best_memory;
  return best_memory;
}

void MeshfreeMapper::default_policy_select_target_processors(
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
    rt->replace_default_mapper(new MeshfreeMapper(rt->get_mapper_runtime(), machine, proc), proc);
  }
}

class LinearShardingFunctor : public ShardingFunctor {
public:
  LinearShardingFunctor();
  LinearShardingFunctor(const LinearShardingFunctor &rhs);
  virtual ~LinearShardingFunctor(void);
public:
  LinearShardingFunctor& operator=(const LinearShardingFunctor &rhs);
public:
  template<int DIM>
  size_t linearize_point(const Realm::IndexSpace<DIM,coord_t> &is,
                         const Realm::Point<DIM,coord_t> &point) const;
public:
  virtual ShardID shard(const DomainPoint &point,
                        const Domain &full_space,
                        const size_t total_shards);
};

//--------------------------------------------------------------------------
LinearShardingFunctor::LinearShardingFunctor()
  : ShardingFunctor()
//--------------------------------------------------------------------------
{
}

//--------------------------------------------------------------------------
LinearShardingFunctor::LinearShardingFunctor(
                                           const LinearShardingFunctor &rhs)
  : ShardingFunctor()
//--------------------------------------------------------------------------
{
  // should never be called
  assert(false);
}

//--------------------------------------------------------------------------
LinearShardingFunctor::~LinearShardingFunctor(void)
//--------------------------------------------------------------------------
{
}

//--------------------------------------------------------------------------
LinearShardingFunctor& LinearShardingFunctor::operator=(
                                           const LinearShardingFunctor &rhs)
//--------------------------------------------------------------------------
{
  // should never be called
  assert(false);
  return *this;
}

//--------------------------------------------------------------------------
template<int DIM>
size_t LinearShardingFunctor::linearize_point(
                               const Realm::IndexSpace<DIM,coord_t> &is,
                               const Realm::Point<DIM,coord_t> &point) const
//--------------------------------------------------------------------------
{
  Realm::AffineLinearizedIndexSpace<DIM,coord_t> linearizer(is);
  return linearizer.linearize(point);
}

//--------------------------------------------------------------------------
ShardID LinearShardingFunctor::shard(const DomainPoint &point,
                                     const Domain &full_space,
                                     const size_t total_shards)
//--------------------------------------------------------------------------
{
#ifdef DEBUG_LEGION
  assert(point.get_dim() == full_space.get_dim());
#endif
  size_t domain_size = full_space.get_volume();
  switch (point.get_dim())
  {
    case 1:
      {
        const DomainT<1,coord_t> is = full_space;
        const Point<1,coord_t> p1 = point;
        return linearize_point<1>(is, p1)  * total_shards / domain_size;
      }
    case 2:
      {
        const DomainT<2,coord_t> is = full_space;
        const Point<2,coord_t> p2 = point;
        return linearize_point<2>(is, p2)  * total_shards / domain_size;
      }
    case 3:
      {
        const DomainT<3,coord_t> is = full_space;
        const Point<3,coord_t> p3 = point;
        return linearize_point<3>(is, p3)  * total_shards / domain_size;
      }
    default:
      assert(false);
  }
  return 0;
}

void register_mappers() {
  LinearShardingFunctor *sharding_functor = new LinearShardingFunctor();

  Runtime::preregister_sharding_functor(1, sharding_functor);
  Runtime::add_registration_callback(create_mappers);
}
