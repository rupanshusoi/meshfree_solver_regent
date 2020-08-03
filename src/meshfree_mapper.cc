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
};

MeshfreeMapper::MeshfreeMapper(MapperRuntime *rt, Machine machine, Processor local)
  : DefaultMapper(rt, machine, local)
{
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
