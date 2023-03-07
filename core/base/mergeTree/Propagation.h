#pragma once

#include <boost/heap/fibonacci_heap.hpp>
#include <vector>

namespace ttk::mt {

  template <typename IT>
  struct Propagation {

    // union find members
    Propagation<IT> *parent{this};

    IT branchId{-1};

    // propagation data
    std::vector<IT> criticalPoints;
    boost::heap::fibonacci_heap<std::pair<IT, IT>> queue;
    bool interrupted{false};

    inline Propagation *find() {
      Propagation *p;
#pragma omp atomic read
      p = this->parent;

      if(p == this)
        return this;
      else {
        p = p->find();
#pragma omp atomic write
        this->parent = p;
        return this->parent;
      }
    }

    static inline Propagation<IT> *unify(Propagation<IT> *p0,
                                         Propagation<IT> *p1,
                                         const IT *orderArray,
                                         const bool mergeHeaps = true) {

      if(!p1)
        return p0;

      auto master = p0->find();
      auto slave = p1->find();

      if(master == slave)
        return master;

      // determine master and slave based on rank
      if(orderArray[master->criticalPoints[0]]
         < orderArray[slave->criticalPoints[0]]) {
        auto temp = master;
        master = slave;
        slave = temp;
      }

// update union find tree
#pragma omp atomic write
      slave->parent = master;

      // merge f. heaps
      if(mergeHeaps)
        master->queue.merge(slave->queue);

      return master;
    }
  };
} // namespace ttk::mt