/// \ingroup baseCode
/// \class ttk::AtomicUF
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa vtkContourForests.cpp %for a usage example.

#ifndef ATOMICUF_H
#define ATOMICUF_H

#include <memory>
#include <vector>

#include "DataTypes.h"
#include "Structures.h"

namespace ttk
{
namespace ftm
{

   class AtomicUF
   {
     private:
      unsigned   rank_;
      AtomicUF * parent_;
      SharedData data_;

     public:
      inline explicit AtomicUF(idVertex extrema = nullVertex) : rank_(0), data_(extrema)
      {
         parent_ = this;
      }

      // heavy recursif
      inline AtomicUF *find()
      {
         if (parent_ == this)
            return this;
         else {
            decltype(parent_) tmp = parent_->find();
#pragma omp atomic write
            parent_ = tmp;

            return parent_;
         }
      }

      // Shared data get/set

      inline idVertex getExtrema(void) const
      {
         return data_.extrema;
      }

      inline CurrentState* getFirstState(void)
      {
#ifndef withKamikaze
          if(data_.states.size() != 1){
             std::cout << "AtomicUF :: getFirstState : nb state 1 !=  " << data_.states.size()
                       << std::endl;
          }
#endif
          return data_.states[0];
      }

      inline AtomicVector<CurrentState *> &getStates(void)
      {
         return data_.states;
      }

      inline size_t getNbStates(void) const
      {
         return data_.states.size();
      }

      inline AtomicVector<idSuperArc> &getOpenedArcs(void)
      {
         return data_.openedArcs;
      }

      inline void clearOpenedArcs(void)
      {
         data_.openedArcs.reset();
      }

      inline void setExtrema(idVertex v)
      {
         data_.extrema = v;
      }

      inline void reserveData(const size_t& s)
      {
          data_.reserve(s);
      }

      inline void addState(CurrentState *s)
      {
         data_.addState(s);
      }

      inline void addArcToClose(idSuperArc a)
      {
         data_.addArc(a);
      }

      inline CurrentState *mergeStates(void)
      {
         CurrentState *s       = data_.states[0];
         const auto &  nbState = data_.states.size();

         for (valence i = 1; i < nbState; ++i) {
            s->merge(*data_.states[i]);
            delete data_.states[i];
         }

         data_.states.reset(1);
         return s;
      }

      // UF get / set

      inline int getRank() const
      {
         return rank_;
      }

      inline void setRank(const int &rank)
      {
#pragma omp atomic write
         rank_ = rank;
      }

      inline void setParent(AtomicUF *parent)
      {
#pragma omp atomic write
         parent_ = parent;
      }

      static inline AtomicUF *makeUnion(AtomicUF *uf0, AtomicUF *uf1)
      {
         uf0 = uf0->find();
         uf1 = uf1->find();

         if (uf0 == uf1) {
            return uf0;
         } else if (uf0->getRank() > uf1->getRank()) {
            uf1->setParent(uf0);
            uf0->data_.merge(uf1->data_);
            return uf0;
         } else if (uf0->getRank() < uf1->getRank()) {
            uf0->setParent(uf1);
            uf1->data_.merge(uf0->data_);
            return uf1;
         } else {
            uf1->setParent(uf0);
            uf0->setRank(uf0->getRank() + 1);
            uf0->data_.merge(uf1->data_);
            return uf0;
         }

         return NULL;
      }

      static inline AtomicUF *makeUnion(std::vector<AtomicUF *> &sets)
      {
         AtomicUF *n = NULL;

         if (!sets.size())
            return NULL;

         if (sets.size() == 1)
            return sets[0];

         for (int i = 0; i < (int)sets.size() - 1; i++)
            n = makeUnion(sets[i], sets[i + 1]);

         return n;
      }

      inline bool operator<(const AtomicUF &other) const
      {
         return rank_ < other.rank_;
      }

      inline bool operator>(const AtomicUF &other) const
      {
         return rank_ > other.rank_;
      }
   };

}
}

#endif /* end of include guard: ATOMICUF_H */

