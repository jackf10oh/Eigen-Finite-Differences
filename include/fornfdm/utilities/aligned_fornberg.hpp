// aligned_fornberg.hpp
//
// same STL interface as fornberg,
// but places the beginning of each order 
// of the weights to begin at fixed strides
//
// JAF 5/8/2026 

#ifndef FORNFDM_UTILS_ALIGNED_FORNBERG_H
#define FORNFDM_UTILS_ALIGNED_FORNBERG_H

#include<cassert>
#include<iterator> 

namespace fornfdm{
namespace utils{

template<class ForwardItIn, class BidItOut>
BidItOut aligned_fornberg(
  ForwardItIn start, ForwardItIn end, 
  const typename std::iterator_traits<ForwardItIn>::value_type& x_bar, 
  std::size_t order,
  std::size_t alignment, 
  BidItOut dest
)
{
  using Scalar = typename std::iterator_traits<ForwardItIn>::value_type; 
  auto num_nodes_used = std::distance(start, end); 

  assert((num_nodes_used > order) && "Fornberg algorithm requires # nodes >= order+1.");
  assert((num_nodes_used <= alignment) && "ALigned Fornberg algorithm requires # nodes < alignment.");

  switch(num_nodes_used)
  {
    case 0:
      // silly edge case of 0 nodes :P 
      return dest;
      break; 
    case 1:
      // silly edge case of 1 node :P 
      *dest = 1.0; 
      return std::next(dest,alignment); 
      break; 

    case 2:
      if(order == 1)
      {
        Scalar x0 = *start; 
        ++start; 
        Scalar x1 = *start; 
        *dest = (x_bar - x1)/(x0-x1); 
        ++dest; 
        *dest = (x_bar - x0)/(x1-x0); 
        std::advance(dest, alignment-1); 
        *dest = 1/(x0-x1); 
        ++dest; 
        *dest = 1/(x1-x0); 
        return std::next(dest,alignment-1); 
      }
      else // order == 0 
      {
        Scalar x0 = *start; 
        ++start; 
        Scalar x1 = *start; 
        *dest = (x_bar - x1)/(x0-x1); 
        ++dest; 
        *dest = (x_bar - x0)/(x1-x0); 
        return std::next(dest,alignment-1); 
      }
      break; 

    default: // num_nodes_used >= 3
      *dest = 1.0; 
      Scalar c1=1.0, c2, c3; 
      std::size_t counter=1;
      // for number of nodes n=2, ..., N (first node was zero index)
      std::advance(dest,alignment); 
      BidItOut write01; 
      auto penultimate_node = start;
      auto node = std::next(start);  
      for(; node != end; ++node)
      {
        c2=1.0; 
        write01 = dest;
        // loop over all previous nodes 0, ..., n-1
        for(auto old_node=start; old_node!=node; ++old_node)
        {
          c3 = *node - *old_node; 
          c2 *= c3; 
          auto write02 = write01; 
          if(counter <= order)
          {
              // Explicitly zero an entry 
              *write02 = 0.0; 
          }
          // if old_n == n-1 we must use the very last old column 
          // to update the newest nodes column
          if(old_node == penultimate_node)
          {
            auto write03 = std::next(write02); 
            auto read01 = write02; 
            auto read02 = std::prev(read01,alignment); 
            for(auto m=std::min(order,counter); m>0; --m)
            {
              // (3.9) -------------------------------------
              *write03 = (c1/c2) * ( *read02 * m - *read01 * (*old_node - x_bar));
              read01 = read02; // jumps read01 back -alignment
              write03 = std::next(read01); // jumps write03 back -alignment 
              std::advance(read02, -alignment); 
            }
            // (3.7) -------------------------------------
            *write03 = *std::prev(write03) * (c1/c2) * (x_bar - *old_node);
          }

          for(auto m=std::min(counter,order); m>0; --m)
          {
            // (3.8) ----------------------------------------
            auto read03 = std::prev(write02, alignment); 
            *write02 = (*write02 * (*node-x_bar) - *read03 * m) / c3;
            write02 = read03; 
          }
          // (3.6) -------------------------------------------------
          *write02 = *write02 * ((*node - x_bar)/(*node - *old_node));
          ++write01; 
        }
        c1 = c2;  
        if(counter<order)
        {
          dest = std::next(write01, alignment - counter); // copy back from write01 uses less increments. 
        } 
        ++counter; 
        ++penultimate_node; 
      }
      return std::next(write01,alignment - num_nodes_used + 1); 
      break;
  }  
}

}
}

#endif 