// Fornberg2.hpp
//
// Cleaner STL like interface to Fornberg Algorithm 
// that takes forward iterators for input
// and bidirectional iterator for  
//
// JAF 5/8/2026 

#ifndef FDM_FORNBERG2_H
#define FDM_FORNBERG2_H 

#include<cassert>
#include<iterator> 

namespace fdm{
namespace utils{

template<class ForwardItIn, class BidItOut>
BidItOut forberg2(
  ForwardItIn start, ForwardItIn end, 
  const typename std::iterator_traits<ForwardItIn>::value_type& x_bar, std::size_t order, 
  BidItOut dest  
)
{
  using Scalar = typename std::iterator_traits<ForwardItIn>::value_type; 
  auto num_nodes_used = std::distance(start, end); 

  assert((num_nodes_used > order) && "Fornberg algorithm requires # nodes >= order+1.");

  switch(num_nodes_used)
  {
    case 0:
      // silly edge case of 0 nodes :P 
      return dest;
      break; 
    case 1:
      // silly edge case of 1 node :P 
      *dest = 1.0; 
      return std::next(dest); 
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
        ++dest; 
        *dest = 1/(x0-x1); 
        ++dest; 
        *dest = 1/(x1-x0); 
        return std::next(dest); 
      }
      else // order == 0 
      {
        Scalar x0 = *start; 
        ++start; 
        Scalar x1 = *start; 
        *dest = (x_bar - x1)/(x0-x1); 
        ++dest; 
        *dest = (x_bar - x0)/(x1-x0); 
        return std::next(dest); 
      }
      break; 

    default: // num_nodes_used >= 3
      // for(auto it=dest; it!=std::next(dest,num_nodes_used*(order+1)); ++it) *it = 0.0; 
      *dest = 1.0; 
      Scalar c1=1.0, c2, c3; 
      std::size_t counter=1;
      // for number of nodes n=2, ..., N (first node was zero index)
      std::advance(dest,num_nodes_used); 
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
          // if old_n == n-1 we must use the very last old column 
          // to update the newest nodes column
          if(old_node == penultimate_node)
          {
            auto write03 = std::next(write02); 
            auto write04 = std::next(write03,num_nodes_used); 
            auto read = write02; 
            for(auto m=std::min(order,counter); m>0; --m)
            {
              // (3.9) -------------------------------------
              *write03 = (c1/c2) * ( *std::prev(read,num_nodes_used) * m - *read * (*old_node - x_bar));
              if(counter+1 < order)
              {
                *write04 = 0.0; 
                write04 = write03; 
              }  
              std::advance(write03, -num_nodes_used); 
              read = std::prev(write03);
            }
            // (3.7) -------------------------------------
            *write03 = *std::prev(write03) * (c1/c2) * (x_bar - *old_node);
            if(counter+1 < order)
            {
              *write04 =0.0;
            } 
          }
          for(auto m=std::min(counter,order); m>0; --m)
          {
            // (3.8) ----------------------------------------
            *write02 = (*write02 * (*node-x_bar) - *std::prev(write02,num_nodes_used) * m) / c3;
            std::advance(write02, -num_nodes_used); 
          }
          // (3.6) -------------------------------------------------
          *write02 = *write02 * ((*node - x_bar)/(*node - *old_node));
          ++write01; 
        }
        c1 = c2; 
        std::advance(dest, (counter<order) ? num_nodes_used : 0); 
        ++counter; 
        ++penultimate_node; 
      }
      return std::next(write01); 
      break;
  }  
}

}
}

#endif 