// fornberg.hpp
//
// Cleaner STL like interface to Fornberg Algorithm 
// that takes forward iterators for input
// and bidirectional iterator for  
//
// JAF 5/8/2026 

#ifndef FORNFDM_UTILS_FORNBERG_H
#define FORNFDM_UTILS_FORNBERG_H

#include<cassert>
#include<iterator> 

namespace fornfdm{
namespace utils{

// ================================================
// Random Access + Random Access 
// ================================================

template<
  class RandIn, 
  class RandDest
>
std::enable_if_t<
    (std::is_base_of_v<std::random_access_iterator_tag, typename std::iterator_traits<RandIn>::iterator_category> &&
    std::is_base_of_v<std::random_access_iterator_tag, typename std::iterator_traits<RandDest>::iterator_category>),
RandDest> fornberg(
  RandIn start, RandIn end, 
  const typename std::iterator_traits<RandIn>::value_type& x_bar, std::size_t order, 
  RandDest dest  
)
{ 
  using Scalar = typename std::iterator_traits<RandIn>::value_type; 
  std::size_t num_nodes_used = std::distance(start, end); 
  assert((num_nodes_used > order) && "Fornberg algorithm requires # nodes >= order+1.");
  assert((order >= 0) && "Fornberg algorithm requires order >= 0.");

  dest[0] = 1.0; 
  Scalar c1=1.0, c2, c3; 
  std::size_t counter=1;
  std::size_t write_idx_01 = num_nodes_used;
  // for number of nodes n=2, ..., N (first node was zero index)      
  for(std::size_t node_idx = 1; node_idx != num_nodes_used; ++node_idx)
  {
    c2=1.0;
    std::size_t penultimate_node_idx = node_idx-1;
    // loop over all previous nodes 0, ..., n-1
    std::size_t old_node_idx=0;
    for(; old_node_idx!=node_idx; ++old_node_idx)
    {
      c3 = start[node_idx] - start[old_node_idx]; 
      c2 *= c3; 
      if(counter <= order)
      {
          // Explicitly zero an entry 
          dest[write_idx_01] = 0.0; 
      }
      // if old_n == n-1 we must use the very last old column 
      // to update the newest nodes column
      if(old_node_idx == penultimate_node_idx)
      {
        std::size_t read_idx_01 = write_idx_01; 
        for(auto m=std::min(order,counter); m>0; --m)
        {
          // (3.9) -------------------------------------
          dest[read_idx_01+1] = (c1/c2) * ( dest[read_idx_01 - num_nodes_used] * m - dest[read_idx_01] * (start[old_node_idx] - x_bar));
          read_idx_01 -= num_nodes_used;
        }
        // (3.7) -------------------------------------
        dest[read_idx_01+1] = dest[read_idx_01] * (c1/c2) * (x_bar - start[old_node_idx]);
      }

      std::size_t write_idx_02 = write_idx_01; 
      for(auto m=std::min(counter,order); m>0; --m)
      {
        // (3.8) ----------------------------------------
        dest[write_idx_02] = (dest[write_idx_02] * (start[node_idx]-x_bar) - dest[write_idx_02 - num_nodes_used] * m) / c3;
        write_idx_02 -= num_nodes_used; 
      }
      // (3.6) -------------------------------------------------
      dest[write_idx_02] = dest[write_idx_02] * ((start[node_idx] - x_bar)/(start[node_idx] - start[old_node_idx]));
      ++write_idx_01;
    }
    c1 = c2; 
    if(counter<order)
    {
      write_idx_01 += (num_nodes_used - counter);
    } 
    else
    {
      write_idx_01 -= counter;
    }
    ++counter; 
  }
  // returns 1 past the end
  return std::next(dest, (order+1)*num_nodes_used); 
}

// ================================================
// Forward Iterator + Bidirectional Iterator
// ================================================

template<
  class ForwardItIn, 
  class BidItOut
>
std::enable_if_t<
    (!std::is_base_of_v<std::random_access_iterator_tag, typename std::iterator_traits<ForwardItIn>::iterator_category> ||
    !std::is_base_of_v<std::random_access_iterator_tag, typename std::iterator_traits<BidItOut>::iterator_category>),
BidItOut> fornberg(
  ForwardItIn start, ForwardItIn end, 
  const typename std::iterator_traits<ForwardItIn>::value_type& x_bar, std::size_t order, 
  BidItOut dest  
)
{
  static_assert(std::is_base_of_v<std::forward_iterator_tag, typename std::iterator_traits<ForwardItIn>::iterator_category>, "input to fornberg algo must be forward iter");
  static_assert(std::is_base_of_v<std::bidirectional_iterator_tag, typename std::iterator_traits<BidItOut>::iterator_category>, "output iter for fornberg algo must be bidirectional");
    
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
            auto read02 = std::prev(read01,num_nodes_used); 
            for(auto m=std::min(order,counter); m>0; --m)
            {
              // (3.9) -------------------------------------
              *write03 = (c1/c2) * ( *read02 * m - *read01 * (*old_node - x_bar));
              read01 = read02; // jumps read01 back -num_nodes_used
              write03 = std::next(read01); // jumps write03 back -num_nodes_used 
              std::advance(read02, -num_nodes_used); 
            }
            // (3.7) -------------------------------------
            *write03 = *std::prev(write03) * (c1/c2) * (x_bar - *old_node);
          }

          for(auto m=std::min(counter,order); m>0; --m)
          {
            auto read03 = std::prev(write02, num_nodes_used); 
            // (3.8) ----------------------------------------
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
          dest = std::next(write01, num_nodes_used - counter); // copy back from write01 uses less increments. 
        } 
        ++counter; 
        ++penultimate_node; 
      }
      return std::next(write01); 
      break;
  }  
}

} // end namespace utils 
} // end namespace fornfdm 

#endif 