// linear_interpolation.hpp 
//
// Utility file with some functions around 
// linear interpolation and locating bounding interals 
//
// JAF 4/12/2026 

#ifndef FORNFDM_UTILS_LINEARINTERPOLATION_H
#define FORNFDM_UTILS_LINEARINTERPOLATION_H 

#include<cassert>
#include<utility> // std::pair
#include<algorithm> // std::lower_bound
#include<iterator>

namespace fornfdm{
  namespace utils{

template<typename Iterator>
auto make_subinterval(
  typename std::iterator_traits<Iterator>::value_type x, 
  Iterator start, Iterator stop)
{
  // runtime checks 
  assert((std::distance(start,stop) >= 2) && "error: size of v < 2"); 
  assert((x >= (*start)) && "error: c < v[0]"); 

  // right side b in [a,b].
  auto after = std::lower_bound(start, stop, x);
  if(after == stop)
  {
    assert((x == *after) && "right bound != v.cend()");
    return std::pair(std::prev(after), after); 
  }

  if(after == start)
  {
    return std::pair(start, std::next(start)); 
  }
  else
  {
    return std::pair(std::prev(after),after); 
  }

  // std::lower_bound should seriously be renamed 
  // to reflect the fact it an the supremum 
  // grumble grumble
}; 

template<typename DomainIterator, typename ValueIterator>
auto linear_interpolation(
  typename std::iterator_traits<DomainIterator>::value_type x, 
  DomainIterator d_start, DomainIterator d_stop, 
  ValueIterator v_start)
{
  auto bounding_interval =  make_subinterval(x, d_start, d_stop); 

  auto offset = std::distance(d_start,bounding_interval.first); 
  auto it = std::next(v_start, offset); 
  auto y1 = *it; 
  it++; 
  auto y2 = *it; 

  return y1 + (y2-y1) * (x - *bounding_interval.first) / (*bounding_interval.second - *bounding_interval.first); 
}

template<typename DomainIterator, typename ValueType>
ValueType linear_interpolation(
  typename std::iterator_traits<DomainIterator>::value_type x, 
  DomainIterator d_start, DomainIterator d_stop,
  ValueType y1, ValueType y2)
{
  auto bounding_interval =  make_subinterval(x, d_start, d_stop);  

  return y1 + (y2-y1) * (x - *bounding_interval.first) / (*bounding_interval.second - *bounding_interval.first); 
}

  } // end namespace utils 
} // end namespace fornfdm 

#endif // linear_interpolation.hpp 