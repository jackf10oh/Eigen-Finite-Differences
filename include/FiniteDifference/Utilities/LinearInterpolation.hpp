// LinearInterpolation.hpp 
//
// Utility file with some functions around 
// linear interpolation and locating bounding interals 
//
// JAF 4/12/2026 

#ifndef FDM_UTILS_LINEARINTERPOLATION_H
#define FDM_UTILS_LINEARINTERPOLATION_H 

#include<utility> // std::pair
#include<algorithm> // std::lower_bound

namespace fdm{
  namespace utils{

template<typename Iterator>
auto make_subinterval(double x, Iterator start, Iterator stop)
{
  // runtime checks 
  if(std::distance(start,stop) < 2) throw std::runtime_error("size of v < 2"); 
  if(x < (*start)) throw std::runtime_error("c < v[0]"); 

  // right side b in [a,b].
  auto after = std::lower_bound(start, stop, x);
  if(after == stop) throw std::runtime_error("right bound == v.cend()"); 

  // if b == v[0] bump it by 1. 
  auto before = (after==stop) ? after++ : std::prev(after); 
  return std::pair(before, after); 

  // std::lower_bound should seriously be renamed 
  // to reflect the fact it an the supremum 
  // grumble grumble
}; 

template<typename DomainIterator, typename ValueIterator>
double linear_interpolation(
  double x, 
  DomainIterator d_start, DomainIterator d_stop, 
  ValueIterator v_start)
{
  auto bounding_interval =  make_subinterval(x, d_start, d_stop); 

  auto offset = std::distance(d_start,bounding_interval.first); 
  auto it = std::next(v_start, offset); 
  double y1 = *it; 
  it++; 
  double y2 = *it; 

  return y1 + (y2-y1) * (x - *bounding_interval.first) / (*bounding_interval.second - *bounding_interval.first); 
}

template<typename DomainIterator>
double linear_interpolation(
  double x, 
  DomainIterator d_start, DomainIterator d_stop,
  double y1, double y2)
{
  auto bounding_interval =  make_subinterval(x, d_start, d_stop);  

  return y1 + (y2-y1) * (x - *bounding_interval.first) / (*bounding_interval.second - *bounding_interval.first); 
}

  } // end namespace utils 
} // end namespace fdm 

#endif // LinearInterpolation.hpp 