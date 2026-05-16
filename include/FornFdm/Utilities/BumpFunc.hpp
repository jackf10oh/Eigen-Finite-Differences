// BumpFunc.hpp
//
// Bump function on [L,R] with center c and height h
//
// JAF 1/13/2026

#ifndef FORNFDM_UTILS_BUMPFUNC_H
#define FORNFDM_UTILS_BUMPFUNC_H 

namespace fornfdm{
  namespace utils{

struct BumpFunc
{
    fornfdm::Scalar L = 0.0;
    fornfdm::Scalar R = 1.0; 
    fornfdm::Scalar c = 0.5;
    fornfdm::Scalar h = 1.0; 
    fornfdm::Real focus = 1.0; 
    
    // Member Funcs 
    fornfdm::Scalar operator()(fornfdm::Scalar x) const 
    {
      if( c<L || c>R) throw std::runtime_error("Error: Bad args to BumpFunc operator(). center outside [L,R]");
      // if(focus < 1.0) throw std::runtime_error("Error: Bar args to BumpFunc operator(). focus must be >= 1"); 
      if(x <= L)
      {
        return 0.0; 
      }
      else if(x < c)
      {
        // cubic spline formula 1 + (x-c)/(c-L) + a*(x-c)^2*(x-L) + b*(x-c)*(x-L)^2
        // a and b can be derived by hand ... 
        fornfdm::Scalar a = 1 / ((L-c)*(L-c)*(L-c)); 
        fornfdm::Scalar b = 1 / ((L-c)*(c-L)*(c-L)); 
        fornfdm::Scalar prefocus = 1 + (x-c)/(c-L) + a*(x-c)*(x-c)*(x-L) + b*(x-c)*(x-L)*(x-L); 
        return h * std::pow(prefocus, focus); 
      }
      else if(x <= R)
      {
        fornfdm::Scalar a = 1 / ((R-c)*(R-c)*(R-c)); 
        fornfdm::Scalar b = 1 / ((R-c)*(c-R)*(c-R)); 
        fornfdm::Scalar prefocus = 1 + (x-c)/(c-R) + a*(x-c)*(x-c)*(x-R) + b*(x-c)*(x-R)*(x-R); 
        return h * std::pow(prefocus, focus); 
      }
      else
      {
        return 0.0; 
      }
    }
};

  } // end namespace utils 
} // end namespace fornfdm 

#endif // BumpFunc.hpp