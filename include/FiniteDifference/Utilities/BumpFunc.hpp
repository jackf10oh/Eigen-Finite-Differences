// BumpFunc.hpp
//
// Bump function on [L,R] with center c and height h
//
// JAF 1/13/2026

#ifndef BUMPFUNC_H
#define BUMPFUNC_H 

struct BumpFunc
{
    double L = 0.0;
    double R = 1.0; 
    double c = 0.5;
    double h = 1.0; 
    double focus = 1.0; 
    
    // Member Funcs 
    double operator()(double x) const 
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
        double a = 1 / ((L-c)*(L-c)*(L-c)); 
        double b = 1 / ((L-c)*(c-L)*(c-L)); 
        double prefocus = 1 + (x-c)/(c-L) + a*(x-c)*(x-c)*(x-L) + b*(x-c)*(x-L)*(x-L); 
        return h * std::pow(prefocus, focus); 
      }
      else if(x <= R)
      {
        double a = 1 / ((R-c)*(R-c)*(R-c)); 
        double b = 1 / ((R-c)*(c-R)*(c-R)); 
        double prefocus = 1 + (x-c)/(c-R) + a*(x-c)*(x-c)*(x-R) + b*(x-c)*(x-R)*(x-R); 
        return h * std::pow(prefocus, focus); 
      }
      else
      {
        return 0.0; 
      }
    }
};

#endif // BumpFunc.hpp