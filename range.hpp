/*
Anthony Caznona 
20887843
Winter 2022
University of Waterloo
ECE 204 Douglas Harder

-------------- Work Cited --------------
Richard Brent,
Algorithms for Minimization without Derivatives,
Dover, 2002,
ISBN: 0-486-41998-3,
LC: QA402.5.B74. 
*/

#include <stdlib.h>
#include <functional>
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;

double range( std::function< double(double) > f,
 double a,
 double b,
 double max_freq );

double golden(
  std::function< double(double) > f, double a, double b,
  double eps_step, double eps_abs,
  unsigned int max_iterations
);

double gldn_min(
  std::function< double(double) > f, double a, double b,
  double eps_step, double eps_abs,
  unsigned int max_iterations
);

double brent( std::function< double(double) > f,
 double a,
 double b,
 double err_tolr);

double brent_max( std::function< double(double) > f,
 double a,
 double b,
 double err_tolr);


double range( std::function< double(double) > f,
 double a,
 double b,
 double max_freq ){

  double err = 1e-17;

  double step { a + 1.0/max_freq };
  double glb_max {brent_max(f, a, b, err)};
  double glb_min {brent(f, a, b, err)};

 while(step <= b){
    double loc_max {brent_max(f, a, step, err)};
    double loc_min {brent(f, a, step, err)};
    if (loc_max > glb_max){glb_max = loc_max;}
    if (loc_min < glb_min){glb_min = loc_min;}
    a = step;
    step += 1/max_freq;
  }

  cout << "max: " << glb_max << endl << "min: "<< glb_min << endl;
  return glb_max - glb_min;
}



double gldn_min(
  std::function< double(double) > f, double a, double b,
  double eps_step, double eps_abs,
  unsigned int max_iterations
) {
assert( a < b );

  double const INV_PHI{ (std::sqrt(5.0) - 1.0)/2.0 };
  double width{ (b - a)*INV_PHI };
  double x1{ b - width };
  double x2{ a + width };

  double f1{ f( x1 ) };
  double f2{ f( x2 ) };

  assert( std::isfinite( f1 ) && std::isfinite( f2 ) );

  for ( unsigned int k{0}; k < max_iterations; ++k ) {
    width *= INV_PHI;

    if ( f1 < f2 ) {
      if ( (width < eps_step) && (std::abs( f1 - f2 ) < eps_abs) ) {
        return f1;
      }

       b = x2;
      x2 = x1;
      f2 = f1;
      x1 = b - width;
      f1 = f( x1 );

      assert( std::isfinite( f1 ) );
    } else if ( f2 < f1 ) {
      if ( (width < eps_step) && (std::abs( f1 - f2 ) < eps_abs) ) {
        return  f2;
      }

       a = x1;
      x1 = x2;
      f1 = f2;
      x2 = a + width;
      f2 = f( x2 );

      assert( std::isfinite( f2 ) );
    } else {
      assert( f1 == f2 );

      if ( width < eps_step ) {
        return f1;
      }

      a = x1;
      b = x2;

      width *= INV_PHI*INV_PHI;

      x1 = b - width;
      x2 = a + width;

      f1 = f( x1 );
      f2 = f( x2 );

      assert( std::isfinite( f1 ) && std::isfinite( f2 ) );
    }
  }

  return NAN;
}



double golden(
  std::function< double(double) > f, double a, double b,
  double eps_step, double eps_abs,
  unsigned int max_iterations
) {
  assert( a < b );

  double const INV_PHI{ (std::sqrt(5.0) - 1.0)/2.0 };
  double width{ (b - a)*INV_PHI };
  double x1{ b - width };
  double x2{ a + width };

  double f1 = f( x1 ) ;
  double f2 = f( x2 ) ;

  assert( std::isfinite( f1 ) && std::isfinite( f2 ) );

  for ( unsigned int k{0}; k < max_iterations; ++k ) {
    width *= INV_PHI;

    if ( f1 > f2 ) {
      if ( (width < eps_step) && (std::abs( f1 - f2 ) < eps_abs) ) {
        return f1;
      }

      b = x2;
      x2 = x1;
      f2 = f1;
      x1 = b - width;
      f1 = f( x1 ) ;

      assert( std::isfinite( f1 ) );
    } else if (f2 > f1) {
      if ( (width < eps_step) && (std::abs( f1 - f2 ) < eps_abs) ) {
        return f2;
      }

       a = x1;
      x1 = x2;
      f1 = f2;
      x2 = a + width;
      f2 = f( x2 ) ;

      assert( std::isfinite( f2 ) );
    } else {
      assert( f1 == f2 );

      if ( width < eps_step ) {
        return f1;
      }

      a = x1;
      b = x2;

      width *= INV_PHI*INV_PHI;

      x1 = b - width;
      x2 = a + width;

      f1 = f( x1 ) ;
      f2 = f( x2 ) ;

      assert( std::isfinite( f1 ) && std::isfinite( f2 ) );
    }
  }

  return NAN;
}



double brent( std::function< double(double) > f,
 double a,
 double b,
 double err_tolr){

  const double gldn_sqr {0.5 * ( 3.0 - sqrt ( 5.0 ) )};
  const double r8_eps {2.220446049250313E-016};

  double eps {sqrt(r8_eps)};

  double strt_a {a}; 
  double strt_b {b};

  double x {strt_a + gldn_sqr * (b - a)};
  double v , w {x}; 

  double fx = f(x);
  double fw, fv {fx};
  double fu {0.0};  
    
  double m {0.0};
  double tol {0.0};
  double tol2 {0.0};
  double e {0.0};
  double d {0.0};
  double u {0.0};
  double iter {0.0};

  while (true){
    m = (strt_a + strt_b) / 2 ;
    tol = eps * abs(x) + err_tolr;
    tol2 = 2.0 * tol;

    if ( abs ( x - m ) <= tol2 - 0.5 * ( strt_b - strt_a ) ){break;}

    double r {0.0};
    double q, p {r};

    if (tol < abs(e)){

    r = ( x - w ) * ( fx - fv );
    q = ( x - v ) * ( fx - fw );
    p = ( x - v ) * q - ( x - w ) * r;
    q = 2.0 * ( q - r );

    if ( 0.0 < q ){ p = - p;}

    q = abs(q);
    r = e;
    e = d;
    }

    if ( abs ( p ) < abs ( 0.5 * q * r ) &&
     q * ( strt_a - x ) < p &&
     p < q * ( strt_b - x ) )
    {
      d = p / q;
      u = x + d;

      if ( ( u - strt_a ) < tol2 || ( strt_b - u ) < tol2 )
      {
        if ( x < m ) { d = tol;}
        else {d = - tol;}
      }
    } 

    else
    {
      if ( x < m ) {e = strt_b - x;}
      else {e = strt_a - x;}
      d = gldn_sqr * e;
    }

    if ( tol <= abs ( d ) ) {u = x + d;}
    else if ( 0.0 < d ) {u = x + tol;}
    else {u = x - tol;}

    fu = f(u);

    if ( fu <= fx )
    {
      if ( u < x ) {strt_b = x;}
      else {strt_a = x;}
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    }
    else
    {
      if ( u < x ) {strt_a = u;}
      else {strt_b = u;}

      if ( fu <= fw || w == x )
      {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if ( fu <= fv || v == x || v == w )
      {
        v = u;
        fv = fu;
      }
    }
  }

  return fx;
}


double brent_max( std::function< double(double) > f,
 double a,
 double b,
 double err_tolr){

  const double gldn_sqr {0.5 * ( 3.0 - sqrt ( 5.0 ) )};
  const double r8_eps {2.220446049250313E-016};

  double eps {sqrt(r8_eps)};

  double strt_a {a}; 
  double strt_b {b};

  double x {strt_a + gldn_sqr * (b - a)};
  double v , w {x}; 

  double fx = f(x);
  double fw, fv {fx};
  double fu {0.0};  
    
  double m {0.0};
  double tol {0.0};
  double tol2 {0.0};
  double e {0.0};
  double d {0.0};
  double u {0.0};
  double iter {0.0};

  while (true){
    m = (strt_a + strt_b) / 2 ;
    tol = eps * abs(x) + err_tolr;
    tol2 = 2.0 * tol;

    if ( abs ( x - m ) <= tol2 - 0.5 * ( strt_b - strt_a ) ){break;}

    double r {0.0};
    double q, p {r};

    if (tol < abs(e)){

    r = ( x - w ) * ( fx - fv );
    q = ( x - v ) * ( fx - fw );
    p = ( x - v ) * q - ( x - w ) * r;
    q = 2.0 * ( q - r );

    if ( 0.0 < q ){ p = - p;}

    q = abs(q);
    r = e;
    e = d;
    }

    if ( abs ( p ) > abs ( 0.5 * q * r ) &&
     q * ( strt_a - x ) < p &&
     p < q * ( strt_b - x ) )
    {
      d = p / q;
      u = x + d;

      if ( ( u - strt_a ) < tol2 || ( strt_b - u ) < tol2 )
      {
        if ( x < m ) { d = tol;}
        else {d = - tol;}
      }
    } 

    else
    {
      if ( x < m ) {e = strt_b - x;}
      else {e = strt_a - x;}
      d = gldn_sqr * e;
    }

    if ( tol <= abs ( d ) ) {u = x + d;}
    else if ( 0.0 < d ) {u = x + tol;}
    else {u = x - tol;}

    fu = f(u);

    if ( fu >= fx )
    {
      if ( u < x ) {strt_b = x;}
      else {strt_a = x;}
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    }
    else
    {
      if ( u < x ) {strt_a = u;}
      else {strt_b = u;}

      if ( fu >= fw || w == x )
      {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if ( fu >= fv || v == x || v == w )
      {
        v = u;
        fv = fu;
      }
    }
  }

  return fx;
}