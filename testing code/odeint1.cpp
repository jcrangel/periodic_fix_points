//Windows compilation: g++ -I  D:\\boost_1_67_0 exampleBoost.cpp

#include <iostream>
#include <boost/numeric/odeint.hpp>


using namespace std;
using namespace boost::numeric::odeint;


/* we solve the simple ODE x' = 3/(2t^2) + x/(2t)
 * with initial condition x(1) = 0.
 * Analytic solution is x(t) = sqrt(t) - 1/t
 */
typedef std::vector< double > state_type;

void rhs( const double x , double &dxdt , const double t )
{
    dxdt = 3.0/(2.0*t*t) + x/(2.0*t);
}

void write_cout( const double &x , const double t )
{
    cout << t << '\t' << x << endl;
}

// state_type = double 
// the steper type
typedef runge_kutta_dopri5< double > stepper_type;

int main()
{
    double x = 0.0;    

    //this version just print the solution
    // //                                  abserr , relerr ,                          
    integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) ,
                        rhs , x , 1.0 , 10.0 , 0.1 , write_cout );
    //     // integrate rhs(x) with from t=1 to t=10 with initial step size = 0.1 


}