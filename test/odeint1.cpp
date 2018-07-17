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

/* The rhs of x' = f(x) */
class harm_osc {

    double m_gam;

public:
    harm_osc( double gam ) : m_gam(gam) { }

    void operator() ( const state_type &x , state_type &dxdt , const double /* t */ )
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - m_gam*x[1];
    }
};
//]
// state_type = double 
// the steper type


int main()
{
 //    state_type x = {0,0};    

 //    //this version just print the solution
 //    // //                                  abserr , relerr ,     
 //    harm_osc harm(0.1);
	// typedef runge_kutta_dopri5< double > stepper_type;
 //    integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) ,
 //                        harm , x , 1.0 , 10.0 , 0.1 );
    //     // integrate rhs(x) with from t=1 to t=10 with initial step size = 0.1 
    //the integrator gives always the last element at the end
    // cout << endl << x[0] <<" , "<< x[1]<< endl;

    for(double i: vector<double> {10,20,30} ){
		cout<< i <<" "  ;  	
    }

}