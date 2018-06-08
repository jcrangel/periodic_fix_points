/*
 * rosenbrock4.cpp
 *
 * Copyright 2009-2012 Karsten Ahnert
 * Copyright 2009-2012 Mario Mulansky
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <iostream>
#include <fstream>
#include <utility>
#include "all.h"
#include <boost/numeric/odeint.hpp>

#include <boost/phoenix/core.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;



//[ stiff_system_definition
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

struct stiff_system
{
    void operator()( const vector_type &x , vector_type &dxdt , double /* t */ )
    {
        dxdt[ 0 ] = -101.0 * x[ 0 ] - 100.0 * x[ 1 ];
        dxdt[ 1 ] = 2*x[ 0 ];
    }
};

struct stiff_system_jacobi
{
    void operator()( const vector_type & /* x */ , matrix_type &J , const double & /* t */ , vector_type &dfdt )
    {
        J( 0 , 0 ) = -101.0;
        J( 0 , 1 ) = -100.0;
        J( 1 , 0 ) = 2.0;
        J( 1 , 1 ) = 0.0;
        dfdt[0] = 0.0;
        dfdt[1] = 0.0;
    }
};
//]

//template <class T>
struct push_back_state_and_time2
{
    std::vector< std::vector<double> >& m_states;
    std::vector<double> & m_times;

    push_back_state_and_time2( std::vector< std::vector<double> > &states, std::vector<double> &times )
        : m_states( states ), m_times( times ) { }

    void operator()( const std::vector<double> &x, double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
	void operator()(const vector_type &x, double t)
	{
		//std::vector<double> v(x.size());
		//std::copy(x.begin(), x.end(), v.begin());
		m_states.push_back(toStdVectorD(x));
		m_times.push_back(t);
	}
};


/*
//[ stiff_system_alternative_definition
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

struct stiff_system
{
    template< class State >
    void operator()( const State &x , State &dxdt , double t )
    {
        ...
    }
};

struct stiff_system_jacobi
{
    template< class State , class Matrix >
    void operator()( const State &x , Matrix &J , const double &t , State &dfdt )
    {
        ...
    }
};
//]
 */



int main( int argc , char **argv )
{
//    typedef rosenbrock4< double > stepper_type;
//    typedef rosenbrock4_controller< stepper_type > controlled_stepper_type;
//    typedef rosenbrock4_dense_output< controlled_stepper_type > dense_output_type;
    //[ integrate_stiff_system
    vector_type x( 2);
    x[0]=7;
    x[1]=-3;
    std::vector<std::vector<double>> u;
    std::vector<double> tt;
    size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-6 , 1.0e-6 ) ,
            make_pair( stiff_system() , stiff_system_jacobi() ) ,
            x , 0.0 , 50.0 , 0.01, push_back_state_and_time2(u, tt));
    //]
    //clog << num_of_steps << end;
    for( size_t i=0; i<= tt.size() - 1 ; i++ )
    {
        clog << tt[i] << '\t' << u[i][0] << '\t' << u[i][1] << '\n';
    }
    clog << x(0) << " " << x(1)<< endl;


	cin.get();
//    typedef runge_kutta_dopri5< vector_type > dopri5_type;
//    typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
//    typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;
    //[ integrate_stiff_system_alternative

   vector_type x2( 2);
   x2[0]=7;
   x2[1]=-3;
   size_t num_of_steps2 = integrate_const( make_dense_output< runge_kutta_dopri5< vector_type > >( 1.0e-6 , 1.0e-6 ) ,
           stiff_system() , x2 , 0.0 , 50.0 , 0.01,cout << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << "\n" );
   //]
   clog << num_of_steps2 << endl;
   clog << x2(0) << " " << x2(1)<< endl;
   cin.get();
    return 0;
}
