#include "all.h"
#include "itikBanks.h"
#include "newtonPoincare.h"
#include "fixPointsInSpace.h"
#include <fstream>
#include <utility>


#include <boost/phoenix/core.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;



//[ stiff_system_definition
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;


int main(int argc, char **argv)
{
	//    typedef rosenbrock4< double > stepper_type;
	//    typedef rosenbrock4_controller< stepper_type > controlled_stepper_type;
	//    typedef rosenbrock4_dense_output< controlled_stepper_type > dense_output_type;
	//[ integrate_stiff_system
	itikBanks fun(PARAMETERS);
	vector_type x(3);
	x[0] = 0;
	x[1] = 0;
	x[2] = 0.06;

	size_t num_of_steps = integrate_const(make_dense_output< rosenbrock4< double > >(1.0e-6, 1.0e-6),
		make_pair(stiff_system(), stiff_system_jacobi()),
		x, 0.0, 50.0, 0.01,
		cout << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << "\n");
	//]
	clog << num_of_steps << endl;
	cout << x[0] << " " << x[1];
	cin.get();

	//    typedef runge_kutta_dopri5< vector_type > dopri5_type;
	//    typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
	//    typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;
	//[ integrate_stiff_system_alternative

	vector_type x2(2, 1.0);

	size_t num_of_steps2 = integrate_const(make_dense_output< runge_kutta_dopri5< vector_type > >(1.0e-6, 1.0e-6),
		stiff_system(), x2, 0.0, 50.0, 0.01,
		cout << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << "\n");
	//]
	clog << num_of_steps2 << endl;

	cin.get();
	return 0;
}