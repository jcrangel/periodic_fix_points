/*
All the functions related to integration of diferential equations
*/

#ifndef INTEGRATIONODE_H
#define INTEGRATIONODE_H


#include "util.h"



//[ integrate_observer
struct push_back_state_and_time
{
	std::vector< stateType>& m_states;
	std::vector<double> & m_times;

	push_back_state_and_time(std::vector< stateType > &states, std::vector<double> &times)
		: m_states(states), m_times(times) { }

	void operator()(const stateType &x, double t)
	{
		m_states.push_back(x);
		m_times.push_back(t);
	}

	void operator()(const vectorBoost &x, double t)
	{
		m_states.push_back(toStdVectorD(x));
		m_times.push_back(t);
	}
};

//Integrate a ODE system and returns the last value of the integration in x
//In the stiff version the system should be a class that have the Jacobian as a functor
template <class T>
void integrateSystem(T &system, stateType initialCondition, stateType &x, double t0, double tf) {

	typedef runge_kutta_dopri5<stateType> error_stepper_type;

	size_t steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6),
		system, initialCondition, t0, tf, 0.001);
	x = initialCondition;


}

//Integrate a ODE system and returns the last value of the integration in x
//In the stiff version the system should be a class that have the Jacobian as a functor
template <class T>
void integrateStiffSystem(T &system, stateType initialCondition, stateType &x, double t0, double tf) {

	vectorBoost v = toBoostVectorD(initialCondition);;

	//std::copy(initialCondition.begin(), initialCondition.end(), v.begin());

	size_t num_of_steps = integrate_const(make_dense_output< rosenbrock4< double > >(1.0e-6, 1.0e-6),
		std::make_pair(system, system),
		v, t0, tf, 0.01);

	//std::vector<double> res(STATE_SIZE);
	std::copy(v.begin(), v.end(), x.begin());

}

//Integrate a ODE system and returns the last value of the integration in x
//In the stiff version the system should be a class that have the Jacobian as a functor
//save the data into u and t
template <class T>
void integrateStiffSystem(T &system, stateType initialCondition,
	double t0, double tf, std::vector<stateType> &u, stateType &t) {

	vectorBoost v = toBoostVectorD(initialCondition);

	//std::copy(initialCondition.begin(), initialCondition.end(), v.begin());

	size_t num_of_steps = integrate_const(make_dense_output< rosenbrock4< double > >(1.0e-6, 1.0e-6),
		std::make_pair(system, system),
		v, t0, tf, 0.01, push_back_state_and_time(u, t));

}





#endif