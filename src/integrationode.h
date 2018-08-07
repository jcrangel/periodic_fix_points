/*
All the functions related to integration of diferential equations
*/

#ifndef INTEGRATIONODE_H
#define INTEGRATIONODE_H


#include "util.h"

using namespace boost::numeric::odeint;


/**********************************************************************************************//**
 * @struct	push_back_state_and_time
 *
 * @brief	Integrator  observer, for saving the states each iteration.
 *
 * @author	Iron
 * @date	7/31/2018
 **************************************************************************************************/

struct push_back_state_and_time
{
	std::vector< StateType>& m_states;
	std::vector<Doub> & m_times;

	push_back_state_and_time(std::vector< StateType > &states, std::vector<Doub> &times)
		: m_states(states), m_times(times) { }

	void operator()(const StateType &x, Doub t)
	{
		m_states.push_back(x);
		m_times.push_back(t);
	}

	void operator()(const VectorBoost &x, Doub t)
	{
		m_states.push_back(toStateType(x));
		m_times.push_back(t);
	}
};


/**********************************************************************************************//**
 * @fn	template <class T> void integrateSystem(T &system, StateType initialCondition, StateType &x, Doub t0, Doub tf)
 *
 * @brief	Integrate a ODE system and returns the last value of the integration in x In the stiff version.
 * 			the system should be a class that have the Jacobian as a functor. 
 *
 * @author	Iron
 * @date	7/31/2018
 *
 * @tparam	T	Generic type parameter.
 * @param [in,out]	system				The system.
 * @param 		  	initialCondition	The initial condition.
 * @param [in,out]	x					A StateType to process.
 * @param 		  	t0					The initial time.
 * @param 		  	tf					The final time.
 **************************************************************************************************/

template <class T>
void integrateSystem(T &system, StateType initialCondition, StateType &x, Doub t0, Doub tf) {

	typedef runge_kutta_dopri5<StateType> error_stepper_type;

	size_t steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6),
		system, initialCondition, t0, tf, 0.001);
	x = initialCondition;


}


/**********************************************************************************************//**
 * @fn	template <class T> void integrateStiffSystem(T &system, StateType initialCondition, StateType &x, Doub t0, Doub tf)
 *
 * @brief	Integrate a ODE system and returns the last value of the integration in x.System is stiff version
 * 			 the system should be a class that have the Jacobian as a functor. 
 *
 * @author	Iron
 * @date	7/31/2018
 *
 * @tparam	T	Generic type parameter.
 * @param [in,out]	system				The system.
 * @param 		  	initialCondition	The initial condition.
 * @param [in,out]	x					A StateType to process.
 * @param 		  	t0					The t0.
 * @param 		  	tf					The tf.
 **************************************************************************************************/

template <class T>
void integrateStiffSystem(T &system, StateType initialCondition, StateType &x, Doub t0, Doub tf) {

	VectorBoost v = toBoostVectorD(initialCondition);;

	//std::copy(initialCondition.begin(), initialCondition.end(), v.begin());

	size_t num_of_steps = integrate_const(make_dense_output< rosenbrock4< Doub > >(1.0e-6, 1.0e-6),
		std::make_pair(system, system ),
		v, t0, tf, 0.01);

	//std::vector<Doub> res(STATE_SIZE);
	std::copy(v.begin(), v.end(), x.begin());

}

//Integrate a ODE system and returns the last value of the integration in x
//In the stiff version the system should be a class that have the Jacobian as a functor
//save the data into u and t
template <class T>
void integrateStiffSystem(T &system, StateType initialCondition,
	Doub t0, Doub tf, std::vector<StateType> &u, StateType &t) {

	VectorBoost v = toBoostVectorD(initialCondition);

	//std::copy(initialCondition.begin(), initialCondition.end(), v.begin());

	size_t num_of_steps = integrate_const(make_dense_output< rosenbrock4< Doub > >(1.0e-6, 1.0e-6),
		std::make_pair(system, system),
		v, t0, tf, 0.01, push_back_state_and_time(u, t));

}





#endif