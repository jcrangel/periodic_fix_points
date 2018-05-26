#ifndef NEWTONPOINCARE_H
#define NEWTONPOINCARE_H

#include "all.h"
//Gets the jacobian following the procedure in eq (6)&(7) of the
//paper Wei 2014
template <class T>
Matrix3d DFode(T &fun,stateType initialCondition,double tau,double d)
{
// xp=[0, tau];
	stateType xp{0,tau};

	initialCondition[2] = initialCondition[2] + d;
	std::vector<stateType> u;
    std::vector<double> tt;
    typedef runge_kutta_cash_karp54<stateType> error_stepper_type;
    size_t steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6),
                   fun, initialCondition, xp[0], xp[1], 0.001,
                   push_back_state_and_time(u, tt));
    // u save the data column wise, the firts colum has the first state data. It would
    //be nice to have in row 
    //This create a transpose of u, with dimentions: STATE_SIZE * steps (3 x steps)
    // with elements equal to zero
    std::vector < std::vector < double>> state(STATE_SIZE, std::vector<double>(u.size()));
    transpose(u,state);
     //Check if its correct
    //  for( size_t i=0; i <= steps; i++ )
    //  {
    //     std::cout << tt[i] << '\t' << state[0][i] << '\t'
    //     	     << state[1][i] << '\t' << state[2][i] << '\n';
    //  }
    std::vector<stateType> V;
    std::vector<double> t;
	//!!!! this param shoul exist only on one place not to be repated
    std::vector<double> param = { 1,2.5,0.5,1.5,4.5,1,0.2,0.5 };
    std::vector<double> identityMatrixVector = {1,0,0,0,1,0,0,0,1};
    itikBanksJacobian Dfu(param,state,tt);
    steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6), 
            Dfu, identityMatrixVector, xp[0], xp[1], 0.001,
            push_back_state_and_time(V, t));
 //   for( size_t i=0; i <= steps; i++ )
	//{
 //       std::cout << V[i][0] << '\t' << V[i][1] << '\t'
 //                 << V[i][2] << '\t' << V[i][3] << '\t' << V[i][4] << '\t'
 //                 << V[i][5] << '\t' << V[i][6] << '\t' << V[i][7] << '\t'
 //                 << V[i][8] << '\n';
 //   }

	// Returns the last matrix since 
	// DF(x) = v(tau).See page 12 research notebook
	return reshapeVectorToMatrix(V.back());
};



template <class T>
fixPoint newtonPoincare(T &fun, stateType initialCondition, double tau, double d) {
    int n = 0;            //%initialize iteration counter
    double eps = 1;          //%initialize error
    std::vector<double> xp{ 0, tau }; //% define the span of the computational domain
    double tol = 1e-6;
    int maxStep = 100;

    Vector3d Xk(initialCondition.data());
    Vector3d Pxk;
    Matrix3d df;
    std::vector<stateType> u;
    std::vector<double> t;
    Matrix3d A;
    Matrix3d I = Matrix3d::Identity();
    Vector3d y;
    Vector3d Xplus;
    size_t steps;
    typedef runge_kutta_cash_karp54<stateType> error_stepper_type;

    //stateType initCond;

    while (n < maxStep) {                   //!@_@
        df = DFode(fun, std::vector<double> {Xk[0], Xk[1], Xk[2]}, tau, d);  //Jacobian
        initialCondition[0] = Xk[0];
        initialCondition[1] = Xk[1];
        initialCondition[2] = Xk[2]+d;

        steps = integrate_adaptive(make_controlled<error_stepper_type>(1.0e-10, 1.0e-6),
            fun,                //!@_@                          
            initialCondition,
            xp[0], xp[1], 0.001, push_back_state_and_time(u, t));
        //[t, u] = ode15s(functionName, xp, [Xk(1), Xk(2), Xk(3) + d], [], [parameters d]);
        
        //!@This should be done in a bette wway
        Pxk[0] = u.back()[0];
        Pxk[1] = u.back()[1];
        Pxk[2] = u.back()[2];

        A = df - I;
        //Solve Ay = (Pxk-Xk)
        y = A.colPivHouseholderQr().solve(Pxk-Xk);
        Xplus = Xk - y;
        
        eps = y.norm();
        if (eps < tol)
            break;
    
        Xk = Xplus;  //%update x
        n++;
        
    }
    // 0 means unstable
    double stability = 0;
    bool convergence = false;
    if (n < maxStep){
        convergence = true;
        
        if (isStable(df))
            stability = 1;

    }
    fixPoint ss(convergence, stability, Xk);
    return ss;
}


#endif