/**********************************************************************************************//**
 * Copyright 2018 Julio C. Rangel
 * @file	murinemodelv2.h
 *
 * @brief	This is the murine model for the system on the Dendritic Immunotherapy Improvement for an Optimal
Control Murine Model. 2017 Hindawi 
 **************************************************************************************************/

#ifndef MURINEMODEL_H
#define MURINEMODEL_H

#include "util.h"


class MurineModelv2 {
private:
	std::vector<Doub> parameters;
	const int systemSize = 8; // number of states
	int controlIndex = 3; // The index of the state where the control happens


};

#endif