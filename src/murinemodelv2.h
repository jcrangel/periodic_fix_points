/**********************************************************************************************/ /**
 * Copyright 2018 Julio C. Rangel
 * @file	murinemodelv2.h
 *
 * @brief	This is the murine model for the system on the Dendritic Immunotherapy Improvement for an Optimal
Control Murine Model. 2017 Hindawi
 **************************************************************************************************/

#ifndef MURINEMODEL_H
#define MURINEMODEL_H

#include "util.h"
#include <cmath>
class MurineModelv2
{
private:
    std::vector<Doub> parameters;
    const int systemSize = 8; // number of states
    int controlIndex = 3;	 // The index of the state where the control happens
    const int numParameters = 28;
    Doub aH, MuH, rH, KH, aC, MuC, rC, KC, KT, MuD, rI, MuIC, MuI;
    Doub rT, aT, eT, hT, aTBeta, eTBeta, rTBeta, MuBeta, aGammaC;
    Doub MuGamma, gMl, aMlGamma, eMlGamma, MuMl;

    void loadParameters()
    {
        aH = parameters[0];
        MuH = parameters[1];
        rH = parameters[2];
        KH = parameters[3];
        aC = parameters[4];
        MuC = parameters[5];
        rC = parameters[6];
        KC = parameters[7];
        KT = parameters[8];
        MuD = parameters[9];
        rI = parameters[10];
        MuIC = parameters[11];
        MuI= parameters[13];
        rT = parameters[14];
        aT = parameters[15];
        eT = parameters[16];
        hT = parameters[17];
        aTBeta = parameters[18];
        eTBeta = parameters[19];
        rTBeta = parameters[20];
        MuBeta = parameters[21];
        aGammaC = parameters[22];
        MuGamma = parameters[23];
        gMl = parameters[24];
        aMlGamma = parameters[25];
        eMlGamma = parameters[26];
        MuMl = parameters[27];
    }

    void changeParameters(std::vector<Doub> newParameters)
    {
        parameters = newParameters;
        loadParameters();
    }
public:
    MurineModelv2(std::vector<Doub> parameters) : parameters(parameters)
    {
        loadParameters();
    }

    void operator()(const StateType &x, StateType &dxdt, Doub /*t*/)
    {
        Doub T = x[0];
        Doub H = x[1];
        Doub CTL = x[2];
        Doub Den = x[3];
        Doub IL2 = x[4];
        Doub FBeta = x[5];
        Doub FGamma = x[6];
        Doub Ml = x[7];

        dxdt[0] = rT * T * std::log(KT / T) - ( aT * T * CTL * Ml / (eT + Ml) ) * ( (aTBeta * FBeta + eTBeta) / (FBeta + eTBeta));
        dxdt[1] = aH - MuH * H + rH * Den*(H*(1 - H / KH));
        dxdt[2] = aC - MuC * CTL + rC * IL2*(CTL*(1 - CTL / KC));
        dxdt[3] = -MuD * Den * CTL;
        dxdt[4] = rI * H * Den - MuIC * CTL * IL2 - MuI * IL2;
        dxdt[5] = rTBeta * T - MuBeta * FBeta;
        dxdt[6] = aGammaC * CTL - MuGamma * FGamma;
        dxdt[7] = gMl + (aMlGamma * FGamma) / (FGamma + eMlGamma) - MuMl * Ml;
    }

    int getSystemSize()
    {
        return systemSize;
    }
    int getControlIndex()
    {
        return controlIndex;
    }
};

#endif
