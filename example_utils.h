#ifndef CT_EXAMPLE_UTILS_H
#define CT_EXAMPLE_UTILS_H

#include "cantera/base/Array.h"
#include "cantera/base/plots.h"

// Save the temperature, density, pressure, and mole fractions at one
// time
template<class G, class A>
void saveSoln(int i, double time, double dTdt, double dTdt_2, double HRR, double dSdT, double TME, double ChE, double ExD, double f_Ex_Des, double f_En_Conv, double f_Fu_Conv, const G& gas, A& soln)
{
    soln(0,i) = time;
	soln(1,i) = gas.temperature();
	soln(2,i) = gas.pressure()/Cantera::OneAtm;
	soln(3,i) = dTdt;
	soln(4,i) = dTdt_2;
	soln(5,i) = HRR;
	soln(6,i) = dSdT;
	soln(7,i) = TME;
	soln(8,i) = ChE;
	soln(9,i) = ExD;
	soln(10,i) = f_Ex_Des;
	soln(11,i) = f_En_Conv;
	soln(12,i) = f_Fu_Conv;
	//soln(13,i) = lambda;
	int nsp = gas.nSpecies();
	//for (int k = 0; k < nsp; k++) {
    //    soln(14+k,i) = u_SenCoe[k];
    //}
    gas.getMoleFractions(&soln(13,i));
}

template<class G, class A>
void saveSoln(double time, double dTdt, double dTdt_2, double HRR, double dSdT, double TME, double ChE, double ExD, double f_Ex_Des, double f_En_Conv, double f_Fu_Conv, const G& gas, A& soln)
{
    soln.resize(soln.nRows(), soln.nColumns() + 1);
    int back = soln.nColumns() - 1;
    soln(0,back) = time;
	soln(1,back) = gas.temperature();
	soln(2,back) = gas.pressure()/Cantera::OneAtm;
	soln(3,back) = dTdt;
	soln(4,back) = dTdt_2;
	soln(5,back) = HRR;
	soln(6,back) = dSdT;
	soln(7,back) = TME;
	soln(8,back) = ChE;
	soln(9,back) = ExD;
	soln(10,back) = f_Ex_Des;
	soln(11,back) = f_En_Conv;
	soln(12,back) = f_Fu_Conv;
	//soln(13,back) = lambda;
    int nsp = gas.nSpecies();
	//for (int k = 0; k < nsp; k++) {
    //    soln(14+k,back) = u_SenCoe[k];
    //}
    for (int k = 0; k < nsp; k++) {
        soln(13+k,back) = gas.moleFraction(k);
    }
}

template<class G, class V>
void makeDataLabels(const G& gas, V& names)
{
    int nsp = gas.nSpecies();
    names.resize(nsp + 13);
    names[0] = "time (s)";
    names[1] = "Temperature (K)";
    names[2] = "Pressure (atm)";
	names[3] = "dTdt (K/s)";
	names[4] = "dTdt_2 (K/s^2)";
	names[5] = "HRR (J/m^3/s)";
	names[6] = "dSdT (J/m^3/K/s)";
	names[7] = "TherMechE (J/m^3)";
	names[8] = "ChemE (J/m^3)";
	names[9] = "ExD (J/m^3)";
	names[10] = "f_Ex_Des (%)";
	names[11] = "f_En_Conv (%)";
	names[12] = "f_Fu_Conv (%)";
	//names[13] = "lambda (-)";

    for (int k = 0; k < nsp; k++) {
        names[13+k] = gas.speciesName(k);
    }
	//for (int k = 0; k < nsp; k++) {
    //    names[nsp+14+k] = gas.speciesName(k);
    //}
}

template<class G, class A>
void plotSoln(const std::string& fname, const std::string& fmt,
              const std::string& title, const G& gas, const A& soln)
{
    std::vector<std::string> names;
    makeDataLabels(gas, names);
    writePlotFile(fname, fmt, title, names, soln);
}

#endif
