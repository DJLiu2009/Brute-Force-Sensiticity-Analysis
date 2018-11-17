/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/
/*
#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "fvcSmooth.H"*/
//#include "SLGThermo.H"//for ESF
//#include "basicSprayCloud.H"//for spray
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/////////////////////////////////////////////////////////////
//
// zero-dimensional kinetics example program
//
/////////////////////////////////////////////////////////////

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.
#include<math.h>
#include "cantera/zerodim.h"
#include "cantera/IdealGasMix.h"
#include "example_utils.h"
#include "scalarMatrices.H"

#define SENS_REAC 1
#define SENS_Cp 2
#define SENS_Hf 3
#define SENS_S0 4

using namespace Cantera;
using std::cout;
using std::endl;

int kinetics1(int np, void* p)
{
    // create an ideal gas mixture that corresponds to GRI-Mech 3.0
	int SENS_type;
	bool ThermalEffect;
	double tend;
	double FuelMole;
	size_t SENS_reaction_index;
	size_t SENS_species_index;

    IdealGasMix gas("mech/DME_chem.cti", "");
	IdealGasMix gas_0("mech/DME_chem.cti", "");
	gas_0.setState_TP(298.15, OneAtm);
	int nsp = gas.nSpecies();
	int nre = gas.nReactions();
	
	ThermalEffect=false;  // [True]: set the thermal effect [on]. [False]:set the thermal effect [off].
	SENS_type=2;         // 1: Reaction, 2: Heat_capacity, 3: Formation Enthalpy, 4: Entropy
	tend = 0.1;          // endtime, s
	std::string Fuel="CH3OCH3";
	doublereal Tinit=700.0;
	doublereal Pinit=30.0;
	doublereal phi=1.0;
	
	vector_fp x(nsp, 0.0);
	doublereal C_atoms = gas.nAtoms(gas.speciesIndex(Fuel),gas.elementIndex("C"));
	doublereal H_atoms = gas.nAtoms(gas.speciesIndex(Fuel),gas.elementIndex("H"));
	doublereal O_atoms = gas.nAtoms(gas.speciesIndex(Fuel),gas.elementIndex("O"));
	doublereal Oxy_stoic = C_atoms + H_atoms/4.0 - O_atoms/2.0;
	doublereal Air_stoic = 0.21 / Oxy_stoic;
	x[gas.speciesIndex(Fuel)] = 1.0;
	x[gas.speciesIndex("O2")] = 0.21 / phi / Air_stoic;
	x[gas.speciesIndex("N2")] = 0.79 / phi/ Air_stoic;
	
	vector_fp C_num(nsp,0);
	vector_fp H_num(nsp,0);
	vector_fp O_num(nsp,0);
	vector_fp N_num(nsp,0);
	vector_fp Ar_num(nsp,0);
	vector_fp He_num(nsp,0);
	vector_fp RS_CKP(nsp,0);
	for (int i=0; i<nsp; i++)
	{
		C_num[i]=gas.nAtoms(i,gas.elementIndex("C"));
		H_num[i]=gas.nAtoms(i,gas.elementIndex("H"));
		O_num[i]=gas.nAtoms(i,gas.elementIndex("O"));
		N_num[i]=gas.nAtoms(i,gas.elementIndex("N"));
		Ar_num[i]=gas.nAtoms(i,gas.elementIndex("Ar"));
		He_num[i]=gas.nAtoms(i,gas.elementIndex("He"));
	}

	IdealGasMix Air_0("mech/Methane_chem.cti", "gas");
	Air_0.setState_TPX(298.15, OneAtm, "O2:0.2035, N2:0.7567, CO2:0.0003, H2O:0.0303,AR:0.0091948, HE:0.0000052");
	vector_fp ParCKP_air(nsp,0);
	Air_0.getChemPotentials(&ParCKP_air[0]);
	double ENCKP_O2 = ParCKP_air[Air_0.speciesIndex("O2")];
	double ENCKP_N2 = ParCKP_air[Air_0.speciesIndex("N2")];
	double ENCKP_CO2 = ParCKP_air[Air_0.speciesIndex("CO2")];
	double ENCKP_H2O = ParCKP_air[Air_0.speciesIndex("H2O")];
	double ENCKP_AR = ParCKP_air[Air_0.speciesIndex("AR")];
	double ENCKP_HE = ParCKP_air[Air_0.speciesIndex("HE")];
	for (int i=0; i<nsp; i++)
	{
		RS_CKP[i]=C_num[i]*ENCKP_CO2+0.5*H_num[i]*ENCKP_H2O+(0.5*O_num[i]-C_num[i]-0.25*H_num[i])*ENCKP_O2+0.5*N_num[i]*ENCKP_N2+Ar_num[i]*ENCKP_AR+He_num[i]*ENCKP_HE;
	}

	Foam::scalarRectangularMatrix StoichCoeff(nsp,nre);
	for (int i=0; i<nsp; i++)
	{
		for (int j=0; j<nre; j++)
		{
			StoichCoeff[i][j]=gas.reactantStoichCoeff(i,j) - gas.productStoichCoeff(i,j);
		}
	}
	
	cout << "Constant-volume ignition of "
         << Fuel<< "/oxygen/nitrogen mixture."
         << endl;
    // create a 2D array to hold the output variables,
    Array2D soln(nsp+13, 1);
	std::ofstream f("Global_Sensitivity.csv");
	f << "SpeciesName" << "," << "Tau_dTdt_max" << "," << "Tau_OH_max" << "," << "dTdt" << "," << "f_Ex_Des"  
	  << "," << "f_Ex_Los" << "," << "f_En_Conv" << "," << "f_Ex_InC" << "," << "f_Ex_ThM" << std::endl;
	  

	
	for (int k=0; k<nsp; k++)
	{
		clock_t t0 = clock(); // save start time
	// set the initial state
		SENS_reaction_index=1;
		SENS_species_index=k;	
		gas.setState_TPX(Tinit, Pinit*OneAtm, x.data());
		FuelMole=gas.molarDensity()*gas.moleFraction(Fuel);

		double tm=0.0; // starttime, s
		vector_fp Time;
		vector_fp dSdt;
		vector_fp dTdt;
		vector_fp TherM;
		vector_fp ChemE;
		vector_fp Etot;
		vector_fp OH_mf;
		
		int nstep=0;	
		double T_dot;
		double f_OH;
		double dTdt_2=0.0;
		double Tau_dTdt_max=0.0;
		double Tau_OH_max=0.0;
		double HRR;
		double ROEG;
		double TME;
		double ChE;
		double LHV;
		double f_Fu_Conv;
		double f_En_Conv;
		double Ex_D=0.0;
		double f_Ex_Des=0.0;
		double f_Ex_Los;
		double f_Ex_InC;
		double f_Ex_ThM;
		
		vector_fp Xmf(nsp);
		vector_fp NetPR(nsp);
		vector_fp NetROP(nre);
		vector_fp ParIntE(nsp);
		vector_fp ParH(nsp);
		vector_fp ParCKP(nsp);
		vector_fp ParCKP_0(nsp);
		vector_fp dS_kdt(nre);

		// create a reactor
		IdealGasReactor r;
		IdealGasReactor r_0;
		r.insert(gas);     // setThermoMgr:,setKineticsMgr: m_chem is true
		r_0.insert(gas_0);

		// create a container object to run the simulation
		// and add the reactor to it
		ReactorNet sim;
		ReactorNet sim_0;
		sim.addReactor(r);
		sim_0.addReactor(r_0);
	
		//Set the Sensitivity object
		if (SENS_type==SENS_REAC) {
			r.addGlobalSensitivityReaction(SENS_reaction_index);
		}
		else if (SENS_type==SENS_Cp) {
			r.addGlobalSensitivitySpeciesHeat_Capacity(SENS_species_index);
			r_0.addGlobalSensitivitySpeciesHeat_Capacity(SENS_species_index);
		}
		else if (SENS_type==SENS_Hf) {
			r.addGlobalSensitivitySpeciesEnthalpy(SENS_species_index);
			r_0.addGlobalSensitivitySpeciesEnthalpy(SENS_species_index);
		}
		else if (SENS_type==SENS_S0) {
			r.addGlobalSensitivitySpeciesEntropy(SENS_species_index);
			r_0.addGlobalSensitivitySpeciesEntropy(SENS_species_index);
		}
		else {
		throw CanteraError("Brute-force Sensitivity Analysis",
							"Unknown Sensitivity object type: {}.", SENS_type);
		}
			
		r.applyGlobalSensitivity();
		r_0.applyGlobalSensitivity();
		
		r.setThermalEffect(ThermalEffect);  
		r.setGlobalSens(false); // True: set the Global Sensitivity [off].False:set the Global Sensitivity [on]. 
								// Do not modify for global sensitivity analysis!!!
		
	// main loop
		while (tm <= tend) 
		{
			Time.push_back(tm);
			HRR=0;
			ROEG=0;
			TME=0;
			ChE=0;
			LHV=0;
			fill(dS_kdt.begin(), dS_kdt.end(), 0.0);
	
			gas.getMoleFractions(&Xmf[0]);
			
			if (ThermalEffect) r.resetGlobalSensitivity();
			gas.getNetProductionRates(&NetPR[0]);
			gas.getNetRatesOfProgress(&NetROP[0]);
			if (ThermalEffect) r.applyGlobalSensitivity();
			
			gas.getPartialMolarIntEnergies(&ParIntE[0]);
			gas.getPartialMolarEnthalpies(&ParH[0]);
			gas.getChemPotentials(&ParCKP[0]);
			
			gas_0.setMoleFractions(&Xmf[0]);
			gas_0.getChemPotentials(&ParCKP_0[0]);
	
			for (int i=0; i<nsp; i++)
			{
				HRR -= NetPR[i]*ParIntE[i];
				TME -= Xmf[i]*ParCKP_0[i];
				ChE -= Xmf[i]*RS_CKP[i];
				LHV += Xmf[i]*(ParH[i]+(C_num[i]+H_num[i]/4-O_num[i]/2)*ParH[gas.speciesIndex("O2")]-(C_num[i]*ParH[gas.speciesIndex("CO2")]+H_num[i]/2*ParH[gas.speciesIndex("H2O")]+N_num[i]/2*ParH[gas.speciesIndex("N2")]));
			}
			T_dot=HRR/gas.density()/gas.cv_mass();
			dTdt.push_back(T_dot);
			OH_mf.push_back(gas.moleFraction("OH"));
			Etot.push_back(gas.molarDensity()*LHV);
			f_En_Conv=100*(1-Etot[nstep]/Etot[0]);
			f_Fu_Conv=100*(1-gas.molarDensity()*gas.moleFraction(Fuel)/FuelMole);
			
			ChE-=TME;
			ChE*=gas.molarDensity();
			TME+=(gas.intEnergy_mole()+OneAtm*gas.molarVolume()-298.15*gas.entropy_mole());
			TME*=gas.molarDensity();
	
			for (int j=0; j<nre; j++)
			{
				for (int i=0; i<nsp; i++)
				{
					dS_kdt[j] += StoichCoeff[i][j]*ParCKP[i]*NetROP[j]/gas.temperature();
				}
				
				ROEG += dS_kdt[j];
			}
			dSdt.push_back(ROEG);
			TherM.push_back(TME);
			ChemE.push_back(ChE);
			f_Ex_Los=100*(TherM[0]+ChemE[0]-TME-ChE)/ChemE[0];
			if (nstep>0)
			{
				dTdt_2 = (dTdt[nstep]-dTdt[nstep-1])/(tm-Time[nstep-1]);
				Ex_D += (298.15*(ROEG+dSdt[nstep-1])*(tm-Time[nstep-1])/2);
				f_Ex_Des=100*Ex_D/ChemE[0];
			}
	
			//cout << "Time = " << tm << " s,"
			//	<< "Temperature = " << gas.temperature() << " K." << endl;
	
			if (nstep%5==0)
			{
				saveSoln(tm, T_dot, dTdt_2, HRR, ROEG, TME, ChE, Ex_D, f_Ex_Des, f_En_Conv, f_Fu_Conv, gas, soln);
			}
	
			tm=sim.step();
			nstep++;
		}
		T_dot=dTdt[0];
		f_OH=OH_mf[0];
		for (int i=0; i<nstep; i++)
		{
			if (T_dot<dTdt[i])
			{
				T_dot=dTdt[i];
				Tau_dTdt_max=Time[i];
			}
			if (f_OH<OH_mf[i])
			{
				f_OH=OH_mf[i];
				Tau_OH_max=Time[i];
			}
		}
		
		f_Ex_InC = 100 * ChemE[nstep-1] / ChemE[0];
		f_Ex_ThM = 100 * (TherM[nstep-1]-TherM[0]) / ChemE[0];
		
		clock_t t1 = clock(); // save end time
		doublereal tmm = 1.0*(t1 - t0)/CLOCKS_PER_SEC;
		
		cout << "\nDisplay the final results: "<< endl;
		cout << " Computing time = " << tmm << " s." << endl;
		cout << "Speciesindex  = " << k+1 << " . " << endl;
		cout << "SpeciesName  = " << gas.speciesName(k) << " . " << endl;
		cout << "Tau_dTdt_max = " << 1000*Tau_dTdt_max << " ms. " << endl;
		cout << "Tau_OH_max   = " << 1000*Tau_OH_max << " ms. " << endl;
		cout << "Temperature  = " << gas.temperature() << " K. " << endl;
		cout << "Pressure     = " << gas.pressure()/OneAtm << " atm. " << endl;
		cout << "HRR          = " << HRR << " J/m^3/s. " << endl;
		cout << "dTdt         = " << dTdt[nstep-1] << " K/s. " << endl;
		cout << "dTdt_2       = " << dTdt_2 << " K/s^2. " << endl;
		cout << "dSdt         = " << dSdt[nstep-1] << " J/m^3/K/s." << endl;
		cout << "TherM        = " << TME << " J/m^3." << endl;
		cout << "ChemE        = " << ChE << " J/m^3." << endl;
		cout << "f_Ex_Des     = " << f_Ex_Des << "%." << endl;
		cout << "f_Ex_Los     = " << f_Ex_Los << "%." << endl;
		cout << "f_En_Conv    = " << f_En_Conv << "%." << endl;
		cout << "f_Fuel_Conv  = " << f_Fu_Conv << "%." << endl;
		cout << "f_Ex_InC     = " << f_Ex_InC << "%." << endl;
		cout << "f_Ex_ThM     = " << f_Ex_ThM << "%." << endl << endl;
		
		f << gas.speciesName(k) << "," << 1000*Tau_dTdt_max << "," << 1000*Tau_OH_max 
		  << "," << dTdt[nstep-1] << "," << f_Ex_Des << "," << f_Ex_Los << "," << f_En_Conv 
		  << "," << f_Ex_InC << "," << f_Ex_ThM << std::endl;
		
		// make a Tecplot data file and an Excel spreadsheet
		std::string plotTitle = "kinetics example 1: constant-volume ignition.";
		plotSoln("kin1.dat", "TEC", plotTitle, gas, soln);
		plotSoln("kin1.csv", "XL", plotTitle, gas, soln);
	
		// print the timing data
		//cout << " number of residual function evaluations = "
		//	<< sim.integrator().nEvals() << "." << endl;
		//cout << " time per evaluation = " << tmm/sim.integrator().nEvals()
		//	<< " s." << endl << endl;
		//cout << "Output files:" << endl
		//	<< "  kin1.csv  (Excel CSV file)" << endl
		//	<< "  kin1.dat  (Tecplot data file)" << endl << endl;
	
		r.resetGlobalSensitivity();
		r_0.resetGlobalSensitivity();
	}
    return 0;
}

int main(int argc, char *argv[])
{    
    try 
    {
        int retn = kinetics1(0, 0);
        appdelete();
        return retn;
    } 
    catch (CanteraError& err) 
    {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        cout << " terminating... " << endl;
        appdelete();
        return -1;
    }
}
// ************************************************************************* //
