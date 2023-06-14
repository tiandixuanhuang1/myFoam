/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: Open Source CFD
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
------------------------------------------------------------------------------- 
License
    This file isn't part of foam-extend nor OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    myRiemannSolverFoam

Description
    Density-based compressible implicit steady-state & transient flow solver
    using Riemann Solver.


Author
    Lei Wu

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "bound.H"
#include "boundMinMax.H"
#include "numericFlux.H"  // 数值格式
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#include "postProcess.H"
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createDynamicFvMesh.H"
	#include "createFields.H"
	#include "createTimeControls.H"
	
//	#include "createRDeltaTau.H"
  Info << "1712" << endl;	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	turbulence->validate();
	// Patch for correct calculation of meshPhi (U needs old value)
	{
		auto dummy = fvc::ddt(U);
	}
	
	Info << "\n Starting time loop\n" << endl;
	
	scalar CoNum = 0.0;
	scalar meanCoNum = 0.0;
	
	// Time advance
	while (runTime.run())
	{
		#include "readTimeControls.H"
		#include "readFieldBounds.H"
		
		// Calculate sound sonic
		surfaceScalarField amaxSf
		(
			"amaxSf",
			mag(fvc::interpolate(U) & mesh.Sf())
		  + mesh.magSf() * fvc::interpolate( sqrt(thermo.Cp()/thermo.Cv()/thermo.psi()) )
		);
		
		#include "compressibleCFLNO.H"
		#include "setDeltaT.H"
		
/* 		if (LTS)
		{
			#include "setRDeltaTau.H"
		} */
		
		Info << "\n Time = " << runTime.value() << nl;
		mesh.update();
		
		// Saves the value of the conserved quantity at the last time.
		rhoOld = rho;
		rhoUOld = rhoU;
		rhoEOld = rhoE;
		
		for (label cycle = 0; cycle < TVDCoeffU.size(); cycle++)
		{
			MRF.correctBoundaryVelocity(U);
			dbnsFlux.computeFlux();
			dimensionedScalar dt = runTime.deltaT();

			rho  = TVDCoeffU[cycle]*rhoOld + (1.0-TVDCoeffU[cycle])*rho 
				 - TVDCoeffR[cycle]*dt*fvc::div(dbnsFlux.rhoFlux());
			rhoU = TVDCoeffU[cycle]*rhoUOld + (1.0-TVDCoeffU[cycle])*rhoU 
				 - TVDCoeffR[cycle]*dt*fvc::div(dbnsFlux.rhoUFlux());
			rhoE = TVDCoeffU[cycle]*rhoEOld + (1.0-TVDCoeffU[cycle])*rhoE 
				 - TVDCoeffR[cycle]*dt*fvc::div(dbnsFlux.rhoEFlux());
		
      
      U.ref() = rhoU()/rho();
      U.correctBoundaryConditions();
      rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

      e = rhoE/rho - 0.5*magSqr(U);
      e.correctBoundaryConditions();
      thermo.correct();
      rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));

      p.ref() = rho()/thermo.psi();
      p.correctBoundaryConditions();
      rho.boundaryFieldRef() == thermo.psi().boundaryField()*p.boundaryField();

		}
		
		runTime++;
		
		runTime.write();
		
		Info << "	ExecutionTime = "
			 << runTime.elapsedCpuTime()
			 << " s\n" << endl;
		 
	}
	
	Info << "\n end \n" << endl;
	
	return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
