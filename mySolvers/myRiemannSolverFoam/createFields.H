#include "createMRF.H"

Info << "Reading thermophysical properties \n" << endl;

autoPtr<psiThermo> pThermo
(
	psiThermo::New(mesh)
);

psiThermo& thermo = pThermo();
thermo.validate(args.executable(),"e");  // e 表示内能

// Primitive variables

volScalarField& e = thermo.he();
volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();


Info << "Reading field rho \n" << endl;
volScalarField rho
(
	IOobject
	(
		"rho",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	thermo.rho()
);


Info << "Reading field U\n" << endl;
volVectorField U
(
	IOobject
	(
			"U",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
	),
	mesh
);


// Conservative variables

volVectorField rhoU
(
	IOobject
	(
		"rhoU",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	rho*U
);

volScalarField rhoE
(
	IOobject
	(
		"rhoE",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	rho*(e + 0.5*magSqr(U))
);


// Create numeric flux
numericFlux dbnsFlux
(
	p,
	U,
	T,
	thermo,
	MRF
);


// Create mass flux alias for easier coupling with other code components
const surfaceScalarField& phi = dbnsFlux.rhoFlux();


Info << "Creating turbulence model \n" << endl;

autoPtr<compressible::turbulenceModel> turbulence
(
	compressible::turbulenceModel::New
	(
		rho,
		U,
		phi,
		thermo
	)
);

// define TVD coefficient
scalarList TVDCoeffU(0);

TVDCoeffU.append(1.0);
TVDCoeffU.append(3.0/4.0);
TVDCoeffU.append(1.0/3.0);

scalarList TVDCoeffR(0);

TVDCoeffR.append(1.0);
TVDCoeffR.append(1.0/4.0);
TVDCoeffR.append(2.0/3.0);


// define rho, rhoU, rhoE at last time.
volScalarField rhoOld
(
	IOobject
	(
		"rhoOld",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("rhoOld", dimensionSet(1,-3,0,0,0,0,0),scalar(0))
);

volVectorField rhoUOld
(
	IOobject
	(
		"rhoUOld",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedVector("rhoUOld", dimensionSet(1,-2,-1,0,0,0,0),vector(0.0,0.0,0.0))
);

volScalarField rhoEOld
(
	IOobject
	(
		"rhoEOld",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("rhoEOld", dimensionSet(1,-1,-2,0,0,0,0),scalar(0))
);
