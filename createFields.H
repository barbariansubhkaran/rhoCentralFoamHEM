#include "createRDeltaT.H"


Info<< "Reading field U\n" << endl;
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

Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
MytwoPhaseMixtureThermo mixture(U);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

volScalarField& e = mixture.thermo1().he();

e = alpha1*mixture.thermo1().he()+alpha2*mixture.thermo2().he();


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
   alpha1*mixture.thermo1().rho()+alpha2*mixture.thermo2().rho()
);

volScalarField rhovSat
(
    IOobject
    (
        "rhovSat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhovSat", dimDensity, 1.0)
);
volScalarField rholSat
(
    IOobject
    (
        "rholSat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rholSat", dimDensity, 1.0)
);


interpolationTable<scalar> tempSeriesRhovSat
(
 runTime.path()/runTime.caseConstant()/"rhovSat.dat"
);
interpolationTable<scalar> tempSeriesRholSat
(
 runTime.path()/runTime.caseConstant()/"rholSat.dat"
);

Info<< "Reading thermophysical properties\n" << endl;

IOdictionary thermophysicalProperties
(
 IOobject
 (
  "thermophysicalProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
 )
);

scalar TmaxSat(readScalar(thermophysicalProperties.lookup("TmaxSat")));
scalar TminSat(readScalar(thermophysicalProperties.lookup("TminSat")));


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

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

surfaceScalarField phi("phi", fvc::flux(rhoU));

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        mixture
    )
);

Info<< "Creating compressibilityModel\n" << endl;
autoPtr<barotropicCompressibilityModel> psiModel =
   barotropicCompressibilityModel::New
   (
	thermophysicalProperties,
	alpha1,
	rhovSat,
	rholSat
   );

const volScalarField& psi = psiModel->psi();

IOdictionary thermophysicalPropertiesLiquid
(
   IOobject
   (
	  "thermophysicalProperties.liquid",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ_IF_MODIFIED,
	IOobject::NO_WRITE
   )
);

dimensionedScalar K(dimPressure,readScalar(thermophysicalPropertiesLiquid.lookup("K")));
