Info << "calculating T for mixture" << endl;
// Code to solve for temperature T of mixture from internal energy e

const volScalarField T0 = T;
const volScalarField f = e;

volScalarField Test = T0;
volScalarField Tnew = T0;
volScalarField eTol = e*1.0e-4;

volScalarField F = mixture.he(p,T0);
volScalarField dFdT = mixture.Cp();

T = Test;

mixture.correctThermo();
mixture.correct();

int iter = 0;
do
{
	Test = T;
	F = mixture.he(p,Test);
	dFdT = mixture.Cp();

	Tnew = T - (F-f)/dFdT;
	T = Tnew;

	#include "calcAlpha.H"
	
	mixture.correctThermo();
	mixture.correct();
	
	if (iter++ > 200)
	{
		FatalErrorInFunction
			<< "Maximum number of iterations exceeded: 200" << abort(FatalError);
	}
} while (max(mag(T - Test)).value() > 250*1e-4);
