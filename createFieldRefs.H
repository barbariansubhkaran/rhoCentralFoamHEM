volScalarField& p = mixture.p();
volScalarField& T = mixture.T();
const volScalarField& mu = alpha1*mixture.thermo1().mu()+alpha2*mixture.thermo2().mu();

bool inviscid(true);
if (max(mu.primitiveField()) > 0.0)
{
    inviscid = false;
}
