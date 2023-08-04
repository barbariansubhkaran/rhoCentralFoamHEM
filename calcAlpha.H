forAll ( alpha1, celli) //loop through cell centres
{
   if(T[celli]>TmaxSat)
   {
	rhovSat[celli] = rho[celli];
	rholSat[celli] = rhovSat[celli]*2;
   }
   else
   {
	if(T[celli]<TminSat)
	{
   	   rholSat[celli] = rho[celli]*2;
	   rhovSat[celli] = rho[celli];
	}
	else
	{
	   rhovSat[celli] = tempSeriesRhovSat(T[celli]);
	   rholSat[celli] = tempSeriesRholSat(T[celli]);
	}
  }
  alpha1[celli]=(rho[celli]-rholSat[celli])/(rhovSat[celli]-rholSat[celli]);
}

forAll ( mesh.boundaryMesh() , ipatch ) {

   if(alpha1.boundaryField()[ipatch].type()!=fixedValueFvPatchScalarField::typeName)
   {
	fvPatchScalarField& faces_alpha = alpha1.boundaryFieldRef()[ipatch];

	forAll(faces_alpha, faceI)  // Loop over each face of the patch
	{ 
	   if(T.boundaryField()[ipatch][faceI]>TmaxSat)
	   {
		rhovSat.boundaryFieldRef()[ipatch][faceI] = rho.boundaryFieldRef()[ipatch][faceI];
		rholSat.boundaryFieldRef()[ipatch][faceI] = rhovSat.boundaryFieldRef()[ipatch][faceI]*2;
	   }
	   else
	   {
		if(T.boundaryField()[ipatch][faceI]<TminSat) 
		{
  		   rholSat.boundaryFieldRef()[ipatch][faceI] = rho.boundaryFieldRef()[ipatch][faceI]*2;
		   rhovSat.boundaryFieldRef()[ipatch][faceI] = rho.boundaryFieldRef()[ipatch][faceI];
		}
		else 
		{
		   rhovSat.boundaryFieldRef()[ipatch][faceI] = tempSeriesRhovSat(T.boundaryField()[ipatch][faceI]);
		   rholSat.boundaryFieldRef()[ipatch][faceI] = tempSeriesRholSat(T.boundaryField()[ipatch][faceI]);
		}
	   }
	   faces_alpha[faceI] = (rho.boundaryField()[ipatch][faceI]-rholSat.boundaryField()[ipatch][faceI])/(rhovSat.boundaryField()[ipatch][faceI] - rholSat.boundaryField()[ipatch][faceI]);
	   faces_alpha[faceI] = min(1.0,max(faces_alpha[faceI],0.0));
	}
   }
}
alpha1 = min(1.0,max(alpha1,0.0));
alpha2 = 1-alpha1;
