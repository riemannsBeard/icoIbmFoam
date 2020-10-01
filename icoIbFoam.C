/*---------------------------------------------------------------------------*\
 =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    icoIbFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "Time.H"
#include <fstream> 
#include <iomanip>
#include <stdio.h>

const fileName getOutputDir(const Time& fvm)
{
    // Create the output folder if it doesn't exist yet
    fileName outputDir;

    if (Pstream::parRun())
    {
        outputDir = fvm.time().path()/".."/"postProcessing";
    }
    else
    {
        outputDir = fvm.time().path()/"postProcessing";
    }

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    return outputDir;
}



vector calculaDistancia
(
	const fvMesh& fvm,
	const DynamicList<label>& vent,
	const label& i,
	const pointField& pf,
	const label& j,
	const vector& v
)
{
	vector tmp;

	tmp.x() = (fvm.C()[vent[i]].component(0) - pf[j].x())/v.x();
	tmp.y() = (fvm.C()[vent[i]].component(1) - pf[j].y())/v.y();
	tmp.z() = (fvm.C()[vent[i]].component(2) - pf[j].z())/v.z();
	
	return tmp;
}


DynamicList<label> setVentana
(
	const fvMesh& fvm,
	const pointIOField& pf
)
{
	const scalar radius = max(pf.component(0));
	const volVectorField& C = fvm.C();

	DynamicList<label> ventana;

	forAll(C, I)
	{
		if (mag(C[I]) <= 1.25*radius && mag(C[I]) >= 0.75*radius)
		{
			ventana.append(I);
		}
	}
	return ventana;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "\nh: " << h << endl;

	scalar deltaH;
	vector delta;
	vector rTmp;

	vectorField Fl(puntosIO.size(), pTraits<vector>::zero);
	vectorField Ul(puntosIO.size(), pTraits<vector>::zero);
	vectorField Ud = Fl;

	const scalar& pi = Foam::constant::mathematical::pi;

	//- Busco los indices de los puntos de la malla que cumplen
	// el criterio de la distancia

	const DynamicList<label>& ventana = setVentana(mesh, puntosIO); 



	//
	
	//- Comienza bucle temporal

	//

	Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNo.H"
  		
		fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );
	
		//- Primera prediccion del campo de velocidades sin tener
		// en cuenta el forcing del IBM
        solve(UEqn == -fvc::grad(p));
		
		//- Busco los indices de los puntos de la malla que cumplen
		// el criterio de la distancia

		forAll(puntosIO, j)
		{
			Ul[j] = vector(0.0, 0.0, 0.0);

			forAll(ventana, i)
			{
				const vector& rTmp = calculaDistancia
				(
					mesh, ventana, i, puntosIO, j, h
				);
				
				if
				(
					(abs(rTmp.x()) < 2.0) &&
					(abs(rTmp.y()) < 2.0) &&
					(abs(rTmp.z()) < 2.0)
				)
				{
					delta.x() =	(1.0/4.0)*(1.0 + Foam::cos(pi*0.5*rTmp.x()));
					delta.y() =	(1.0/4.0)*(1.0 + Foam::cos(pi*0.5*rTmp.y()));
					delta.z() =	(1.0/4.0)*(1.0 + Foam::cos(pi*0.5*rTmp.z()));

					deltaH = (1.0/(h.x()*h.y()*h.z()))*
						delta.x()*delta.y()*delta.z();
				}
					
				else
				{
					deltaH = 0.0;
				}
				Ul[j] += U[ventana[i]]*deltaH*(h.x()*h.y()*h.z());
			}
		}
		Fl = (Ud - Ul)/runTime.deltaTValue();


		forAll(ventana, i)
		{
			f[ventana[i]] = vector(0.0, 0.0, 0.0);

			forAll(puntosIO, j)
			{
				const vector& rTmp = calculaDistancia
				(
					mesh, ventana, i, puntosIO, j, h
				);
				
				if
				(
					(abs(rTmp.x()) < 2.0) &&
					(abs(rTmp.y()) < 2.0) &&
					(abs(rTmp.z()) < 2.0)
				)
				{
					delta.x() =	(1.0/4.0)*(1.0 + Foam::cos(pi*0.5*rTmp.x()));
					delta.y() =	(1.0/4.0)*(1.0 + Foam::cos(pi*0.5*rTmp.y()));
					delta.z() =	(1.0/4.0)*(1.0 + Foam::cos(pi*0.5*rTmp.z()));

					deltaH = (1.0/(h.x()*h.y()*h.z()))*
						delta.x()*delta.y()*delta.z();
				}

				else
				{
					deltaH = 0.0;
				}

				f[ventana[i]] += Fl[j]*deltaH*(h.x()*h.y()*h.z());
			}
		}

		//- Segunda prediccion del campo de velocidades teniendo
		// en cuenta el forcing del IBM
		UEqn -= f;
		solve(UEqn == -fvc::grad(p));

		// --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA("HbyA", U);

//			surfaceVectorField rAUf("rAUf", fvc::interpolate(f*rAU));

			HbyA = rAU*UEqn.H();

//			surfaceScalarField phif(rAUf & mesh.Sf());

            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                (fvc::interpolate(HbyA) & mesh.Sf())
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
			);

			adjustPhi(phiHbyA, U, p);	

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
					fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
				);

				pEqn.solve();
				
                if (nonOrth == nNonOrthCorr)
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"
			
			U = HbyA - rAU*fvc::grad(p);// + runTime.deltaT()*f;
            U.correctBoundaryConditions();
		}


		// Integro y adimensionalizo por unidad de profundidad 
		vector Cf(-2.0*gSum(f)*h.x()*h.y());
		Info << "Cf: " << Cf << endl;

//		getOutputDir(runTime);
		#include "escribeFuerzas.H"

        runTime.write();

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
