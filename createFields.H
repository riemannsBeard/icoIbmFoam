    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    IOdictionary spacing
    (
        IOobject
        (
            "spacing",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

	vector h = vector(spacing.lookup("h"));


	pointIOField puntosIO
    (
        IOobject
        (
            "puntosIO",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );



//- Diccionario y entradas para la cinematica del IBM
/*
	IOdictionary IBMvelocities
    (
        IOobject
        (
            "IBMvelocities",
            runTime.constant(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    dimensionedVector translationalU
    (
        IBMvelocities.lookup("translationalU")
    );

    dimensionedVector rotationalU
    (
        IBMvelocities.lookup("rotationalU")
    );

    dimensionedVector rotationalCenter
    (
        transportProperties.lookup("rotationCenter")
    );
*/
//----------------------------------------------------

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


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


    volVectorField f
    (
        IOobject
        (
            "f",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedVector
		(
			"0",
			dimLength/(dimTime*dimTime),
			pTraits<vector>::zero
		)
 	);


#   include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
