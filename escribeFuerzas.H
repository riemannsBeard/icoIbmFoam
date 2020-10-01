if (Pstream::master())
{
	fileName forcesDir;
	if (Pstream::parRun())
	{
		forcesDir = runTime.time().path()/"..";
	}
	else
	{
		forcesDir = runTime.time().path();
	}

	word forcesFileName = "forces.dat";
	fileName forces = forcesDir/forcesFileName;

	std::ofstream os
	(
		forces.c_str(),
		ios_base::app
	);

	os << runTime.time().value() << tab << tab <<
		"(" << Cf.x() << " " << Cf.y() << " " << Cf.z() << ")" << std::endl;
}

