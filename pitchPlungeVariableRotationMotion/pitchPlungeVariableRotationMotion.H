/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::solidBodyMotionFunctions::pitchPlungeVariableRotationMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Oscillating displacement.

SourceFiles
    pitchPlungeVariableRotationMotion.C

\*---------------------------------------------------------------------------*/

#ifndef pitchPlungeVariableRotationMotion_H
#define pitchPlungeVariableRotationMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class pitchPlungeVariableRotationMotion Declaration
\*---------------------------------------------------------------------------*/

class pitchPlungeVariableRotationMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Amplitude
        vector linearAmplitude_, rotationalAmplitude_;

        //- Radial velocity
        scalar linearOmega_, rotationalOmega_;

	scalar thetaPhase_, theta0_, linearPhase_, linear0_;

	// - Start parameter
	scalar ldecay_, rdecay_;

 	//- Point of rotation(varies with time possibly)
 	vector origin_;

	//- Stroke angle
	scalar strokeAngle_;

    // Private Member Functions

        //- No copy construct
        pitchPlungeVariableRotationMotion(const pitchPlungeVariableRotationMotion&) = delete;

        //- No copy assignment
        void operator=(const pitchPlungeVariableRotationMotion&) = delete;


public:

    //- Runtime type information
    TypeName("pitchPlungeVariableRotationMotion");


    // Constructors

        //- Construct from components
        pitchPlungeVariableRotationMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new pitchPlungeVariableRotationMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~pitchPlungeVariableRotationMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
