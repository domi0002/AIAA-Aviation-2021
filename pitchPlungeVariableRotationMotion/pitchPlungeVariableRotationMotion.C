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

\*---------------------------------------------------------------------------*/

#include "pitchPlungeVariableRotationMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(pitchPlungeVariableRotationMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        pitchPlungeVariableRotationMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::pitchPlungeVariableRotationMotion::pitchPlungeVariableRotationMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::pitchPlungeVariableRotationMotion::
~pitchPlungeVariableRotationMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::pitchPlungeVariableRotationMotion::transformation() const
{
    scalar t = time_.value();

    scalar displacementy   = (1-exp(-ldecay_*t))*(linear0_ + linearAmplitude_.y()*sin(linearOmega_*t + linearPhase_));

    scalar displacementx   = -displacementy/tan(strokeAngle_);

    vector displacement(displacementx, displacementy, 0);
    
    vector eulerAngles    = theta0_*vector(0,0,1) +rotationalAmplitude_*sin(rotationalOmega_*t + thetaPhase_);
	   eulerAngles   *= (1-exp(-rdecay_*t));

    vector angleDeg = eulerAngles;

    eulerAngles *= degToRad();

  
    // Origin moves along the stroke plane
    point movingOrigin_(origin_ + displacement);


    Info << " LDecay = " << ldecay_ << nl
         << " RDecay = " << rdecay_ << nl
	 << " Displ  = " << displacement << endl
	 << " Rot    = " << angleDeg << endl;
  
    quaternion R(quaternion::XYZ, eulerAngles);
    septernion TR(septernion(-origin_ + -displacement)*R*septernion(origin_));


    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::pitchPlungeVariableRotationMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.readEntry("linearAmplitude", linearAmplitude_);
    SBMFCoeffs_.readEntry("rotationalAmplitude", rotationalAmplitude_);
    SBMFCoeffs_.readEntry("linearOmega", linearOmega_);
    SBMFCoeffs_.readEntry("rotationalOmega", rotationalOmega_);
    SBMFCoeffs_.readEntry("linearPhase", linearPhase_);
    SBMFCoeffs_.readEntry("linear0", linear0_);
    SBMFCoeffs_.readEntry("thetaPhase", thetaPhase_);
    SBMFCoeffs_.readEntry("theta0", theta0_);

    SBMFCoeffs_.readEntry("strokeAngle", strokeAngle_);
    SBMFCoeffs_.readEntry("origin", origin_);
    SBMFCoeffs_.readEntry("lineardecay", ldecay_);
    SBMFCoeffs_.readEntry("rotationaldecay", rdecay_);

    return true;
}


// ************************************************************************* //
