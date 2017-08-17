#ifndef VORTICITY_DISTRIBUTION_H
#define VORTICITY_DISTRIBUTION_H

#include "Core\Math\vector3D.h"
#include "vorton.h"
#include <vector>
#include <math.h>

/*
	interface for defining a vorticity distribution
*/

class IVorticityDistribution
{
public:
	virtual Vector3D	GetDomainSize() const = 0;
	virtual void AssignVorticity(Vector3D & vorticity, const Vector3D & position, const Vector3D & vCenter) const = 0;
};

/*
	VortexRing, vorticity in this ring is not guaranteed to be solenoidal
	vCenter in AssignVorticity is the center of the ring,
	mThickness is the radius of the ring
	mRadius is the radius of the circle(vCenter and ring's center)
*/

class VortexRing : public IVorticityDistribution
{
public:
	VortexRing( const float & fRaduis, const float & fThickness, const Vector3D & vDirection )
		: mRaduis(fRaduis)
		, mThickness(fThickness)
		, mDirection(vDirection)
	{
		//vDirection must be normalized
	}

	virtual Vector3D	GetDomainSize() const
	{
		const float boxSideLength = 2.0f * ( mRaduis + mThickness );
		return Vector3D(1.0f, 1.0f, 1.0f) * boxSideLength;
	}

	virtual void		AssignVorticity(Vector3D & vorticity, const Vector3D & position, const Vector3D & vCenter) const
	{
		const Vector3D	vFromCenter		= position - vCenter;
		const float		fdistAlongDir	= dot( mDirection,  vFromCenter );
		const Vector3D	vPtOnLine		= vCenter + mDirection * fdistAlongDir; //Project point on the center-dir line
		const Vector3D  vRho			= position - vPtOnLine;
		const float		fRho			= vRho.norm();
		const float		fRadCore = sqrtf( POW2(fRho - mRaduis)
										+ POW2(fdistAlongDir ) );
		//cosf : The nearer the position to ring core, The bigger the vorticity is.
		if (fRadCore < mThickness)
		{
			Vector3D	vRhoDir		= vRho;
			vRhoDir.normalize();
			Vector3D	vVortDir	= cross( mDirection, vRhoDir );
			const float fVortSize	= 0.5f * ( cosf( PI * fRadCore / mThickness ) + 1.0f);
			vorticity				= fVortSize * vVortDir;
		}
		else
		{
			vorticity = Vector3D( 0.0f, 0.0f, 0.0f );
		}
	}

	float		mRaduis;
	float		mThickness;
	Vector3D	mDirection;
};

/*
	Vorticity in the shape of a vortex ring.
	The vorticity specified by this class derives from taking the curl of
	a localized jet.  The vorticity is therefore guaranteed to be solenoidal,
	to within the accuracy the discretization affords.
*/

class JetRing : public IVorticityDistribution
{
public:

	JetRing(const float & fRadiusSlug, const float & fThickness, const Vector3D & vDirection)
		: mRadiusSlug( fRadiusSlug )
		, mThickness( fThickness )
		, mRadiusOuter( fRadiusSlug + fThickness)
		, mDirection( vDirection )
	{

	}

	virtual Vector3D	GetDomainSize() const
	{
		const float		boxSideLength	= 2.0f * mRadiusOuter;
		return Vector3D( 1.0f, 1.0f, 1.0f ) * boxSideLength;
	}

	virtual void		AssignVorticity(Vector3D & vorticity, const Vector3D & position, const Vector3D & vCenter) const
	{
		const Vector3D	vFromCenter		= position - vCenter;
		const float		fdistAlongDir	= dot(vFromCenter, mDirection);
		const Vector3D	vPtOnLine		= vCenter + fdistAlongDir * mDirection;
		const Vector3D	vRho			= position - vPtOnLine;
		const float		fRho			= vRho.norm();

		if ( (fRho < mRadiusOuter) && (fRho > mRadiusSlug) )
		{	
			const float		streamwiseProfile	= (fabs(fdistAlongDir) < mRadiusSlug) ? 0.5f * (cosf(PI * fdistAlongDir / mRadiusSlug) + 1.0f) : 0.0f;
			const float		radialProfile		= sinf( PI * (fRho - mRadiusSlug) / mThickness );
			const float		fVortSize			= streamwiseProfile * radialProfile * PI / mThickness;
			
			Vector3D		vRhoDir				= vRho;
			vRhoDir.normalize();
			Vector3D vVortDir = cross( mDirection, vRhoDir );
			vorticity = vVortDir * fVortSize;

		}
		else
		{
			vorticity = Vector3D(0.0f, 0.0f, 0.0f);
		}


	}

	float		mRadiusSlug;
	float		mThickness;
	float		mRadiusOuter;
	Vector3D	mDirection;
};

class VortexTube : public IVorticityDistribution
{
public:
	
	VortexTube(const float & fDiameter, const float & fVariation, const float & fWidth, const int & iPeriods, const int & iLocation)
		: mRadius(0.5f * fDiameter)
		, mVariation(fVariation)
		, mWidth(fWidth)
		, mWavenumber(float(iPeriods))
		, mLocation(iLocation)
	{
	}

	virtual Vector3D	GetDomainSize() const
	{
		return Vector3D( 8.0f * mRadius, mWidth, 8.0f * mRadius );
	}

	virtual void		AssignVorticity( Vector3D & vorticity, const Vector3D & position, const Vector3D & vCenter) const
	{
		if (0 == mLocation)
		{
			const Vector3D	posRel			= position - vCenter;
			const float		fRho			= sqrtf( POW2( posRel.x ) + POW2( posRel.z ) );
			const float		fModulation		= 1.0f - mVariation * ( cosf( TWO_PI *  mWavenumber * posRel.y / mWidth ) - 1.0f );
			const float		fRadiusLocal	= mRadius * fModulation;
			if ( fRho < fRadiusLocal )
			{
				const float vortY	= 0.5f * ( cosf( PI * fRho / fRadiusLocal ) + 1);
				vorticity			= Vector3D(0.0f, vortY, 0.0f);
			}
			else
			{
				vorticity			= Vector3D(0.0f, 0.0f, 0.0f);
			}
		}
		else if ( 1 == mLocation )
		{	//center is at vCenter + (0, 0, 1) 
			const Vector3D	posRel			= position - vCenter - Vector3D( 0.0f, 0.0f, 1.0f * mRadius );
			const float		fRho			= sqrtf(POW2(posRel.x) + POW2(posRel.z));
			const float		fModulation		= 1.0f - mVariation * (cosf(TWO_PI *  mWavenumber * posRel.y / mWidth) - 1.0f);
			const float		fRadiusLocal	= mRadius * fModulation;
			if (fRho < fRadiusLocal)
			{
				const float vortY	= 0.5f * (cosf(PI * fRho / fRadiusLocal) + 1);
				vorticity			= Vector3D(0.0f, vortY, 0.0f);
			}
			else
			{
				vorticity			= Vector3D(0.0f, 0.0f, 0.0f);
			}
		}
		else if (-1 == mLocation )
		{	//center is at vCenter - (0, 0, 1)
			const Vector3D	posRel			= position - vCenter - Vector3D(0.0f, 0.0f, -1.0f * mRadius);;
			const float		fRho			= sqrtf(POW2(posRel.y) + POW2(posRel.z));
			const float		fModulation		= 1.0f - mVariation * (cosf(TWO_PI *  mWavenumber * posRel.x / mWidth) - 1.0f);
			const float		fRadiusLocal	= mRadius * fModulation;
			if (fRho < fRadiusLocal)
			{
				const float vortX	= 0.5f * ( cosf(PI * fRho / fRadiusLocal) + 1 );
				vorticity			= Vector3D(vortX, 0.0f, 0.0f);
			}
			else
			{
				vorticity			= Vector3D(0.0f, 0.0f, 0.0f);
			}
		}
	}


	float   mRadius;   ///< Maximum radius of vortex tube
	float   mVariation;   ///< Amplitude of radius variation
	float   mWidth;   ///< Spanwise width of domain. (one period length)
	float   mWavenumber;   ///< Number of full periods of spanwise variation to fit in domain
	int     mLocation;   ///< HACK: one of a few hard-coded locations of tube
};

class VortexSheet : public IVorticityDistribution
{
public:

	VortexSheet(const float & fThickness, const float & fVariation, const float & fWidth)
		: mThickness(fThickness)
		, mVariation(fVariation)
		, mWidth(fWidth)
	{
	}

	virtual Vector3D	GetDomainSize() const
	{
		return Vector3D( 14.0f * mThickness, mWidth, 14.0f * mThickness );
	}

	virtual void		AssignVorticity( Vector3D & vorticity, const Vector3D & position, const Vector3D & vCenter )
	{
		const float yOverWidth = position.y / mWidth;
		const float d = 1.0f - 0.5f * mVariation * (cosf(TWO_PI * yOverWidth) - 1.0f);
		const float zOverD = position.z / d;
		vorticity.x = 0.0f;
		const float s = sechf(zOverD);
		vorticity.y = s * s / d;
		const float t = tanhf(zOverD);
		vorticity.z = t * t * PI * mVariation * zOverD / (mWidth * d) * sinf(TWO_PI * yOverWidth);
		if (vorticity.norm2() < 0.01f)
		{   // When vorticity is small, force it to zero, to keep number of vortons down.
			vorticity = Vector3D(0.0f, 0.0f, 0.0f);
		}
	}

	float   mThickness;
	float   mVariation;
	float   mWidth;
};

class VortexNoise : public IVorticityDistribution
{
public:
	/*! \brief Initialize parameters for vortex noise

	\param shape - dimensions of box with noisy vorticity

	*/
	VortexNoise(const Vector3D & vBox)
		: mBox(vBox)
		, mAmplitude(1.0f, 1.0f, 1.0f)
	{
		if (0.0f == vBox.z)
		{   // Domain is 2D (in XY plane).
			// Make vorticity purely vertical.
			mAmplitude = Vector3D(0.0f, 0.0f, 1.0f);
		}
	}

	virtual Vector3D GetDomainSize(void) const
	{
		return mBox;
	}

	virtual void AssignVorticity(Vector3D & vorticity, const Vector3D & position, const Vector3D & vCenter) const
	{
		vorticity = RandomSpread(mAmplitude);
	}

	Vector3D    mBox;
	Vector3D    mAmplitude;
};

extern void AddCornerVortons( std::vector<Vorton> & vortons, const Vector3D & vMin, const Vector3D & vMax);
extern void AssignVorticity( std::vector<Vorton> & vortons, float fMagnitude, unsigned numVortonsMax, const IVorticityDistribution & vorticityDistribution);
extern void VortDistributionTest(float fMagnitude, unsigned numVortonsMax, const IVorticityDistribution & vorticityDistribution);

#endif