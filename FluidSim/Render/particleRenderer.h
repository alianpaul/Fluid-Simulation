
class Matrix4x4;

class ParticleRenderer
{
public:
	struct ParticleIndex
	{
		int		mPcl;		///< index into particle arrays
		float	mDepth;		///< distance along view forward direction
	};
	 
	ParticleRenderer( const char * pParticleData, size_t stride, size_t offsetToAngVel, size_t offsetToSize );
	~ParticleRenderer();

	void Render( double timeNow, float timeChange, size_t numParticles );

	void SetParticleData(const char * pParticleData) { mParticleData = pParticleData; }

private:

	void FillVertexBufferSlice(const double & timeNow, const class Matrix4x4 & viewMatrix, size_t iPclStart, size_t iPclEnd);

	const char *                mParticleData;   ///< Dynamic array of particles
	size_t                      mStride;   ///< Number of bytes between particles
	size_t                      mOffsetToAngVel;   ///< Number of bytes to angular velocity
	size_t                      mOffsetToSize;   ///< Number of bytes to size
	unsigned char *             mVertexBuffer;   ///< Address of vertex buffer
	size_t                      mVertexBufferCapacity;   ///< number of vertices this buffer can hold
	ParticleIndex *				mIndices;   ///< buffer used to sort particles
	size_t                      mIndicesCapacity;   ///< number of elements in mIndices
};