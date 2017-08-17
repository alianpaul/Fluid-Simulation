#include "particleRenderer.h"
#include "Core\Math\matrix4x4.h"
#include <stdlib.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <thread>
#include <vector>

using namespace std;

struct VertexFormatPos3Tex2
{
	float tu, tv;
	float px, py, pz;
};

#define VertexFormatFlags_Pos3Tex2  (GL_T2F_V3F)

//For debugging
struct VertexFormatPos3
{
	float px, py, pz;
};
#define VertexFormatFlags_Pos3 (GL_V3F)

ParticleRenderer::ParticleRenderer( const char * pParticleData, size_t stride, size_t offsetToAngVel, size_t offsetToSize )
	:mParticleData(pParticleData)
	, mStride(stride)
	, mOffsetToAngVel(offsetToAngVel)
	, mOffsetToSize(offsetToSize)
	, mVertexBuffer(0)
	, mVertexBufferCapacity(0)
	, mIndices(0)
	, mIndicesCapacity(0)
{
}

ParticleRenderer::~ParticleRenderer()
{
	delete[] mVertexBuffer;
	mVertexBuffer = 0;
	mVertexBufferCapacity = 0;
	if (mIndices != 0)
	{
		free(mIndices);
		mIndices = 0;
		mIndicesCapacity = 0;
	}
}

void ParticleRenderer::FillVertexBufferSlice(const double & timeNow, const class Matrix4x4 & viewMatrix, size_t iPclStart, size_t iPclEnd)
{
	//DEBUG PARTICLES
	VertexFormatPos3 *	pVertices	= (VertexFormatPos3 *) mVertexBuffer;
	for (unsigned iPcl = iPclStart; iPcl < iPclEnd; ++iPcl)
	{
		const char *		pPos		= mParticleData + iPcl * mStride;
		const Vector3D &	rPclPos		= *((Vector3D *)pPos);

		pVertices[iPcl].px = rPclPos.x;
		pVertices[iPcl].py = rPclPos.y;
		pVertices[iPcl].pz = rPclPos.z;
	}
	
}

void ParticleRenderer::Render(double timeNow, float timeChange, size_t numParticles)
{
	if (0 == numParticles)
	{
		return;
	}

	unsigned vertexFormatFlags	= VertexFormatFlags_Pos3 ;
	unsigned vertexFormatSize = sizeof(VertexFormatPos3);

	if ((0 != mVertexBuffer)
		&& (mVertexBufferCapacity < numParticles))
	{
		delete[] mVertexBuffer;
		mVertexBuffer			= 0;
		mVertexBufferCapacity	= 0;
	}

	//TODO: The MODELVIEW Matrix for the Material

	if (0 == mVertexBuffer)
	{
		mVertexBuffer			= new unsigned char[vertexFormatSize * numParticles];
		mVertexBufferCapacity	= numParticles;
	}

	//TODO: For the Fancy Particles

	thread* threadsWorker[THREAD_NUM];
	const unsigned uStride = numParticles / THREAD_NUM;
	for (size_t iThread = 0; iThread < THREAD_NUM; iThread++)
	{
		unsigned iStart		= iThread * uStride;
		unsigned iEnd       = (iThread == THREAD_NUM - 1) ? numParticles : iStart + uStride;
		threadsWorker[iThread] = new thread(&ParticleRenderer::FillVertexBufferSlice, this, timeNow, Matrix4x4(), iStart, iEnd);
	}

	for (unsigned iThread = 0; iThread < THREAD_NUM; ++iThread)
	{
		threadsWorker[iThread]->join();
		delete threadsWorker[iThread];
	}


	glInterleavedArrays( vertexFormatFlags, vertexFormatSize, mVertexBuffer );
	glDrawArrays( GL_POINTS, 0, (GLsizei)numParticles );
}