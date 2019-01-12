//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#ifndef DETOUROBSTACLEAVOIDANCE_H
#define DETOUROBSTACLEAVOIDANCE_H

#include "DetourCommon.h"

struct dtObstacleCircle
{
	dtFloat p[3];				///< Position of the obstacle
	dtFloat vel[3];			///< Velocity of the obstacle
	dtFloat dvel[3];			///< Velocity of the obstacle
	dtFloat rad;				///< Radius of the obstacle
	dtFloat dp[3], np[3];		///< Use for side selection during sampling.
};

struct dtObstacleSegment
{
	dtFloat p[3], q[3];		///< End points of the obstacle segment
	bool touch;
};


class dtObstacleAvoidanceDebugData
{
public:
	dtObstacleAvoidanceDebugData();
	~dtObstacleAvoidanceDebugData();
	
	bool init(const int maxSamples);
	void reset();
	void addSample(const dtFloat* vel, const dtFloat ssize, const dtFloat pen,
				   const dtFloat vpen, const dtFloat vcpen, const dtFloat spen, const dtFloat tpen);
	
	void normalizeSamples();
	
	inline int getSampleCount() const { return m_nsamples; }
	inline const dtFloat* getSampleVelocity(const int i) const { return &m_vel[i*3]; }
	inline dtFloat getSampleSize(const int i) const { return m_ssize[i]; }
	inline dtFloat getSamplePenalty(const int i) const { return m_pen[i]; }
	inline dtFloat getSampleDesiredVelocityPenalty(const int i) const { return m_vpen[i]; }
	inline dtFloat getSampleCurrentVelocityPenalty(const int i) const { return m_vcpen[i]; }
	inline dtFloat getSamplePreferredSidePenalty(const int i) const { return m_spen[i]; }
	inline dtFloat getSampleCollisionTimePenalty(const int i) const { return m_tpen[i]; }

private:
	// Explicitly disabled copy constructor and copy assignment operator.
	dtObstacleAvoidanceDebugData(const dtObstacleAvoidanceDebugData&);
	dtObstacleAvoidanceDebugData& operator=(const dtObstacleAvoidanceDebugData&);

	int m_nsamples;
	int m_maxSamples;
	dtFloat* m_vel;
	dtFloat* m_ssize;
	dtFloat* m_pen;
	dtFloat* m_vpen;
	dtFloat* m_vcpen;
	dtFloat* m_spen;
	dtFloat* m_tpen;
};

dtObstacleAvoidanceDebugData* dtAllocObstacleAvoidanceDebugData();
void dtFreeObstacleAvoidanceDebugData(dtObstacleAvoidanceDebugData* ptr);


static const int DT_MAX_PATTERN_DIVS = 32;	///< Max numver of adaptive divs.
static const int DT_MAX_PATTERN_RINGS = 4;	///< Max number of adaptive rings.

struct dtObstacleAvoidanceParams
{
	dtFloat velBias;
	dtFloat weightDesVel;
	dtFloat weightCurVel;
	dtFloat weightSide;
	dtFloat weightToi;
	dtFloat horizTime;
	unsigned char gridSize;	///< grid
	unsigned char adaptiveDivs;	///< adaptive
	unsigned char adaptiveRings;	///< adaptive
	unsigned char adaptiveDepth;	///< adaptive
};

class dtObstacleAvoidanceQuery
{
public:
	dtObstacleAvoidanceQuery();
	~dtObstacleAvoidanceQuery();
	
	bool init(const int maxCircles, const int maxSegments);
	
	void reset();

	void addCircle(const dtFloat* pos, const dtFloat rad,
				   const dtFloat* vel, const dtFloat* dvel);
				   
	void addSegment(const dtFloat* p, const dtFloat* q);

	int sampleVelocityGrid(const dtFloat* pos, const dtFloat rad, const dtFloat vmax,
						   const dtFloat* vel, const dtFloat* dvel, dtFloat* nvel,
						   const dtObstacleAvoidanceParams* params,
						   dtObstacleAvoidanceDebugData* debug = 0);

	int sampleVelocityAdaptive(const dtFloat* pos, const dtFloat rad, const dtFloat vmax,
							   const dtFloat* vel, const dtFloat* dvel, dtFloat* nvel,
							   const dtObstacleAvoidanceParams* params, 
							   dtObstacleAvoidanceDebugData* debug = 0);
	
	inline int getObstacleCircleCount() const { return m_ncircles; }
	const dtObstacleCircle* getObstacleCircle(const int i) { return &m_circles[i]; }

	inline int getObstacleSegmentCount() const { return m_nsegments; }
	const dtObstacleSegment* getObstacleSegment(const int i) { return &m_segments[i]; }

private:
	// Explicitly disabled copy constructor and copy assignment operator.
	dtObstacleAvoidanceQuery(const dtObstacleAvoidanceQuery&);
	dtObstacleAvoidanceQuery& operator=(const dtObstacleAvoidanceQuery&);

	void prepare(const dtFloat* pos, const dtFloat* dvel);

	dtFloat processSample(const dtFloat* vcand, const dtFloat cs,
						const dtFloat* pos, const dtFloat rad,
						const dtFloat* vel, const dtFloat* dvel,
						const dtFloat minPenalty,
						dtObstacleAvoidanceDebugData* debug);

	dtObstacleAvoidanceParams m_params;
	dtFloat m_invHorizTime;
	dtFloat m_vmax;
	dtFloat m_invVmax;

	int m_maxCircles;
	dtObstacleCircle* m_circles;
	int m_ncircles;

	int m_maxSegments;
	dtObstacleSegment* m_segments;
	int m_nsegments;
};

dtObstacleAvoidanceQuery* dtAllocObstacleAvoidanceQuery();
void dtFreeObstacleAvoidanceQuery(dtObstacleAvoidanceQuery* ptr);


#endif // DETOUROBSTACLEAVOIDANCE_H
