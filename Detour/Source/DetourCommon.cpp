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

#include "DetourCommon.h"
#include "DetourMath.h"

//////////////////////////////////////////////////////////////////////////////////////////

void dtClosestPtPointTriangle(dtFloat* closest, const dtFloat* p,
							  const dtFloat* a, const dtFloat* b, const dtFloat* c)
{
	// Check if P in vertex region outside A
	dtFloat ab[3], ac[3], ap[3];
	dtVsub(ab, b, a);
	dtVsub(ac, c, a);
	dtVsub(ap, p, a);
	dtFloat d1 = dtVdot(ab, ap);
	dtFloat d2 = dtVdot(ac, ap);
	if (d1 <= 0.0f && d2 <= 0.0f)
	{
		// barycentric coordinates (1,0,0)
		dtVcopy(closest, a);
		return;
	}
	
	// Check if P in vertex region outside B
	dtFloat bp[3];
	dtVsub(bp, p, b);
	dtFloat d3 = dtVdot(ab, bp);
	dtFloat d4 = dtVdot(ac, bp);
	if (d3 >= 0.0f && d4 <= d3)
	{
		// barycentric coordinates (0,1,0)
		dtVcopy(closest, b);
		return;
	}
	
	// Check if P in edge region of AB, if so return projection of P onto AB
	dtFloat vc = d1*d4 - d3*d2;
	if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
	{
		// barycentric coordinates (1-v,v,0)
		dtFloat v = d1 / (d1 - d3);
		closest[0] = a[0] + v * ab[0];
		closest[1] = a[1] + v * ab[1];
		closest[2] = a[2] + v * ab[2];
		return;
	}
	
	// Check if P in vertex region outside C
	dtFloat cp[3];
	dtVsub(cp, p, c);
	dtFloat d5 = dtVdot(ab, cp);
	dtFloat d6 = dtVdot(ac, cp);
	if (d6 >= 0.0f && d5 <= d6)
	{
		// barycentric coordinates (0,0,1)
		dtVcopy(closest, c);
		return;
	}
	
	// Check if P in edge region of AC, if so return projection of P onto AC
	dtFloat vb = d5*d2 - d1*d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
	{
		// barycentric coordinates (1-w,0,w)
		dtFloat w = d2 / (d2 - d6);
		closest[0] = a[0] + w * ac[0];
		closest[1] = a[1] + w * ac[1];
		closest[2] = a[2] + w * ac[2];
		return;
	}
	
	// Check if P in edge region of BC, if so return projection of P onto BC
	dtFloat va = d3*d6 - d5*d4;
	if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f)
	{
		// barycentric coordinates (0,1-w,w)
		dtFloat w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		closest[0] = b[0] + w * (c[0] - b[0]);
		closest[1] = b[1] + w * (c[1] - b[1]);
		closest[2] = b[2] + w * (c[2] - b[2]);
		return;
	}
	
	// P inside face region. Compute Q through its barycentric coordinates (u,v,w)
	dtFloat denom = 1.0f / (va + vb + vc);
	dtFloat v = vb * denom;
	dtFloat w = vc * denom;
	closest[0] = a[0] + ab[0] * v + ac[0] * w;
	closest[1] = a[1] + ab[1] * v + ac[1] * w;
	closest[2] = a[2] + ab[2] * v + ac[2] * w;
}

bool dtIntersectSegmentPoly2D(const dtFloat* p0, const dtFloat* p1,
							  const dtFloat* verts, int nverts,
							  dtFloat& tmin, dtFloat& tmax,
							  int& segMin, int& segMax)
{
	static const dtFloat EPS = DT_FLT_EPSILON;
	
	tmin = 0;
	tmax = 1;
	segMin = -1;
	segMax = -1;
	
	dtFloat dir[3];
	dtVsub(dir, p1, p0);
	
	for (int i = 0, j = nverts-1; i < nverts; j=i++)
	{
		dtFloat edge[3], diff[3];
		dtVsub(edge, &verts[i*3], &verts[j*3]);
		dtVsub(diff, p0, &verts[j*3]);
		const dtFloat n = dtVperp2D(edge, diff);
		const dtFloat d = dtVperp2D(dir, edge);
		if (dtMathFabsf(d) < EPS)
		{
			// S is nearly parallel to this edge
			if (n < 0)
				return false;
			else
				continue;
		}
		const dtFloat t = n / d;
		if (d < 0)
		{
			// segment S is entering across this edge
			if (t > tmin)
			{
				tmin = t;
				segMin = j;
				// S enters after leaving polygon
				if (tmin > tmax)
					return false;
			}
		}
		else
		{
			// segment S is leaving across this edge
			if (t < tmax)
			{
				tmax = t;
				segMax = j;
				// S leaves before entering polygon
				if (tmax < tmin)
					return false;
			}
		}
	}
	
	return true;
}

dtFloat dtDistancePtSegSqr2D(const dtFloat* pt, const dtFloat* p, const dtFloat* q, dtFloat& t)
{
	dtFloat pqx = q[0] - p[0];
	dtFloat pqz = q[2] - p[2];
	dtFloat dx = pt[0] - p[0];
	dtFloat dz = pt[2] - p[2];
	dtFloat d = pqx*pqx + pqz*pqz;
	t = pqx*dx + pqz*dz;
	if (d > 0) t /= d;
	if (t < 0) t = 0;
	else if (t > 1) t = 1;
	dx = p[0] + t*pqx - pt[0];
	dz = p[2] + t*pqz - pt[2];
	return dx*dx + dz*dz;
}

void dtCalcPolyCenter(dtFloat* tc, const unsigned short* idx, int nidx, const dtFloat* verts)
{
	tc[0] = 0.0f;
	tc[1] = 0.0f;
	tc[2] = 0.0f;
	for (int j = 0; j < nidx; ++j)
	{
		const dtFloat* v = &verts[idx[j]*3];
		tc[0] += v[0];
		tc[1] += v[1];
		tc[2] += v[2];
	}
	const dtFloat s = 1.0f / nidx;
	tc[0] *= s;
	tc[1] *= s;
	tc[2] *= s;
}

bool dtClosestHeightPointTriangle(const dtFloat* p, const dtFloat* a, const dtFloat* b, const dtFloat* c, dtFloat& h)
{
	dtFloat v0[3], v1[3], v2[3];
	dtVsub(v0, c,a);
	dtVsub(v1, b,a);
	dtVsub(v2, p,a);
	
	const dtFloat dot00 = dtVdot2D(v0, v0);
	const dtFloat dot01 = dtVdot2D(v0, v1);
	const dtFloat dot02 = dtVdot2D(v0, v2);
	const dtFloat dot11 = dtVdot2D(v1, v1);
	const dtFloat dot12 = dtVdot2D(v1, v2);
	
	// Compute barycentric coordinates
	// const dtFloat invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
	const dtFloat denom = (dot00 * dot11 - dot01 * dot01);
	const dtFloat u = (dot11 * dot02 - dot01 * dot12) / denom;
	const dtFloat v = (dot00 * dot12 - dot01 * dot02) / denom;

	// The (sloppy) epsilon is needed to allow to get height of points which
	// are interpolated along the edges of the triangles.
	static const dtFloat EPS = 1e-4f;
	
	// If point lies inside the triangle, return interpolated ycoord.
	if (u >= -EPS && v >= -EPS && (u+v) <= 1+EPS)
	{
		h = a[1] + v0[1]*u + v1[1]*v;
		return true;
	}
	
	return false;
}

/// @par
///
/// All points are projected onto the xz-plane, so the y-values are ignored.
bool dtPointInPolygon(const dtFloat* pt, const dtFloat* verts, const int nverts)
{
	// TODO: Replace pnpoly with triArea2D tests?
	int i, j;
	bool c = false;
	for (i = 0, j = nverts-1; i < nverts; j = i++)
	{
		const dtFloat* vi = &verts[i*3];
		const dtFloat* vj = &verts[j*3];
		if (((vi[2] > pt[2]) != (vj[2] > pt[2])) &&
			(pt[0] < (vj[0]-vi[0]) * (pt[2]-vi[2]) / (vj[2]-vi[2]) + vi[0]) )
			c = !c;
	}
	return c;
}

bool dtDistancePtPolyEdgesSqr(const dtFloat* pt, const dtFloat* verts, const int nverts,
							  dtFloat* ed, dtFloat* et)
{
	// TODO: Replace pnpoly with triArea2D tests?
	int i, j;
	bool c = false;
	for (i = 0, j = nverts-1; i < nverts; j = i++)
	{
		const dtFloat* vi = &verts[i*3];
		const dtFloat* vj = &verts[j*3];
		if (((vi[2] > pt[2]) != (vj[2] > pt[2])) &&
			(pt[0] < (vj[0]-vi[0]) * (pt[2]-vi[2]) / (vj[2]-vi[2]) + vi[0]) )
			c = !c;
		ed[j] = dtDistancePtSegSqr2D(pt, vj, vi, et[j]);
	}
	return c;
}

static void projectPoly(const dtFloat* axis, const dtFloat* poly, const int npoly,
						dtFloat& rmin, dtFloat& rmax)
{
	rmin = rmax = dtVdot2D(axis, &poly[0]);
	for (int i = 1; i < npoly; ++i)
	{
		const dtFloat d = dtVdot2D(axis, &poly[i*3]);
		rmin = dtMin(rmin, d);
		rmax = dtMax(rmax, d);
	}
}

inline bool overlapRange(const dtFloat amin, const dtFloat amax,
						 const dtFloat bmin, const dtFloat bmax,
						 const dtFloat eps)
{
	return ((amin+eps) > bmax || (amax-eps) < bmin) ? false : true;
}

/// @par
///
/// All vertices are projected onto the xz-plane, so the y-values are ignored.
bool dtOverlapPolyPoly2D(const dtFloat* polya, const int npolya,
						 const dtFloat* polyb, const int npolyb)
{
	const dtFloat eps = 1e-4f;
	
	for (int i = 0, j = npolya-1; i < npolya; j=i++)
	{
		const dtFloat* va = &polya[j*3];
		const dtFloat* vb = &polya[i*3];
		const dtFloat n[3] = { vb[2]-va[2], 0, -(vb[0]-va[0]) };
		dtFloat amin,amax,bmin,bmax;
		projectPoly(n, polya, npolya, amin,amax);
		projectPoly(n, polyb, npolyb, bmin,bmax);
		if (!overlapRange(amin,amax, bmin,bmax, eps))
		{
			// Found separating axis
			return false;
		}
	}
	for (int i = 0, j = npolyb-1; i < npolyb; j=i++)
	{
		const dtFloat* va = &polyb[j*3];
		const dtFloat* vb = &polyb[i*3];
		const dtFloat n[3] = { vb[2]-va[2], 0, -(vb[0]-va[0]) };
		dtFloat amin,amax,bmin,bmax;
		projectPoly(n, polya, npolya, amin,amax);
		projectPoly(n, polyb, npolyb, bmin,bmax);
		if (!overlapRange(amin,amax, bmin,bmax, eps))
		{
			// Found separating axis
			return false;
		}
	}
	return true;
}

// Returns a random point in a convex polygon.
// Adapted from Graphics Gems article.
void dtRandomPointInConvexPoly(const dtFloat* pts, const int npts, dtFloat* areas,
							   const dtFloat s, const dtFloat t, dtFloat* out)
{
	// Calc triangle araes
	dtFloat areasum = 0.0f;
	for (int i = 2; i < npts; i++) {
		areas[i] = dtTriArea2D(&pts[0], &pts[(i-1)*3], &pts[i*3]);
		areasum += dtMax(0.001f, areas[i]);
	}
	// Find sub triangle weighted by area.
	const dtFloat thr = s*areasum;
	dtFloat acc = 0.0f;
	dtFloat u = 1.0f;
	int tri = npts - 1;
	for (int i = 2; i < npts; i++) {
		const dtFloat dacc = areas[i];
		if (thr >= acc && thr < (acc+dacc))
		{
			u = (thr - acc) / dacc;
			tri = i;
			break;
		}
		acc += dacc;
	}
	
	dtFloat v = dtMathSqrtf(t);
	
	const dtFloat a = 1 - v;
	const dtFloat b = (1 - u) * v;
	const dtFloat c = u * v;
	const dtFloat* pa = &pts[0];
	const dtFloat* pb = &pts[(tri-1)*3];
	const dtFloat* pc = &pts[tri*3];
	
	out[0] = a*pa[0] + b*pb[0] + c*pc[0];
	out[1] = a*pa[1] + b*pb[1] + c*pc[1];
	out[2] = a*pa[2] + b*pb[2] + c*pc[2];
}

inline dtFloat vperpXZ(const dtFloat* a, const dtFloat* b) { return a[0]*b[2] - a[2]*b[0]; }

bool dtIntersectSegSeg2D(const dtFloat* ap, const dtFloat* aq,
						 const dtFloat* bp, const dtFloat* bq,
						 dtFloat& s, dtFloat& t)
{
	dtFloat u[3], v[3], w[3];
	dtVsub(u,aq,ap);
	dtVsub(v,bq,bp);
	dtVsub(w,ap,bp);
	dtFloat d = vperpXZ(u,v);
	if (dtMathFabsf(d) < DT_FLT_EPSILON) return false;
	s = vperpXZ(v,w) / d;
	t = vperpXZ(u,w) / d;
	return true;
}

