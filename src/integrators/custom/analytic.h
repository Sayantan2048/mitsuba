#pragma once
#if !defined(__MITSUBA_PLUGIN_ANALYTIC_H_)
#define __MITSUBA_PLUGIN_ANALYTIC_H_
MTS_NAMESPACE_BEGIN

class Analytic {
public:
    // Assuming shading point at (0,0,0) with shading normal (0,0,1), and v1, v2 projected on the norm-ball.
    static Float integrateEdge(const Vector &v1, const Vector &v2) {
        Float cosTheta = math::clamp(dot(v1, v2), -1.0f, 1.0f);
        Float theta = std::acos(cosTheta);
        Float d = cross(v1,v2).z;
        Float t = theta/std::sin(theta);
        t = std::isnan(t) ? 1.0f : t;
        return d * t;
    }

    static Float integrate(const Vector &a, const Vector &b, const Vector &c) {
        Vector quad[6];
        
        int index = 0;
        
        // If all vertices on the light source are visible from shade point
        if (a.z >= 0 && b.z >= 0 && c.z >= 0) {
            Vector e0 = normalize(a);
            Vector e1 = normalize(b);
            Vector e2 = normalize(c);
            return integrateEdge(e0, e1) + integrateEdge(e1, e2) + integrateEdge(e2, e0);
        }
        // If none of the verices are visible from shade point
        else if (a.z <= 0 && b.z <=0 && c.z <= 0)
            return 0.0f;
        
        // Some but not all verices are visible from shade point.
        // It is very important to clip the light triangle before normalizing. 
        // i.e instead of clipping the triangle on surface of ball, clip it in full space and then project on ball. This is because when you project on unit himisphere, the distances are scaled non-linearly. 
        
        Vector intersection;
        if (a.z >= 0)
            quad[index++] = normalize(a);
        if (linePlaneIntersect(a, b, intersection))
            quad[index++] = normalize(intersection);
        if (b.z >= 0)
            quad[index++] = normalize(b);
        if (linePlaneIntersect(b, c, intersection))
            quad[index++] = normalize(intersection);
        if (c.z >= 0)
            quad[index++] = normalize(c);
        if (linePlaneIntersect(c, a, intersection))
            quad[index++] = normalize(intersection);
        
        Float result = 0.0f;
        for (int i = 0; i < index; i++)
            result += integrateEdge(quad[i], quad[(i+1)%index]);

        return result;
    }
    
    // Find the intersection between line segment AB and plane with normal (0,0,1).
    // There are various possibilites: either A or B or both is on the palne, AB does not intersect internally.
    static bool linePlaneIntersect(const Vector &A,  const Vector &B, Vector &intersection) {
        Float eps = 1e-15;
       
        Float ABz = A.z * B.z;
        // either A or B or both is on the plane
        // or A and B both are on same side of plane.
        if (std::abs(ABz) <= eps || ABz > 0)
            return false;

        Float t = -A.z/(A.z - B.z);
        intersection.x = A.x + (A.x - B.x) * t;
        intersection.y = A.y + (A.y - B.y) * t;
        intersection.z = 0.0f;
        
        return true;
    }

};
MTS_NAMESPACE_END
#endif