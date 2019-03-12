#pragma once
#if !defined(__MITSUBA_PLUGIN_ANALYTIC_H_)
#define __MITSUBA_PLUGIN_ANALYTIC_H_
MTS_NAMESPACE_BEGIN

#define MAX(a, b) ((a) > (b) ? (a) : (b))
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

     // w is expected to be unit length, in world coordinates
     // We want w in our new XY plane. i.e the plane containing normal and w is XY plane.
    static void getRotMat(const Intersection &its, const Vector &w, Float &cosTheta, Matrix3x3 &rotMat) {
        // First we need a coordinate system in which w has phi = 0(wrt new X axis) and normal is (0, 0, 1)
        cosTheta = dot(its.shFrame.n, w);
        if (cosTheta < 0) {
            rotMat.setZero();
            return;
        }

        // what happens when w == shFrame.n?
        // well in that case, the required specifications for wi are already met and we can use whatever local basis we have. 
        else if (cosTheta > 1 - 1e-10) {
            Matrix3x3 m(its.shFrame.s, its.shFrame.t, its.shFrame.n);
            m.transpose(rotMat);
            return;
        }

        Vector newX = normalize(w - its.shFrame.n * cosTheta); // row 0
        Vector newY = cross(its.shFrame.n, newX); // row 1
        //rotMat = Matrix3x3(newX, newY, its.shFrame.n);
        Matrix3x3 m(newX, newY, its.shFrame.n); // insert as columns
        m.transpose(rotMat); // make them rows
        //Log(EInfo, "%f and %f", rotMat.row(0).y, m.col(0).y);
    }

    // computes brdf * cosine
    static Spectrum approxBrdfEval(const BSDFSamplingRecord &bRec, const Matrix3x3 &mInv, const Float mInvDet, const Float &amplitude, const Spectrum &specular, const Spectrum diffuse) {
        if (Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0)
            return Spectrum(0.0f);
        
        Vector w_ = mInv * bRec.wo;
        Float length = w_.length();
        
        // Also note that if mInv is scaled by a scalar, it is still okay as length ^ 3 is proportional to mInv.det(). 
        Float jacobian = mInvDet / (length * length * length);

        return (specular * MAX(0, w_.z / length) * jacobian * amplitude + diffuse * Frame::cosTheta(bRec.wo)) * INV_PI;
    }

};
MTS_NAMESPACE_END
#endif