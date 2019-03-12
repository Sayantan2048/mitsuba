/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include "analytic.h"
MTS_NAMESPACE_BEGIN

/* This is a direct illumination with approximate specular + exact diffuse renderer. Also does not cast any shadows.
 * Also when comparing with regular microfacet distributions, one should disable the fresnel term or set it to 1. Otherwise there
 * may be issues with the intensity. However the specular highlights seem very close with or without disabling fresnel term.
 */

class LtcIntegrator : public SamplingIntegrator {
public:
    LtcIntegrator(const Properties &props) : SamplingIntegrator(props) {
        m_hideEmitters = props.getBoolean("hideEmitters", false);
    }

    /// Unserialize from a binary data stream
    LtcIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        m_hideEmitters = stream->readBool();
        m_subIntegrator = static_cast<SamplingIntegrator *>(manager->getInstance(stream));
        configure(); // Inherited/defined in ConfigurableObject::configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        stream->writeBool(m_hideEmitters);
        manager->serialize(stream, m_subIntegrator.get());
    }

    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);
        m_subIntegrator->configureSampler(scene, sampler);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        const Class *cClass = child->getClass();

        if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
            if (!cClass->derivesFrom(MTS_CLASS(SamplingIntegrator)))
                Log(EError, "The sub-integrator must be derived from the class SamplingIntegrator");
            m_subIntegrator = static_cast<SamplingIntegrator *>(child);
        } else {
            Integrator::addChild(name, child);
        }
    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
        const RenderJob *job, int sceneResID, int sensorResID,
        int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID,
                sensorResID, samplerResID);
        if (!m_subIntegrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
            return false;

         auto emitters = scene->getEmitters();

        for (auto emitter : emitters) {
            if (!emitter->isOnSurface())
                Log(EInfo, "Ignoring light sources other than area light.");
            else {
                if (emitter->getShape() == NULL)
                    Log(EInfo, "Ignoring emitter with no shape.");
                else if (typeid(*(emitter->getShape())) != typeid(TriMesh))
                    Log(EInfo, "Ignoring emitter geometry other than TriMesh. RectMesh is possible but not yet supported.");
            }
        }
        Log(EInfo, "Running LTC integrator.");
        return true;
    }

    Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
        Intersection &its = rRec.its;
        const Scene *scene = rRec.scene;
       
        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        if (!rRec.rayIntersect(ray)) {
            /* If no intersection could be found, possibly return
               radiance from a background emitter or return zero.*/
            return Spectrum(0.0f);
        }

         // We do not shade a light source.
        if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
            return its.Le(-ray.d);

        Matrix3x3 rotMat;
        Float cosThetaIncident;
        Analytic::getRotMat(its, -ray.d, cosThetaIncident, rotMat);
       
        if (cosThetaIncident < 0)
            return Spectrum(0.0f);
        
        Float thetaIncident = std::acos(cosThetaIncident);
        Float amplitude;
        Matrix3x3 mInv;

        const BSDF *bsdf = its.getBSDF(ray);
        Spectrum diffuseLi = m_subIntegrator->Li(ray, rRec);

        try {
            bsdf->transform(its, thetaIncident, mInv, amplitude);
        }
        catch(...) {
            Log(EInfo, "Caught except");
            return diffuseLi;
        }
       
        Spectrum sum(0.0f);
        Matrix3x3 mInv_rot = mInv * rotMat;
        for (auto emitter : scene->getEmitters()) {
            if (emitter->isOnSurface() && 
                emitter->getShape() != NULL && 
                typeid(*(emitter->getShape())) == typeid(TriMesh)) {
                
                const TriMesh *triMesh = static_cast<const TriMesh *>(emitter->getShape());
                const Triangle *triangles = triMesh->getTriangles();
                const Point *vertexPositions = triMesh->getVertexPositions();
        
                for (size_t i = 0; i < triMesh->getTriangleCount(); i++) {
                    Float result;

                    if (true) {
                    // Vector between shade point to the vertices on the light source
                    // in specially designed local coordinate of shade point.
                        Vector e0 = (rotMat * (vertexPositions[triangles[i].idx[2]] - its.p));
                        Vector e1 = (rotMat * (vertexPositions[triangles[i].idx[1]] - its.p));
                        Vector e2 = (rotMat * (vertexPositions[triangles[i].idx[0]] - its.p));

                        //Log(EInfo, "%f", dot(normalize(e0), normalize(mInv*e0)));

                        // transform to bsdf space. Note that uniform scaling of mInv does not affect the integration result.
                        // Also clipping of polygons shouldn't be affected by uniform scaling as points above horizon always stays above horizon, no matter how far or close 
                        // they are from the shade point. In other words, solid angle formed by the polygon on hemisphere does not change with uniform scaling, hence result
                        // stays same.
                        e0 = mInv * e0;
                        e1 = mInv * e1;
                        e2 = mInv * e2;

                        result = Analytic::integrate(e0, e1, e2);
                    }
                    // same as other branch but more optimized and less readable.
                    else {
                        Vector e0 = (mInv_rot * (vertexPositions[triangles[i].idx[2]] - its.p));
                        Vector e1 = (mInv_rot * (vertexPositions[triangles[i].idx[1]] - its.p));
                        Vector e2 = (mInv_rot * (vertexPositions[triangles[i].idx[0]] - its.p));

                        result = Analytic::integrate(e0, e1, e2);
                    }

                    // Note that each triangle is considered a light source, hence we apply single sided or double sided processing here.
                    if (true) // One sided light source
                        result = result > 0.0f ? result : 0.0f;
                    else // double sided light source
                        result = std::abs(result);
                    
                    sum += result * emitter->getRadiance();          
                }
            }
        }
        
        return diffuseLi + bsdf->getSpecularReflectance(its) * sum * amplitude * 0.5f * INV_PI;    
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "LtcIntegrator" << endl;
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ref<SamplingIntegrator> m_subIntegrator;
    bool m_hideEmitters;
};

MTS_IMPLEMENT_CLASS_S(LtcIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(LtcIntegrator, "Linearly transformed cosines integrator");
MTS_NAMESPACE_END
