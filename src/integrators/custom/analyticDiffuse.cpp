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

/* This is a direct illumination with diffuse component only renderer. Also does not cast any shadows.
 *
 */

class AnalyticDiffuseIntegrator : public SamplingIntegrator {
public:
    AnalyticDiffuseIntegrator(const Properties &props) : SamplingIntegrator(props) { 
        m_hideEmitters = props.getBoolean("hideEmitters", false);
    }

    /// Unserialize from a binary data stream
    AnalyticDiffuseIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        m_hideEmitters = stream->readBool();
        configure(); // Inherited/defined in ConfigurableObject::configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        stream->writeBool(m_hideEmitters);
    }

    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);
    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
        const RenderJob *job, int sceneResID, int sensorResID,
        int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID,
                sensorResID, samplerResID);
        
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

        return true;
    }
    
    Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        Spectrum Li(0.0f);
        Intersection &its = rRec.its;
        
        const Scene *scene = rRec.scene;
       
        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        if (!rRec.rayIntersect(ray)) {
            /* If no intersection could be found, possibly return
               radiance from a background emitter */
            return Spectrum(0.0f);
        }

        // We do not shade a light source.
        if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
            return its.Le(-ray.d);

        Spectrum sum(0.0f);
        
        for (auto emitter : scene->getEmitters()) {
            if (emitter->isOnSurface() && 
                emitter->getShape() != NULL && 
                typeid(*(emitter->getShape())) == typeid(TriMesh)) {
                
                const TriMesh *triMesh = static_cast<const TriMesh *>(emitter->getShape());
                const Triangle *triangles = triMesh->getTriangles();
                const Point *vertexPositions = triMesh->getVertexPositions();
        
                for (size_t i = 0; i < triMesh->getTriangleCount(); i++) {
                    // Vector between shade point to the vertices on the light source
                    Vector e0 = its.toLocal(vertexPositions[triangles[i].idx[2]] - its.p);
                    Vector e1 = its.toLocal(vertexPositions[triangles[i].idx[1]] - its.p);
                    Vector e2 = its.toLocal(vertexPositions[triangles[i].idx[0]] - its.p);

                    Float result = Analytic::integrate(e0, e1, e2);

                    // Note that each triangle is considered a light source, hence we apply single sided or double sided processing here.
                    if (true) // One sided light source
                        result = result > 0.0f ? result : 0.0f;
                    else // double sided light source
                        result = std::abs(result);
                                               
                    sum += result * emitter->getRadiance();          
                }
            }
        }
              
        const BSDF *bsdf = its.getBSDF(ray);
        return bsdf->getDiffuseReflectance(its) * sum * 0.5f * INV_PI;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "AnalyticIntegrator" << endl;
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    bool m_hideEmitters;
};

MTS_IMPLEMENT_CLASS_S(AnalyticDiffuseIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(AnalyticDiffuseIntegrator, "Analytic diffuse integrator");
MTS_NAMESPACE_END
