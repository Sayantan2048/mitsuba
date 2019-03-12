#include <mitsuba/render/scene.h>
#include "analytic.h"
MTS_NAMESPACE_BEGIN

class DirectCvIntegrator : public SamplingIntegrator {
    public:
    DirectCvIntegrator(const Properties &props) : SamplingIntegrator(props) {
        /* Number of samples to take using the emitter sampling technique */
        m_emitterSamples = props.getSize("emitterSamples", 1);
        m_hideEmitters = props.getBoolean("hideEmitters", false);

        m_brdfSamples = m_emitterSamples;
        m_emitterSamples = 0;
        //m_brdfSamples = 0;
        Assert(m_emitterSamples + m_brdfSamples> 0);
    }

    /// Unserialize from a binary data stream
    DirectCvIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        m_emitterSamples = stream->readSize();
        m_hideEmitters = stream->readBool();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        stream->writeSize(m_emitterSamples);
        stream->writeBool(m_hideEmitters);
    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
        const RenderJob *job, int sceneResID, int sensorResID,
        int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID,
                sensorResID, samplerResID);
        if (!m_subIntegrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
            return false;

        sceneSize = scene->getBSphere().radius;
        
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
        Log(EInfo, "Running LTC control variate integrator.");
        return true;
    }

    void configure() {
        SamplingIntegrator::configure();
    }

    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);
        if (m_emitterSamples > 1)
            sampler->request2DArray(m_emitterSamples);
        if (m_brdfSamples > 1)
            sampler->request2DArray(m_brdfSamples);
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

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        Spectrum Li(0.0f);
        Point2 sample;

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
        bsdf->transform(its, thetaIncident, mInv, amplitude);
        const Spectrum specularReflectance = bsdf->getSpecularReflectance(its);
        const Spectrum diffuseReflectance = bsdf->getDiffuseReflectance(its);
        Float mInvDet = mInv.det();
        // Probably one might want to change its.shFrame to rotMat. Be careful, one must change all the local frame vectors to the new frame before setting the its.shFrame.
        its.wi = rotMat * its.toWorld(its.wi);
        its.shFrame.s = rotMat.row(0);
        its.shFrame.t = rotMat.row(1);
        its.shFrame.n = Normal(rotMat.row(2));
        
        Point2 *sampleArray;
        size_t numDirectSamples = m_emitterSamples;

        if (numDirectSamples > 1) {
            sampleArray = rRec.sampler->next2DArray(numDirectSamples);
        } else {
            sample = rRec.nextSample2D(); sampleArray = &sample;
        }

        DirectSamplingRecord dRec(its);
        Intersection hitLoc;
        if (bsdf->getType() & BSDF::ESmooth) {
            /* Only use direct illumination sampling when the surface's
               BSDF has smooth (i.e. non-Dirac delta) component */
            for (size_t i=0; i < numDirectSamples; ++i) {
                /* Estimate the direct illumination if this is requested */
                // Do not test for visibility yet.
                Spectrum valueUnhindered = scene->sampleEmitterDirect(dRec, sampleArray[i], false);
                const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

                if (dRec.pdf > 0 && // emitter is NULL when dRec.pdf is zero. 
                    emitter->isOnSurface() && // The next three condition essentially checks if the emitter is a mesh light.
                    emitter->getShape() != NULL && 
                    typeid(*(emitter->getShape())) == typeid(TriMesh)) {
                    
                    // Test for visibility
                    Ray shadowRay(its.p, dRec.d, ray.time);
                    Float notInShadow = 0;
                    if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter() && hitLoc.shape->getEmitter() == emitter)
                        notInShadow = 1;
                                
                    if (!valueUnhindered.isZero()) {
                        /* Allocate a record for querying the BSDF */
                        /* Evaluate BSDF * cos(theta) */
                        BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
                        Spectrum bsdfVal(0.0f);
                        if (notInShadow > 0.0f) {
                            bsdfVal = bsdf->eval(bRec);
                        }
                        const Spectrum bsdfValApprox = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);
                        Li += valueUnhindered * ( bsdfVal ) / (Float)m_emitterSamples;
                    }
                }
            }
        }

        size_t numBRDFSamples = m_brdfSamples;
        if (numBRDFSamples > 1) {
            sampleArray = rRec.sampler->next2DArray(numBRDFSamples);
        } else {
            sample = rRec.nextSample2D(); sampleArray = &sample;
        }

        for (size_t i=0; i < numBRDFSamples; ++i) {
            /* Sample BSDF * cos(theta) and also request the local density */
            Float bsdfPdf;

            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, sampleArray[i]);

            if (bsdfPdf > 0) {
                const Vector wo = its.toWorld(bRec.wo);

                // Trace a ray in this direction
                Ray shadowRay(its.p, wo, ray.time);
                Float notInShadow = 0;
                Spectrum value(0.0f);
                if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter()) {
                    const Emitter *tempEmitter = hitLoc.shape->getEmitter();
                    
                    if (tempEmitter != NULL &&
                        tempEmitter->isOnSurface() && // The next three condition essentially checks if the emitter is a mesh light.
                        tempEmitter->getShape() != NULL && 
                        typeid(*(tempEmitter->getShape())) == typeid(TriMesh)) {
                            // If the shadowRay hits a light source
                            notInShadow = 1;
                            value = hitLoc.Le(-shadowRay.d);
                        }
                }
                // In case the shadow ray didn't hit a light source.
                //if (notInShadow < 1)
                    value = intersectEmitter(scene, shadowRay);

                //dRec.setQuery(shadowRay, hitLoc);
                const Spectrum bsdfValApprox = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);
             
                Li += value * bsdfVal / (Float) numBRDFSamples; 
            }
        }

        //Float diff = Li.abs().average() - m_subIntegrator->Li(ray, rRec).abs().average();
        //Li += m_subIntegrator->Li(ray, rRec);
        //if (diff > 0)
            //Li.fromLinearRGB(0, 1, 0);
        //Li.clampNegative();
        return Li;
    }

    Spectrum intersectEmitter(const Scene *scene, const Ray &shadowRay) const {
        Float distance = sceneSize;
        Intersection its;
        Spectrum Li(0.0f);
        
        // Brute force ray-emitter intersection.
        // Find the closest emitter.
        for (auto emitter : scene->getEmitters()) {
            if (emitter->isOnSurface() && 
                emitter->getShape() != NULL && 
                typeid(*(emitter->getShape())) == typeid(TriMesh)) {

                const TriMesh *triMesh = static_cast<const TriMesh *>(emitter->getShape());
                const Triangle *triangles = triMesh->getTriangles();
                const Point *vertexPositions = triMesh->getVertexPositions();

                for (size_t i = 0; i < triMesh->getTriangleCount(); i++) {
                    Float u,v;// Unused barycentric coords
                    Float t; 
                    Vector n = cross(vertexPositions[triangles[i].idx[1]] - vertexPositions[triangles[i].idx[0]], 
                                vertexPositions[triangles[i].idx[2]] - vertexPositions[triangles[i].idx[1]]);

                    if (Triangle::rayIntersect(vertexPositions[triangles[i].idx[0]], vertexPositions[triangles[i].idx[1]], vertexPositions[triangles[i].idx[2]], shadowRay, u, v, t)
                        && t < distance && dot(n, -shadowRay.d) > 0) {
                        Li = emitter->getRadiance();
                    } 
                }
            }
        }

        return Li;
    }

    
    std::string toString() const {
        std::ostringstream oss;
        oss << "DirectCvIntegrator[" << endl
            << "  emitterSamples = " << m_emitterSamples << "," << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    size_t m_emitterSamples;
    size_t m_brdfSamples;
    ref<SamplingIntegrator> m_subIntegrator;
    bool m_hideEmitters;
    Float sceneSize;
};

MTS_IMPLEMENT_CLASS_S(DirectCvIntegrator, false, SamplingIntegrator);
MTS_EXPORT_PLUGIN(DirectCvIntegrator, "Ltc control variate direct illumination integrator");
MTS_NAMESPACE_END