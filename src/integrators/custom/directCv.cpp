#include <mitsuba/render/scene.h>
#include "analytic.h"
MTS_NAMESPACE_BEGIN

class DirectCvIntegrator : public SamplingIntegrator {
    public:
    DirectCvIntegrator(const Properties &props) : SamplingIntegrator(props) {
        /* Number of samples to take using the emitter sampling technique */
        m_emitterSamples = props.getSize("emitterSamples", 1);
        m_brdfSamples = props.getSize("brdfSamples", 0);
        m_approxBrdfSamples = props.getSize("approxBrdfSamples", 0);
        m_uniformSamples = props.getSize("uniformSamples", 0);
        m_cosineSamples = props.getSize("cosineSamples", 0);
        m_hideEmitters = props.getBoolean("hideEmitters", false);
       
        //m_brdfSamples = 0;
        Assert(m_emitterSamples + m_brdfSamples + m_approxBrdfSamples + m_uniformSamples + m_cosineSamples > 0);
        // Either do MIS or approxBrdfSamples
        Assert((m_emitterSamples + m_brdfSamples) * m_approxBrdfSamples * m_uniformSamples * m_cosineSamples == 0);
    }

    /// Unserialize from a binary data stream
    DirectCvIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        m_emitterSamples = stream->readSize();
        m_brdfSamples = stream->readSize();
        m_approxBrdfSamples = stream->readSize();
        m_uniformSamples = stream->readSize();
        m_cosineSamples = stream->readSize();
        m_hideEmitters = stream->readBool();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        stream->writeSize(m_emitterSamples);
        stream->writeSize(m_brdfSamples);
        stream->writeSize(m_approxBrdfSamples);
        stream->writeSize(m_uniformSamples);
        stream->writeSize(m_cosineSamples);
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

        size_t sum = m_emitterSamples + m_brdfSamples;
        m_weightBRDF = 1 / (Float) m_brdfSamples;
        m_weightEmitter = 1 / (Float) m_emitterSamples;
        m_fracBRDF = m_brdfSamples / (Float) sum;
        m_fracEmitter = m_emitterSamples / (Float) sum;
    }

    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);
        if (m_emitterSamples > 1)
            sampler->request2DArray(m_emitterSamples);
        if (m_brdfSamples > 1)
            sampler->request2DArray(m_brdfSamples);
        if (m_approxBrdfSamples > 1)
            sampler->request2DArray(m_approxBrdfSamples);
        if (m_uniformSamples > 1)
            sampler->request2DArray(m_uniformSamples);
        if (m_cosineSamples > 1)
            sampler->request2DArray(m_cosineSamples);
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

        const BSDF *brdf = its.getBSDF(ray);
        brdf->transform(its, thetaIncident, mInv, amplitude);
        const Spectrum specularReflectance = brdf->getSpecularReflectance(its);
        const Spectrum diffuseReflectance = brdf->getDiffuseReflectance(its);
        Float mInvDet = mInv.det();
        // Probably one might want to change its.shFrame to rotMat. Be careful, one must change all the local frame vectors to the new frame before setting the its.shFrame.
        its.wi = rotMat * its.toWorld(its.wi);
        its.shFrame.s = rotMat.row(0);
        its.shFrame.t = rotMat.row(1);
        its.shFrame.n = Normal(rotMat.row(2));
        
        Point2 *sampleArray;
        
        if (m_emitterSamples > 1) {
            sampleArray = rRec.sampler->next2DArray(m_emitterSamples);
        } else {
            sample = rRec.nextSample2D(); sampleArray = &sample;
        }

        DirectSamplingRecord dRec(its);
        Intersection hitLoc;
        if (brdf->getType() & BSDF::ESmooth) {
            /* Only use direct illumination sampling when the surface's
               BSDF has smooth (i.e. non-Dirac delta) component */
            for (size_t i=0; i < m_emitterSamples; ++i) {
                /* Estimate the direct illumination if this is requested */
                // Do not test for visibility yet.
                Spectrum valueUnhindered = scene->sampleEmitterDirect(dRec, sampleArray[i], false); // returns radiance / pdfEmitter
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
                        const Spectrum brdfVal = brdf->eval(bRec);
                        
                        const Spectrum brdfValApprox = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);

                        Float brdfPdf = m_brdfSamples > 0 ? brdf->pdf(bRec) : 0;

                        /* Weight using the power heuristic */
                        const Float weight = miWeight(dRec.pdf * m_fracEmitter,
                                brdfPdf *  m_fracBRDF) * m_weightEmitter;

                        //Log(EInfo, "Weight %f", weight);
                        Li += valueUnhindered * (brdfVal * notInShadow - brdfValApprox) * weight;
                    }
                }
            }
        }

        if (m_brdfSamples > 1) {
            sampleArray = rRec.sampler->next2DArray(m_brdfSamples);
        } else {
            sample = rRec.nextSample2D(); sampleArray = &sample;
        }

        for (size_t i=0; i < m_brdfSamples; ++i) {
            /* Sample BSDF * cos(theta) and also request the local density */
            Float brdfPdf;

            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum brdfVal = brdf->sample(bRec, brdfPdf, sampleArray[i]);

            if (brdfPdf > 0) {
                const Vector wo = its.toWorld(bRec.wo);

                // Trace a ray in this direction
                Ray shadowRay(its.p, wo, ray.time);
                Float notInShadow = 0;
                Spectrum valueUnhindered(0.0f);
                if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter()) {
                    const Emitter *tempEmitter = hitLoc.shape->getEmitter();
                    
                    if (tempEmitter != NULL &&
                        tempEmitter->isOnSurface() && // The next three condition essentially checks if the emitter is a mesh light.
                        tempEmitter->getShape() != NULL && 
                        typeid(*(tempEmitter->getShape())) == typeid(TriMesh)) {
                            // If the shadowRay hits a light source
                            notInShadow = 1;
                            valueUnhindered = hitLoc.Le(-shadowRay.d);
                            dRec.setQuery(shadowRay, hitLoc);
                        }
                }
                // In case the shadow ray didn't hit a light source.
                if (notInShadow < 1)
                    valueUnhindered = intersectEmitter(scene, shadowRay, dRec);

                if (valueUnhindered.isZero())
                    continue;

                const Spectrum brdfValApprox = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance) / brdfPdf;

                const Float emitterPdf = m_emitterSamples > 0 ? scene->pdfEmitterDirect(dRec) : 0;
                
                const Float weight = miWeight(brdfPdf * m_fracBRDF,
                    emitterPdf * m_fracEmitter) * m_weightBRDF;

                Li += valueUnhindered * ( brdfVal * notInShadow - brdfValApprox) * weight;
            }
        }

        if (m_approxBrdfSamples > 1) {
            sampleArray = rRec.sampler->next2DArray(m_approxBrdfSamples);
        } else {
            sample = rRec.nextSample2D(); sampleArray = &sample;
        }

        Matrix3x3 m(0.0f);
        if (m_approxBrdfSamples > 0 && mInv.invert(m))
        for (size_t i=0; i < m_approxBrdfSamples; ++i) {
            /* Sample BSDF * cos(theta) and also request the local density */
            Float approxBrdfPdf;

            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum brdfValApprox = Analytic::sample(bRec, approxBrdfPdf, sampleArray[i], m, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);

            if (approxBrdfPdf > 0) {
                const Vector wo = its.toWorld(bRec.wo);

                // Trace a ray in this direction
                Ray shadowRay(its.p, wo, ray.time);
                Float notInShadow = 0;
                Spectrum valueUnhindered(0.0f);
                if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter()) {
                    const Emitter *tempEmitter = hitLoc.shape->getEmitter();
                    
                    if (tempEmitter != NULL &&
                        tempEmitter->isOnSurface() && // The next three condition essentially checks if the emitter is a mesh light.
                        tempEmitter->getShape() != NULL && 
                        typeid(*(tempEmitter->getShape())) == typeid(TriMesh)) {
                            // If the shadowRay hits a light source
                            notInShadow = 1;
                            valueUnhindered = hitLoc.Le(-shadowRay.d);
                        }
                }
                // In case the shadow ray didn't hit a light source.
                if (notInShadow < 1)
                    valueUnhindered = intersectEmitter(scene, shadowRay, dRec);

                if (valueUnhindered.isZero())
                    continue;

                const Spectrum brdfVal = brdf->eval(bRec) / approxBrdfPdf;
             
                Li += valueUnhindered * (brdfVal * notInShadow - brdfValApprox) / (Float) m_approxBrdfSamples;
            }
        }

        
        if (m_uniformSamples > 1) {
            sampleArray = rRec.sampler->next2DArray(m_uniformSamples);
        } else {
            sample = rRec.nextSample2D(); sampleArray = &sample;
        }

        for (size_t i=0; i < m_uniformSamples; ++i) {
           
            Float uniformPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            bRec.wo = warp::squareToUniformHemisphere(sampleArray[i]);
            uniformPdf = warp::squareToUniformHemispherePdf();

            if (uniformPdf > 0) {
                const Vector wo = its.toWorld(bRec.wo);

                // Trace a ray in this direction
                Ray shadowRay(its.p, wo, ray.time);
                Float notInShadow = 0;
                Spectrum valueUnhindered(0.0f);
                if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter()) {
                    const Emitter *tempEmitter = hitLoc.shape->getEmitter();
                    
                    if (tempEmitter != NULL &&
                        tempEmitter->isOnSurface() && // The next three condition essentially checks if the emitter is a mesh light.
                        tempEmitter->getShape() != NULL && 
                        typeid(*(tempEmitter->getShape())) == typeid(TriMesh)) {
                            // If the shadowRay hits a light source
                            notInShadow = 1;
                            valueUnhindered = hitLoc.Le(-shadowRay.d);
                        }
                }
                // In case the shadow ray didn't hit a light source.
                if (notInShadow < 1)
                    valueUnhindered = intersectEmitter(scene, shadowRay, dRec);

                if (valueUnhindered.isZero())
                    continue;

                const Spectrum brdfVal = brdf->eval(bRec) / uniformPdf;
                const Spectrum brdfValApprox = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance) / uniformPdf;

                Li += valueUnhindered * (brdfVal * notInShadow - brdfValApprox) / (Float) m_uniformSamples;
            }
        }

        if (m_cosineSamples > 1) {
            sampleArray = rRec.sampler->next2DArray(m_cosineSamples);
        } else {
            sample = rRec.nextSample2D(); sampleArray = &sample;
        }

        for (size_t i=0; i < m_cosineSamples; ++i) {
           
            Float cosinePdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            bRec.wo = warp::squareToCosineHemisphere(sampleArray[i]);
            cosinePdf = warp::squareToCosineHemispherePdf(bRec.wo);
            bRec.wo /= bRec.wo.length();

            if (cosinePdf > 0) {
                const Vector wo = its.toWorld(bRec.wo);

                // Trace a ray in this direction
                Ray shadowRay(its.p, wo, ray.time);
                Float notInShadow = 0;
                Spectrum valueUnhindered(0.0f);
                if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter()) {
                    const Emitter *tempEmitter = hitLoc.shape->getEmitter();
                    
                    if (tempEmitter != NULL &&
                        tempEmitter->isOnSurface() && // The next three condition essentially checks if the emitter is a mesh light.
                        tempEmitter->getShape() != NULL && 
                        typeid(*(tempEmitter->getShape())) == typeid(TriMesh)) {
                            // If the shadowRay hits a light source
                            notInShadow = 1;
                            valueUnhindered = hitLoc.Le(-shadowRay.d);
                        }
                }
                // In case the shadow ray didn't hit a light source.
                if (notInShadow < 1)
                    valueUnhindered = intersectEmitter(scene, shadowRay, dRec);

                if (valueUnhindered.isZero())
                    continue;

                const Spectrum brdfVal = brdf->eval(bRec) / cosinePdf;
                const Spectrum brdfValApprox = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance) / cosinePdf;

                Li += valueUnhindered * (brdfVal * notInShadow - brdfValApprox) / (Float) m_cosineSamples;
            }
        }

        Li += m_subIntegrator->Li(ray, rRec);
        
        //Li = Li.abs();
        Li.clampNegative();
        return Li;
    }

    Spectrum intersectEmitter(const Scene *scene, const Ray &shadowRay, DirectSamplingRecord &dRec) const {
        Float distance = sceneSize;
        Intersection its;
        Spectrum Li(0.0f);

        dRec.object = NULL;
                
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

                        // Fill dRec
                        dRec.object = emitter;
                        dRec.d = shadowRay.d;
                        dRec.n = normalize(n);
                        dRec.dist = t;
                        // Not required for MeshLights.
                        //dRec.measure = ESolidAngle;
                        //dRec.p = u * vertexPositions[triangles[i].idx[0]] + v * vertexPositions[triangles[i].idx[1]] + (1 - u - v) * vertexPositions[triangles[i].idx[2]];
                        //dRec.uv = Point2(u, v);

                        distance = t;
                    } 
                }
            }
        }

        return Li;
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        //pdfA *= pdfA; pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }
    
    std::string toString() const {
        std::ostringstream oss;
        oss << "DirectCvIntegrator[" << endl
            << "  emitterSamples = " << m_emitterSamples << "," << endl
            << "  brdfSamples = " << m_brdfSamples << "," << endl
            << "  approxBrdfSamples = " << m_approxBrdfSamples << "," << endl
            << "  cosineSamples = " << m_cosineSamples << "," << endl
            << "  uniformSamples = " << m_uniformSamples << "," << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    size_t m_emitterSamples;
    size_t m_brdfSamples;
    size_t m_approxBrdfSamples;
    size_t m_uniformSamples;
    size_t m_cosineSamples;
    Float m_fracBRDF, m_fracEmitter;
    Float m_weightBRDF, m_weightEmitter;
    ref<SamplingIntegrator> m_subIntegrator;
    bool m_hideEmitters;
    Float sceneSize;
};

MTS_IMPLEMENT_CLASS_S(DirectCvIntegrator, false, SamplingIntegrator);
MTS_EXPORT_PLUGIN(DirectCvIntegrator, "Ltc control variate direct illumination integrator");
MTS_NAMESPACE_END