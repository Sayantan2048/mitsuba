#include <mitsuba/render/scene.h>
#include "analytic.h"
MTS_NAMESPACE_BEGIN

class PathCv : public SamplingIntegrator {
public:
    PathCv(const Properties &props) : SamplingIntegrator(props) {
        m_explicitConnect = props.getBoolean("explicitConnect", true); // Request a non-recursive ray
        m_explicitSamples = props.getSize("explicitSamples", 1);
        m_useApproxBrdf = props.getBoolean("useApproxBrdf", false);
        m_sampleApproxBrdf = props.getBoolean("sampleApproxBrdf", false);

        Assert(!m_explicitConnect || (m_explicitConnect && m_explicitSamples > 0));
    }

    /// Unserialize from a binary data stream
    PathCv(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        m_explicitConnect = stream->readBool();
        m_explicitSamples = stream->readSize();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        stream->writeBool(m_explicitConnect);
        stream->writeSize(m_explicitSamples);
    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
        const RenderJob *job, int sceneResID, int sensorResID,
        int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID,
                sensorResID, samplerResID);
        if (!m_subIntegrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
            return false;

        sceneSize = scene->getBSphere().radius * 50;
        
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
        Log(EInfo, "Running LTC control variate path integrator.");
        return true;
    }

    void configure() {
        SamplingIntegrator::configure();
        size_t explicitSamples = m_explicitConnect ? m_explicitSamples : 0;
        size_t sum = 1 + explicitSamples;
        m_weightExplicit = 1 / (Float) explicitSamples;
        m_weightImplicit = 1;
        m_fracExplicit = explicitSamples / (Float) sum;
        m_fracImplicit = 1 / (Float) sum;
    }

    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);
        
        if (m_explicitSamples > 1)
            sampler->request2DArray(m_explicitSamples);
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

    /*
    void renderBlock(const Scene *scene,
        const Sensor *sensor, Sampler *sampler, ImageBlock *block,
        const bool &stop, const std::vector< TPoint2<uint8_t> > &points) const {

        Float diffScaleFactor = 1.0f /
            std::sqrt((Float) sampler->getSampleCount());

        bool needsApertureSample = sensor->needsApertureSample();
        bool needsTimeSample = sensor->needsTimeSample();

        RadianceQueryRecord rRec(scene, sampler);
        Point2 apertureSample(0.5f);
        Float timeSample = 0.5f;
        RayDifferential sensorRay;

        block->clear();

        Spectrum *liArray = new Spectrum[sampler->getSampleCount()];
        Point2 *samplePos = new Point2[sampler->getSampleCount()];

        uint32_t queryType = RadianceQueryRecord::ESensorRay;

        if (!sensor->getFilm()->hasAlpha()) // Don't compute an alpha channel if we don't have to
            queryType &= ~RadianceQueryRecord::EOpacity;

        for (size_t i = 0; i<points.size(); ++i) {
            Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
            if (stop)
                break;

            rRec.pixelPosition = offset;
            rRec.directLtc = Spectrum(0.0f);
            rRec.stochasticLtc = Spectrum(0.0f);
            rRec.notInShadow = 0; 
            
            sampler->generate(offset);

            for (size_t j = 0; j<sampler->getSampleCount(); j++) {
                rRec.newQuery(queryType, sensor->getMedium());
                samplePos[j] = Point2(Point2(offset) + Vector2(rRec.nextSample2D()));
                rRec.sampleIndex = (int)j;
                if (needsApertureSample)
                    apertureSample = rRec.nextSample2D();
                if (needsTimeSample)
                    timeSample = rRec.nextSample1D();

                Spectrum spec = sensor->sampleRayDifferential(
                    sensorRay, samplePos[j], apertureSample, timeSample);

                sensorRay.scaleDifferential(diffScaleFactor);
                Li(sensorRay, rRec);
                //spec *= rRec.notInShadow;
                liArray[j] = spec;
                sampler->advance();
            }

            //Log(EInfo, "%f", rRec.notInShadow / sampler->getSampleCount());
            for (size_t j = 0; j<sampler->getSampleCount(); j++) {
                block->put(samplePos[j],  0.01 * Spectrum( rRec.notInShadow / (sampler->getSampleCount() * (1 + m_explicitConnect ? m_explicitSamples : 0))), rRec.alpha);
                
            }
        }
    }*/

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);

        // return immediately if the camera ray doesn't intersect the scene
        if (!rRec.rayIntersect(ray))
            return Spectrum(0.0f);

        // return Le if camera ray intersects a light source
         if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance))
             return its.Le(-ray.d);

        Spectrum throughput(1.0f);
        Spectrum accumulate(0.0f);
        Float terminationProbability = 0.05f;
        DirectSamplingRecord dRec(its);
        Point2 sample;
        Intersection hitLoc;
        int bounce = 0;

        while(true) {
            // Do some setup
            Matrix3x3 rotMat;
            Float cosThetaIncident;
            Analytic::getRotMat(its, -ray.d, cosThetaIncident, rotMat);

            if (cosThetaIncident < 0)
                break;
    
            Float thetaIncident = std::acos(cosThetaIncident);
            Float amplitude;
            Matrix3x3 mInv;

            const BSDF *brdf = its.getBSDF(ray);
            brdf->transform(its, thetaIncident, mInv, amplitude);
            const Spectrum specularReflectance = brdf->getSpecularReflectance(its);
            const Spectrum diffuseReflectance = brdf->getDiffuseReflectance(its);
            Float mInvDet = mInv.det();
            Matrix3x3 m(0.0f);

            // No point in proceeding furthur if cannot sample implicit ray according to approx-brdf.
            if (!mInv.invert(m))
                break;
            
            // Probably one might want to change its.shFrame to rotMat. Be careful, one must change all the local frame vectors to the new frame before setting the its.shFrame.
            its.wi = rotMat * its.toWorld(its.wi);
            its.shFrame.s = rotMat.row(0);
            its.shFrame.t = rotMat.row(1);
            its.shFrame.n = Normal(rotMat.row(2));

            if (!(brdf->getType() & BSDF::ESmooth))
                Log(EError, "Delta brdfs are not supported.");

            // Emitter sampling
if (m_explicitConnect) {            
            Point2 *sampleArray;
        
            if (m_explicitSamples > 1) {
                sampleArray = rRec.sampler->next2DArray(m_explicitSamples);
            } else {
                sample = rRec.nextSample2D(); sampleArray = &sample;
            }

            for (size_t i=0; i < m_explicitSamples; ++i) {            
                const Spectrum valueUnhinderedFirstHit = scene->sampleEmitterDirect(dRec, sampleArray[i], false);
                const Emitter *emitter = static_cast<const Emitter *>(dRec.object);
                if (dRec.pdf >= Epsilon && // emitter is NULL when dRec.pdf is zero. 
                    emitter->isOnSurface() && // The next three condition essentially checks if the emitter is a mesh light.
                    emitter->getShape() != NULL && 
                    typeid(*(emitter->getShape())) == typeid(TriMesh)) {
                    
                    // Test for visibility
                    Ray shadowRay(its.p, dRec.d, ray.time);
                    Spectrum valueUnhinderedAll(0.0f);
                    // sum of radiance from all light sources igonoring all blockers and occlusion by light source itself.
                    intersectEmitter(scene, shadowRay, valueUnhinderedAll);
                    valueUnhinderedAll /= dRec.pdf;

                    if (valueUnhinderedAll.isZero())
                        continue;

                    // Check if the shade point is in shadow.
                    Float notInShadow = 0;
                    if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter() && hitLoc.shape->getEmitter() == emitter)
                        notInShadow = 1;
                    
                    if (bounce == 0)
                        rRec.notInShadow += notInShadow * (valueUnhinderedAll.getLuminance() > Epsilon ? 1.0f: 0.0f);
                                   
                    /* Allocate a record for querying the BSDF */
                    /* Evaluate BSDF * cos(theta) */
                    BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
                    const Spectrum brdfVal = m_useApproxBrdf ? Spectrum(0.0f) : brdf->eval(bRec);
                    const Spectrum brdfValApprox = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);
                    const Float brdfPdf = m_sampleApproxBrdf ? Analytic::pdf(bRec, mInv, mInvDet, specularReflectance, diffuseReflectance) : brdf->pdf(bRec);

                    const Float explicitMisWeight = miWeight(dRec.pdf * m_fracExplicit,
                            brdfPdf *  m_fracImplicit) * m_weightExplicit;

                    accumulate += throughput * valueUnhinderedFirstHit * notInShadow * (m_useApproxBrdf ? brdfValApprox : brdfVal) * explicitMisWeight;

                    /*
                    if (bounce > 0)
                        accumulate += -throughput * valueUnhinderedAll * brdfValApprox * explicitMisWeight;
                    else {
                        rRec.stochasticLtc += valueUnhinderedAll * brdfValApprox * explicitMisWeight;
                        rRec.notInShadow += notInShadow * (valueUnhinderedAll.getLuminance() > Epsilon ? 1.0f: 0.0f);
                    }*/
                }
            }
}
            // Get a sample for recursive connection
            sample = rRec.nextSample2D();

            // Brdf sampling    
            Float brdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum brdfValSampled = m_sampleApproxBrdf ? Analytic::sample(bRec, brdfPdf, sample, m, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance) :
                brdf->sample(bRec, brdfPdf, sample); // brdfeval/pdf
            
            if (brdfPdf <= Epsilon)
                break;

            // Set recursive ray
            ray = Ray(its.p, its.toWorld(bRec.wo) , ray.time);
            if (!(rRec.type & RadianceQueryRecord::EIntersection))
                rRec.type ^= RadianceQueryRecord::EIntersection;
            
            if (!rRec.rayIntersect(ray))
                break;
            
            dRec = DirectSamplingRecord(its);
            
            //if (bounce > 0)
                //accumulate += throughput * m_subIntegrator->Li(ray, rRec);
            
            const Spectrum brdfVal = m_useApproxBrdf ?  Spectrum(0.0f) : (m_sampleApproxBrdf ? brdf->eval(bRec) / brdfPdf : brdfValSampled);
            const Spectrum brdfValApprox = m_sampleApproxBrdf ? brdfValSampled : Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance) / brdfPdf;

            throughput *= (m_useApproxBrdf ? brdfValApprox : brdfVal);

            Spectrum valueUnhinderedAll(0.0f);
            // Note that this intersection with the light source is with shadowRay(same as recursive ray), so do not put the intersection before setting recursive ray.
            intersectEmitter(scene, ray, dRec, valueUnhinderedAll);

            if (!valueUnhinderedAll.isZero()) {
                Spectrum valueDirect(0.0f);
                if (its.isEmitter()) {
                    const Emitter *tempEmitter = its.shape->getEmitter();
                    
                    if (tempEmitter != NULL &&
                        tempEmitter->isOnSurface() && // The next three condition essentially checks if the emitter is a mesh light.
                        tempEmitter->getShape() != NULL && 
                        typeid(*(tempEmitter->getShape())) == typeid(TriMesh)) {
                            // If the shadowRay hits a light source
                            valueDirect = its.Le(-ray.d);
                            dRec.setQuery(ray, its);
                        }
                }

                Float explicitPdf = m_explicitConnect && dRec.object != NULL ? scene->pdfEmitterDirect(dRec) : 0;
                Float implicitMisWeight = miWeight(brdfPdf * m_fracImplicit,
                    explicitPdf * m_fracExplicit) * m_weightImplicit;

                accumulate += throughput * valueDirect * (m_useApproxBrdf ? brdfValApprox : brdfVal) * implicitMisWeight;
            }

            if (its.isEmitter())
                break;
/*
            
            

            if (bounce > 0)
                accumulate += -throughput * (valueUnhinderedAll * brdfValApprox * implicitMisWeight);
            else
                rRec.stochasticLtc += valueUnhinderedAll * brdfValApprox * implicitMisWeight;
*/

            
             // Russian Roulette path termination
            if (rRec.nextSample1D() < terminationProbability)
				break;
           
            throughput /= (1 - terminationProbability);

            if (m_explicitSamples > 1)
                 rRec.sampler->request2DArray(m_explicitSamples);
           
            bounce++;
       }
    
       //accumulate.clampNegative();
       return accumulate;
    }
# if 0
    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        Spectrum throughput(1.0f);
        Spectrum accumulate(0.0f);
        Float terminationProbability = 0.05f;
        Point2 sample;
        Intersection hitLoc;
        int bounce = 0;

        while(true) {
            /* Perform the first ray intersection (or ignore if the
            intersection has already been provided). */
            if (!rRec.rayIntersect(ray))
                break;

             // We do not shade a light source.
            if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)) {
                accumulate += throughput * its.Le(-ray.d);
                break;
            }

            // Do some setup
            Matrix3x3 rotMat;
            Float cosThetaIncident;
            Analytic::getRotMat(its, -ray.d, cosThetaIncident, rotMat);

            if (cosThetaIncident < 0)
                break;
    
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

            DirectSamplingRecord dRec(its);

            // Get a sample for recursive path
            sample = rRec.nextSample2D();

            // Russian Roulette path termination
            if (sample.y < terminationProbability)
				break;
            else {
                throughput /= (1 - terminationProbability);
                sample.y = (sample.y - terminationProbability) / (1 - terminationProbability);
            }

            if (!(brdf->getType() & BSDF::ESmooth))
                Log(EError, "Delta brdfs are not supported.");

            // Brdf sampling    
            Float brdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            const Spectrum brdfVal = brdf->sample(bRec, brdfPdf, sample); // brdfeval/pdf
                      
            if (brdfPdf <= Epsilon)
                break;
            
            accumulate += throughput * m_subIntegrator->Li(ray, rRec);
            
            const Spectrum brdfValApprox = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance) / brdfPdf;  
            // Set recursive ray
            ray = Ray(its.p, its.toWorld(bRec.wo) , ray.time);
            if (!(rRec.type & RadianceQueryRecord::EIntersection))
                rRec.type ^= RadianceQueryRecord::EIntersection;
            
            Spectrum valueUnhinderedAll(0.0f);
            // Note that this intersection with the light source is with shadowRay(same as recursive ray), so do not put the intersection before setting recursive ray.
            intersectEmitter(scene, ray, dRec, valueUnhinderedAll);
            
            Float explicitPdf = m_explicitConnect && dRec.object != NULL ? scene->pdfEmitterDirect(dRec) : 0;
           
            Float implicitMisWeight = miWeight(brdfPdf * m_fracImplicit,
                    explicitPdf * m_fracExplicit) * m_weightImplicit;

            accumulate += -throughput * (valueUnhinderedAll * brdfValApprox * implicitMisWeight);

            // Emitter sampling
if (m_explicitConnect) {            
            sample = rRec.nextSample2D();
            
            const Spectrum valueUnhinderedFirstHit = scene->sampleEmitterDirect(dRec, sample, false);
            const Emitter *emitter = static_cast<const Emitter *>(dRec.object);
            if (dRec.pdf >= Epsilon && // emitter is NULL when dRec.pdf is zero. 
                emitter->isOnSurface() && // The next three condition essentially checks if the emitter is a mesh light.
                emitter->getShape() != NULL && 
                typeid(*(emitter->getShape())) == typeid(TriMesh)) {
                    
                // Test for visibility
                Ray shadowRay(its.p, dRec.d, ray.time);
                Spectrum valueUnhinderedAll(0.0f);
                // sum of radiance from all light sources igonoring all blockers and occlusion by light source itself.
                intersectEmitter(scene, shadowRay, valueUnhinderedAll);
                valueUnhinderedAll /= dRec.pdf;

                // Check if the shade point is in shadow.
                Float notInShadow = 0;
                if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter() && hitLoc.shape->getEmitter() == emitter)
                    notInShadow = 1;
                                
                /* Allocate a record for querying the BSDF */
                /* Evaluate BSDF * cos(theta) */
                BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
                const Spectrum brdfVal = brdf->eval(bRec);
                const Spectrum brdfValApprox = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);

                Float brdfPdf = brdf->pdf(bRec);

                const Float explicitMisWeight = miWeight(dRec.pdf * m_fracExplicit,
                            brdfPdf *  m_fracImplicit) * m_weightExplicit;

                accumulate += throughput * valueUnhinderedFirstHit * notInShadow * brdfVal * explicitMisWeight;
                accumulate += -throughput * valueUnhinderedAll * brdfValApprox * explicitMisWeight; 

            }
}
            throughput *= brdfVal * implicitMisWeight;
            bounce++;
       }
    
       //accumulate.clampNegative();
       return accumulate;
    }
#endif
    Spectrum intersectEmitter(const Scene *scene, const Ray &shadowRay, DirectSamplingRecord &dRec, Spectrum &LiAllHit) const {
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
                        && dot(n, -shadowRay.d) > 0) {
                        
                        LiAllHit += emitter->getRadiance();
                        if (t < distance) {
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
        }

        return Li;
    }

    void intersectEmitter(const Scene *scene, const Ray &shadowRay, Spectrum &LiAllHit) const {
        Intersection its;
                        
        // Brute force ray-emitter intersection.
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
                        && dot(n, -shadowRay.d) > 0) {
                        
                        LiAllHit += emitter->getRadiance();
                       
                    } 
                }
            }
        }
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        //pdfA *= pdfA; pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }
    
    std::string toString() const {
        std::ostringstream oss;
        oss << "PathCv[" << endl
            << "  explicitConnect = " << m_explicitConnect << "," << endl
            << "  explicitSamples = " << m_explicitSamples << "," << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    bool m_explicitConnect;
    bool m_useApproxBrdf, m_sampleApproxBrdf;
    size_t m_explicitSamples;
    Float m_fracExplicit, m_fracImplicit;
    Float m_weightExplicit, m_weightImplicit;
    ref<SamplingIntegrator> m_subIntegrator;
    Float sceneSize;
};

MTS_IMPLEMENT_CLASS_S(PathCv, false, SamplingIntegrator);
MTS_EXPORT_PLUGIN(PathCv, "Ltc control variate path tracing integrator");
MTS_NAMESPACE_END