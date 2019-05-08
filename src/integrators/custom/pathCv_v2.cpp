#include <mitsuba/render/scene.h>
#include "analytic.h"
MTS_NAMESPACE_BEGIN

// The v2 version only uses approx-brdf and approx-brdf sampling.
// There are three #if's, first #if is vanilla explicit path tracing, second #if uses emitter-brdf joint sampling as explicit connection, last #if uses brdf-sampling as explicit connection.
// The inner #if is to switch between control variate and ratio estimator. 

class PathCv_v2 : public SamplingIntegrator {
public:
    PathCv_v2(const Properties &props) : SamplingIntegrator(props) {
        m_explicitConnect = props.getBoolean("explicitConnect", true); // Request a non-recursive ray
        m_collectAll = props.getBoolean("collectAll", true);
        m_whichBounce = props.getSize("whichBounce", 0);
        m_explicitSamples = props.getSize("explicitSamples", 1);
        m_emitterSamples = m_explicitSamples;
        Assert(!m_explicitConnect || (m_explicitConnect && m_explicitSamples > 0));
    }

    /// Unserialize from a binary data stream
    PathCv_v2(Stream *stream, InstanceManager *manager)
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

        m_triEmitters = 0; 
        for (auto emitter : scene->getEmitters())
            if (emitter->isOnSurface() && 
                emitter->getShape() != NULL && 
                typeid(*(emitter->getShape())) == typeid(TriMesh)) {
                
                const TriMesh *triMesh = static_cast<const TriMesh *>(emitter->getShape());
                size_t triEmitters = triMesh->getTriangleCount();
                m_triEmitters += triEmitters;

                for (size_t i = 0; i < triEmitters; i++)
                    m_triEmitterList.push_back(emitter);
            }
        

        Log(EInfo, "Running LTC control variate path integrator v2.");
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
    }

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

        rRec.triEmitterAreaList = new Float[2 * m_triEmitters];
        rRec.triEmitterAreaLumList = new Float[2 * m_triEmitters];
        rRec.triEmitterSurfaceNormalBuffer = new Vector[2 * m_triEmitters];
        rRec.triEmitterVertexBuffer = new Vector[2 * m_triEmitters * 3];
        rRec.triEmitterRadianceBuffer = new Spectrum[m_triEmitters];
        
        rRec.emitterSampleValues = new Spectrum[m_emitterSamples];
        rRec.emitterSampleDirections = new Vector[m_emitterSamples];
        rRec.emitterPdf = new Float[m_emitterSamples];
        rRec.emitterIndices = new size_t[m_emitterSamples];

        for (size_t i = 0; i < m_triEmitters; i++) {
            rRec.triEmitterAreaList[i] = 0;
            rRec.triEmitterAreaLumList[i] = 0;
            rRec.triEmitterSurfaceNormalBuffer[i] = Vector(0.0f);
            rRec.triEmitterSurfaceNormalBuffer[m_triEmitters + i] = Vector(0.0f);
            rRec.triEmitterVertexBuffer[3 * i] = Vector(0.0f);
            rRec.triEmitterVertexBuffer[3 * i + 1] = Vector(0.0f);
            rRec.triEmitterVertexBuffer[3 * i + 2] = Vector(0.0f);
            rRec.triEmitterVertexBuffer[3 * m_triEmitters + 3 * i] = Vector(0.0f);
            rRec.triEmitterVertexBuffer[3 * m_triEmitters + 3 * i + 1] = Vector(0.0f);
            rRec.triEmitterVertexBuffer[3 * m_triEmitters + 3 * i + 2] = Vector(0.0f);
        }

        for (size_t i = 0; i < m_emitterSamples; i++) {
            rRec.emitterSampleValues[i] = Spectrum(0.0f);
            rRec.emitterSampleDirections[i] = Vector(0.0f);
            rRec.emitterPdf[i] = 0;
            rRec.emitterIndices[i] = -1;
        }

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
                spec *= Li(sensorRay, rRec);
                liArray[j] = spec;
                sampler->advance();
            }

            //Log(EInfo, "%f", rRec.notInShadow / sampler->getSampleCount());
            for (size_t j = 0; j<sampler->getSampleCount(); j++) {
                block->put(samplePos[j],  liArray[j], rRec.alpha);
            }
        }


        delete []liArray;
        delete []samplePos;

        delete []rRec.triEmitterAreaList;
        delete []rRec.triEmitterAreaLumList;
        delete []rRec.triEmitterSurfaceNormalBuffer;
        delete []rRec.triEmitterVertexBuffer;
        delete []rRec.triEmitterRadianceBuffer;

        delete []rRec.emitterSampleValues;
        delete []rRec.emitterSampleDirections;
        delete []rRec.emitterPdf;
        delete []rRec.emitterIndices;
    }

//#define RATIO

#if 1 //Standard explicit path tracing with n bounces of emitter sample and one recursive sample using approx brdf as ltc and evaluation.
#ifndef RATIO
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
        Spectrum throughputLtc(1.0f);
        Spectrum accumulate(0.0f);
        Spectrum accumulateLtc(0.0f);
        Float terminationProbability = 0.05f;
        Intersection hitLoc;
        size_t bounce = 0;
        //size_t nEmitterSamples = m_explicitSamples;

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
            const Float specularLum = specularReflectance.getLuminance();
            const Spectrum diffuseReflectance = brdf->getDiffuseReflectance(its);
            const Float diffuseLum = diffuseReflectance.getLuminance();
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

            DirectSamplingRecord dRec(its);
            
            if (m_collectAll || bounce == m_whichBounce)
                accumulateLtc += throughputLtc * Analytic::ltcIntegrate(scene, its.p, rotMat, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);

            // Emitter sampling
if (m_explicitConnect) {            
            for (size_t i=0; i < m_explicitSamples; ++i) {            
                const Spectrum valueUnhinderedFirstHit = scene->sampleEmitterDirect(dRec, rRec.nextSample2D(), false);
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
                    
                    //if (bounce == 0)
                    //    rRec.notInShadow += notInShadow * (valueUnhinderedAll.getLuminance() > Epsilon ? 1.0f: 0.0f);
                                   
                    /* Allocate a record for querying the BSDF */
                    /* Evaluate BSDF * cos(theta) */
                    BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
                    const Spectrum brdfVal = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);
                    const Float brdfPdf = Analytic::pdf(bRec, mInv, mInvDet, specularLum, diffuseLum);

                    const Float explicitMisWeight = miWeight(dRec.pdf * m_fracExplicit,
                            brdfPdf *  m_fracImplicit) * m_weightExplicit;
                    
                    if ( m_collectAll || bounce == m_whichBounce) {
                        accumulate += throughput * valueUnhinderedFirstHit * notInShadow * brdfVal * explicitMisWeight;
                        accumulateLtc -= throughputLtc * valueUnhinderedAll * brdfVal * explicitMisWeight;
                    }
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

            // Brdf sampling    
            Float brdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum brdfValSampled = Analytic::sample(bRec, brdfPdf, rRec.nextSample2D(), m, mInv, mInvDet, amplitude, specularReflectance, specularLum, diffuseReflectance, diffuseLum); // brdfeval/pdf
            
            if (brdfValSampled.isZero() || brdfPdf < Epsilon)
                break;

            // Set recursive ray
            ray = Ray(its.p, its.toWorld(bRec.wo) , ray.time);
            if (!(rRec.type & RadianceQueryRecord::EIntersection))
                rRec.type ^= RadianceQueryRecord::EIntersection;
            
            if (!rRec.rayIntersect(ray))
                break;
            
            //if (bounce > 0)
                //accumulate += throughput * m_subIntegrator->Li(ray, rRec);

            throughput *= brdfValSampled;
            throughputLtc *= brdfValSampled;

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

                Float explicitPdf = m_explicitConnect && dRec.object != NULL ? scene->pdfEmitterDirect(dRec): 0;
                Float implicitMisWeight = miWeight(brdfPdf * m_fracImplicit,
                    explicitPdf * m_fracExplicit) * m_weightImplicit;
                if (m_collectAll || bounce == m_whichBounce) {
                    accumulate += throughput * valueDirect * implicitMisWeight;
                    //accumulateLtc -= throughputLtc * valueUnhinderedAll * implicitMisWeight;
                }
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
            if (rRec.nextSample1D() <= terminationProbability)
                break;
            
            throughput /= (1 - terminationProbability);
            throughputLtc /= (1 - terminationProbability);

            bounce++;
       }
     
        return accumulate;
    }
#else
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
        Spectrum throughputLtc(1.0f);
        Spectrum accumulate(0.0f);
        Float terminationProbability = 0.05f;
        Intersection hitLoc;
        int bounce = 0;
        size_t nEmitterSamples = m_explicitSamples;

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
            const Float specularLum = specularReflectance.getLuminance();
            const Spectrum diffuseReflectance = brdf->getDiffuseReflectance(its);
            const Float diffuseLum = diffuseReflectance.getLuminance();
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

            DirectSamplingRecord dRec(its);
            
            Spectrum LiLtc = throughputLtc * Analytic::ltcIntegrate(scene, its.p, rotMat, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);
            Spectrum LiWithVisibility(0.0f);
            Spectrum LiWithoutVisibility(0.0f);

            // Emitter sampling
if (m_explicitConnect) {            
            for (size_t i=0; i < m_explicitSamples; ++i) {            
                const Spectrum valueUnhinderedFirstHit = scene->sampleEmitterDirect(dRec, rRec.nextSample2D(), false);
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
                    const Spectrum brdfVal = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);
                    const Float brdfPdf = Analytic::pdf(bRec, mInv, mInvDet, specularLum, diffuseLum);

                    const Float explicitMisWeight = miWeight(dRec.pdf * m_fracExplicit,
                            brdfPdf *  m_fracImplicit) * m_weightExplicit;
                    
                    LiWithVisibility += throughput * valueUnhinderedFirstHit * notInShadow * brdfVal * explicitMisWeight;
                    LiWithoutVisibility += throughputLtc * valueUnhinderedAll * brdfVal * explicitMisWeight;
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

            // Brdf sampling    
            Float brdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum brdfValSampled = Analytic::sample(bRec, brdfPdf, rRec.nextSample2D(), m, mInv, mInvDet, amplitude, specularReflectance, specularLum, diffuseReflectance, diffuseLum); // brdfeval/pdf
            
            if (brdfValSampled.isZero() || brdfPdf < Epsilon) {
                if (m_collectAll || bounce == m_whichBounce)
                    accumulate += getRatio(LiLtc, LiWithVisibility, LiWithoutVisibility);
                break;
            }

            // Set recursive ray
            ray = Ray(its.p, its.toWorld(bRec.wo) , ray.time);
            if (!(rRec.type & RadianceQueryRecord::EIntersection))
                rRec.type ^= RadianceQueryRecord::EIntersection;
            
            if (!rRec.rayIntersect(ray)) {
                 if (m_collectAll || bounce == m_whichBounce)
                    accumulate += getRatio(LiLtc, LiWithVisibility, LiWithoutVisibility);
                break;
            }
            
            //if (bounce > 0)
                //accumulate += throughput * m_subIntegrator->Li(ray, rRec);

            throughput *= brdfValSampled;
            throughputLtc *= brdfValSampled;

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

                Float explicitPdf = m_explicitConnect && dRec.object != NULL ? scene->pdfEmitterDirect(dRec): 0;
                Float implicitMisWeight = miWeight(brdfPdf * m_fracImplicit,
                    explicitPdf * m_fracExplicit) * m_weightImplicit;
               
                LiWithVisibility += throughput * valueDirect * implicitMisWeight;
                LiWithoutVisibility += throughputLtc * valueUnhinderedAll * implicitMisWeight;
            }

            if (m_collectAll || bounce == m_whichBounce)
                accumulate += getRatio(LiLtc, LiWithVisibility, LiWithoutVisibility);

            if (its.isEmitter())
                break;
            
/*
            
            

            if (bounce > 0)
                accumulate += -throughput * (valueUnhinderedAll * brdfValApprox * implicitMisWeight);
            else
                rRec.stochasticLtc += valueUnhinderedAll * brdfValApprox * implicitMisWeight;
*/

            // Russian Roulette path termination
            if (rRec.nextSample1D() <= terminationProbability)
                break;
            
            throughput /= (1 - terminationProbability);
            throughputLtc /= (1 - terminationProbability);

            bounce++;
       }

        return accumulate;
    }
#endif
#endif

# if 0 // Path tracing with emitter-brdf sampling as explicit connection and brdf sampling as recursive part!
#ifndef RATIO
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
        Spectrum throughputLtc(1.0f);
        Spectrum accumulate(0.0f);
        Spectrum accumulateLtc(0.0f);
        Float terminationProbability = 0.05f;
        Intersection hitLoc;
        int bounce = 0;
        size_t nEmitterSamples = m_explicitSamples;

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
            const Float specularLum = specularReflectance.getLuminance();
            const Spectrum diffuseReflectance = brdf->getDiffuseReflectance(its);
            const Float diffuseLum = diffuseReflectance.getLuminance();
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

            DirectSamplingRecord dRec(its);
            Spectrum ltcEval = Analytic::ltcIntegrateAndSample(scene, its.p, rotMat, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance, 
                m_triEmitters, rRec.triEmitterAreaList, rRec.triEmitterAreaLumList, rRec.areaNormSpecular, rRec.areaNormLumSpecular, rRec.omegaSpecular, 
                rRec.areaNormDiffuse, rRec.areaNormLumDiffuse, rRec.omegaDiffuse,
                rRec.triEmitterSurfaceNormalBuffer, rRec.triEmitterVertexBuffer, rRec.triEmitterRadianceBuffer);
            
            if (m_collectAll || bounce == m_whichBounce)
                accumulateLtc += throughputLtc * ltcEval;

            // Emitter sampling
if (m_explicitConnect) {
            Analytic::getEmitterSamples(m_triEmitters, rRec.triEmitterAreaList, rRec.triEmitterAreaLumList, rRec.areaNormSpecular, rRec.areaNormLumSpecular, rRec.omegaSpecular,
                rRec.areaNormDiffuse, rRec.areaNormLumDiffuse, rRec.omegaDiffuse, cosThetaIncident, rRec.triEmitterSurfaceNormalBuffer, rRec.triEmitterVertexBuffer, rRec.triEmitterRadianceBuffer, // some buffers to work with
                m, mInv, mInvDet, amplitude, // domain transformation
                specularLum, diffuseLum, // material properties for sampling
                rRec, // sampler
                m_emitterSamples, rRec.emitterSampleValues, rRec.emitterSampleDirections, rRec.emitterPdf, rRec.emitterIndices); // samples
 
            for (size_t i=0; i < m_explicitSamples; ++i) {            
                if (rRec.emitterPdf[i] >= Epsilon && rRec.emitterIndices[i] >= 0) {
                    
                    // Test for visibility
                    Ray shadowRay(its.p, its.toWorld(rRec.emitterSampleDirections[i]), ray.time);
                    Spectrum valueUnhinderedAll(0.0f);
                    // sum of radiance from all light sources igonoring all blockers and occlusion by light source itself.
                    intersectEmitter(scene, shadowRay, valueUnhinderedAll);
                    valueUnhinderedAll /= rRec.emitterPdf[i];

                    if (valueUnhinderedAll.isZero())
                        continue;

                    // Check if the shade point is in shadow.
                    Float notInShadow = 0;
                    if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter() && hitLoc.shape->getEmitter() == m_triEmitterList[rRec.emitterIndices[i]].get())
                        notInShadow = 1;
                    
                    if (bounce == 0)
                        rRec.notInShadow += notInShadow * (valueUnhinderedAll.getLuminance() > Epsilon ? 1.0f: 0.0f);
                                   
                    /* Allocate a record for querying the BSDF */
                    /* Evaluate BSDF * cos(theta) */
                    BSDFSamplingRecord bRec(its, rRec.emitterSampleDirections[i]);
                    const Spectrum brdfVal = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);
                    const Float brdfPdf = Analytic::pdf(bRec, mInv, mInvDet, specularLum, diffuseLum);

                    const Float explicitMisWeight = miWeight(rRec.emitterPdf[i] * m_fracExplicit,
                            brdfPdf *  m_fracImplicit) * m_weightExplicit;
                    
                    //if (std::abs(rRec.emitterSampleValues[i].getLuminance() - valueUnhinderedAll.getLuminance()) > Epsilon)
                        //Log(EInfo, "%f %f", rRec.emitterSampleValues[i].getLuminance(), valueUnhinderedAll.getLuminance());
                    
                    //Log(EInfo, "%f %f", rRec.emitterPdf[i], rRec.emitterSampleValues[i].getLuminance());
                    if (m_collectAll || bounce == m_whichBounce) {
                        accumulate += throughput * rRec.emitterSampleValues[i] * brdfVal * notInShadow * explicitMisWeight;
                        accumulateLtc -= throughputLtc * valueUnhinderedAll * brdfVal * explicitMisWeight;
                    }
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

            // Brdf sampling    
            Float brdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum brdfValSampled = Analytic::sample(bRec, brdfPdf, rRec.nextSample2D(), m, mInv, mInvDet, amplitude, specularReflectance, specularLum, diffuseReflectance, diffuseLum); // brdfeval/pdf
           
            if (brdfValSampled.isZero() || brdfPdf < Epsilon)
                break;

            // Set recursive ray
            ray = Ray(its.p, its.toWorld(bRec.wo) , ray.time);
            if (!(rRec.type & RadianceQueryRecord::EIntersection))
                rRec.type ^= RadianceQueryRecord::EIntersection;
            
            if (!rRec.rayIntersect(ray))
                break;
            
            throughput *= brdfValSampled;
            throughputLtc *= brdfValSampled;

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
                        typeid(*(tempEmitter->getShape())) == typeid(TriMesh))
                            // If the shadowRay hits a light source
                            valueDirect = its.Le(-ray.d);
                        
                }

                Float explicitPdf = m_explicitConnect && dRec.object != NULL ? Analytic::emitterPdf(bRec.wo, m_triEmitters,
                    rRec.triEmitterAreaList,  rRec.triEmitterAreaLumList, rRec.areaNormSpecular, rRec.areaNormLumSpecular, rRec.omegaSpecular, 
                    rRec.areaNormDiffuse, rRec.areaNormLumDiffuse, rRec.omegaDiffuse, cosThetaIncident, rRec.triEmitterSurfaceNormalBuffer, rRec.triEmitterVertexBuffer,
                    mInv, mInvDet, specularLum, diffuseLum): 0;
                
                Float implicitMisWeight = miWeight(brdfPdf * m_fracImplicit,
                    explicitPdf * m_fracExplicit) * m_weightImplicit;
                
                if (m_collectAll || bounce == m_whichBounce) {
                    accumulate += throughput * valueDirect * implicitMisWeight;
                    accumulateLtc -= throughputLtc * valueUnhinderedAll * implicitMisWeight;
                }
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
            if (rRec.nextSample1D() <= terminationProbability)
                break;
            
            throughput /= (1 - terminationProbability);
            throughputLtc /= (1 - terminationProbability);

           
            bounce++;
       }

        return accumulate + accumulateLtc;
    }
#else
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
        Spectrum throughputLtc(1.0f);
        Spectrum accumulate(0.0f);
        Float terminationProbability = 0.05f;
        Intersection hitLoc;
        int bounce = 0;
        size_t nEmitterSamples = m_explicitSamples;

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
            const Float specularLum = specularReflectance.getLuminance();
            const Spectrum diffuseReflectance = brdf->getDiffuseReflectance(its);
            const Float diffuseLum = diffuseReflectance.getLuminance();
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

            DirectSamplingRecord dRec(its);
            Spectrum LiLtc = throughputLtc * Analytic::ltcIntegrateAndSample(scene, its.p, rotMat, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance, 
                m_triEmitters, rRec.triEmitterAreaList, rRec.triEmitterAreaLumList, rRec.areaNormSpecular, rRec.areaNormLumSpecular, rRec.omegaSpecular, 
                rRec.areaNormDiffuse, rRec.areaNormLumDiffuse, rRec.omegaDiffuse,
                rRec.triEmitterSurfaceNormalBuffer, rRec.triEmitterVertexBuffer, rRec.triEmitterRadianceBuffer);
            
            Spectrum LiWithVisibility(0.0f);
            Spectrum LiWithoutVisibility(0.0f);

            // Emitter sampling
if (m_explicitConnect) {
            Analytic::getEmitterSamples(m_triEmitters, rRec.triEmitterAreaList, rRec.triEmitterAreaLumList, rRec.areaNormSpecular, rRec.areaNormLumSpecular, rRec.omegaSpecular,
            rRec.areaNormDiffuse, rRec.areaNormLumDiffuse, rRec.omegaDiffuse, cosThetaIncident, rRec.triEmitterSurfaceNormalBuffer, rRec.triEmitterVertexBuffer, rRec.triEmitterRadianceBuffer, // some buffers to work with
                m, mInv, mInvDet, amplitude, // domain transformation
                specularLum, diffuseLum, // material properties for sampling
                rRec, // sampler
                m_emitterSamples, rRec.emitterSampleValues, rRec.emitterSampleDirections, rRec.emitterPdf, rRec.emitterIndices); // samples
 
            for (size_t i=0; i < m_explicitSamples; ++i) {            
                if (rRec.emitterPdf[i] >= Epsilon && rRec.emitterIndices[i] >= 0) {
                    
                    // Test for visibility
                    Ray shadowRay(its.p, its.toWorld(rRec.emitterSampleDirections[i]), ray.time);
                    Spectrum valueUnhinderedAll(0.0f);
                    // sum of radiance from all light sources igonoring all blockers and occlusion by light source itself.
                    intersectEmitter(scene, shadowRay, valueUnhinderedAll);
                    valueUnhinderedAll /= rRec.emitterPdf[i];

                    if (valueUnhinderedAll.isZero())
                        continue;

                    // Check if the shade point is in shadow.
                    Float notInShadow = 0;
                    if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter() && hitLoc.shape->getEmitter() == m_triEmitterList[rRec.emitterIndices[i]].get())
                        notInShadow = 1;
                    
                    if (bounce == 0)
                        rRec.notInShadow += notInShadow * (valueUnhinderedAll.getLuminance() > Epsilon ? 1.0f: 0.0f);
                                   
                    /* Allocate a record for querying the BSDF */
                    /* Evaluate BSDF * cos(theta) */
                    BSDFSamplingRecord bRec(its, rRec.emitterSampleDirections[i]);
                    const Spectrum brdfVal = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);
                    const Float brdfPdf = Analytic::pdf(bRec, mInv, mInvDet, specularLum, diffuseLum);

                    const Float explicitMisWeight = miWeight(rRec.emitterPdf[i] * m_fracExplicit,
                            brdfPdf *  m_fracImplicit) * m_weightExplicit;
                    
                    //if (std::abs(rRec.emitterSampleValues[i].getLuminance() - valueUnhinderedAll.getLuminance()) > Epsilon)
                        //Log(EInfo, "%f %f", rRec.emitterSampleValues[i].getLuminance(), valueUnhinderedAll.getLuminance());
                    
                    //Log(EInfo, "%f %f", rRec.emitterPdf[i], rRec.emitterSampleValues[i].getLuminance());
                    LiWithVisibility += throughput * rRec.emitterSampleValues[i] * brdfVal * notInShadow * explicitMisWeight;
                    LiWithoutVisibility += throughputLtc * valueUnhinderedAll * brdfVal * explicitMisWeight;
                   
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

            // Brdf sampling    
            Float brdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum brdfValSampled = Analytic::sample(bRec, brdfPdf, rRec.nextSample2D(), m, mInv, mInvDet, amplitude, specularReflectance, specularLum, diffuseReflectance, diffuseLum); // brdfeval/pdf
           
            if (brdfValSampled.isZero() || brdfPdf < Epsilon) {
                if (m_collectAll || bounce == m_whichBounce)
                    accumulate += getRatio(LiLtc, LiWithVisibility, LiWithoutVisibility);
                break;
            }

            // Set recursive ray
            ray = Ray(its.p, its.toWorld(bRec.wo) , ray.time);
            if (!(rRec.type & RadianceQueryRecord::EIntersection))
                rRec.type ^= RadianceQueryRecord::EIntersection;
            
            if (!rRec.rayIntersect(ray)) {
                if (m_collectAll || bounce == m_whichBounce)
                    accumulate += getRatio(LiLtc, LiWithVisibility, LiWithoutVisibility);
                break;
            }
            
            throughput *= brdfValSampled;
            throughputLtc *= brdfValSampled;

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
                        typeid(*(tempEmitter->getShape())) == typeid(TriMesh))
                            // If the shadowRay hits a light source
                            valueDirect = its.Le(-ray.d);
                        
                }

                Float explicitPdf = m_explicitConnect && dRec.object != NULL ? Analytic::emitterPdf(bRec.wo, m_triEmitters,
                    rRec.triEmitterAreaList,  rRec.triEmitterAreaLumList, rRec.areaNormSpecular, rRec.areaNormLumSpecular, rRec.omegaSpecular, 
                    rRec.areaNormDiffuse, rRec.areaNormLumDiffuse, rRec.omegaDiffuse, cosThetaIncident, rRec.triEmitterSurfaceNormalBuffer, rRec.triEmitterVertexBuffer,
                    mInv, mInvDet, specularLum, diffuseLum): 0;
                
                Float implicitMisWeight = miWeight(brdfPdf * m_fracImplicit,
                    explicitPdf * m_fracExplicit) * m_weightImplicit;
                
                LiWithVisibility += throughput * valueDirect * implicitMisWeight;
                LiWithoutVisibility += throughputLtc * valueUnhinderedAll * implicitMisWeight;
                
            }
            
            if (m_collectAll || bounce == m_whichBounce)
                accumulate += getRatio(LiLtc, LiWithVisibility, LiWithoutVisibility);

            if (its.isEmitter())
                break;
/*
            
            

            if (bounce > 0)
                accumulate += -throughput * (valueUnhinderedAll * brdfValApprox * implicitMisWeight);
            else
                rRec.stochasticLtc += valueUnhinderedAll * brdfValApprox * implicitMisWeight;
*/

            // Russian Roulette path termination
            if (rRec.nextSample1D() <= terminationProbability)
                break;
            
            throughput /= (1 - terminationProbability);
            throughputLtc /= (1 - terminationProbability);

           
            bounce++;
       }

        return accumulate;
    }
#endif
#endif

#if 0 // Path tracing with brdf-sampling as explicit sample and one recursive sample using brdf-sampling.
#ifndef RATIO    
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
        Spectrum throughputLtc(1.0f);
        Spectrum accumulate(0.0f);
        Spectrum accumulateLtc(0.0f);
        Float terminationProbability = 0.05f;
        Intersection hitLoc;
        int bounce = 0;
        size_t nEmitterSamples = m_explicitSamples;

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
            const Float specularLum = specularReflectance.getLuminance();
            const Spectrum diffuseReflectance = brdf->getDiffuseReflectance(its);
            const Float diffuseLum = diffuseReflectance.getLuminance();
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

            DirectSamplingRecord dRec(its);
            
            if (m_collectAll || bounce == m_whichBounce)
                accumulateLtc += throughputLtc * Analytic::ltcIntegrate(scene, its.p, rotMat, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);

            // Brdf sampling
if (m_explicitConnect) {            
            for (size_t i=0; i < m_explicitSamples; ++i) {
                Float brdfPdf;
                BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
                Spectrum brdfValSampled = Analytic::sample(bRec, brdfPdf, rRec.nextSample2D(), m, mInv, mInvDet, amplitude, specularReflectance, specularLum, diffuseReflectance, diffuseLum); // brdfeval/pdf            
                if (brdfValSampled.isZero() || brdfPdf < Epsilon)
                    continue;
                
                Ray shadowRay(its.p, its.toWorld(bRec.wo) , ray.time);
                Spectrum valueUnhinderedAll(0.0f);
                intersectEmitter(scene, shadowRay, valueUnhinderedAll);

                if (valueUnhinderedAll.isZero())
                    continue;
                
                Spectrum valueDirect(0.0f);
                if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter()) {
                    const Emitter *tempEmitter = hitLoc.shape->getEmitter();
                    
                    if (tempEmitter != NULL &&
                        tempEmitter->isOnSurface() && // The next three condition essentially checks if the emitter is a mesh light.
                        tempEmitter->getShape() != NULL && 
                        typeid(*(tempEmitter->getShape())) == typeid(TriMesh))
                            // If the shadowRay hits a light source
                            valueDirect = hitLoc.Le(-shadowRay.d);    
                }
                if (m_collectAll || bounce == m_whichBounce) {
                    accumulate += throughput * valueDirect * brdfValSampled * m_fracExplicit;
                    accumulateLtc -= throughputLtc * valueUnhinderedAll * brdfValSampled * m_fracExplicit;
                }
            }
}
            // Brdf sampling    
            Float brdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum brdfValSampled = Analytic::sample(bRec, brdfPdf, rRec.nextSample2D(), m, mInv, mInvDet, amplitude, specularReflectance, specularLum, diffuseReflectance, diffuseLum); // brdfeval/pdf
            
            if (brdfValSampled.isZero() || brdfPdf < Epsilon)
                break;

            // Set recursive ray
            ray = Ray(its.p, its.toWorld(bRec.wo) , ray.time);
            if (!(rRec.type & RadianceQueryRecord::EIntersection))
                rRec.type ^= RadianceQueryRecord::EIntersection;
            
            if (!rRec.rayIntersect(ray))
                break;

            throughput *= brdfValSampled;
            throughputLtc *= brdfValSampled;

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

                if (m_collectAll || bounce == m_whichBounce) {
                    accumulate += throughput * valueDirect * m_fracImplicit;
                    accumulateLtc -= throughputLtc * valueUnhinderedAll * m_fracImplicit;
                }
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
            if (rRec.nextSample1D() <= terminationProbability)
                break;
            
            throughput /= (1 - terminationProbability);
            throughputLtc /= (1 - terminationProbability);

            bounce++;
       }
       
       return accumulate + accumulateLtc;
    }
#else
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
        Spectrum throughputLtc(1.0f);
        Spectrum accumulate(0.0f);
        Float terminationProbability = 0.05f;
        Intersection hitLoc;
        int bounce = 0;
        size_t nEmitterSamples = m_explicitSamples;

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
            const Float specularLum = specularReflectance.getLuminance();
            const Spectrum diffuseReflectance = brdf->getDiffuseReflectance(its);
            const Float diffuseLum = diffuseReflectance.getLuminance();
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

            DirectSamplingRecord dRec(its);
                      
            Spectrum LiLtc = throughputLtc * Analytic::ltcIntegrate(scene, its.p, rotMat, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);
            Spectrum LiWithVisibility(0.0f);
            Spectrum LiWithoutVisibility(0.0f);

            // Brdf sampling
if (m_explicitConnect) {            
            for (size_t i=0; i < m_explicitSamples; ++i) {
                Float brdfPdf;
                BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
                Spectrum brdfValSampled = Analytic::sample(bRec, brdfPdf, rRec.nextSample2D(), m, mInv, mInvDet, amplitude, specularReflectance, specularLum, diffuseReflectance, diffuseLum); // brdfeval/pdf            
                if (brdfValSampled.isZero() || brdfPdf < Epsilon)
                    continue;
                
                Ray shadowRay(its.p, its.toWorld(bRec.wo) , ray.time);
                Spectrum valueUnhinderedAll(0.0f);
                intersectEmitter(scene, shadowRay, valueUnhinderedAll);

                if (valueUnhinderedAll.isZero())
                    continue;
                
                Spectrum valueDirect(0.0f);
                if (scene->rayIntersect(shadowRay, hitLoc) && hitLoc.isEmitter()) {
                    const Emitter *tempEmitter = hitLoc.shape->getEmitter();
                    
                    if (tempEmitter != NULL &&
                        tempEmitter->isOnSurface() && // The next three condition essentially checks if the emitter is a mesh light.
                        tempEmitter->getShape() != NULL && 
                        typeid(*(tempEmitter->getShape())) == typeid(TriMesh))
                            // If the shadowRay hits a light source
                            valueDirect = hitLoc.Le(-shadowRay.d);    
                }
                LiWithVisibility += throughput * valueDirect * brdfValSampled * m_fracExplicit;
                LiWithoutVisibility += throughputLtc * valueUnhinderedAll * brdfValSampled * m_fracExplicit;
                
            }
}
            // Brdf sampling    
            Float brdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum brdfValSampled = Analytic::sample(bRec, brdfPdf, rRec.nextSample2D(), m, mInv, mInvDet, amplitude, specularReflectance, specularLum, diffuseReflectance, diffuseLum); // brdfeval/pdf
            
            if (brdfValSampled.isZero() || brdfPdf < Epsilon) {
                if (m_collectAll || bounce == m_whichBounce)
                    accumulate += getRatio(LiLtc, LiWithVisibility, LiWithoutVisibility);
                break;
            }

            // Set recursive ray
            ray = Ray(its.p, its.toWorld(bRec.wo) , ray.time);
            if (!(rRec.type & RadianceQueryRecord::EIntersection))
                rRec.type ^= RadianceQueryRecord::EIntersection;
            
            if (!rRec.rayIntersect(ray)) {
                if (m_collectAll || bounce == m_whichBounce)
                    accumulate += getRatio(LiLtc, LiWithVisibility, LiWithoutVisibility);
                break;
            }

            throughput *= brdfValSampled;
            throughputLtc *= brdfValSampled;

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

                LiWithVisibility += throughput * valueDirect * m_fracImplicit;
                LiWithoutVisibility += throughputLtc * valueUnhinderedAll * m_fracImplicit;
                
            }

            if (m_collectAll || bounce == m_whichBounce)
                accumulate += getRatio(LiLtc, LiWithVisibility, LiWithoutVisibility);

            if (its.isEmitter())
                break;

/*
            
            

            if (bounce > 0)
                accumulate += -throughput * (valueUnhinderedAll * brdfValApprox * implicitMisWeight);
            else
                rRec.stochasticLtc += valueUnhinderedAll * brdfValApprox * implicitMisWeight;
*/

            // Russian Roulette path termination
            if (rRec.nextSample1D() <= terminationProbability)
                break;
            
            throughput /= (1 - terminationProbability);
            throughputLtc /= (1 - terminationProbability);

            bounce++;
       }
       
       return accumulate;
    }
#endif
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
        pdfA *= pdfA; pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    inline Spectrum getRatio(const Spectrum &LiLtc, const Spectrum &LiWithVisibility, const Spectrum &LiWithoutVisibility) const {
        Spectrum Li = LiWithVisibility * LiLtc;

        if (LiWithoutVisibility[0] <= 0)
            Li[0] = Li[0] > 0.0f ? 100000000000000000000000000000.0f : 0.0f;
        else
            Li[0] /= LiWithoutVisibility[0];

        
        if (LiWithoutVisibility[1] <= 0)
            Li[1] = Li[1] > 0 ? 100000000000000000000000000000.0f : 0.0f;
        else
            Li[1] /= LiWithoutVisibility[1];


        if (LiWithoutVisibility[2] <= 0)
            Li[2] = Li[2] > 0 ? 100000000000000000000000000000.0f : 0.0f;
        else
            Li[2] /= LiWithoutVisibility[2];
        
        return Li;
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
    bool m_collectAll;
    size_t m_whichBounce;
    size_t m_explicitSamples, m_emitterSamples;
    Float m_fracExplicit, m_fracImplicit;
    Float m_weightExplicit, m_weightImplicit;
    size_t m_triEmitters;
    Float sceneSize;

    ref_vector<Emitter> m_triEmitterList; // A list of triEmitters
};

MTS_IMPLEMENT_CLASS_S(PathCv_v2, false, SamplingIntegrator);
MTS_EXPORT_PLUGIN(PathCv_v2, "Ltc control variate path tracing v2 integrator");
MTS_NAMESPACE_END