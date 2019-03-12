#include <mitsuba/render/scene.h>
#include "analytic.h"
MTS_NAMESPACE_BEGIN

class DirectRatioIntegrator : public SamplingIntegrator {
    public:
    DirectRatioIntegrator(const Properties &props) : SamplingIntegrator(props) {
        /* Number of samples to take using the emitter sampling technique */
        m_emitterSamples = props.getSize("emitterSamples", 1);
        m_hideEmitters = props.getBoolean("hideEmitters", false);

        Assert(m_emitterSamples > 0);
    }

    /// Unserialize from a binary data stream
    DirectRatioIntegrator(Stream *stream, InstanceManager *manager)
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

    void configure() {
        SamplingIntegrator::configure();
    }

    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);
        if (m_emitterSamples > 1)
            sampler->request2DArray(m_emitterSamples);
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
        Log(EInfo, "Running LTC ratio integrator.");
        return true;
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
        // Probably one might want to change its.shFrame to rotMat. Be careful, one must change all the local variables to the new frame before setting the its.shFrame.
        its.wi = rotMat * its.toWorld(its.wi);

        Point2 *sampleArray;
        size_t numDirectSamples = m_emitterSamples;

        if (numDirectSamples > 1) {
            sampleArray = rRec.sampler->next2DArray(numDirectSamples);
        } else {
            sample = rRec.nextSample2D(); sampleArray = &sample;
        }

        DirectSamplingRecord dRec(its);
        Intersection hitLoc;

        Spectrum LiWithVisibility(0.0f);
        Spectrum LiWithoutVisibility(0.0f);

        if (bsdf->getType() & BSDF::ESmooth) {
            /* Only use direct illumination sampling when the surface's
               BSDF has smooth (i.e. non-Dirac delta) component */
            for (size_t i=0; i<numDirectSamples; ++i) {
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
                        BSDFSamplingRecord bRec(its, rotMat * dRec.d);
                        Spectrum bsdfVal(0.0f);
                        if (notInShadow > 0.0f) {
                            bsdfVal = bsdf->eval(bRec);
                        }
                        const Spectrum bsdfValApprox = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance);
                        LiWithVisibility += valueUnhindered * bsdfValApprox * notInShadow / (Float)m_emitterSamples;
                        LiWithoutVisibility += valueUnhindered * bsdfValApprox / (Float)m_emitterSamples;
                    }
                }
            }
        }

        
        Spectrum LiWithoutVisibilityInv = Spectrum(1.0f) / LiWithoutVisibility;
        if (!LiWithoutVisibilityInv.isValid()) {
            LiWithoutVisibilityInv[0] = std::isnan(LiWithoutVisibilityInv[0]) || std::isinf(LiWithoutVisibilityInv[0]) ? 0 : LiWithoutVisibilityInv[0];
            LiWithoutVisibilityInv[1] = std::isnan(LiWithoutVisibilityInv[1]) || std::isinf(LiWithoutVisibilityInv[1]) ? 0 : LiWithoutVisibilityInv[1];
            LiWithoutVisibilityInv[2] = std::isnan(LiWithoutVisibilityInv[2]) || std::isinf(LiWithoutVisibilityInv[2]) ? 0 : LiWithoutVisibilityInv[2];
        }

        Spectrum Li = LiWithoutVisibilityInv * LiWithVisibility * m_subIntegrator->Li(ray, rRec);

        //LiWithVisibility[0];// m_subIntegrator->Li(ray, rRec);
        //if (diff > 0)
            //Li.fromLinearRGB(0, 1, 0);
        //Li.clampNegative();

        return Li;
    }

    
    std::string toString() const {
        std::ostringstream oss;
        oss << "DirectRatioIntegrator[" << endl
            << "  emitterSamples = " << m_emitterSamples << "," << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    size_t m_emitterSamples;
    ref<SamplingIntegrator> m_subIntegrator;
    bool m_hideEmitters;
};

MTS_IMPLEMENT_CLASS_S(DirectRatioIntegrator, false, SamplingIntegrator);
MTS_EXPORT_PLUGIN(DirectRatioIntegrator, "Ratio estimator with LTC for direct illumination");
MTS_NAMESPACE_END