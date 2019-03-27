#include <mitsuba/render/scene.h>
#include "analytic.h"
MTS_NAMESPACE_BEGIN

class PathCv : public SamplingIntegrator {
    public:
    PathCv(const Properties &props) : SamplingIntegrator(props) {
        m_explicitConnect = props.getBoolean("explicit", true); // Request a non-recursive ray
        m_explicitType = props.getSize("explicitType", 0); // The type of explicit connection e.g. emitter, brdf, approx brdf, cosine, uniform.
        m_explicitSamples = props.getSize("explicitSamples", 1); // Number of explicit samples.

        m_explicitConnect = false;
        m_explicitSamples = 0;
        Assert(!m_explicitConnect || (m_explicitConnect && m_explicitSamples > 0));
        Assert(m_explicitType <= 4);
    }

    /// Unserialize from a binary data stream
    PathCv(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        m_explicitConnect = stream->readBool();
        m_explicitType = stream->readSize();
        m_explicitSamples = stream->readSize();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        stream->writeBool(m_explicitConnect);
        stream->writeSize(m_explicitType);
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

        size_t sum = m_explicitSamples + 1;
        m_weightExplicit = 1 / (Float) m_explicitSamples;
        m_weightImplicit = 1;
        m_fracExplicit = m_explicitSamples / (Float) sum;
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

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        Spectrum throughput(1.0f);
        Spectrum accumulate(0.0f);
        Float terminationProbability = 0.05f;
        Point2 sample;
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
            
            Float brdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            const Spectrum brdfVal = brdf->sample(bRec, brdfPdf, sample); // brdfeval/pdf
                      
            if (brdfPdf <= Epsilon)
                break;
            
            accumulate += throughput * m_subIntegrator->Li(ray, rRec);
            
            const Spectrum brdfValApprox = Analytic::approxBrdfEval(bRec, mInv, mInvDet, amplitude, specularReflectance, diffuseReflectance) / brdfPdf;  
           
            Float explicitPdf = 0;
           
            Float implicitMisWeight = miWeight(brdfPdf * m_fracImplicit,
                    explicitPdf * m_fracExplicit) * m_weightImplicit;

            // Set recursive ray
            ray = Ray(its.p, its.toWorld(bRec.wo) , ray.time);
            if (!(rRec.type & RadianceQueryRecord::EIntersection))
                rRec.type ^= RadianceQueryRecord::EIntersection;
            
            
                Spectrum valueUnhinderedAll(0.0f);
                // Note that this intersection with the light source is with shadowRay(same as recursive ray), so do not put the intersection before setting recursive ray.
                intersectEmitter(scene, ray, valueUnhinderedAll);
                accumulate += -throughput * (valueUnhinderedAll * brdfValApprox * implicitMisWeight);
                if (!brdfValApprox.isValid())
                    Log(EInfo, "Invalid %f", brdfPdf);
            

            throughput *= brdfVal * implicitMisWeight;
            bounce++;
       }
    
       //accumulate.clampNegative();
       return accumulate;
    }

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
            << "  explicitType = " << m_explicitType << "," << endl
            << "  explicitSamples = " << m_explicitSamples << "," << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    bool m_explicitConnect;
    size_t m_explicitType;
    size_t m_explicitSamples;
    Float m_fracExplicit, m_fracImplicit;
    Float m_weightExplicit, m_weightImplicit;
    ref<SamplingIntegrator> m_subIntegrator;
    Float sceneSize;
};

MTS_IMPLEMENT_CLASS_S(PathCv, false, SamplingIntegrator);
MTS_EXPORT_PLUGIN(PathCv, "Ltc control variate path tracing integrator");
MTS_NAMESPACE_END