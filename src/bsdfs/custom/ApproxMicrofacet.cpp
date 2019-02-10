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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

// Contains the LTC matrices for given(discreet) alpha and theta.
#include "ggx.dat"
#include "beckmann.dat"

class ApproxMicrofacet : public BSDF {
public:
    enum EType {
        /// Beckmann distribution derived from Gaussian random surfaces
        EBeckmann         = 0,

        /// GGX: Long-tailed distribution for very rough surfaces (aka. Trowbridge-Reitz distr.)
        EGGX              = 1,
    };

    ApproxMicrofacet(const Properties &props) : BSDF(props) {
        m_specularReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("specularReflectance", Spectrum(1.0f)));
        m_diffuseReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));
        m_alpha = new ConstantFloatTexture(props.getFloat("alpha", 0.2f));

        if (props.hasProperty("distribution")) {
            std::string distr = boost::to_lower_copy(props.getString("distribution"));
            if (distr == "beckmann")
                m_type = EBeckmann;
            else if (distr == "ggx")
                m_type = EGGX;
            else
                SLog(EError, "Specified an invalid distribution \"%s\", must be "
                    "\"beckmann\" or \"ggx\" \"as\"!", distr.c_str());
        }
        else
            m_type = EGGX;

    }

    ApproxMicrofacet(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_alpha = static_cast<Texture *>(manager->getInstance(stream));
        m_type = (EType) stream->readUInt();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        manager->serialize(stream, m_specularReflectance.get());
        manager->serialize(stream, m_diffuseReflectance.get());
        manager->serialize(stream, m_alpha.get());
        stream->writeUInt(m_type);
    }

    void configure() {
        /* Verify the input parameter and fix them if necessary */
        /* Verify the input parameters and fix them if necessary */
        m_specularReflectance = ensureEnergyConservation(m_specularReflectance, "specularReflectance", 1.0f);
        m_diffuseReflectance = ensureEnergyConservation(m_diffuseReflectance, "diffuseReflectance", 1.0f);

        bool isAlphaConstant = m_alpha->isConstant();

        m_components.clear();

        m_components.push_back(EGlossyReflection | EFrontSide
            | ((isAlphaConstant && m_specularReflectance->isConstant())
                ? 0 : ESpatiallyVarying));
        m_components.push_back(EDiffuseReflection | EFrontSide
            | ((isAlphaConstant && m_diffuseReflectance->isConstant())
                ? 0 : ESpatiallyVarying));

        ltcDataMatSize = m_type == EGGX ? ggxSize : beckmannSize;
        ltcDataMat = m_type == EGGX ? ggxDataMat : beckmannDataMat;
        ltcDataMatInv = m_type == EGGX ? ggxDataMatInv : beckmannDataMatInv;
        ltcDataMatAmp = m_type == EGGX ? ggxDataMatAmp : beckmannDataMatAmp;

        BSDF::configure();
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        return m_diffuseReflectance->eval(its);
    }

    Spectrum getSpecularReflectance(const Intersection &its) const {
        return m_specularReflectance->eval(its);
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        return Spectrum(0.0);
    }
    
    void transform(const Intersection &its, const Float &theta, const Matrix3x3 &rotate, Matrix3x3 &mInv, Float &amplitude) const {
        Float alpha = m_alpha->eval(its).average();
        Float blendeAlpha, blendeTheta;
        Vector4i lerp;
        getLerp(theta, alpha, lerp, blendeTheta, blendeAlpha);
        //Log(EInfo, "%d %d %d %d %f %f", lerp.x, lerp.y, lerp.z, lerp.w, blendeTheta, blendeAlpha);
        
        mInv = getMatrixInv(lerp, blendeTheta , blendeAlpha);
        amplitude = getAmp(lerp, blendeTheta , blendeAlpha);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        /*if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        return warp::squareToCosineHemispherePdf(bRec.wo);*/
        return 0.0f;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        /*if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        bRec.wo = warp::squareToCosineHemisphere(sample);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;
        return m_reflectance->eval(bRec.its);*/
        return Spectrum(0.0f);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        /*if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        bRec.wo = warp::squareToCosineHemisphere(sample);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;
        pdf = warp::squareToCosineHemispherePdf(bRec.wo);
        return m_reflectance->eval(bRec.its);*/
        return Spectrum(0.0f);
    }

    // This is required for putting another brdf inside this one.
    // May come handy later on when I test with LTC with CV
    void addChild(const std::string &name, ConfigurableObject *child) {
        /*
        Log(EInfo, "%f", m_alpha);
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
            if (name == "specularReflectance")
                m_specularReflectance = static_cast<Texture *>(child);
            else if (name == "diffuseReflectance")
                m_diffuseReflectance = static_cast<Texture *>(child);
            else
                BSDF::addChild(name, child);
        } else if (name == "alpha")
            Log(EInfo, "I don't know what to do yet");
            //m_alpha = static_cast<Float>(child);
        else
            BSDF::addChild(name, child); */
    }

    Float getRoughness(const Intersection &its, int component) const {
         Assert(component == 0 || component == 1);

        if (component == 0)
            return m_alpha->eval(its).average();
        else
            return std::numeric_limits<Float>::infinity();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "ApproxMicrofacet[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  specularReflectance = " << indent(m_specularReflectance->toString()) << endl
            << "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << endl
            << "  alpha = " << indent(m_alpha->toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ref<Texture> m_diffuseReflectance;
    ref<Texture> m_specularReflectance;
    ref<Texture> m_alpha;
    EType m_type;

    const float (*ltcDataMat)[9];
    const float (*ltcDataMatInv)[9];
    const float *ltcDataMatAmp;
    int ltcDataMatSize;
    
    void getLerp(const Float theta, const Float alpha, Vector4i &lerp, Float &blendeTheta, Float &blendeAlpha) const {
        Float theta_ = theta / (0.5f*M_PI) * ltcDataMatSize;
        Float alpha_ = sqrtf(alpha) * ltcDataMatSize;
        int thetaLo = MAX(0, MIN(ltcDataMatSize-1, (int)floor(theta_)));
        int thetaHi = MAX(0, MIN(ltcDataMatSize-1, (int)ceil(theta_)));;
        int alphaLo = MAX(0, MIN(ltcDataMatSize-1, (int)floor(alpha_)));
        int alphaHi = MAX(0, MIN(ltcDataMatSize-1, (int)ceil(alpha_)));

        if (theta_ < thetaLo)
            theta_ = (Float)thetaLo;
        else if (theta_ > thetaHi)
            theta_ = (Float)thetaHi;

        if (alpha_ < alphaLo)
            alpha_ = (Float)alphaLo;
        else if (alpha_ > alphaHi)
            alpha_ = (Float)alphaHi;

        /*
        Float cosThetaLo = cos(thetaLo * 0.5f * M_PI / (ltcDataMatSize - 1));
        Float cosThetaHi = cos(thetaHi * 0.5f * M_PI / (ltcDataMatSize - 1));
        Float cosTheta = cos(theta);
        */
                    
        blendeTheta = (thetaHi == thetaLo) ? (Float)1.0 : (theta_ - thetaLo) / (thetaHi - thetaLo);
        //blendeTheta = (thetaHi == thetaLo) ? (Float)1.0 : (cosThetaLo - cosTheta) / (cosThetaLo - cosThetaHi);
        blendeAlpha = (alphaLo == alphaHi) ? (Float)1.0 : (alpha_ - alphaLo)/ (Float)(alphaHi - alphaLo);
        
        lerp.x = alphaLo;
        lerp.y = thetaLo;
        lerp.z = alphaHi;
        lerp.w = thetaHi;
    }

    Matrix3x3 getMatrixInv(const Vector4i &lerp, const Float blendeTheta , const Float blendeAlpha) const {
        int alphaLo = lerp.x;
        int thetaLo = lerp.y;
        int alphaHi = lerp.z;
        int thetaHi = lerp.w;
        Matrix3x3 ret;
        
        // Matrix3x3 constructor assumes the input array in row major order.
        // The data however is column major(standard in GLSL). Hence a transpose is required. 
        ((1.0f - blendeTheta) * ((1.0f - blendeAlpha) * Matrix3x3(ltcDataMatInv[alphaLo + thetaLo * ltcDataMatSize]) + 
                blendeAlpha * Matrix3x3(ltcDataMatInv[alphaHi + thetaLo * ltcDataMatSize])) +
                blendeTheta * ((1.0f - blendeAlpha) * Matrix3x3(ltcDataMatInv[alphaLo + thetaHi * ltcDataMatSize]) + 
                blendeAlpha * Matrix3x3(ltcDataMatInv[alphaHi + thetaHi * ltcDataMatSize]))).transpose(ret);
        
        return ret;
    }

    Matrix3x3 getMatrix(const Vector4i &lerp, const Float blendeTheta , const Float blendeAlpha) const {
        int alphaLo = lerp.x;
        int thetaLo = lerp.y;
        int alphaHi = lerp.z;
        int thetaHi = lerp.w;
        
        Matrix3x3 ret;
        ((1.0f - blendeTheta) * ((1.0f - blendeAlpha) * Matrix3x3(ltcDataMat[alphaLo + thetaLo * ltcDataMatSize]) + 
                blendeAlpha * Matrix3x3(ltcDataMat[alphaHi + thetaLo * ltcDataMatSize])) +
                blendeTheta * ((1.0f - blendeAlpha) * Matrix3x3(ltcDataMat[alphaLo + thetaHi * ltcDataMatSize]) + 
                blendeAlpha * Matrix3x3(ltcDataMat[alphaHi + thetaHi * ltcDataMatSize]))).transpose(ret);
        
        return ret;
    }

     Float getAmp(const Vector4i &lerp, const Float blendeTheta , const Float blendeAlpha) const {
        int alphaLo = lerp.x;
        int thetaLo = lerp.y;
        int alphaHi = lerp.z;
        int thetaHi = lerp.w;
        
        return (1.0f - blendeTheta) * ((1.0f - blendeAlpha) * ltcDataMatAmp[alphaLo + thetaLo * ltcDataMatSize] + 
                blendeAlpha * ltcDataMatAmp[alphaHi + thetaLo * ltcDataMatSize]) +
                blendeTheta * ((1.0f - blendeAlpha) * ltcDataMatAmp[alphaLo + thetaHi * ltcDataMatSize] + 
                blendeAlpha * ltcDataMatAmp[alphaHi + thetaHi * ltcDataMatSize]) ;
    }
};

MTS_IMPLEMENT_CLASS_S(ApproxMicrofacet, false, BSDF)
MTS_EXPORT_PLUGIN(ApproxMicrofacet, "Approx Microfacet BRDF")
MTS_NAMESPACE_END
