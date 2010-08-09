#include <mitsuba/render/bsdf.h>

MTS_NAMESPACE_BEGIN

/**
 * Perfect specular BRDF (i.e. an ideal mirror)
 */
class Mirror : public BSDF {
public:
	Mirror(const Properties &props) 
		: BSDF(props) {
		m_reflectance = props.getSpectrum("reflectance", Spectrum(0.8f));
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDeltaReflection;
		m_combinedType = m_type[0];
		m_usesRayDifferentials = false;
	}

	Mirror(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_reflectance = Spectrum(stream);
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDeltaReflection;
		m_combinedType = m_type[0];
		m_usesRayDifferentials = false;
	}

	virtual ~Mirror() {
		delete[] m_type;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		m_reflectance.serialize(stream);
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return Spectrum(0.0f);
	}
	
	Spectrum f(const BSDFQueryRecord &bRec) const {
		return Spectrum(0.0f);
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		return 0.0f;
	}

	inline void reflect(const Vector &wi, Vector &wo) const {
		wo = Vector(-wi.x, -wi.y, wi.z);
	}

	inline Spectrum sample(BSDFQueryRecord &bRec) const {
		Float pdf;
		return Mirror::sample(bRec, pdf);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf) const {
		if (!(bRec.typeMask & m_combinedType))
			return Spectrum(0.0f);
		reflect(bRec.wi, bRec.wo);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDeltaReflection;
		pdf = std::abs(Frame::cosTheta(bRec.wo));
		return m_reflectance;
	}

	Float pdfDelta(const BSDFQueryRecord &bRec) const {
		return std::abs(Frame::cosTheta(bRec.wo));
	}
	
	Spectrum fDelta(const BSDFQueryRecord &bRec) const {
		return m_reflectance;
	}


	MTS_DECLARE_CLASS()
protected:
	Spectrum m_reflectance;
};


MTS_IMPLEMENT_CLASS_S(Mirror, false, BSDF)
MTS_EXPORT_PLUGIN(Mirror, "Mirror BRDF");
MTS_NAMESPACE_END
