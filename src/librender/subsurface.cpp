#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

Subsurface::Subsurface(const Properties &props)
 : NetworkedObject(props) {
	/* Skim milk data from "A Practical Model for Subsurface scattering" (Jensen et al.) */
	Spectrum defaultSigmaS, defaultSigmaA;

	defaultSigmaA.fromLinearRGB(0.0014f, 0.0025f, 0.0142f);
	defaultSigmaS.fromLinearRGB(0.7f, 1.22f, 1.9f);

	m_sizeMultiplier = props.getFloat("sizeMultiplier", 1);
	/* Scattering coefficient */
	m_sigmaS = props.getSpectrum("sigmaS", defaultSigmaS);
	/* Absorption coefficient */
	m_sigmaA = props.getSpectrum("sigmaA", defaultSigmaA);
	/* Refractive index of the object */
	m_eta = props.getFloat("eta", 1.5f);
		
	m_sigmaS *= m_sizeMultiplier;
	m_sigmaA *= m_sizeMultiplier;
	m_sigmaT = m_sigmaS + m_sigmaA;
}

Subsurface::Subsurface(Stream *stream, InstanceManager *manager) :
	NetworkedObject(stream, manager) {
	m_sigmaS = Spectrum(stream);
	m_sigmaA = Spectrum(stream);
	m_eta = stream->readFloat();
	m_sizeMultiplier = stream->readFloat();
	unsigned int shapeCount = stream->readUInt();
	for (unsigned int i=0; i<shapeCount; ++i)
		m_shapes.push_back(static_cast<Shape *>(manager->getInstance(stream)));
	m_sigmaT = m_sigmaS + m_sigmaA;
}
	
void Subsurface::setParent(ConfigurableObject *parent) {
	if (parent->getClass()->derivesFrom(Shape::m_theClass)) {
		Shape *shape = static_cast<Shape *>(parent);
		m_shapes.push_back(shape);
		m_configured = false;
	} else {
		Log(EError, "IsotropicDipole: Invalid child node!");
	}
}

void Subsurface::serialize(Stream *stream, InstanceManager *manager) const {
	NetworkedObject::serialize(stream, manager);

	m_sigmaS.serialize(stream);
	m_sigmaA.serialize(stream);
	stream->writeFloat(m_eta);
	stream->writeFloat(m_sizeMultiplier);
	stream->writeUInt(m_shapes.size());
	for (unsigned int i=0; i<m_shapes.size(); ++i)
		manager->serialize(stream, m_shapes[i]);
}

Subsurface::~Subsurface() {
}

MTS_IMPLEMENT_CLASS(Subsurface, true, NetworkedObject)
MTS_NAMESPACE_END
