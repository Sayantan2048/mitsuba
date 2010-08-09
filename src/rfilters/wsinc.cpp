#include <mitsuba/render/rfilter.h>

MTS_NAMESPACE_BEGIN

/**
 * Windowed version of the ideal low-pass filter (with a Lanczos
 * envelope) -- may produce ringing artifacts near sharp edges. 
 */
class WindowedSincFilter : public ReconstructionFilter {
public:
	WindowedSincFilter(const Properties &props) 
		: ReconstructionFilter(props) {
		/* Half filter size in pixels */
		Float halfSize = props.getFloat("halfSize", 3.0f);
		/* Number of cycles, after which the sinc function should be truncated */
		m_cycles = props.getFloat("cycles", 3.0f);
		m_size = Vector2(halfSize, halfSize);
	}

	WindowedSincFilter(Stream *stream, InstanceManager *manager) 
		: ReconstructionFilter(stream, manager) {
		Float halfSize = stream->readFloat();
		m_size = Vector2(halfSize, halfSize);
		m_cycles = stream->readFloat();
	}
 
	void serialize(Stream *stream, InstanceManager *manager) const {
		ReconstructionFilter::serialize(stream, manager);
		stream->writeFloat(m_size.x);
		stream->writeFloat(m_cycles);
	}

	Float evaluate(Float x, Float y) const {
		return lanczosSinc(x / m_size.x, m_cycles)
			 * lanczosSinc(y / m_size.y, m_cycles);
	}

	MTS_DECLARE_CLASS()
protected:
	Float m_cycles;
};

MTS_IMPLEMENT_CLASS_S(WindowedSincFilter, false, ReconstructionFilter);
MTS_EXPORT_PLUGIN(WindowedSincFilter, "Windowed sinc filter");
MTS_NAMESPACE_END
