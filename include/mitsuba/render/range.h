#if !defined(__RANGE_WU_H)
#define __RANGE_WU_H

#include <mitsuba/core/sched.h>

MTS_NAMESPACE_BEGIN

/**
 * A work unit specifying e.g. a range of particles to be traced
 * by a worker.
 */
class MTS_EXPORT_RENDER RangeWorkUnit : public WorkUnit {
public:
	inline void set(const WorkUnit *wu) {
		const RangeWorkUnit *other = static_cast<const RangeWorkUnit *>(wu);
		m_rangeStart = other->m_rangeStart;
		m_rangeEnd = other->m_rangeEnd;
	}
	
	inline void load(Stream *stream) {
		m_rangeStart = (size_t) stream->readULong();
		m_rangeEnd = (size_t) stream->readULong();
	}

	inline void save(Stream *stream) const {
		stream->writeULong(m_rangeStart);
		stream->writeULong(m_rangeEnd);
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "RangeWorkUnit[rangeStart=" << m_rangeStart
			<< ", rangeEnd=" << m_rangeEnd << "]" << endl;
		return oss.str();
	}

	inline void setRange(size_t start, size_t end) {
		m_rangeStart = start;
		m_rangeEnd = end;
	}

	inline size_t getRangeStart() const {
		return m_rangeStart;
	}

	inline size_t getRangeEnd() const {
		return m_rangeEnd;
	}

	inline size_t getSize() const {
		return m_rangeEnd - m_rangeStart + 1;
	}

	MTS_DECLARE_CLASS()
private:
	size_t m_rangeStart, m_rangeEnd;
};

MTS_NAMESPACE_END

#endif /* __RANGE_WU */
