#if !defined(__SAMPLER_H)
#define __SAMPLER_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/random.h>

MTS_NAMESPACE_BEGIN

/**
 * Abstract sample generator. Generates a high-dimensional sample
 * and provides 1D and 2D pieces of it to the integrator as requested. 
 * An implementation will ideally provide a sequence that is well 
 * laid-out in sample space using either stratification or QMC techniques.
 */
class MTS_EXPORT_RENDER Sampler : public ConfigurableObject {
public:
	/**
	 * Create a clone of this sampler. The clone is allowed to be different
	 * to some extent, e.g. a pseudorandom generator should be based on a
	 * different random seed compared to the original. All other parameters,
	 * are copied exactly.
	 */
	virtual ref<Sampler> clone() = 0;

	/**
	 * Generate new samples - called initially and every time the generated 
	 * samples have been exhausted. When used in conjunction with a 
	 * SampleIntegrator, this will be called before starting to render 
	 * each pixel.
	 */
	virtual void generate();

	/// Advance to the next sample
	virtual void advance();

	/// Manually set the current sample index
	virtual void setSampleIndex(uint64_t sampleIndex);

	/// Retrieve the next component value from the current sample
	virtual Float next1D() = 0;

	/// Retrieve the next two component values from the current sample
	virtual Point2 next2D() = 0;

	/**
	 * Retrieve the next 2D array of values from the current sample.
	 * Note that this is different from just calling <tt>next2D()</tt>
	 * repeatedly - this function will generally return a set of 2D vectors,
	 * which are not only well-laid out over all samples at the current pixel,
	 * but also with respect to each other. Note that this 2D array has to be
	 * requested initially using <tt>request2DArray</tt> and later, they have 
	 * to be retrieved in the same same order and size configuration as the 
	 * requests. An exception is thrown when a mismatch is detected.
	 */
	Point2 *next2DArray(unsigned int size);

	/// Same as above, but 1D
	Float *next1DArray(unsigned int size);

	/**
	 * Request that a 2D array will be made available for 
	 * later consumption by next2DArray(). This must be called
	 * before generate(). See 'next2DArray' for a description
	 * of this feature.
	 */
	virtual void request2DArray(unsigned int size);

	/// Same as above, but 1D
	virtual void request1DArray(unsigned int size);

	/**
	 * Return an uniformly distributed number on [0, 1).
	 * Throws an error when the underlying implementation is 
	 * fully deterministic (e.g. QMC).
	 */
	virtual Float independent1D() = 0;

	/**
	 * Return an uniformly distributed 2D vector on [0, 1)x[0, 1).
	 * Throws an error when the underlying implementation is 
	 * fully deterministic (e.g. QMC).
	 */
	virtual Point2 independent2D() = 0;

	/// Return total number of samples
	inline uint64_t getSampleCount() const { return m_sampleCount; }
	
	/// Return the current sample index
	inline uint64_t getSampleIndex() const { return m_sampleIndex; }

	/// Serialize this sampler to disk
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Return the properties of this sampler
	inline const Properties &getProperties() const { return m_properties; }

	/// Return a string description
	virtual std::string toString() const = 0;

	MTS_DECLARE_CLASS()
protected:
	/// Construct a new sampler
	Sampler(const Properties &props);

	/// Unserialize a sampler
	Sampler(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Sampler();
protected:
	uint64_t m_sampleCount;
	uint64_t m_sampleIndex;
	std::vector<unsigned int> m_req1D, m_req2D;
	std::vector<Float *> m_sampleArrays1D;
	std::vector<Point2 *> m_sampleArrays2D;
	int m_sampleDepth1DArray, m_sampleDepth2DArray;
	Properties m_properties;
};

MTS_NAMESPACE_END

#endif /* __SAMPLER_H */
