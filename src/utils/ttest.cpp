#include <mitsuba/mitsuba.h>
#include <boost/math/distributions/students_t.hpp>
#include <fstream>

using namespace mitsuba;

struct Sample {
	Float value;
	Float variance;
	int nSamples;
};

static std::vector<Sample> parseMFile(std::ifstream &is) {
	std::string line;
	std::vector<Sample> result;
	while (!is.eof() && !is.fail()) {
		std::getline(is, line);
		std::vector<std::string> tokens = tokenize(line, " \t;,[]");
		if ((tokens.size() % 3) != 0) {
			cerr << "The number of entries is not divisible by 3 - does this file "
				 << "contain variance information?" << endl;
			exit(-1);
		}

		for (size_t i=0; i<tokens.size(); ) {
			Sample sample;
			char *end_ptr = NULL;
			sample.value = (Float) strtod(tokens[i++].c_str(), &end_ptr);
			if (*end_ptr != '\0') { cerr << "Could not parse file!" << endl; exit(-1); }
			sample.variance = (Float) strtod(tokens[i++].c_str(), &end_ptr);
			if (*end_ptr != '\0') { cerr << "Could not parse file!" << endl; exit(-1); }
			sample.nSamples = (int) strtol(tokens[i++].c_str(), &end_ptr, 10);
			if (*end_ptr != '\0') { cerr << "Could not parse file!" << endl; exit(-1); }
			result.push_back(sample);
		}
	}
	return result;
}

	
int main(int argc, char **argv) {
	using namespace boost::math;

	if (argc < 3) {
		cout << "This utility runs a t-test on the equality of means for two M-files" << endl;
		cout << "developed using the 'mfilm' film implementation. Both files must" << endl;
		cout << "contain variance information." << endl << endl;
		cout << "Syntax: ttest <file1.m> <file2.m>" << endl;
		return 0;
	}
	std::ifstream is1(argv[1]);
	std::ifstream is2(argv[2]);
	if (is1.fail()) {
		cerr << "Could not open " << argv[1] << "!" << endl;
		return -1;
	}
	if (is2.fail()) {
		cerr << "Could not open " << argv[2] << "!" << endl;
		return -1;
	}
	std::vector<Sample> input1 = parseMFile(is1);
	std::vector<Sample> input2 = parseMFile(is2);

	if (input1.size() != input2.size()) {
		cerr << "Error: the inputs are of different size!" << endl;
		return -1;
	}

	for (size_t i=0; i<input1.size(); ++i) {
		/* Studentized statistic for testing the equality of two means */
		Float Tstat = (input1[i].value - input2[i].value) 
			/ std::sqrt(input1[i].variance / input1[i].nSamples 
				+ input2[i].variance / input2[i].nSamples);

		/* Satterthwaites' approximation */
		Float num = std::pow(input1[i].variance/input1[i].nSamples
			+input2[i].variance/input2[i].nSamples,2);
		Float denom = (Float) 1 / (input1[i].nSamples-1) * 
			std::pow(input1[i].variance/input1[i].nSamples, 2)
			+ (Float) 1 / (input2[i].nSamples-1) * 
			std::pow(input2[i].variance/input2[i].nSamples, 2);
		Float df = num/denom;
		students_t dist(df);
		Float pval = (Float) (2*cdf(complement(dist, std::abs(Tstat))));
		cout << "  Entry " << i << ": t-statistic=" << Tstat << ", df=" << df << ", p-value=" << pval<< endl;
		if (pval < 0.01)
			cout << "     => The t-test REJECTS!" << endl;
		else
			cout << "     => The t-test accepts." << endl;
	}
	is1.close();
	is2.close();
	return 0;
}

