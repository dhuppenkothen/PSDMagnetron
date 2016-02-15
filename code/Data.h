#ifndef _Data_
#define _Data_

#include <vector>

class Data
{
	private:
		std::vector<double> f, y, m;
		std::vector<double> f_left, f_right;

		// Some useful summaries
		double f_min, f_max, f_range, df;
		void compute_summaries();

	public:
		Data();
		void load(const char* filename);

		// Getters
		const std::vector<double>& get_f() const { return f; }
		const std::vector<double>& get_f_left() const { return f_left; }
		const std::vector<double>& get_f_right() const { return f_right; }
		const std::vector<double>& get_y() const { return y; }
		const std::vector<double>& get_m() const { return m; }
		double get_f_min() const { return f_min; }
		double get_f_max() const { return f_max; }
		double get_f_range() const { return f_range; }
		double get_df() const { return df; }

	// Singleton
	private:
		static Data instance;
	public:
		static Data& get_instance() { return instance; }
};

#endif

