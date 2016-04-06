#ifndef _MyConditionalPrior_
#define _MyConditionalPrior_

#include "DNest4/code/RJObject/ConditionalPriors/ConditionalPrior.h"

class MyNarrowConditionalPrior:public DNest4::ConditionalPrior
{
	private:
		// Limits
		double x_min, x_max;

		// Mean of amplitudes and widths
		double mu_loga;

		// Uniform for log-skews
		double mu_logq, sigma_logq; // Midpoint and half-width

		double perturb_hyperparameters(DNest4::RNG& rng);

	public:
		MyNarrowConditionalPrior(double x_min, double x_max);

		void from_prior(DNest4::RNG& rng);

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;

		static const int weight_parameter = 1;
};

class MyWideConditionalPrior:public DNest4::ConditionalPrior
{
        private:
                // Limits
                double x_min, x_max;

                // Mean of amplitudes and widths
                double mu_amp;

                // Uniform for log-skews
                double mu_logq, sigma_logq; // Midpoint and half-width

                double perturb_hyperparameters(DNest4::RNG& rng);

        public:
                MyWideConditionalPrior(double x_min, double x_max);

                void from_prior(DNest4::RNG& rng);

                double log_pdf(const std::vector<double>& vec) const;
                void from_uniform(std::vector<double>& vec) const;
                void to_uniform(std::vector<double>& vec) const;

                void print(std::ostream& out) const;

                static const int weight_parameter = 1;
};

#endif

