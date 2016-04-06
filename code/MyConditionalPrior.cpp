#include "MyConditionalPrior.h"
#include "DNest4/code/Utils.h"
#include "Data.h"
#include <cmath>

using namespace DNest4;

//MyWideConditionalPrior::MyWideConditionalPrior(double x_min, double x_max)
//:x_min(x_min)
//,x_max(x_max)
//{

//}

void MyWideConditionalPrior::from_prior(RNG& rng)
{
	// Cauchy prior on the mean of the exponential amplitude prior
	mu_amp = tan(M_PI*(0.97*rng.rand() - 0.485));
	mu_amp = exp(mu_amp);

	// Uniform prior on the mean of the Laplacian logq prior:
	mu_logq = (log(2.)-log(1E-5))*rng.rand() + log(1E-5);
	// Uniform prior on the width of the Laplacian logq prior:
	sigma_logq = 2.*rng.rand(); 

//	mu_widths = exp(log(1E-3*(x_max - x_min)) + log(1E3)*rng.rand());

//	a = -10. + 20.*rng.rand();
//	b = 2.*rng.rand();
}

double MyWideConditionalPrior::perturb_hyperparameters(RNG& rng)
{
	double logH = 0.;

	int which = rng.rand_int(3);

	if(which == 0)
	{
		mu_amp = log(mu_amp);
		mu_amp = (atan(mu_amp)/M_PI + 0.485)/0.97;
		mu_amp += rng.randh();
 		wrap(mu_amp, 0., 1.);
		mu_amp = tan(M_PI*(0.97*mu_amp - 0.485));
		mu_amp = exp(mu_amp);
	}
	if(which == 1)
	{
		// check this!
		mu_logq += rng.randh()*(log(2.)-log(1E-5)); //log(100)*pow(10., log(2.) - log(100.)*rng.rand())*rng.randn();
		wrap(mu_logq, log(1E-5), log(2.));
		//mu_logq = mod(mu_logq - log(.2), log(100.)) + log(2.);
	}
	if(which == 2)
	{
		sigma_logq += 2.*rng.randh();
		wrap(sigma_logq, 0., 2.);
	}
	return logH;
}

double MyWideConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
	if(vec[0] < x_min || vec[0] > x_max || vec[1] < 0.0) 
		return -1E300;

	return -log(mu_amp) - vec[1]/mu_amp - log(2.*sigma_logq) - 
			std::abs(vec[2]-mu_logq)/sigma_logq;
	//return -log(mu) - vec[1]/mu - log(mu_widths)
	//		- (vec[2] - min_width)/mu_widths - log(2.*b*vec[3]);
}

void MyWideConditionalPrior::from_uniform(std::vector<double>& vec) const
{
	vec[0] = x_min + (x_max - x_min)*vec[0];
	vec[1] = -mu_amp*log(1. - vec[1]);
	if (vec[2] < 0.5)
		vec[2] = mu_logq + sigma_logq*log(2.*vec[2]);
	else
		vec[2] = mu_logq - sigma_logq*log(2. - 2.*vec[2]);

}

void MyWideConditionalPrior::to_uniform(std::vector<double>& vec) const
{
	vec[0] = (vec[0] - x_min)/(x_max - x_min);
	vec[1] = 1. - exp(-vec[1]/mu_amp);
	if (vec[2] < mu_logq)
		vec[2] = 0.5*exp((vec[2] - mu_logq)/sigma_logq);
	else
		vec[2] = 1.0 - 0.5*exp((mu_logq - vec[2])/sigma_logq);

}

void MyWideConditionalPrior::print(std::ostream& out) const
{
	out<<mu_amp<<' '<<mu_logq<<' '<<sigma_logq<<' ';
}


MyNarrowConditionalPrior::MyNarrowConditionalPrior(double x_min, double x_max)
:x_min(x_min)
,x_max(x_max)
{

}

void MyNarrowConditionalPrior::from_prior(RNG& rng)
{
	// Cauchy prior on the mean of the exponential amplitude prior
	mu_amp = tan(M_PI*(0.97*rng.rand() - 0.485));
	mu_amp = exp(mu_amp);

	// Uniform prior on the mean of the Laplacian logq prior:
	mu_logq = (log(100.)-log(1E-5))*rng.rand() + log(1E-5);
	// Uniform prior on the width of the Laplacian logq prior:
	sigma_logq = 2.*rng.rand(); 

//	mu_widths = exp(log(1E-3*(x_max - x_min)) + log(1E3)*rng.rand());

//	a = -10. + 20.*rng.rand();
//	b = 2.*rng.rand();
}

double MyNarrowConditionalPrior::perturb_hyperparameters(RNG& rng)
{
	double logH = 0.;

	int which = rng.rand_int(3);

	if(which == 0)
	{
		mu_amp = log(mu_amp);
		mu_amp = (atan(mu_amp)/M_PI + 0.485)/0.97;
		mu_amp += rng.randh();
 		wrap(mu_amp, 0., 1.);
		mu_amp = tan(M_PI*(0.97*mu_amp - 0.485));
		mu_amp = exp(mu_amp);
	}
	if(which == 1)
	{
		// check this!
		mu_logq += rng.randh()*(log(100.)-log(1E-5)); //log(100)*pow(10., log(2.) - log(100.)*rng.rand())*rng.randn();
		wrap(mu_logq, log(1E-5), log(100));
		//mu_logq = mod(mu_logq - log(.2), log(100.)) + log(2.);
	}
	if(which == 2)
	{
		sigma_logq += 2.*rng.randh();
		wrap(sigma_logq, 0., 2.);
	}
	return logH;
}

double MyNarrowConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
	if(vec[0] < x_min || vec[0] > x_max || vec[1] < 0.0) 
		return -1E300;

	return -log(mu_amp) - vec[1]/mu_amp - log(2.*sigma_logq) - 
			std::abs(vec[2]-mu_logq)/sigma_logq;
	//return -log(mu) - vec[1]/mu - log(mu_widths)
	//		- (vec[2] - min_width)/mu_widths - log(2.*b*vec[3]);
}

void MyNarrowConditionalPrior::from_uniform(std::vector<double>& vec) const
{
	vec[0] = x_min + (x_max - x_min)*vec[0];
	vec[1] = -mu_amp*log(1. - vec[1]);
	if (vec[2] < 0.5)
		vec[2] = mu_logq + sigma_logq*log(2.*vec[2]);
	else
		vec[2] = mu_logq + sigma_logq*log(2. - 2.*vec[2]);

}

void MyNarrowConditionalPrior::to_uniform(std::vector<double>& vec) const
{
	vec[0] = (vec[0] - x_min)/(x_max - x_min);
	vec[1] = 1. - exp(-vec[1]/mu_amp);
	if (vec[2] < mu_logq)
		vec[2] = 0.5*exp((vec[2] - mu_logq)/sigma_logq);
	else
		vec[2] = 1.0 - 0.5*exp((mu_logq - vec[2])/sigma_logq);

}

void MyNarrowConditionalPrior::print(std::ostream& out) const
{
	out<<mu_amp<<' '<<mu_logq<<' '<<sigma_logq<<' ';
}


