#include "MyModel.h"
#include "DNest4/code/Utils.h"
#include "Data.h"
#include <cmath>

using namespace std;
using namespace DNest4;

const Data& MyModel::data = Data::get_instance();
#include <iostream>

MyModel::MyModel()
:narrowlorentzians(3, 100, false, MyNarrowConditionalPrior(data.get_f_min(), data.get_f_max()))
,widelorentzians(3, 100, false, MyWideConditionalPrior(0.0, data.get_f_max()))
,mu(data.get_f().size())
{
}

void MyModel::calculate_mu()
{
        const vector<double>& f_left = data.get_f_left();
        const vector<double>& f_right = data.get_f_right();

        // Update or from scratch?
//        bool update_narrow = (narrowlorentzians.get_added().size() < narrowlorentzians.get_components().size());
//	bool update_wide = (widelorentzians.get_added().size() < widelorentzians.get_components().size());

        // Get the components
//        const vector< vector<double> >& narrowcomponents = (update_narrow)?(narrowlorentzians.get_added()):
//                                (narrowlorentzians.get_components());

//	const vector< vector<double> >& widecomponents = (update_wide)?(widelorentzians.get_added()):
//				(widelorentzians.get_components());

        // Set the background level
//        if(!(update_narrow || update_wide))
//                mu.assign(mu.size(), background);
	mu.assign(mu.size(), background);
	const vector< vector<double> >& narrowcomponents = narrowlorentzians.get_components();
	const vector< vector<double> >& widecomponents = widelorentzians.get_components();


	double f0, amplitude, q;
	double gamma, fac;
        for(size_t j=0; j<widecomponents.size(); j++)
        {
                f0 = widecomponents[j][0];
                amplitude = widecomponents[j][1];
                q = exp(widecomponents[j][2]);
		
		gamma = f0/q;
		
		fac = 0.5*amplitude*gamma/M_PI;

                for(size_t i=0; i<mu.size(); i++)
                {
			// Integral over the Lorentzian distribution
			mu[i] += -fac*(2./gamma)*atan((2.*(f0-f_right[i]))/gamma) +
				  fac*(2./gamma)*atan((2.*(f0-f_left[i]))/gamma); 
                } 

	// this is a OU process; we're currently not going to use that, but it might come 
	// in handy later, so we'll leave it commented out at the appropriate places
//        vector<double> y(mu.size());
//        double alpha = exp(-1./noise_L);
 
//        for(size_t i=0; i<mu.size(); i++)
//        {
//                if(i==0)
//                        y[i] = noise_sigma/sqrt(1. - alpha*alpha)*noise_normals[i];
//                else
//                        y[i] = alpha*y[i-1] + noise_sigma*noise_normals[i];
//                mu[i] *= exp(y[i]);
//        }


        }
        for(size_t j=0; j<narrowcomponents.size(); j++)
        {
                f0 = narrowcomponents[j][0];
                amplitude = narrowcomponents[j][1];
                q = exp(narrowcomponents[j][2]);

		gamma = f0/q;
		
                fac = 0.5*amplitude*gamma/M_PI;

                for(size_t i=0; i<mu.size(); i++)
                {
                        // Integral over the Lorentzian distribution
                        mu[i] += -fac*(2./gamma)*atan((2.*(f0-f_right[i]))/gamma) +
                                  fac*(2./gamma)*atan((2.*(f0-f_left[i]))/gamma);
                }

        // this is a OU process; we're currently not going to use that, but it might come 
        // in handy later, so we'll leave it commented out at the appropriate places
//        vector<double> y(mu.size());
//        double alpha = exp(-1./noise_L);

//        for(size_t i=0; i<mu.size(); i++)
//        {
//                if(i==0)
//                        y[i] = noise_sigma/sqrt(1. - alpha*alpha)*noise_normals[i];
//                else
//                        y[i] = alpha*y[i-1] + noise_sigma*noise_normals[i];
//                mu[i] *= exp(y[i]);
//        }


        }

}

void MyModel::from_prior(RNG& rng)
{
	background = tan(M_PI*(0.97*rng.rand() - 0.485));
	background = exp(background);
	narrowlorentzians.from_prior(rng);
	widelorentzians.from_prior(rng);
 
	// this, too belongs to the noise process we're not using 
//        noise_sigma = exp(log(1E-3) + log(1E3)*rng.rand());
//        noise_L = exp(log(1E-2*Data::get_instance().get_t_range())
//                        + log(1E3)*rng.rand());

        calculate_mu();

}

double MyModel::perturb(RNG& rng)
{
	double logH = 0.;

        if(rng.rand() <= 0.2)
        {
                for(size_t i=0; i<mu.size(); i++)
                        mu[i] -= background;

                background = log(background);
                background = (atan(background)/M_PI + 0.485)/0.97;
                background += pow(10., 1.5 - 6.*rng.rand())*rng.randn();
                background = mod(background, 1.);
                background = tan(M_PI*(0.97*background - 0.485));
                background = exp(background);

                for(size_t i=0; i<mu.size(); i++)
                        mu[i] += background;
        }
        else if(rng.rand() <= 0.7)
        {
                logH += narrowlorentzians.perturb(rng);
//              spikes.consolidate_diff();
                calculate_mu();
        }

	else
	{
		logH += widelorentzians.perturb(rng);
		calculate_mu();
	}
//        else if(rng.rand() <= 0.5)
//        {
//                noise_sigma = log(noise_sigma);
//                noise_sigma += log(1E3)*rng.randh();
//                wrap(noise_sigma, log(1E-3), log(1.));
//                noise_sigma = exp(noise_sigma);

//                noise_L = log(noise_L);
//                noise_L += log(1E3)*rng.randh();
//                wrap(noise_L, log(1E-2*Data::get_instance().get_t_range()), log(10.*Data::get_instance().get_t_range()));
//                noise_L = exp(noise_L);

//                calculate_mu();
//        }
//        else
//        {
//                int num = exp(log((double)noise_normals.size())*rng.rand());
//                for(int i=0; i<num; i++)
//                {
//                        int k = rng.rand_int(noise_normals.size());
//                        noise_normals[k] = rng.randn();
//                }
//                calculate_mu();
//        }


	return logH;
}

double MyModel::log_likelihood() const
{
        const vector<double>& f = data.get_f();
        const vector<double>& y = data.get_y();
		const vector<double>& m = data.get_m();

		// NOTE: This log likelihood is missing a constant factor (not dependent
		
        double logl = 0.;
	    for(size_t i=0; i<f.size(); i++)
			logl += -m[i]*y[i]/mu[i] - 0.5*log(mu[i]) - 0.5*(1./m[i] - 1.)*log(y[i]);

	return logl;
}

void MyModel::print(std::ostream& out) const
{
        out<<background<<' ';
        narrowlorentzians.print(out);
        widelorentzians.print(out);

	for(size_t i=0; i<mu.size(); i++)
                out<<mu[i]<<' ';

}

string MyModel::description() const
{
	return string("");
}

