#include "Data.h"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

Data Data::instance;

Data::Data()
{

}

void Data::load(const char* filename)
{
	fstream fin(filename, ios::in);
	if(!fin)
	{
		cerr<<"# Failed to open file "<<filename<<"."<<endl;
		return;
	}

	f.clear();
	y.clear();
	m.clear();

	double temp1, temp2, temp3;
	while(fin>>temp1 && fin>>temp2 && fin>>temp3)
	{
		f.push_back(temp1);
		y.push_back(temp2);
		m.push_back(temp3);
	}

	fin.close();
	cout<<"# Found "<<t.size()<<" points in file "<<filename<<"."<<endl;

	compute_summaries();
}

void Data::compute_summaries()
{
	f_min = *min_element(f.begin(), f.end());
	f_max = *max_element(f.begin(), f.end());
	f_range = f_max - f_min;
	df = f[1] - f[0];

	// Left and right edges of the data bins
	f_left.assign(f.size(), 0.);
	f_right.assign(f.size(), 0.);
	for(size_f i=0; i<f.size(); i++)
	{
		f_left[i] = f[i] - 0.5*df;
		f_right[i] = f[i] + 0.5*df;
	}
}

