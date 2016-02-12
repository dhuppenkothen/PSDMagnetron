#include <iostream>
#include "DNest4/code/Start.h"
#include "MyModel.h"
#include "Data.h"

using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{
	CommandLineOptions data_file(argc, argv);
	Data::get_instance().load(data_file.get_data_file().c_str());
	Sampler<MyModel> sampler = setup<MyModel>(argc, argv);
	sampler.run();
	return 0;
}

