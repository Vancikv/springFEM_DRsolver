// CDEM_solver.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Domain.h"


int main(int argc, char** argv)
{
	//dom.load_from_file("c:/Users/Kapsak/Dropbox/CDEM/Cpp_data/3pb/3pb.txt");
	//dom.solve(1., 2., 4000, "i:/dokumenty/Cpp_CDEM/3pb/3pbout",1);
	//std::cout << dom.nodes[58].v_disp << std::endl;
	//std::cin.get();
	std::string f, o;
	double t1, t2;
	int ns, of, is;
	if (argc == 1)
	{
		std::cout << "This is a simple solver for FEM domain of elements connected via springs." << std::endl 
			<< "Following is a list of arguments, all of which must be given" << std::endl << std::endl;
		std::cout << "-f\tSpecifies input file path." << std::endl;
		std::cout << "-o\tSpecifies output file path." << std::endl;
		std::cout << "-t1\tTime - load function reaches 1.0." << std::endl;
		std::cout << "-t2\tTime - maximum time." << std::endl;
		std::cout << "-ns\tNumber of steps." << std::endl;
		std::cout << "-of\tOutput frequency." << std::endl;
		std::cout << "-is\tNumber of inner steps." << std::endl;
	}else if (argc == 15){
		for (int i = 1; i < 15; i += 2)
		{
			if (std::string(argv[i]) == "-f") {
				f = argv[i + 1];
			}
			else if (std::string(argv[i]) == "-o") {
				o = argv[i + 1];
			}
			else if (std::string(argv[i]) == "-t1") {
				t1 = std::stod(argv[i + 1]);
			}
			else if (std::string(argv[i]) == "-t2") {
				t2 = std::stod(argv[i + 1]);
			}
			else if (std::string(argv[i]) == "-ns") {
				ns = std::stoi(argv[i + 1]);
			}
			else if (std::string(argv[i]) == "-of") {
				of = std::stoi(argv[i + 1]);
			}
			else if (std::string(argv[i]) == "-is") {
				is = std::stoi(argv[i + 1]);
			}
			else {
				std::cout << "Unknown argument " << argv[i] << ", calculation cancelled" << std::endl;
				return 1;
			}
		}
		Domain dom = Domain();
		dom.load_from_file(f);
		dom.solve(t1, t2, ns, o,of,is);
	}else {
		std::cout << "Wrong number of arguments.";
	}
    return 0;
}