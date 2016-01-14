// CDEM_solver.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Domain.h"


int main()
{
	Domain dom = Domain();
	dom.load_from_file("c:/Users/Kapsak/Dropbox/CDEM/Cpp_data/3pb/3pb.txt");
	dom.solve(0.1, 0.2, 50, "c:/Users/Kapsak/Dropbox/CDEM/Cpp_data/3pb/3pbout",1);
	std::cout << dom.nodes[58].v_disp << std::endl;
	/*dom.elements[0].set_matrices();
	std::cout << dom.elements[0].K_loc << std::endl;
	dom.elements[0].print_self();
	std::cout << dom.elements[0].M_loc << std::endl;
	std::cout << dom.elements[0].C_loc << std::endl;*/
	std::cin.get();
    return 0;
}