// CDEM_solver.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Domain.h"
using namespace std;

int main()
{
	Domain dom = Domain();
	dom.load_from_file("c:/Users/Kapsak/Dropbox/CDEM/Cpp_data/test/test.TXT");
	dom.solve(0.01, 0.02, 20);
	std::cout << dom.nodes[3].v_disp << endl;
	/*dom.elements[0].set_matrices();
	std::cout << dom.elements[0].K_loc << endl;
	dom.elements[0].print_self();*/
	std::cin.get();
    return 0;
}