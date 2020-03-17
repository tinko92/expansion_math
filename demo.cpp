#include <iostream>
#include <chrono>
#include "dynamic_exact.h"

int main() {
	using A = expansion<double>;
	auto R = A(1e100) + A(1e-20) + A(1e-20) + A(1e-40) - A(1e100) - A(2.0) * A(1e-20);
	std::cout << R.estimate() << "\n"; //completely unaffected by loss of significance
	


	std::cout.precision(15);
	auto start = std::chrono::system_clock::now();
	A f(0.0);
	for(double i = 1.0; i < 2e7; i += 1.0)
		f = f + A(1.0/(i * i));
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "elapsed time: " << elapsed.count() << "ms\n"; //acceptable speed
	std::cout << "result: " << f.estimate() << "\n";	    //precise results
	


	auto det1 = (A(1e-20) - A(1e20)) *
                   (A(2.0) - A(2e20)) -
                   (A(2e-20) - A(2e20)) *
                   (A(1.0 + 1e-10) - A(1e20));
	auto det2 = (A(1e-20) - A(1e20)) *
                   (A(2.0) - A(2e20)) -
                   (A(2e-20) - A(2e20)) *
                   (A(1.0) - A(1e20));
	auto det3 = (A(1e-20) - A(1e20)) *
                   (A(2.0) - A(2e20)) -
                   (A(2e-20) - A(2e20)) *
                   (A(1.0 - 1e-10) - A(1e20));
	std::cout << (det1.estimate() > 0 ? 1 : (det1.estimate() == 0 ? 0 : -1)) << "\n"
		  << (det2.estimate() > 0 ? 1 : (det2.estimate() == 0 ? 0 : -1)) << "\n"
		  << (det3.estimate() > 0 ? 1 : (det3.estimate() == 0 ? 0 : -1)) << "\n";
	//suitable to build robust geometric predicates.
	return 0;
}
