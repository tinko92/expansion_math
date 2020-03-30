#include <iostream>
#include <chrono>
#include "dynamic_exact.h"
#include "static_exact.h"

int main() {
	using A = expansion<double>;
	using B = float_wrapper<double>;
	auto x = A(1e100) + A(1e-20) + A(1e-20) + A(1e-40) - A(1e100) - A(2.0) * A(1e-20);
	std::cout << x.estimate() << "\n"; 
	//output: 1e-40
	//completely unaffected by loss of significance


	std::cout.precision(15);
	auto start = std::chrono::system_clock::now();
	A s(0.0);
	for(double i = 1.0; i < 2e7; i += 1.0)
		s = s + A(1.0/(i * i));
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "elapsed time: " << elapsed.count() << "ms\n";
	//output: elapsed time: 698ms
	//acceptable speed
	std::cout << "result: " << s.estimate() << "\n";
	//result: 1.64493401684823
	//precise results


	auto det1 = ((B(1e-20) + B(-1e20)) *
                   (B(2.0) + B(-2e20)) +
                   (B(-2e-20) + B(2e20)) *
                   (B(1.0 + 1e-10) + B(-1e20))).sign();
	auto det2 = ((B(1e-20) + B(-1e20)) *
                   (B(2.0) + B(-2e20)) +
                   (B(-2e-20) + B(2e20)) *
                   (B(1.0) + B(-1e20))).sign();
	auto det3 = ((B(1e-20) + B(-1e20)) *
                   (B(2.0) + B(-2e20)) +
                   (B(-2e-20) + B(2e20)) *
                   (B(1.0 - 1e-10) + B(-1e20))).sign();
	std::cout << (det1 > 0 ? 1 : (det1 == 0 ? 0 : -1)) << "\n"
		  << (det2 > 0 ? 1 : (det2 == 0 ? 0 : -1)) << "\n"
		  << (det3 > 0 ? 1 : (det3 == 0 ? 0 : -1)) << "\n";
	//output:
	//1
	//0
	//-1
	//suitable to build robust geometric predicates.
	return 0;
}
