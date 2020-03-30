#include <vector>
#include <random>
#include <iostream>
#include <chrono>
#include <boost/geometry/extensions/triangulation/strategies/cartesian/detail/precise_math.hpp>
#include "dynamic_exact.h"
#include "static_exact.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;

const int samples = 10'000'000;

struct problem { std::array<double, 2> a,b,c; };

int main()
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, 1.0);
	int mod = 1000;
	std::cout << "1 in " << mod << " is guaranteed to require higher precision.\n";
	std::vector<problem> problems;
	problems.reserve(samples);
	for(int i=0; i < samples; ++i) {
		double bfac = 1.0, cfac = 1.0;
		if(i % mod == 0) {
			bfac = 1e20;
			cfac = 1e40;
		}
		problems.emplace_back(problem{dis(gen), dis(gen), dis(gen) * bfac, dis(gen) * bfac, dis(gen) * cfac, dis(gen) * cfac});
	}
	auto start = std::chrono::system_clock::now();
	int sum1 = 0;
	for(int i = 0; i < samples; ++i) {
		auto det = boost::geometry::detail::precise_math::orient2d<double>(problems[i].a, problems[i].b, problems[i].c);
		auto r = det > 0 ? 1 : (det < 0 ? -1 : 0);
		sum1 += r;
	}
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Boost robust orient2:\t\t" << sum1 << "\t" << elapsed.count() << "ms\n";
	
	start = std::chrono::system_clock::now();
	int sum2 = 0;
	for(int i = 0; i < samples; ++i) {
	        Point_2 p(problems[i].a[0], problems[i].a[1]), q(problems[i].b[0], problems[i].b[1]), r(problems[i].c[0], problems[i].c[1]);
	        auto det = CGAL::orientation(p, q, r);
	        sum2 += det > 0 ? 1 : (det < 0 ? -1 : 0);
	}
	end = std::chrono::system_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "CGAL exact predicate:\t\t" << sum2 << "\t" << elapsed.count() << "ms\n";

	start = std::chrono::system_clock::now();
	int sum3 = 0;
	for(int i = 0; i < samples; ++i) {
	        auto det = (problems[i].a[0] - problems[i].c[0]) * (problems[i].b[1] - problems[i].c[1])
			- (problems[i].a[1] - problems[i].c[1]) * (problems[i].b[0] - problems[i].c[0]);
	        sum3 += det > 0 ? 1 : (det < 0 ? -1 : 0);
	}
	end = std::chrono::system_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "naive double computation:\t" << sum3 << "\t" << elapsed.count() << "ms\n";

	using A = float_wrapper<double>;
        start = std::chrono::system_clock::now();
        int sum4 = 0;
        for(int i = 0; i < samples; ++i) {
                auto det = ((A(problems[i].a[0]) + A(-problems[i].c[0])) * (A(problems[i].b[1]) + A(-problems[i].c[1]))
                        + (A(-problems[i].a[1]) + A(problems[i].c[1])) * (A(problems[i].b[0]) + A(-problems[i].c[0]))).sign();
                sum4 += det > 0 ? 1 : (det < 0 ? -1 : 0);
        }
        end = std::chrono::system_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "automatically generated exact:\t" << sum4 << "\t" << elapsed.count() << "ms\n";
	return 0;
}
