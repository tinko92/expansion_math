#include <array>
#include <cmath>
#include <limits>
#include <iostream>
#include <cassert>
#include <numeric>
#include "dynamic_exact.h"

template<typename F, typename Left, typename Right>
class sum;

template<typename F, typename Left, typename Right, bool try_skipping_subtrees = false >
class product;

template<typename F>
class float_wrapper {
public:
	using float_type = F;
	static constexpr int error_coefficients_size = 0;
	static constexpr std::array<std::size_t, error_coefficients_size> relative_error_coefficients() {
		std::array<std::size_t, error_coefficients_size> c;
		return c;
	}
	F error_bound() { return F(0.0); }
	static constexpr std::size_t size = 1;
	//dconst float_type& operator[](std::size_t idx) const { assert(idx == 0); return m_val; }
	float_type& operator[](std::size_t idx) { assert(idx == 0); return m_val; }
	float_wrapper(const F& a) : m_val(a) {}
	float_type estimate() const { return m_val; }
	static constexpr bool approx_exact = true;
	static constexpr bool sign_exact = true;
	void eval() {}
	void approx() {}
	F sign() { return m_val; }
	template<typename Rhs>
	sum<F, float_wrapper, Rhs> operator+(Rhs r) {
		return sum<F, float_wrapper, Rhs>(*this, r);
	}
	template<typename Rhs>
	product<F, float_wrapper, Rhs> operator*(Rhs r) {
                return product<F, float_wrapper, Rhs>(*this,r);
        }
private:
	F m_val;
};

template<typename F>
inline void two_sum_tail(const F a, const F b, const F x, F& y) {
        F b_virtual = x - a;
        F a_virtual = x - b_virtual;
        F b_roundoff = b - b_virtual;
        F a_roundoff = a - a_virtual;
        y = a_roundoff + b_roundoff;
}

constexpr std::size_t next_pow2(std::size_t n) {
	std::size_t r = 2;
	++n;
	while(n > r) r *= 2;
	return r;
}

template<typename F, typename Left, typename Right>
class sum {
public:
	template<typename Rhs>
        sum<F, sum, Rhs> operator+(Rhs r) {
                return sum<F, sum, Rhs>(*this,r);
        }       
        template<typename Rhs>
        product<F, sum, Rhs> operator*(Rhs r) {
                return product<F, sum, Rhs>(*this,r);
        }
	using float_type = F;
	static constexpr std::size_t size = Left::size + Right::size;
	static constexpr int error_coefficients_size = 1;
	static constexpr std::array<std::size_t, error_coefficients_size> relative_error_coefficients() {
		std::array<std::size_t, error_coefficients_size> c{1};
		return c;
	}
	float_type& operator[](std::size_t idx) { return m_components[idx]; }
	static constexpr bool approx_exact = false;
	static constexpr bool sign_exact = Left::approx_exact && Right::approx_exact;
	sum(Left& l, Right& r) : m_l(l), m_r(r) {}
	inline void approx() {
		m_l.approx();
		m_r.approx();
		m_components[size - 1] = m_l[Left::size - 1] + m_r[Right::size - 1];
	}
	inline void eval() {
		m_l.eval();
		m_r.eval();
		if constexpr(Left::size == 1 && Right::size == 1) {
			two_sum_tail<float_type>(m_l[0], m_r[0], m_components[size - 1], m_components[0]);
		} else {
        		F Q;
        		std::size_t i_l = 0;
        		std::size_t i_r = 0;
        		std::size_t i = 0;
        		if(std::abs(m_l[i_l]) > std::abs(m_r[i_r]))
                		Q = m_r[i_r++];
        		else
                		Q = m_l[i_l++];
        		if ((i_r < Right::size) && (i_l < Left::size)) {
                		if (std::abs(m_l[i_l]) > std::abs(m_r[i_r]))
                        		fast_two_sum(m_r[i_r++], Q, Q, m_components[i++]);
                		else
                        		fast_two_sum(m_l[i_l++], Q, Q, m_components[i++]);
                		while ((i_r < Right::size) && (i_l < Left::size))
                        		if (std::abs(m_l[i_l]) > std::abs(m_r[i_r]))
                                		two_sum(Q, m_r[i_r++], Q, m_components[i++]);
                        		else
                                		two_sum(Q, m_l[i_l++], Q, m_components[i++]);
        		}
        		while (i_r < Right::size)
                		two_sum(Q, m_r[i_r++], Q, m_components[i++]);
        		while (i_l < Left::size)
				two_sum(Q, m_l[i_l++], Q, m_components[i++]);
			m_components[i] = Q;
		}
	}

	static constexpr int largest_sub_error_coefficient_size = 
		Left::error_coefficients_size > Right::error_coefficients_size ?    
                        Left::error_coefficients_size : Right::error_coefficients_size;
	static constexpr std::array< std::size_t, largest_sub_error_coefficient_size>
		sign_relative_error_coefficients() {
		auto l = Left::relative_error_coefficients();
		auto r = Right::relative_error_coefficients();
		//for the purpose of this demonstration, we only cover the 
		//symmetrical case case of l = r;

		//(1-eps)|A| > ... <=> |A| >= ...
		for(std::size_t i = 1; i < Left::error_coefficients_size; ++i) {
			l[i] += l[i - 1];
		}
		l.back()++; //rounding up
		// * (1 + eps)^2
		for(int i = 0; i < 2; ++i) {
			for(int j = Left::error_coefficients_size - 1; j > 0; --j) {
				l[j] += l[j - 1];
			}
		}
		l.back()++;
		return l;
	}

	static constexpr F relative_error_bound() {
		constexpr auto c = sign_relative_error_coefficients();
		constexpr std::size_t m = next_pow2(c[0]);
		constexpr F eps = std::numeric_limits<F>::epsilon() / 2;
		return (c[0] + eps * ( (c[1] / m + 1) * m ) ) * eps;
	}

	F error_bound() {
		//for this demonstration, we only cover the case in which
		//Left and Right have the same structur
		return    (std::abs(m_l[Left::size - 1]) + std::abs(m_r[Right::size - 1]))
			* relative_error_bound();
	}

	F sign() {
		approx();
		if constexpr(sign_exact) {
			return m_components.back();
		}
		if constexpr(Left::sign_exact && Right::sign_exact) {
			if(    (m_l[Left::size - 1] >= 0 && m_r[Right::size - 1] >= 0)
			    || (m_l[Left::size - 1] <= 0 && m_r[Right::size - 1] <= 0)) {
				return m_components.back();
			}
		}
		if( std::abs(m_components.back()) >= error_bound() ) {
			return m_components.back();
		}
		if constexpr( !Left::sign_exact || !Right::sign_exact) {
			F l = m_l.sign();
			F r = m_r.sign();
			if( (l >= 0 && r >= 0) || (l <= 0 && r <= 0) )
				return l + r;
		}
		eval();
		return std::accumulate(m_components.begin(), m_components.end(), F(0));
	}
private:
	std::array<F, Left::size + Right::size> m_components;
	Left m_l;
	Right m_r;
};

template<typename F>
inline void two_product_tail(const F a, const F b, const F x, F& y) {
        F a_hi, a_lo, b_hi, b_lo;
        split(a, a_hi, a_lo);
        split(b, b_hi, b_lo);
        F err_1 = x - (a_hi * b_hi);
        F err_2 = err_1 - (a_lo * b_hi);
        F err_3 = err_2 - (a_hi * b_lo);
        y = (a_lo * b_lo) - err_3;
}

template<typename F, typename Left, typename Right, bool try_skipping_subtrees>
class product {
public:
	using float_type = F;
	template<typename Rhs>
        sum<F, product, Rhs> operator+(Rhs r) {
                return sum<F, product, Rhs>(*this,r);
        }       
        template<typename Rhs>
        product<F, product, Rhs> operator*(Rhs r) {
                return product<F, product, Rhs>(*this,r);
        }
	float_type& operator[](std::size_t idx) { return m_components[idx]; }
	static constexpr std::size_t size = 2 * Left::size * Right::size;
	static constexpr bool approx_exact = false;
	static constexpr bool sign_exact = Left::sign_exact && Right::sign_exact;
	static constexpr int error_coefficients_size =
		  Left::error_coefficients_size
		+ Right::error_coefficients_size
		+ 1;
        static constexpr std::array<std::size_t, error_coefficients_size> relative_error_coefficients() {
                std::array<std::size_t, error_coefficients_size> c{};
                std::array<std::size_t, Left::error_coefficients_size> l =
			Left::relative_error_coefficients();
		std::array<std::size_t, Left::error_coefficients_size> r =
                        Right::relative_error_coefficients();
		c[0] = 1 + l[0] + r[0];
		for(int i = 1; i < error_coefficients_size; ++i) {
			c[i] =    (i < Left::error_coefficients_size ? l[i] : 0)
				+ (i - 1 < Left::error_coefficients_size ? l[i - 1] : 0)
				+ (i < Right::error_coefficients_size ? r[i] : 0)
                                + (i - 1 < Right::error_coefficients_size ? r[i - 1] : 0);
			for(int j = 0; j <= i - 1; ++j)
				c[i] +=   ( j < Left::error_coefficients_size ? l[j] : 0)
					* ( i - j - 1 < Right::error_coefficients_size ? r[i - j - 1] : 0);
			for(int j = 0; j <= i - 2; ++j)
                                c[i] +=   ( j < Left::error_coefficients_size ? l[j] : 0)
					* ( i - j - 2 < Right::error_coefficients_size ? r[i - j - 2] : 0);
		}
                return c;
        }
	product(Left& l, Right& r) : m_l(l), m_r(r) {}
	inline void approx() {
		m_l.approx();
		if constexpr(try_skipping_subtrees) {
			if(m_l[Left::size - 1] == F(0.0)) {
				m_components[size - 1] = F(0.0);
				return;
			}
		}
		m_r.approx();
		m_components[size - 1] = m_l[Left::size - 1] * m_r[Right::size - 1];
	}
	inline void eval() {
		m_l.eval();
		m_r.eval();
		if constexpr(Left::size == 1 && Right::size == 1) {
			two_product_tail<float_type>(m_l[0], m_r[0], m_components[1], m_components[0]);
		} else if constexpr(Left::size == 2 && Right::size == 2){
			//the following is adapted from two-two-product by Shewchuk
			F l0_hi, l0_lo, r0_hi, r0_lo, l1_hi, l1_lo, r1_hi, r1_lo, tmp1, 
                          tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
                        split(m_l[0], l0_hi, l0_lo);
                        split(m_r[0], r0_hi, r0_lo);
                        two_product(m_l[0], l0_hi, l0_lo, m_r[0], r0_hi, r0_lo, tmp1, m_components[0]);  
                        split(m_l[1], l1_hi, l1_lo);  
                        two_product(m_l[1], l1_hi, l1_lo, m_r[0], r0_hi, r0_lo, tmp2, tmp3);  
                        two_sum(tmp1, tmp3, tmp4, tmp5);  
                        fast_two_sum(tmp2, tmp4, tmp6, tmp7);  
                        split(m_r[1], r1_hi, r1_lo);  
                        two_product(m_l[0], l0_hi, l0_lo, m_r[1], r1_hi, r1_lo, tmp1, tmp3);  
                        two_sum(tmp5, tmp3, tmp4, m_components[1]);  
                        two_sum(tmp7, tmp4, tmp2, tmp1);                
                        two_sum(tmp6, tmp2, tmp8, tmp7);  
                        two_product(m_l[1], l1_hi, l1_lo, m_r[1], r1_hi, r1_lo, tmp2, tmp3);  
                        two_sum(tmp1, tmp3, tmp9, tmp3);
                        two_sum(tmp5, tmp3, tmp1, m_components[2]);
                        two_sum(tmp7, tmp1, tmp4, tmp5);
                        two_sum(tmp8, tmp4, tmp6, tmp7);
                        two_sum(tmp2, tmp9, tmp4, tmp3);
                        two_sum(tmp5, tmp3, tmp2, m_components[3]);
                        two_sum(tmp7, tmp2, tmp1, tmp5);
                        two_sum(tmp6, tmp1, tmp8, tmp7);
                        two_sum(tmp5, tmp4, tmp1, m_components[4]);
                        two_sum(tmp7, tmp1, tmp4, m_components[5]);
                        two_sum(tmp8, tmp4, m_components[7], m_components[6]);
		}//General case is ommitted in this demo
	}
	F sign() {
		F l = m_l.sign();
		if constexpr(try_skipping_subtrees)
			return l != F(0.0) ? l * m_r.sign() : F(0.0);
		else
			return l * m_r.sign();
	}
private:
	std::array<F, 2 * Left::size * Right::size> m_components;
	Left m_l;
	Right m_r;
};
