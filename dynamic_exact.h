#include<vector>
#include<cmath>
#include<limits>
#include<iostream>

template<typename F>
class expansion {
public:
	using float_type = F;
	std::size_t size() const { return m_components.size(); }
	const float_type& operator[](std::size_t idx) const { return m_components[idx]; }
	float_type& operator[](std::size_t idx) { return m_components[idx]; }
	expansion(F a) : m_components(1, a) {}
	expansion() {}
	expansion(std::size_t l) : m_components(l) {}
	void swap(expansion& o) { m_components.swap(o.m_components); }
	float_type estimate() const { float_type r(0); for(const auto& c : m_components) r+=c; return r; }
	void resize(std::size_t i) { m_components.resize(i); }
private:
	std::vector<float_type> m_components;
};

template<typename F>
void print_expansion(const expansion<F>& e) {
	for(std::size_t i = 0; i < e.size(); ++i) {
		std::cout << i << ": " << e[i] << "\n";
	}
}

template<typename F>
inline void fast_two_sum(const F a, const F b, F& x, F& y) {
	x = a + b;
	F b_virtual = x - a;
	y = b - b_virtual;
}

template<typename F>
inline void fast_two_diff(const F a, const F b, F& x, F& y) {
        x = a - b;
        F b_virtual = a - x;
        y = b_virtual - b;
}

template<typename F>
inline void two_sum(const F a, const F b, F& x, F& y) {
	x = a + b;
	F b_virtual = x - a;
	F a_virtual = x - b_virtual;
	F b_roundoff = b - b_virtual;
	F a_roundoff = a - a_virtual;
	y = a_roundoff + b_roundoff;
}

template<typename F>
inline void two_diff(const F a, const F b, F& x, F& y) {
        x = a - b;
        F b_virtual = a - x;
        F a_virtual = x + b_virtual;
        F b_roundoff = b_virtual - b;
        F a_roundoff = a - a_virtual;     
        y = a_roundoff + b_roundoff;
}

template<typename F>
inline void grow_expansion(const expansion<F>& e, F b, expansion<F>& h) {
	F Q = b;
	for(std::size_t i = 0; i < e.size(); ++i)
		two_sum(Q, e[i], Q, h[i]);
	h[e.size()] = Q;
}

template<typename F>
inline void grow_expansion_diff(const expansion<F>& e, F b, expansion<F>& h) {
        F Q = -b;    
        for(std::size_t i = 0; i < e.size(); ++i)
                two_sum(Q, e[i], Q, h[i]);
        h[e.size()] = Q;
}

template<typename F>                                                               
inline void grow_expansion_zero_elimination(const expansion<F>& e, F b, expansion<F>& h) {
	std::size_t h_index = 0;
        F Q = b;
        for(std::size_t i = 0; i < e.size(); ++i) {
		F h_next;
                two_sum(Q, e[i], Q, h_next);
		if(h_next != 0.0)
			h[h_index++] = h_next;
        }
	if(Q != 0.0 || h_index == 0)
        	h[h_index++] = Q;
	h.resize(h_index);
}

template<typename F>                                                               
inline void grow_expansion_diff_zero_elimination(const expansion<F>& e, F b, expansion<F>& h) {
        std::size_t h_index = 0;
        F Q = -b;
        for(std::size_t i = 0; i < e.size(); ++i) {
                F h_next;
                two_sum(Q, e[i], Q, h_next);
                if(h_next != 0.0)
                        h[h_index++] = h_next;
        }
        if(Q != 0.0 || h_index == 0)
                h[h_index++] = Q;
        h.resize(h_index);
}

template<typename F>
void expansion_sum(const expansion<F>& e, const expansion<F>& f, expansion<F>& h) {
  	F Q = f[0];
	std::size_t i = 0;
	for (; i < e.size(); ++i)
		two_sum(Q, e[i], Q, h[i]);
	h[i] = Q;
	std::size_t hlast = i;
	for (std::size_t j = 1; j < f.size(); ++j) {
		Q = f[j];
		for (i = j; i <= hlast; ++i)
			two_sum(Q, h[i], Q, h[i]);
		h[++hlast] = Q;
  	}
}

template<typename F>
void expansion_diff(const expansion<F>& e, const expansion<F>& f, expansion<F>& h) {
        F Q = -f[0];
        std::size_t i = 0;
        for (; i < e.size(); ++i)
                two_sum(Q, e[i], Q, h[i]);
        h[i] = Q;
        std::size_t hlast = i;
        for (std::size_t j = 1; j < f.size(); ++j) {
                Q = -f[j];
                for (i = j; i <= hlast; ++i)
                        two_sum(Q, h[i], Q, h[i]);
                h[++hlast] = Q;
        }
}

template<typename F>
void expansion_sum_zero_elimination(const expansion<F>& e, const expansion<F>& f, expansion<F>& h) {
        h = e;
        for(std::size_t i = 0; i < f.size(); ++i) {
		expansion<F> h2(h.size() + 1);
                grow_expansion_zero_elimination(h, f[i], h2);
		h.swap(h2);
	}
}

template<typename F>
void fast_expansion_sum_zero_elimination(const expansion<F>& e, const expansion<F>& f, expansion<F>& h) {
	std::size_t i_e = 0;
	std::size_t i_f = 0;
	F Q, h_next;
	if (std::abs(f[i_f]) > std::abs(e[i_e]))
		Q = e[i_e++];
	else
		Q = f[i_f++];
	std::size_t i_h = 0;
	if ((i_e < e.size()) && (i_f < f.size())) {
		if (std::abs(f[i_f]) > std::abs(e[i_e]))
			fast_two_sum(e[i_e++], Q, Q, h_next);
		else
			fast_two_sum(f[i_f++], Q, Q, h_next);
		if (h_next != 0.0)
			h[i_h++] = h_next;
		while ((i_e < e.size()) && (i_f < f.size())) {
			if (std::abs(f[i_f]) > std::abs(e[i_e]))
				two_sum(Q, e[i_e++], Q, h_next);
			else
				two_sum(Q, f[i_f++], Q, h_next);
			if (h_next != 0.0)
        			h[i_h++] = h_next;
		}
  	}
	while (i_e < e.size()) {
		two_sum(Q, e[i_e++], Q, h_next);
    		if (h_next != 0.0)
      			h[i_h++] = h_next;
  	}
  	while (i_f < f.size()) {
		two_sum(Q, f[i_f++], Q, h_next);
		if (h_next != 0.0)
			h[i_h++] = h_next;
  	}
  	if ((Q != 0.0) || (i_h == 0)) {
    		h[i_h++] = Q;
	}
	h.resize(i_h);
}

template<typename F>
void fast_expansion_sum(const expansion<F>& e, const expansion<F>& f, expansion<F>& h) {
	F Q;
	std::size_t i_e = 0;
	std::size_t i_f = 0;
	std::size_t i_h = 0;
	if(std::abs(f[i_f]) > std::abs(e[i_e]))
		Q = e[i_e++];
	else
		Q = f[i_f++];
	if ((i_e < e.size()) && (i_f < f.size())) {
		if (std::abs(f[i_f]) > std::abs(e[i_e]))
			fast_two_sum(e[i_e++], Q, Q, h[i_h++]);
		else
			fast_two_sum(f[i_f++], Q, Q, h[i_h++]);
		while ((i_e < e.size()) && (i_f < f.size()))
			if (std::abs(f[i_f]) > std::abs(e[i_e]))
				two_sum(Q, e[i_e++], Q, h[i_h++]);
      			else
				two_sum(Q, f[i_f++], Q, h[i_h++]);
	}
	while (i_e < e.size())
		two_sum(Q, e[i_e++], Q, h[i_h++]);
	while (i_f < f.size())
		two_sum(Q, f[i_f++], Q, h[i_h++]);
	h[i_h] = Q;
}

template<typename F>
void fast_expansion_diff(const expansion<F>& e, const expansion<F>& f, expansion<F>& h) {
        F Q;
        std::size_t i_e = 0;
        std::size_t i_f = 0;
        std::size_t i_h = 0;
        if(std::abs(f[i_f]) > std::abs(e[i_e]))
                Q = e[i_e++];
        else
                Q = -f[i_f++];
        if ((i_e < e.size()) && (i_f < f.size())) {
                if (std::abs(f[i_f]) > std::abs(e[i_e]))
                        fast_two_sum(e[i_e++], Q, Q, h[i_h++]);
                else
                        fast_two_sum(-f[i_f++], Q, Q, h[i_h++]);            
                while ((i_e < e.size()) && (i_f < f.size()))
                        if (std::abs(f[i_f]) > std::abs(e[i_e]))
                                two_sum(Q, e[i_e++], Q, h[i_h++]);
                        else
                                two_sum(Q, -f[i_f++], Q, h[i_h++]);
        }
        while (i_e < e.size())
                two_sum(Q, e[i_e++], Q, h[i_h++]);
        while (i_f < f.size())
                two_sum(Q, -f[i_f++], Q, h[i_h++]);
        h[i_h] = Q;
}

template<typename F>
constexpr F pow(F base, std::size_t exp) {
	F r(1);
	for(std::size_t i = 0; i < exp; ++i)
		r *= base;
	return r;
}

template<typename F, std::size_t s = std::numeric_limits<F>::digits/2>
inline void split(F a, F& a_hi, F& a_lo) {
	F c = (pow(2.0, s) + 1.0) * a;
	F a_big = c - a;
	a_hi = c - a_big;
	a_lo = a - a_hi;
}

template<typename F>
void two_product(F a, F b, F& x, F& y) {
	x = a * b;
	F a_hi, a_lo, b_hi, b_lo;
	split(a, a_hi, a_lo);
	split(b, b_hi, b_lo);
	F err_1 = x - (a_hi * b_hi);
	F err_2 = err_1 - (a_lo * b_hi);
	F err_3 = err_2 - (a_hi * b_lo);
	y = (a_lo * b_lo) - err_3;
}

template<typename F>
void two_product(F a, F b, F b_hi, F b_lo, F& x, F& y) {
	x = a * b;
        F a_hi, a_lo;
        split(a, a_hi, a_lo);
        F err_1 = x - (a_hi * b_hi);
        F err_2 = err_1 - (a_lo * b_hi);
        F err_3 = err_2 - (a_hi * b_lo);
        y = (a_lo * b_lo) - err_3;
}

template<typename F>
void two_product(F a, F a_hi, F a_lo, F b, F b_hi, F b_lo, F& x, F& y) {         
        x = a * b;
        F err_1 = x - (a_hi * b_hi);
        F err_2 = err_1 - (a_lo * b_hi);
        F err_3 = err_2 - (a_hi * b_lo);
        y = (a_lo * b_lo) - err_3;
}

template<typename F>
void scale_expansion(const expansion<F>& e, F b, expansion<F>& h) {
	F Q0, Q1, b_hi, b_lo;
  	split(b, b_hi, b_lo);
	two_product(e[0], b, b_hi, b_lo, Q0, h[0]);
	for (std::size_t i = 1; i < e.size(); i++) {
		F T, t, sum;
		two_product(e[i], b, b_hi, b_lo, T, t);
		two_sum(Q0, t, Q1, h[2 * i - 1]);
		two_sum(T, Q1, Q0, h[2 * i]);
	}
	h[2 * e.size() - 1] = Q0;
}

template<typename F>
void compress(const expansion<F>& e, expansion<F>& h) {
	int bottom = e.size() - 1;
	F Q = e[e.size() - 1], q;
	for (int i = bottom - 2; i >= 0; --i) {
    		fast_two_sum(Q, e[i], Q, q);
    		if (q != 0) {
      			h[bottom--] = Q;
      			Q = q;
    		}
	}
	int top = 0;
	for (int i = bottom + 1; i < e.size(); ++i) {
    		fast_two_sum(h[i], Q, Q, q);
		if (q != 0)
      			h[top++] = q;
  	}
	h[top] = Q;
	h.resize(top + 1);
}

template<typename F>
void expansion_product(const expansion<F>& e, const expansion<F>& f, expansion<F>& h) {
	std::vector<expansion<F>> expansions, expansions_next;
	expansions.reserve(f.size());
	expansions_next.reserve( (f.size() + 1) / 2);
	for(std::size_t i = 0; i < f.size(); ++i) {
		expansion<F> next(e.size() * 2);
		scale_expansion(e, f[i], next);
		expansions.push_back(next);
	}
	while(expansions.size() > 1) {
		for(std::size_t i = 0; i * 2 + 1 < expansions.size(); ++i) {
			expansion<F> next(expansions[i * 2].size() + expansions[i * 2 + 1].size());
			fast_expansion_sum_zero_elimination(expansions[i * 2], expansions[i * 2 + 1], next);
			expansions_next.push_back(next);
		}
		expansions.swap(expansions_next);
	}
	h.swap(expansions[0]);
}

template<typename F>
expansion<F> operator+(const expansion<F>& lhs, const expansion<F>& rhs) {
	expansion<F> h(lhs.size() + rhs.size());
	fast_expansion_sum_zero_elimination(lhs, rhs, h);
	return h;
}

template<typename F>
expansion<F> operator-(const expansion<F>& lhs, const expansion<F>& rhs) {
        expansion<F> h(lhs.size() + rhs.size());
        fast_expansion_diff(lhs, rhs, h);
        return h;
}

template<typename F>
expansion<F> operator*(const expansion<F>& lhs, const expansion<F>& rhs) {        
        expansion<F> h(lhs.size() * rhs.size());
        expansion_product(lhs, rhs, h);
        return h;
}

template<typename F>
expansion<F> operator*(const expansion<F>& lhs, const F& rhs) {
        expansion<F> h(lhs.size() * 2);
        scale_expansion(lhs, rhs, h);
        return h;
}

template<typename F>
expansion<F> operator+(const expansion<F>& lhs, const F& rhs) {
        expansion<F> h(lhs.size() + 1);
        grow_expansion_zero_elimination(lhs, rhs, h);
        return h;
}

template<typename F>
expansion<F> operator+(const F& lhs, const expansion<F>& rhs) {
        expansion<F> h(lhs.size() + 1);         
        grow_expansion_zero_elimination(rhs, lhs, h);
        return h;
}

template<typename F>
expansion<F> operator-(const expansion<F>& lhs, const F& rhs) {
        expansion<F> h(lhs.size() + 1);
        grow_expansion_diff_zero_elimination(lhs, rhs, h);
        return h;
}

template<typename F>
expansion<F> operator-(const F& lhs, const expansion<F>& rhs) {
        expansion<F> h(lhs.size() + 1);
	expansion<F> lhs_ex(lhs);
        fast_expansion_diff(lhs_ex, rhs, h);
        return h;
}

template<typename F>
expansion<F> operator-(const expansion<F>& e) {
        expansion<F> h(e.size());
	for(std::size_t i = 0; i < e.size(); ++i)
		h[i] = -e[i];
        return h;
}
