#include "Vector.h"

vec::vec(uint8_t n) {
	assert(dim<=MAX_LEN_VEC);
	dim = n;
}

// vec Addition
vec vec::operator+(const vec& m) {
	vec b(dim);
	for(uint8_t i=0; i<dim; ++i) {
		b.data[i] = data[i] + m.data[i];
	}
	return b;
}

void vec::operator+=(const vec& m) {
	assert(dim==m.dim);
	for(uint8_t i=0; i<dim; ++i) {
		data[i] += m.data[i];
	}
}

// vec Subtraction
vec vec::operator-(const vec& m) {
	vec b(dim);
	for(uint8_t i=0; i<dim; ++i) {
		b.data[i] = data[i] - m.data[i];
	}
	return b;
}

void vec::operator-=(const vec& m) {
	assert(dim==m.dim);
	for(uint8_t i=0; i<dim; ++i) {
		data[i] -= m.data[i];
	}
}

// Element-by-element multiplication
vec vec::operator*(const vec& m) {
	vec b(dim);
	for(uint8_t i=0; i<dim; ++i) {
		b.data[i] = data[i] * m.data[i];
	}
	return b;
}

void vec::operator*=(const vec& m) {
	assert(dim==m.dim);
	for(uint8_t i=0; i<dim; ++i) {
		data[i] *= m.data[i];
	}
}

// Element-by-element division
vec vec::operator/(const vec& m) {
	vec b(dim);
	for(uint8_t i=0; i<dim; ++i) {
		b.data[i] = data[i] / m.data[i];
	}
	return b;
}

void vec::operator/=(const vec& m) {
	assert(dim==m.dim);
	for(uint8_t i=0; i<dim; ++i) {
		data[i] /= m.data[i];
	}
}

vec vec::operator*(const double k) {
	vec b(dim);
	for(uint8_t i=0; i<dim; ++i) {
		b.data[i] = data[i] * k;
	}
	return b;
}

void vec::operator*=(const double k) {
	for(uint8_t i=0; i<dim; ++i) {
		data[i] *= k;
	}
}

vec vec::operator/(const double k) {
	vec b(dim);
	for(uint8_t i=0; i<dim; ++i) {
		b.data[i] = data[i] / k;
	}
	return b;
}

void vec::operator/=(const double k) {
	for(uint8_t i=0; i<dim; ++i) {
		data[i] /= k;
	}	
}


// Dot product
double vec::dot(const vec& v) {
	double val;
	assert(dim==v.dim);
	for(uint8_t i=0; i<dim; ++i) {
		val += data[i]*v.data[i];
	}
	return val;
}

// Cross product
vec vec::cross(const vec& v) {
	vec b(3);
	b.data[0] = data[1]*v.data[2] - data[2]*v.data[1];
	b.data[1] = data[2]*v.data[0] - data[0]*v.data[2];
	b.data[2] = data[0]*v.data[1] - data[1]*v.data[0];
	return b;
}

// Vector access
double& vec::operator()(uint8_t m) {
	assert(m<dim);
    return data[m];
}


// Vector Initializers
vec vec::zeros(uint8_t n) {
	vec b(n);
	for(int i=0; i<b.dim; ++i) {
		b.data[i] = 0;
	}
	return b;
}

vec vec::ones(uint8_t n) {
	vec b(n);
	for(uint8_t i=0; i<b.dim; ++i) {
		b.data[i] = 1;
	}
	return b;
}

void vec::print(String & str) {
	str="[";
    for (uint8_t i=0;i<dim;++i) {
        if (data[i] >=0.0f){
          str+=" ";
        }
        str+=String(data[i],6);
        str+=" ";
    }
	str+="]\n";

}

