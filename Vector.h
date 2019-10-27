#ifndef VECTOR_H
#define VECTOR_H

#define __ASSERT_USE_STDERR
#include <assert.h>
#include <Arduino.h>

#define MAX_LEN_VEC 6

class vec {
	friend class vec;
	
	public:
		vec(uint8_t n); // Constructor
		
			// mat operations
		vec operator+(const vec& m);		// Addition
		void operator+=(const vec& m);

		vec operator-(const vec& m);		// Subtraction
		void operator-=(const vec& m);

		vec operator/(const vec& m);		// Element-by-el**ement division
		void operator/=(const vec& m);

		vec operator/(const double k);		// scalar division
		void operator/=(const double k);
		
		vec operator*(const vec& m);		// vec multiplication
		void operator*=(const vec& m);

		vec operator*(const double k);		// scalar multiplication
		void operator*=(const double k);

		
		// vec specific Operations
		double dot(const vec& v);		// Dot product
		vec cross(const vec& v);		// Cross product
		
		// Vector access
		double& operator()(uint8_t m);
		
		// Vector Initializers
		static vec zeros(uint8_t n);
		static vec ones(uint8_t n);
		
		// Debug stuff
         void print(String &);			// Dump vec to string
	public:
	
		uint8_t dim;
		double data[MAX_LEN_VEC];
};

#endif
