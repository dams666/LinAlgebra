#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"
#define MAX_LEN_MAT 6

class mat {
	friend class mat;
	friend class vec;
	
	public:
		// Constructor
		mat();
		mat(uint8_t m, uint8_t n);

		// mat operations
		mat operator+(const mat& m);		// Addition
		void operator+=(const mat& m);

		mat operator-(const mat& m);		// Subtraction
		void operator-=(const mat& m);

		mat operator^(const mat& m);		// Element-by-element multiplication
		void operator^=(const mat& m);

		mat operator/(const mat& m);		// Element-by-el**ement division
		void operator/=(const mat& m);

		mat operator*(const mat& m);		// mat multiplication
		void operator*=(const mat& m);

		mat operator*(const double k);		// scalar multiplication
		void operator*=(const double k);

		mat operator/(const double k);		// scalar division
		void operator/=(const double k);

		// mat functions
		double det(); 					// Determinant
		mat t();						// Transpose
		mat inv();						// Inverse
		mat inv2();						// Inverse
		void resize(uint8_t m, uint8_t n);		// Resize
		void copy(const mat& m);				// Copy
		mat cols_cat(const mat& m);			// Horizontal concatenation
		mat rows_cat(const mat& m);			// Vertical concatenation

		// Access functions
		vec get_row(uint8_t m);				// Get row
		vec get_col(uint8_t n);				// Get column
		void set_row(uint8_t m, const vec& row);	// Set row
		void set_col(uint8_t n, const vec& col);	// Set col
		mat get_subm(uint8_t m1, uint8_t m2, uint8_t n1, uint8_t n2);	// Get matrix subset
		void set_subm(uint8_t m, uint8_t n, const mat& subm);	// Set matrix subset
		
		void get_subm(uint8_t m1, uint8_t m2, uint8_t n1, uint8_t n2, mat& b);
		void get_row(uint8_t m, vec& b);
			
		mat cofactor(uint8_t m, uint8_t n);	// Get cofactor matrix
		
		double& operator()(uint8_t m, uint8_t n);

		// Initializers
		static mat zeros(uint8_t m, uint8_t n);
		static mat ones(uint8_t m, uint8_t n);
		static mat identity(uint8_t n);

		// Debug stuff
         void print(String &);			// Dump matrix to string

		// TO-DO: Add relational operators
		
		uint8_t get_rows() const { return rows;}
		uint8_t get_cols() const { return cols;}
	protected:
	
		uint8_t rows;
		uint8_t cols;
		double data[MAX_LEN_MAT][MAX_LEN_MAT];
	
		void row_substract(uint8_t m, const vec& row);
		void row_divide(uint8_t m, const double v );
		
		mat get_mat();					// Return self
};

#endif
