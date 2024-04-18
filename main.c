#include "return_codes.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TOLERANCE 1e-9
#define BALANCE_RADIX 2.0

#define HOUSEHOLDER_VECTOR_NOT_EXIST 101

typedef unsigned int uint;

typedef struct
{
	double real;
	double imaginary;
} ComplexDouble;

char *error_message_table[] = {
	[SUCCESS] = "",
	[ERROR_CANNOT_OPEN_FILE] = "File cannot be opened",
	[ERROR_OUT_OF_MEMORY] = "Not enough memory, memory allocation failed",
	[ERROR_DATA_INVALID] = "The data is invalid",
	[ERROR_PARAMETER_INVALID] = "The parameter or number of parameters is incorrect",

	[ERROR_UNSUPPORTED] = "Unsupported functionality",
};

double VectorNorm(uint size, uint step, const double (*vector)[step])
{
	double result = 0;
	for (uint i = 0; i < size; ++i)
	{
		result += vector[i][0] * vector[i][0];
	}
	return sqrt(result);
}

void FindEigenvalues2x2(uint step, const double (*matrix)[step], ComplexDouble *first_eigenvalue, ComplexDouble *second_eigenvalue)
{
	double diag_mean = (matrix[0][0] + matrix[1][1]) / 2;
	double det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

	double under_square_value = diag_mean * diag_mean - det;

	if (under_square_value < 0)
	{
		*first_eigenvalue = (ComplexDouble){ diag_mean, sqrt(-under_square_value) };
		*second_eigenvalue = (ComplexDouble){ diag_mean, -sqrt(-under_square_value) };
	}
	else
	{
		*first_eigenvalue = (ComplexDouble){ diag_mean + sqrt(under_square_value), 0 };
		*second_eigenvalue = (ComplexDouble){ diag_mean - sqrt(under_square_value), 0 };
	}
}

int GetHouseholderVector2(uint step, const double (*vector)[step], double *householder_vector)
{
	double norm = sqrt(vector[0][0] * vector[0][0] + vector[1][0] * vector[1][0]);
	short p = vector[0][0] < 0 ? -1 : 1;

	if (norm == 0)
	{
		return HOUSEHOLDER_VECTOR_NOT_EXIST;
	}

	householder_vector[0] = vector[0][0] + p * norm;
	householder_vector[1] = vector[1][0];

	double resultNorm = sqrt(householder_vector[0] * householder_vector[0] + householder_vector[1] * householder_vector[1]);

	householder_vector[0] /= resultNorm;
	householder_vector[1] /= resultNorm;

	return 0;
}

int GetHouseholderVector3(uint step, const double (*vector)[step], double *householder_vector)
{
	double norm = sqrt(vector[0][0] * vector[0][0] + vector[1][0] * vector[1][0] + vector[2][0] * vector[2][0]);
	short p = vector[0][0] < 0 ? -1 : 1;

	if (norm == 0)
	{
		return HOUSEHOLDER_VECTOR_NOT_EXIST;
	}

	householder_vector[0] = vector[0][0] + p * norm;
	householder_vector[1] = vector[1][0];
	householder_vector[2] = vector[2][0];

	double resultNorm =
		sqrt(householder_vector[0] * householder_vector[0] + householder_vector[1] * householder_vector[1] +
			 householder_vector[2] * householder_vector[2]);

	householder_vector[0] /= resultNorm;
	householder_vector[1] /= resultNorm;
	householder_vector[2] /= resultNorm;

	return 0;
}

void ApplyHouseholderTransformationLeft2x2(uint width, uint step, double (*matrix)[step], const double *householder_vector)
{
	for (uint j = 0; j < width; ++j)
	{
		double t = 2 * (householder_vector[0] * matrix[0][j] + householder_vector[1] * matrix[1][j]);

		double change_values[] = {
			householder_vector[0] * t,
			householder_vector[1] * t,
		};

		matrix[0][j] -= change_values[0];
		matrix[1][j] -= change_values[1];
	}
}

void ApplyHouseholderTransformationRight2x2(uint height, uint step, double (*matrix)[step], const double *householder_vector)
{
	for (uint i = 0; i < height; ++i)
	{
		double t = 2 * (householder_vector[0] * matrix[i][0] + householder_vector[1] * matrix[i][1]);

		double change_values[] = {
			householder_vector[0] * t,
			householder_vector[1] * t,
		};

		matrix[i][0] -= change_values[0];
		matrix[i][1] -= change_values[1];
	}
}

void ApplyHouseholderTransformationLeft3x3(uint width, uint step, double (*matrix)[step], const double *householder_vector)
{
	for (uint j = 0; j < width; ++j)
	{
		double t =
			2 * (householder_vector[0] * matrix[0][j] + householder_vector[1] * matrix[1][j] + householder_vector[2] * matrix[2][j]);

		double change_values[] = {
			householder_vector[0] * t,
			householder_vector[1] * t,
			householder_vector[2] * t,
		};

		matrix[0][j] -= change_values[0];
		matrix[1][j] -= change_values[1];
		matrix[2][j] -= change_values[2];
	}
}

void ApplyHouseholderTransformationRight3x3(uint height, uint step, double (*matrix)[step], const double *householder_vector)
{
	for (uint i = 0; i < height; ++i)
	{
		double t =
			2 * (householder_vector[0] * matrix[i][0] + householder_vector[1] * matrix[i][1] + householder_vector[2] * matrix[i][2]);

		double change_values[] = {
			householder_vector[0] * t,
			householder_vector[1] * t,
			householder_vector[2] * t,
		};

		matrix[i][0] -= change_values[0];
		matrix[i][1] -= change_values[1];
		matrix[i][2] -= change_values[2];
	}
}

void ComputeCoefficientsForStep(uint step, const double (*submatrix)[step], double *b, double *c)
{
	double trace = submatrix[0][0] + submatrix[1][1];
	double det = submatrix[0][0] * submatrix[1][1] - submatrix[0][1] * submatrix[1][0];

	if (trace * trace > det * 4)
	{
		double s = sqrt(trace * trace - det * 4);
		double l1 = (trace + s) / 2;
		double l2 = (trace - s) / 2;

		if (fabs(l1 - submatrix[1][1]) < fabs(l2 - submatrix[1][1]))
		{
			l2 = l1;
		}
		else
		{
			l1 = l2;
		}

		*b = -l1 - l2;
		*c = l1 * l2;
	}
	else
	{
		*b = -trace;
		*c = det;
	}
}

void BalanceMatrix(uint size, uint step, double (*matrix)[step])
{
	double radix_sqr = BALANCE_RADIX * BALANCE_RADIX;
	uint l = 0;

	while (l == 0)
	{
		l = 1;
		for (int i = 0; i < size; ++i)
		{
			double r = 0;
			double c = 0;

			for (int j = 0; j < size; ++j)
			{
				if (j != i)
				{
					c += fabs(matrix[j][i]);
					r += fabs(matrix[i][j]);
				}
			}

			if (c != 0 && r != 0)
			{
				double g = r / BALANCE_RADIX;
				double f = 1;
				double s = c + r;

				while (c < g)
				{
					f *= BALANCE_RADIX;
					c *= radix_sqr;
				}

				g = r * BALANCE_RADIX;
				while (c > g)
				{
					f /= BALANCE_RADIX;
					c /= radix_sqr;
				}

				if ((c + r) / f < 0.95 * s)
				{
					l = 0;
					g = 1 / f;

					for (int j = 0; j < size; ++j)
					{
						matrix[i][j] *= g;
					}

					for (int j = 0; j < size; ++j)
					{
						matrix[j][i] *= f;
					}
				}
			}
		}
	}
}

void ToHessenbergForm(uint size, uint step, double (*matrix)[step])
{
	for (int r = 0; r < size - 2; ++r)
	{
		double x = 0;
		uint k = r;

		for (int j = r + 2; j < size; ++j)
		{
			if (fabs(matrix[j][r]) > fabs(x))
			{
				x = matrix[j][r];
				k = j;
			}
		}

		if (x)
		{
			for (int j = 0; j < size; ++j)
			{
				double s = matrix[k][j];
				matrix[k][j] = matrix[r + 1][j];
				matrix[r + 1][j] = s;
			}

			for (int j = 0; j < size; ++j)
			{
				double s = matrix[j][k];
				matrix[j][k] = matrix[j][r + 1];
				matrix[j][r + 1] = s;
			}

			for (int i = r + 2; i < size; ++i)
			{
				double p = matrix[i][r] / matrix[r + 1][r];

				for (int j = 0; j < size; ++j)
				{
					matrix[i][j] -= p * matrix[r + 1][j];
				}

				for (int j = 0; j < size; ++j)
				{
					matrix[j][r + 1] += p * matrix[j][i];
				}
			}
		}
	}
}

void QrDoubleStepAlgorithmStep(uint size, uint step, double (*matrix)[step])
{
	double b, c;
	ComputeCoefficientsForStep(step, (const double(*)[]) & matrix[size - 2][size - 2], &b, &c);

	double v[] = {
		matrix[0][0] * matrix[0][0] + matrix[1][0] * matrix[0][1] + b * matrix[0][0] + c,
		matrix[0][0] * matrix[1][0] + matrix[1][0] * matrix[1][1] + b * matrix[1][0],
		matrix[0][0] * matrix[2][0] + matrix[1][0] * matrix[2][1],
	};

	double householder_vector[3];
	int code;

	code = GetHouseholderVector3(1, (const double(*)[])v, (double *)householder_vector);

	if (code != HOUSEHOLDER_VECTOR_NOT_EXIST)
	{
		ApplyHouseholderTransformationLeft3x3(size, step, matrix, householder_vector);
		ApplyHouseholderTransformationRight3x3(size, step, matrix, householder_vector);
	}

	for (uint i = 1; i < size - 1; ++i)
	{
		if (i > size - 3)
		{
			code = GetHouseholderVector2(step, (const double(*)[]) & matrix[i][i - 1], householder_vector);

			if (code == HOUSEHOLDER_VECTOR_NOT_EXIST)
			{
				continue;
			}

			ApplyHouseholderTransformationLeft2x2(size, step, (double(*)[]) & matrix[i][0], householder_vector);
			ApplyHouseholderTransformationRight2x2(size, step, (double(*)[]) & matrix[0][i], householder_vector);

			matrix[size - 1][i - 1] = 0;
		}
		else
		{
			code = GetHouseholderVector3(step, (const double(*)[]) & matrix[i][i - 1], (double *)householder_vector);

			if (code == HOUSEHOLDER_VECTOR_NOT_EXIST)
			{
				continue;
			}

			ApplyHouseholderTransformationLeft3x3(size, step, (double(*)[]) & matrix[i][0], (double *)householder_vector);
			ApplyHouseholderTransformationRight3x3(size, step, (double(*)[]) & matrix[0][i], (double *)householder_vector);

			matrix[i + 2][i - 1] = 0;
		}
	}
}

void QrDoubleStepAlgorithm(uint size, uint step, double (*matrix)[step])
{
	double eps = VectorNorm(size * size, 1, matrix) * TOLERANCE;

	uint start = 0, end = size;
	while (end - start > 2)
	{
		if (fabs(matrix[end - 1][end - 2]) < eps)
		{
			matrix[end - 1][end - 2] = 0;
			end -= 1;
		}
		else if (fabs(matrix[end - 2][end - 3]) < eps)
		{
			matrix[end - 2][end - 3] = 0;
			end -= 2;
		}
		else if (fabs(matrix[start + 1][start]) < eps)
		{
			matrix[start + 1][start] = 0;
			start += 1;
		}
		else if (fabs(matrix[start + 2][start + 1]) < eps)
		{
			matrix[start + 2][start + 1] = 0;
			start += 2;
		}
		else
		{
			QrDoubleStepAlgorithmStep(end - start, size, (double(*)[]) & matrix[start][start]);
		}
	}
}

void FindEigenvalues(uint size, uint step, double (*matrix)[size], ComplexDouble eigenvalues[size])
{
	if (size > 2)
	{
		BalanceMatrix(size, step, matrix);

		ToHessenbergForm(size, step, matrix);

		QrDoubleStepAlgorithm(size, step, matrix);
	}

	uint i = 0;
	while (i < size)
	{
		if (i == size - 1 || matrix[i + 1][i] == 0)
		{
			eigenvalues[i] = (ComplexDouble){ matrix[i][i], 0 };
			i += 1;
		}
		else
		{
			FindEigenvalues2x2(step, (const double(*)[]) & matrix[i][i], &(eigenvalues[i]), &(eigenvalues[i + 1]));
			i += 2;
		}
	}
}

int ScanMatrix(FILE *input_file, double **matrix, uint *size)
{
	if (fscanf(input_file, "%d", size) != 1 || *size == 0)
	{
		return ERROR_DATA_INVALID;
	}

	*matrix = malloc(*size * *size * sizeof(double));

	if (matrix == NULL)
	{
		return ERROR_OUT_OF_MEMORY;
	}

	for (uint i = 0; i < *size * *size; ++i)
	{
		if (fscanf(input_file, "%lf", &(*matrix)[i]) != 1)
		{
			free(*matrix);
			*matrix = NULL;
			return ERROR_DATA_INVALID;
		}
	}

	return SUCCESS;
}

int ReleaseEndCode(int code)
{
	if (code != SUCCESS)
	{
		fprintf(stderr, "%s\n", error_message_table[code]);
	}
	return code;
}

int main(int argc, char **argv)
{
	if (argc != 3)
	{
		return ReleaseEndCode(ERROR_PARAMETER_INVALID);
	}

	FILE *input_file = fopen(argv[1], "r");

	if (input_file == NULL)
	{
		return ReleaseEndCode(ERROR_CANNOT_OPEN_FILE);
	}

	double *matrix_input;
	uint matrix_size;

	int scan_matrix_code = ScanMatrix(input_file, &matrix_input, &matrix_size);

	fclose(input_file);

	if (scan_matrix_code != SUCCESS)
	{
		return ReleaseEndCode(scan_matrix_code);
	}

	ComplexDouble *eigenvalues = malloc(sizeof(ComplexDouble) * matrix_size);

	if (eigenvalues == NULL)
	{
		free(matrix_input);
		return ReleaseEndCode(ERROR_OUT_OF_MEMORY);
	}

	FindEigenvalues(matrix_size, matrix_size, (double(*)[])matrix_input, eigenvalues);

	free(matrix_input);

	FILE *output_file = fopen(argv[2], "w");

	if (output_file == NULL)
	{
		free(eigenvalues);
		return ReleaseEndCode(ERROR_CANNOT_OPEN_FILE);
	}

	for (uint i = 0; i < matrix_size; ++i)
	{
		if (eigenvalues[i].imaginary == 0)
		{
			fprintf(output_file, "%g\n", eigenvalues[i].real);
		}
		else
		{
			fprintf(output_file, "%g %+gi\n", eigenvalues[i].real, eigenvalues[i].imaginary);
		}
	}

	fclose(output_file);
	free(eigenvalues);

	return ReleaseEndCode(SUCCESS);
}
