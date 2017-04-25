#include <stdio.h>
#include <math.h>

int gaussianElimination(double** A, double* b, int dim);

int gaussianElimination(double** A, double* b, int dim)
{
	int i, col, row;
	int pvtRow, pvtCol;
	//A[row][col] = a
	double factor, intermediate;

	for (i = 0; i < dim; i++)
	{

		for(row = i; row < dim; row++)
		{

			if (A[i][i] == 0)	// need to pivot?
			{
				if (row+1 < dim)
				if (A[row+1][i] != 0)
				{
					pvtRow = i;
					for(pvtCol = i; pvtCol < dim; pvtCol++)
					{
						double temp = A[pvtRow][pvtCol];
						A[pvtRow][pvtCol] = A[pvtRow+1][pvtCol];
						A[pvtRow+1][pvtCol] = temp;

						temp = b[pvtRow];
						b[pvtRow] = b[pvtRow+1];
						b[pvtRow+1] = temp;
					}
				}
			}

			if(A[row][i] == 0) continue;

			factor = A[row][i];
			b[row] /= factor;
			if (row != i) b[row] -= b[i];

			
			for(col = i; col < dim; col++)
			{

				A[row][col] /= factor;

				if (row == i) continue;

				A[row][col] -= A[i][col];
			}
			
		}
	}


	for (i = dim-1; i > -1; i--)
	{
		for(col = i; col > i-1; col--)
		{
			for(row = col-1; row > -1; row--)
			{
				b[row] -= A[row][col] * b[i];
				A[row][col] = 0;
			}
		}
	}

	return 1;
}




