#include <vector>
#include <cmath>
#define ll vector<double>
typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));
#define ll4_t vector<double4_t>

using namespace std;

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result)
{
  int next_mul_of_4 = nx + (4 - nx % 4) % 4;
  int nnx = next_mul_of_4 / 4;
  ll4_t normalized(ny * nnx);
#pragma omp parallel for
  for (int j = 0; j < ny; j++)
  {
    double mean = 0.0;
    double magnitude = 0.0;

    for (int i = 0; i < nx; i++)
      mean += data[i + j * nx] / nx;

    for (int i = 0; i < nx; i++)
    {
      int ii = i % 4;
      double v = (data[i + j * nx] - mean);
      normalized[i / 4 + j * nnx][ii] = v;
      magnitude += v * v;
    }
    magnitude = sqrtl(magnitude);

    for (int i = 0; i < nnx; i++)
    {
      normalized[i + j * nnx] /= magnitude;
    }
  }

#pragma omp parallel for
  for (int j = 0; j < ny; j++)
  {
    for (int i = j; i < ny; i++)
    {
      // ll4_t sum(4);
      double4_t sum0 = {0.0};
      double4_t sum1 = {0.0};
      double4_t sum2 = {0.0};
      double4_t sum3 = {0.0};
      int k = 0;
      while (k < nnx - 4)
      {
        // for (int kk = 0; kk < 4; kk++)
        //   sum[kk] = normalized[k + kk + i * nnx] * normalized[k + kk + j * nnx];
        sum0 += normalized[k + 0 + i * nnx] * normalized[k + 0 + j * nnx];
        sum1 += normalized[k + 1 + i * nnx] * normalized[k + 1 + j * nnx];
        sum2 += normalized[k + 2 + i * nnx] * normalized[k + 2 + j * nnx];
        sum3 += normalized[k + 3 + i * nnx] * normalized[k + 3 + j * nnx];
        k += 4;
      }
      while (k < nnx)
      {
        sum0 += normalized[k + 0 + i * nnx] * normalized[k + 0 + j * nnx];
        // sum[0] += normalized[k + i * nnx] * normalized[k + j * nnx];
        k++;
      }

      double cumsum0 = 0.0;
      double cumsum1 = 0.0;
      double cumsum2 = 0.0;
      double cumsum3 = 0.0;
      for (int ii = 0; ii < 4; ii++)
      {
        // cumsum0 += sum[0][ii];
        // cumsum1 += sum[1][ii];
        // cumsum2 += sum[2][ii];
        // cumsum3 += sum[3][ii];
        cumsum0 += sum0[ii];
        cumsum1 += sum1[ii];
        cumsum2 += sum2[ii];
        cumsum3 += sum3[ii];
      }
      result[i + j * ny] = cumsum0 + cumsum1 + cumsum2 + cumsum3;
    }
  }
}
