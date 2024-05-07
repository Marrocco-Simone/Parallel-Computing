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
      ll4_t sum(4);
      int k = 0;
      while (k < nnx - 4)
      {
        for (int kk = 0; kk < 4; kk++)
          sum[kk] += normalized[k + kk + i * nnx] * normalized[k + kk + j * nnx];
        k += 4;
      }
      while (k < nnx)
      {
        sum[0] += normalized[k + i * nnx] * normalized[k + j * nnx];
        k++;
      }

      double4_t cumsum = {0.0};
      for (int ii = 0; ii < 4; ii++)
      {
        cumsum += sum[ii];
      }
      result[i + j * ny] = cumsum[0] + cumsum[1] + cumsum[2] + cumsum[3];
    }
  }
}
