#include <vector>
#include <cmath>
#define ll vector<double>
/** parameter for vectoring operations - dont touch, dependeant on double implementation */
#define N 4
/** parameter for instruction parallelism */
#define S 4
typedef double double4_t __attribute__((vector_size(N * sizeof(double))));
#define ll4_t vector<double4_t>

using namespace std;

/** parallelizable by j */
void normalize_row(int j, int nx, int nnx, const float *data, ll4_t &normalized)
{
  double mean = 0.0;
  double magnitude = 0.0;

  for (int i = 0; i < nx; i++)
    mean += data[i + j * nx] / nx;

  for (int i = 0; i < nx; i++)
  {
    int ii = i % N;
    double v = (data[i + j * nx] - mean);
    normalized[i / N + j * nnx][ii] = v;
    magnitude += v * v;
  }
  magnitude = sqrtl(magnitude);

  for (int i = 0; i < nnx; i++)
  {
    normalized[i + j * nnx] /= magnitude;
  }
}

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
  int next_mul_of_N = nx + (N - nx % N) % N;
  int nnx = next_mul_of_N / N;
  ll4_t normalized(ny * nnx);
#pragma omp parallel for
  for (int j = 0; j < ny; j++)
    normalize_row(j, nx, nnx, data, normalized);

#pragma omp parallel for
  for (int j = 0; j < ny; j++)
  {
    for (int i = j; i < ny; i++)
    {
      ll4_t sum(S);
      int k = 0;
      for (; k < nnx - S; k += S)
      {
        for (int kk = 0; kk < S; kk++)
          sum[kk] += normalized[k + kk + i * nnx] * normalized[k + kk + j * nnx];
      }
      for (; k < nnx; k++)
      {
        sum[0] += normalized[k + i * nnx] * normalized[k + j * nnx];
      }

      double4_t cumsum = {0.0};
      for (int ii = 0; ii < S; ii++)
      {
        cumsum += sum[ii];
      }
      result[i + j * ny] = cumsum[0] + cumsum[1] + cumsum[2] + cumsum[3];
    }
  }
}
