#include <vector>
#include <cmath>
#define ll vector<double>
/** parameter for vectoring operations - dont touch, dependeant on double implementation */
#define N 4
/** parameter for multiple calculations between rows */
#define S 16
typedef double double4_t __attribute__((vector_size(N * sizeof(double))));
#define ll4_t vector<double4_t>

using namespace std;

void normalize_rows(int ny, int nx, int nnx, const float *data, ll4_t &normalized)
{
#pragma omp parallel for
  for (int j = 0; j < ny; j++)
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
}

void set_result_to_zero(int ny, float *result)
{
#pragma omp parallel for
  for (int j = 0; j < ny; j++)
    for (int i = j; i < ny; i++)
      result[i + j * ny] = 0;
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
  normalize_rows(ny, nx, nnx, data, normalized);
  set_result_to_zero(ny, result);

#pragma omp parallel for
  for (int j = 0; j < ny; j++)
    for (int i = j; i < ny; i++)
    {
      // * if calculated from a previous cycle, skip
      if (result[i + j * ny])
        continue;
      if (i - j < S || j > ny - S || i > ny - S)
      {
        // * edge cases, do the simpler version
        double4_t sum = {0.0};
        for (int k = 0; k < nnx; k++)
          sum += normalized[k + i * nnx] * normalized[k + j * nnx];

        double cumsum = 0.0;
        for (int h = 0; h < N; h++)
          cumsum += sum[h];
        result[i + j * ny] = cumsum;
      }
      else
      {
        ll4_t sum(S * S);
        for (int k = 0; k < nnx; k++)
        {
          for (int s1 = 0; s1 < S; s1++)
            for (int s2 = 0; s2 < S; s2++)
              sum[s2 + s1 * S] += normalized[k + (i + s2) * nnx] * normalized[k + (j + s1) * nnx];
        }

        ll cumsum(S * S, 0.0);
        for (int s = 0; s < S * S; s++)
          for (int h = 0; h < N; h++)
            cumsum[s] += sum[s][h];

        for (int s1 = 0; s1 < S; s1++)
          for (int s2 = 0; s2 < S; s2++)
            result[i + s2 + (j + s1) * ny] = cumsum[s2 + s1 * S];
        i += S - 1;
      }
    }
}
