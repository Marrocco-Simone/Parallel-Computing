#include <vector>
#include <cmath>
#define ll vector<double>
/** parameter for vectoring operations - dont touch, dependeant on double implementation */
#define N 4
#define S 2
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

#pragma omp parallel for
  for (int j = 0; j < ny; j++)
    for (int i = j; i < ny; i++)
      result[i + j * ny] = 0;

#pragma omp parallel for
  for (int j = 0; j < ny; j++)
    for (int i = j; i < ny; i++)
    {
      if (result[i + j * ny])
        continue;
      if (i - j < S || j > ny - S || i > ny - S)
      {
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
        ll4_t sum(4);
        for (int k = 0; k < nnx; k++)
        {
          sum[0] += normalized[k + i * nnx] * normalized[k + j * nnx];
          sum[1] += normalized[k + (i + 1) * nnx] * normalized[k + j * nnx];
          sum[2] += normalized[k + i * nnx] * normalized[k + (j + 1) * nnx];
          sum[3] += normalized[k + (i + 1) * nnx] * normalized[k + (j + 1) * nnx];
        }

        double4_t cumsum = {0.0};
        for (int g = 0; g < 4; g++)
          for (int h = 0; h < N; h++)
            cumsum[g] += sum[g][h];

        result[i + j * ny] = cumsum[0];
        result[(i + 1) + j * ny] = cumsum[1];
        result[i + (j + 1) * ny] = cumsum[2];
        result[(i + 1) + (j + 1) * ny] = cumsum[3];
        i++;
      }
    }
}
