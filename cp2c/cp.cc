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

  for (int j = 0; j < ny; j++)
  {
    for (int i = j; i < ny; i++)
    {
      // * instruction level parallelism
      // double sum1 = 0.0;
      // double sum2 = 0.0;
      // double sum3 = 0.0;
      // double sum4 = 0.0;
      double4_t sum = {0.0};
      for (int k = 0; k < nnx; k++)
      {
        // sum1 += normalized[k + 0 + i * nx] * normalized[k + 0 + j * nx];
        // sum2 += normalized[k + 1 + i * nx] * normalized[k + 1 + j * nx];
        // sum3 += normalized[k + 2 + i * nx] * normalized[k + 2 + j * nx];
        // sum4 += normalized[k + 3 + i * nx] * normalized[k + 3 + j * nx];
        sum += normalized[k + i * nnx] * normalized[k + j * nnx];
      }
      // result[i + j * ny] = sum1 + sum2 + sum3 + sum4;
      result[i + j * ny] = sum[0] + sum[1] + sum[2] + sum[3];
    }
  }
}
