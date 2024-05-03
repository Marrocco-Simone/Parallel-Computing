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
  ll normalized(nx * ny, 0);
  for (int j = 0; j < ny; j++)
  {
    double mean = 0.0;
    double magnitude = 0.0;

    for (int i = 0; i < nx; i++)
      mean += data[i + j * nx] / nx;

    for (int i = 0; i < nx; i++)
    {
      normalized[i + j * nx] = (data[i + j * nx] - mean);
      magnitude += normalized[i + j * nx] * normalized[i + j * nx];
    }
    magnitude = sqrtl(magnitude);

    for (int i = 0; i < nx; i++)
    {
      normalized[i + j * nx] /= magnitude;
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
      double4_t sum = {0.0, 0.0, 0.0, 0.0};
      int k = 0;
      while (k < nx - 4)
      {
        // sum1 += normalized[k + 0 + i * nx] * normalized[k + 0 + j * nx];
        // sum2 += normalized[k + 1 + i * nx] * normalized[k + 1 + j * nx];
        // sum3 += normalized[k + 2 + i * nx] * normalized[k + 2 + j * nx];
        // sum4 += normalized[k + 3 + i * nx] * normalized[k + 3 + j * nx];
        double4_t ni = {
            normalized[k + 0 + i * nx],
            normalized[k + 1 + i * nx],
            normalized[k + 2 + i * nx],
            normalized[k + 3 + i * nx],
        };
        double4_t nj = {
            normalized[k + 0 + j * nx],
            normalized[k + 1 + j * nx],
            normalized[k + 2 + j * nx],
            normalized[k + 3 + j * nx],
        };
        sum += ni * nj;
        k += 4;
      }
      while (k < nx)
      {
        // sum1 += normalized[k + 0 + i * nx] * normalized[k + 0 + j * nx];
        sum[0] += normalized[k + 0 + i * nx] * normalized[k + 0 + j * nx];
        k++;
      }
      // result[i + j * ny] = sum1 + sum2 + sum3 + sum4;
      result[i + j * ny] = sum[0] + sum[1] + sum[2] + sum[3];
    }
  }
}
