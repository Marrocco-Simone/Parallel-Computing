#include <vector>
#include <cmath>
#define ll vector<long double>
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
  // ll transposed(nx * ny, 0);
  for (int j = 0; j < ny; j++)
  {
    long double mean = 0.0;
    long double magnitude = 0.0;

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
      // transposed[j + i * ny] = normalized[i + j * nx];
    }
  }

  for (int j = 0; j < ny; j++)
  {
    for (int i = j; i < ny; i++)
    {
      long double sum = 0.0;
      for (int k = 0; k < nx; ++k)
      {
        // sum += normalized[k + i * nx] * transposed[j + k * ny];
        sum += normalized[k + i * nx] * normalized[k + j * nx];
      }
      result[i + j * ny] = sum;
    }
  }
}
