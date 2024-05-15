#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cuda_runtime.h>
// #define ll vector<float>
using namespace std;

static inline void check(cudaError_t err, const char *context)
{
  if (err != cudaSuccess)
  {
    cerr << "CUDA error: " << context << ": "
         << cudaGetErrorString(err) << endl;
    exit(EXIT_FAILURE);
  }
}

#define CHECK(x) check(x, #x)

void set_result_to_zero(int ny, float *result)
{
  for (int j = 0; j < ny; j++)
    for (int i = j; i < ny; i++)
      result[i + j * ny] = 0.0;
}

void normalize_rows(int ny, int nx, const float *data, float *normalized)
{
  for (int j = 0; j < ny; j++)
  {
    float mean = 0.0;
    float magnitude = 0.0;

    for (int i = 0; i < nx; i++)
      mean += data[i + j * nx] / nx;

    for (int i = 0; i < nx; i++)
    {
      normalized[i + j * nx] = (data[i + j * nx] - mean);
      magnitude += normalized[i + j * nx] * normalized[i + j * nx];
    }
    magnitude = sqrtf(magnitude);

    for (int i = 0; i < nx; i++)
    {
      normalized[i + j * nx] /= magnitude;
    }
  }
}

__global__ void calculate_result(int nx, int ny, float *result, float *normalized)
{
  int i = threadIdx.x;
  int j = threadIdx.x;

  float sum = 0.0;
  for (int k = 0; k < nx; ++k)
    sum += normalized[k + i * nx] * normalized[k + j * nx];

  result[i + j * ny] = sum;
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
  set_result_to_zero(ny, result);
  // ll normalized(nx * ny, 0);
  float *normalized = (float *)malloc(nx * ny * sizeof(float));
  normalize_rows(ny, nx, data, normalized);

  // Allocate memory & copy data to GPU
  float *normalizedGPU = NULL;
  CHECK(cudaMalloc((void **)&normalizedGPU, nx * ny * sizeof(float)));
  float *resultGPU = NULL;
  CHECK(cudaMalloc((void **)&resultGPU, ny * ny * sizeof(float)));
  CHECK(cudaMemcpy(normalizedGPU, normalized, nx * ny * sizeof(float), cudaMemcpyHostToDevice));

  calculate_result<<<ny, ny>>>(nx, ny, resultGPU, normalizedGPU);
  CHECK(cudaGetLastError());

  CHECK(cudaMemcpy(result, resultGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
  CHECK(cudaFree(normalizedGPU));
  CHECK(cudaFree(resultGPU));
}
