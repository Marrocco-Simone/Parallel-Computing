#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cuda_runtime.h>
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
  int i = blockIdx.x;
  int j = threadIdx.x;

  if (i >= ny || j >= ny)
    return;

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
  constexpr int steps = 8;
  int nnx = nx + steps - nx % steps;
  int nny = ny + steps - ny % steps;
  // float *normalized = (float *)malloc(nnx * nny * sizeof(float));
  float *normalized = (float *)calloc(nnx * nny, sizeof(float));
  normalize_rows(ny, nx, data, normalized);

  // Allocate memory & copy data to GPU
  float *normalizedGPU = NULL;
  CHECK(cudaMalloc((void **)&normalizedGPU, nnx * nny * sizeof(float)));
  float *resultGPU = NULL;
  CHECK(cudaMalloc((void **)&resultGPU, ny * ny * sizeof(float)));
  CHECK(cudaMemcpy(normalizedGPU, normalized, nnx * nny * sizeof(float), cudaMemcpyHostToDevice));

  calculate_result<<<nny, nny>>>(nx, ny, resultGPU, normalizedGPU);
  CHECK(cudaGetLastError());

  CHECK(cudaMemcpy(result, resultGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
  CHECK(cudaFree(normalizedGPU));
  CHECK(cudaFree(resultGPU));
}
