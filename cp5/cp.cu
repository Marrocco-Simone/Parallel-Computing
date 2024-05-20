#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cuda_runtime.h>
#define STEPS 16
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

__global__ void normalize_rows(int ny, int nx, float *data, float *normalized)
{
  int j = blockIdx.x * STEPS + threadIdx.x;
  if (j >= ny)
    return;

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

__global__ void calculate_result(int nx, int ny, float *result, float *normalized)
{
  int js = blockIdx.x * STEPS;
  int is = threadIdx.x * STEPS;
  int je = min(js + STEPS, ny);
  int ie = min(is + STEPS, ny);

  if (is >= ny || js >= ny || ie < js)
    return;

  for (int j = js; j < je; j++)
    for (int i = is; i < ie; i++)
    {
      float sum = 0.0;
      for (int k = 0; k < nx; ++k)
        sum += normalized[k + i * nx] * normalized[k + j * nx];

      result[i + j * ny] = sum;
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
  int nnx = nx + STEPS - nx % STEPS;
  int nny = ny + STEPS - ny % STEPS;

  float *dataGPU = NULL;
  CHECK(cudaMalloc((void **)&dataGPU, nx * ny * sizeof(float)));
  CHECK(cudaMemcpy(dataGPU, data, nx * ny * sizeof(float), cudaMemcpyHostToDevice));

  float *normalizedGPU = NULL;
  CHECK(cudaMalloc((void **)&normalizedGPU, nnx * nny * sizeof(float)));
  normalize_rows<<<nny / STEPS, STEPS>>>(ny, nx, dataGPU, normalizedGPU);

  float *resultGPU = NULL;
  CHECK(cudaMalloc((void **)&resultGPU, ny * ny * sizeof(float)));

  calculate_result<<<nny / STEPS, nny / STEPS>>>(nx, ny, resultGPU, normalizedGPU);
  CHECK(cudaGetLastError());

  CHECK(cudaMemcpy(result, resultGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
  CHECK(cudaFree(dataGPU));
  CHECK(cudaFree(normalizedGPU));
  CHECK(cudaFree(resultGPU));
}
