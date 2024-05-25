#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cuda_runtime.h>
#include <chrono>
#define STEP 12
#define PADDING 144
using namespace std;

/*
  ! test command
  ./grading test-plain; ./grading test-asan; ./grading test-memcheck-initcheck --remote; ./grading benchmark-cache --remote benchmarks/4b.txt
*/

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

__global__ void normalize_rows(int ny, int nx, int nny, float *data, float *normalizedTransposed)
{
  int j = blockIdx.x * PADDING + threadIdx.x;
  if (j >= ny)
    return;

  float mean = 0.0;
  float magnitude = 0.0;

  for (int k = 0; k < nx; k++)
    mean += data[k + j * nx] / nx;

  for (int k = 0; k < nx; k++)
  {
    normalizedTransposed[j + k * nny] = (data[k + j * nx] - mean);
    magnitude += normalizedTransposed[j + k * nny] * normalizedTransposed[j + k * nny];
  }
  magnitude = sqrtf(magnitude);

  for (int k = 0; k < nx; k++)
  {
    normalizedTransposed[j + k * nny] /= magnitude;
  }
}

__global__ void calculate_result(int nx, int nny, float *result, float *normalizedTransposed)
{
  int is = (threadIdx.x + blockIdx.x * blockDim.x) * STEP;
  int js = (threadIdx.y + blockIdx.y * blockDim.y) * STEP;
  // printf("blockIdx.x: %d, blockIdx.y: %d, threadIdx.x: %d, threadIdx.y: %d, js: %d, is: %d\n", blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y, js, is);

  if (is < js)
    return;

  float sums[STEP][STEP] = {0.0};
  for (int k = 0; k < nx; ++k)
  {
    for (int j = js; j < js + STEP; j++)
      for (int i = is; i < is + STEP; i++)
        sums[i - is][j - js] += normalizedTransposed[i + k * nny] * normalizedTransposed[j + k * nny];
  }

  for (int j = js; j < js + STEP; j++)
    for (int i = is; i < is + STEP; i++)
      result[i + j * nny] = sums[i - is][j - js];
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
  using std::chrono::duration;
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;
  using std::chrono::milliseconds;

  auto t1 = high_resolution_clock::now();

  int nnx = nx + PADDING - nx % PADDING;
  int nny = ny + PADDING - ny % PADDING;

  float *dataGPU = NULL;
  CHECK(cudaMalloc((void **)&dataGPU, nx * ny * sizeof(float)));
  CHECK(cudaMemcpy(dataGPU, data, nx * ny * sizeof(float), cudaMemcpyHostToDevice));

  float *normalizedTransposedGPU = NULL;
  CHECK(cudaMalloc((void **)&normalizedTransposedGPU, nnx * nny * sizeof(float)));
  CHECK(cudaMemset(normalizedTransposedGPU, 0, nnx * nny * sizeof(float)));
  normalize_rows<<<nny / PADDING, PADDING>>>(ny, nx, nny, dataGPU, normalizedTransposedGPU);

  float *resultGPU = NULL;
  CHECK(cudaMalloc((void **)&resultGPU, nny * nny * sizeof(float)));
  CHECK(cudaMemset(resultGPU, 0, nny * nny * sizeof(float)));

  auto t2 = high_resolution_clock::now();

  dim3 dimBlock(STEP, STEP);
  dim3 dimGrid(nny / PADDING, nny / PADDING);
  calculate_result<<<dimGrid, dimBlock>>>(nx, nny, resultGPU, normalizedTransposedGPU);
  CHECK(cudaGetLastError());

  float *result_padded = (float *)malloc(nny * nny * sizeof(float));
  CHECK(cudaMemcpy(result_padded, resultGPU, nny * nny * sizeof(float), cudaMemcpyDeviceToHost));
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < ny; i++)
      result[i + j * ny] = result_padded[i + j * nny];

  free(result_padded);
  CHECK(cudaFree(dataGPU));
  CHECK(cudaFree(normalizedTransposedGPU));
  CHECK(cudaFree(resultGPU));

  auto t3 = high_resolution_clock::now();

  printf("Initialization: %ld ms\n", duration_cast<milliseconds>(t2 - t1).count());
  printf("Main loop: %ld ms\n", duration_cast<milliseconds>(t3 - t2).count());
}
