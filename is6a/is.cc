#include <vector>
#include <chrono>
#include <cstdio>
#define T 8
typedef float float8_t __attribute__((vector_size(T * sizeof(float))));
using namespace std;

constexpr int C = 3;

struct Result
{
    int y0;
    int x0;
    int y1;
    int x1;
    float outer[C];
    float inner[C];
};

int id4_t(int x, int y, int nx)
{
    return x + nx * y;
}

float color(int x, int y, int nx, const float *data)
{
    return data[C * x + C * nx * y];
}

/** O(n^2)*/
float calculate_total_sum(int ny, int nx, const float *data)
{
    float total_sum = 0.0;
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
            total_sum += color(x, y, nx, data);
    return total_sum;
}

/** O(n^2) */
void calculate_sum_from_zero(int ny, int nx, const float *data, vector<float> &sum_from_zero)
{
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
        {
            float point_block = color(x, y, nx, data);
            float prev_left_block = sum_from_zero[id4_t(x, y + 1, nx + 1)];
            float prev_up_block = sum_from_zero[id4_t(x + 1, y, nx + 1)];
            float prev_left_up_block = sum_from_zero[id4_t(x, y, nx + 1)];

            sum_from_zero[id4_t(x + 1, y + 1, nx + 1)] = (prev_left_block + prev_up_block - prev_left_up_block + point_block);
        }
}

/** O(1) - access all combinations of x1 / y1 / x0-1 / y0-1 */
float calculate_avg_in_color(int x0, int x1, int y0, int y1, int nx, vector<float> const &sum_from_zero)
{
    float point_block = sum_from_zero[id4_t(x1, y1, nx)];
    float prev_left_block = sum_from_zero[id4_t(x0, y1, nx)];
    float prev_up_block = sum_from_zero[id4_t(x1, y0, nx)];
    float prev_left_up_block = sum_from_zero[id4_t(x0, y0, nx)];

    float inner = prev_left_up_block + point_block - prev_up_block - prev_left_block;
    return inner;
}

/** O(1) */
void set_result(int x0, int x1, int y0, int y1, int nx, int ny, const std::vector<float> &sum_from_zero, float total_sum, Result &result)
{
    int in_points = (y1 - y0) * (x1 - x0);
    int out_points = nx * ny - in_points;

    float inner = calculate_avg_in_color(x0, x1, y0, y1, nx + 1, sum_from_zero);
    float outer = total_sum - inner;

    result.x0 = x0;
    result.x1 = x1;
    result.y0 = y0;
    result.y1 = y1;
    for (int c = 0; c < C; c++)
    {
        result.inner[c] = inner / in_points;
        result.outer[c] = outer / out_points;
    }
}

/*
This is the function you need to implement. Quick reference:
- x coordinates: 0 <= x < nx
- y coordinates: 0 <= y < ny
- color components: 0 <= c < 3
- input: data[c + 3 * x + 3 * nx * y]
*/
Result segment(int ny, int nx, const float *data)
{
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    auto t1 = high_resolution_clock::now();

    Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};
    float total_sum = calculate_total_sum(ny, nx, data);
    vector<float> sum_from_zero((nx + 1) * (ny + 1));
    calculate_sum_from_zero(ny, nx, data, sum_from_zero);

    vector<float> best_solutions(nx * ny, 0.0);
    int full_points = nx * ny;

    auto t2 = high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic, 2)
    for (int ydim = 1; ydim <= ny; ydim++)
        for (int xdim = 1; xdim <= nx; xdim++)
        {
            int i = (ydim - 1) * nx + xdim - 1;
            float max_err = 0.0;

            int in_points = ydim * xdim;
            int out_points = full_points - in_points;
            float inv_in_points = 1.0 / in_points;
            float inv_out_points = 1.0 / out_points;

            for (int y0 = 0; y0 < ny - ydim + 1; y0++)
            {
                for (int x0 = 0; x0 < nx - xdim + 1; x0++)
                {
                    int y1 = y0 + ydim;
                    int x1 = x0 + xdim;

                    float inner = calculate_avg_in_color(x0, x1, y0, y1, nx + 1, sum_from_zero);
                    float outer = total_sum - inner;

                    float out_err = outer * outer;
                    float in_err = inner * inner;
                    float inv_sq_err = in_err * inv_in_points + out_err * inv_out_points;

                    max_err = max(max_err, inv_sq_err);
                }
            }

            best_solutions[i] = max_err;
        }

    auto t3 = high_resolution_clock::now();

    int min_i = 0;
    for (int i = 0; i < nx * ny; i++)
        if (best_solutions[i] > best_solutions[min_i])
            min_i = i;

    int ydim = min_i / nx + 1;
    int xdim = min_i % nx + 1;

    int y0f = 0;
    int x0f = 0;
    int y1f = 0;
    int x1f = 0;
    float max_err = 0.0;
    int in_points = ydim * xdim;
    int out_points = full_points - in_points;
    float inv_in_points = 1.0 / in_points;
    float inv_out_points = 1.0 / out_points;
    for (int y0 = 0; y0 < ny - ydim + 1; y0++)
    {
        for (int x0 = 0; x0 < nx - xdim + 1; x0++)
        {
            int y1 = y0 + ydim;
            int x1 = x0 + xdim;

            float inner = calculate_avg_in_color(x0, x1, y0, y1, nx + 1, sum_from_zero);
            float outer = total_sum - inner;

            float out_err = outer * outer;
            float in_err = inner * inner;
            float inv_sq_err = in_err * inv_in_points + out_err * inv_out_points;

            if (inv_sq_err > max_err)
            {
                max_err = inv_sq_err;
                y1f = y1;
                x1f = x1;
                y0f = y0;
                x0f = x0;
            }
        }
    }
    set_result(x0f, x1f, y0f, y1f, nx, ny, sum_from_zero, total_sum, result);

    auto t4 = high_resolution_clock::now();

    printf("Initialization: %ld ms\n", duration_cast<milliseconds>(t2 - t1).count());
    printf("Main loop: %ld ms\n", duration_cast<milliseconds>(t3 - t2).count());
    printf("Finalization: %ld ms\n", duration_cast<milliseconds>(t4 - t3).count());

    return result;
}