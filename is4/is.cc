#include <vector>
#include <chrono>
#include <cstdio>
#define T 4
typedef double double4_t __attribute__((vector_size(T * sizeof(double))));
using namespace std;

constexpr double infty = 1.7e308;
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

double color(int x, int y, int c, int nx, const float *data)
{
    return data[c + C * x + C * nx * y];
}

/** O(n^2)*/
void calculate_total_sum(int ny, int nx, const float *data, double4_t &total_sum)
{
#pragma omp parallel for
    for (int c = 0; c < C; c++)
        for (int x = 0; x < nx; x++)
            for (int y = 0; y < ny; y++)
                total_sum[c] += color(x, y, c, nx, data);
}

/** O(n^2) */
void calculate_sum_from_zero(int ny, int nx, const float *data, vector<double4_t> &sum_from_zero)
{
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
        {
            double4_t point_block = {color(x, y, 0, nx, data), color(x, y, 1, nx, data), color(x, y, 2, nx, data), 0.0};
            double4_t prev_left_block = sum_from_zero[id4_t(x, y + 1, nx + 1)];
            double4_t prev_up_block = sum_from_zero[id4_t(x + 1, y, nx + 1)];
            double4_t prev_left_up_block = sum_from_zero[id4_t(x, y, nx + 1)];

            sum_from_zero[id4_t(x + 1, y + 1, nx + 1)] = (prev_left_block + prev_up_block - prev_left_up_block + point_block);
        }
}

/** O(1) - access all combinations of x1 / y1 / x0-1 / y0-1 */
void calculate_avg_in_color(int x0, int x1, int y0, int y1, int nx, vector<double4_t> const &sum_from_zero, double4_t &inner)
{
    double4_t point_block = sum_from_zero[id4_t(x1, y1, nx)];
    double4_t prev_left_block = sum_from_zero[id4_t(x0, y1, nx)];
    double4_t prev_up_block = sum_from_zero[id4_t(x1, y0, nx)];
    double4_t prev_left_up_block = sum_from_zero[id4_t(x0, y0, nx)];

    inner = prev_left_up_block + point_block - prev_up_block - prev_left_block;
}

/** O(1) */
void set_result(int x0, int x1, int y0, int y1, int nx, int ny, const std::vector<double4_t> &sum_from_zero, double4_t &total_sum, Result &result)
{
    int in_points = (y1 - y0) * (x1 - x0);
    int out_points = nx * ny - in_points;

    double4_t inner;
    calculate_avg_in_color(x0, x1, y0, y1, nx + 1, sum_from_zero, inner);
    double4_t outer = total_sum - inner;

    result.x0 = x0;
    result.x1 = x1;
    result.y0 = y0;
    result.y1 = y1;
    for (int c = 0; c < C; c++)
    {
        result.inner[c] = inner[c] / in_points;
        result.outer[c] = outer[c] / out_points;
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
    double4_t total_sum = {0.0};
    calculate_total_sum(ny, nx, data, total_sum);
    vector<double4_t> sum_from_zero((nx + 1) * (ny + 1));
    calculate_sum_from_zero(ny, nx, data, sum_from_zero);

    vector<double> best_solutions(nx * ny, 0.0);
    int full_points = nx * ny;

    auto t2 = high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic, 2)
    for (int ydim = 1; ydim <= ny; ydim++)
        for (int xdim = 1; xdim <= nx; xdim++)
        {
            int i = (ydim - 1) * nx + xdim - 1;
            double max_err = 0.0;

            int in_points = ydim * xdim;
            int out_points = full_points - in_points;
            double inv_in_points = 1.0 / in_points;
            double inv_out_points = 1.0 / out_points;

            for (int y0 = 0; y0 < ny - ydim + 1; y0++)
            {
                for (int x0 = 0; x0 < nx - xdim + 1; x0++)
                {
                    int y1 = y0 + ydim;
                    int x1 = x0 + xdim;

                    double4_t inner;
                    calculate_avg_in_color(x0, x1, y0, y1, nx + 1, sum_from_zero, inner);
                    double4_t outer = total_sum - inner;

                    double out_err = (outer[0] * outer[0] + outer[1] * outer[1] + outer[2] * outer[2]);
                    double in_err = (inner[0] * inner[0] + inner[1] * inner[1] + inner[2] * inner[2]);
                    double inv_sq_err = in_err * inv_in_points + out_err * inv_out_points;

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
    double max_err = 0.0;
    int in_points = ydim * xdim;
    int out_points = full_points - in_points;
    double inv_in_points = 1.0 / in_points;
    double inv_out_points = 1.0 / out_points;
    for (int y0 = 0; y0 < ny - ydim + 1; y0++)
    {
        for (int x0 = 0; x0 < nx - xdim + 1; x0++)
        {
            int y1 = y0 + ydim;
            int x1 = x0 + xdim;

            int in_points = (y1 - y0) * (x1 - x0);
            int out_points = full_points - in_points;

            double4_t inner;
            calculate_avg_in_color(x0, x1, y0, y1, nx + 1, sum_from_zero, inner);
            double4_t outer = total_sum - inner;

            double out_err = (outer[0] * outer[0] + outer[1] * outer[1] + outer[2] * outer[2]);
            double in_err = (inner[0] * inner[0] + inner[1] * inner[1] + inner[2] * inner[2]);
            double inv_sq_err = in_err * inv_in_points + out_err * inv_out_points;

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