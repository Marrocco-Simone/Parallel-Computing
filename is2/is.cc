// #include <limits>
// TODO delete
#include <cstdio>
#include <vector>
using namespace std;

// constexpr double infty = numeric_limits<double>::infinity();
constexpr double infty = 1.7e308;
constexpr int C = 3;
// TODO delete
constexpr double epsilon = 1e-10;

struct Result
{
    int y0;
    int x0;
    int y1;
    int x1;
    float outer[C];
    float inner[C];
};

int id(int x, int y, int c, int nx)
{
    return c + C * x + C * nx * y;
}

double color(int x, int y, int c, int nx, const float *data)
{
    return data[id(x, y, c, nx)];
}

double error(int x, int y, int c, int nx, const float *data, double *area)
{
    return area[c] - color(x, y, c, nx, data);
}

bool is_inside(int x, int y, int x0, int x1, int y0, int y1)
{
    return x >= x0 && x < x1 && y >= y0 && y < y1;
}

double squared_error(int x0, int x1, int y0, int y1, int nx, int ny, const float *data, double *outer, double *inner)
{
    double sq_err = 0.0;
    for (int x = 0; x < nx; x++)
    {
        for (int y = 0; y < ny; y++)
        {
            double *area;
            if (is_inside(x, y, x0, x1, y0, y1))
                area = inner;
            else
                area = outer;
            for (int c = 0; c < C; c++)
            {
                double err = error(x, y, c, nx, data, area);
                sq_err += err * err;
            }
        }
    }
    return sq_err;
}

void avg_colors(int x0, int x1, int y0, int y1, int nx, int ny, const float *data, double *outer, double *inner)
{
    for (int c = 0; c < C; c++)
        outer[c] = inner[c] = 0;

    int inner_points = (x1 - x0) * (y1 - y0);
    int outer_points = nx * ny - inner_points;
    // printf("inner points: %d * %d = %d\n", x1 - x0, y1 - y0, inner_points);

    for (int x = 0; x < nx; x++)
    {
        for (int y = 0; y < ny; y++)
        {
            bool inside = is_inside(x, y, x0, x1, y0, y1);
            // printf("(%d,%d) is inside: %d\n", x, y, inside);
            for (int c = 0; c < C; c++)
            {
                // printf("color %d: %f,\t", c, color(x, y, c, nx, data));
                if (inside)
                    inner[c] += color(x, y, c, nx, data) / inner_points;
                else
                    outer[c] += color(x, y, c, nx, data) / outer_points;
            }
            // printf("\n");
        }
    }
}

/** O(n^2)*/
void calculate_total_avg(int ny, int nx, const float *data, double *total_avg)
{
    int tot = ny * nx;
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
            for (int c = 0; c < C; c++)
                total_avg[c] += color(x, y, c, nx, data) / tot;
}

/** O(n^2) */
void calculate_avg_from_zero(int ny, int nx, const float *data, vector<double> &avg_from_zero)
{
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
            for (int c = 0; c < C; c++)
            {
                printf("zero: %d %d %d\n", x, y, c);
                if (y == 0)
                {
                    if (x == 0)
                    {
                        // * set first element
                        avg_from_zero[id(x, y, c, nx)] = color(x, y, c, nx, data);
                    }
                    else
                    {
                        // * set first row
                        double prev_left_block = avg_from_zero[id(x - 1, y, c, nx)] * ((y + 1) * x);
                        avg_from_zero[id(x, y, c, nx)] = (prev_left_block + color(x, y, c, nx, data)) / ((x + 1) * (y + 1));
                    }
                }
                else if (x == 0)
                {
                    // * set first column
                    double prev_up_block = avg_from_zero[id(x, y - 1, c, nx)] * (y * (x + 1));
                    avg_from_zero[id(x, y, c, nx)] = (prev_up_block + color(x, y, c, nx, data)) / ((x + 1) * (y + 1));
                }
                else
                {
                    double prev_left_block = avg_from_zero[id(x - 1, y, c, nx)] * ((y + 1) * x);
                    double prev_up_block = avg_from_zero[id(x, y - 1, c, nx)] * (y * (x + 1));
                    double prev_left_up_block = avg_from_zero[id(x - 1, y - 1, c, nx)] * (y * x);

                    avg_from_zero[id(x, y, c, nx)] = (prev_left_block + prev_up_block - prev_left_up_block + color(x, y, c, nx, data)) / ((x + 1) * (y + 1));
                }
            }
}

/** O(1) */
void calculate_avg_in_color(int x0, int x1, int y0, int y1, int nx, vector<double> const &avg_from_zero, double *in)
{
    int block_points = (y1 - y0 + 1) * (x1 - x0 + 1);
    for (int c = 0; c < C; c++)
    {
        if (y0 == 0)
        {
            if (x0 == 0)
            {
                in[c] = avg_from_zero[id(x1, y1, c, nx)];
            }
            else
            {
                double point_block = avg_from_zero[id(x1, y1, c, nx)] * ((y1 + 1) * (x1 + 1));
                double prev_left_block = avg_from_zero[id(x0 - 1, y1, c, nx)] * ((y1 + 1) * x0);
                in[c] = (point_block - prev_left_block) / block_points;
            }
        }
        else if (x0 == 0)
        {
            double point_block = avg_from_zero[id(x1, y1, c, nx)] * ((y1 + 1) * (x1 + 1));
            double prev_up_block = avg_from_zero[id(x1, y0 - 1, c, nx)] * (y0 * (x1 + 1));
            in[c] = (point_block - prev_up_block) / block_points;
        }
        else if (y1 == 0)
        {
            // if (x1 == 0)
            // {
            //     // * set first element
            //     in[c] = avg_from_zero[id(x1, y1, c, nx)];
            // }
            // else
            {
                // * set first row
                printf("ids: %d\t%d\n", id(x0 - 1, y1, c, nx), id(x1, y1, c, nx));
                double x0y1_block = avg_from_zero[id(x0 - 1, y1, c, nx)] * x0;
                double x1y1_block = avg_from_zero[id(x1, y1, c, nx)] * (x1 + 1);
                in[c] = (x1y1_block - x0y1_block) / block_points;
            }
        }
        else if (x1 == 0)
        {
            // * set first column
            double x1y0_block = avg_from_zero[id(x1, y0 - 1, c, nx)] * y0;
            double x1y1_block = avg_from_zero[id(x1, y1, c, nx)] * (y1 + 1);
            in[c] = (x1y1_block - x1y0_block) / block_points;
        }
        else
        {
            double x0y0_block = avg_from_zero[id(x0 - 1, y0 - 1, c, nx)] * (y0 * x0);
            double x1y0_block = avg_from_zero[id(x1, y0 - 1, c, nx)] * (y0 * (x1 + 1));
            double x0y1_block = avg_from_zero[id(x0 - 1, y1, c, nx)] * ((y1 + 1) * x0);
            double x1y1_block = avg_from_zero[id(x1, y1, c, nx)] * ((y1 + 1) * (x1 + 1));
            in[c] = (x0y0_block + x1y1_block - x1y0_block - x0y1_block) / block_points;
        }
    }
}

/** O(1) */
void calculate_avg_out_color(int x0, int x1, int y0, int y1, int nx, int ny, double *in, double *total_avg, double *out)
{
    int i_points = (y1 - y0 + 1) * (x1 - x0 + 1);
    int full_points = (nx * ny);
    if (i_points == full_points)
    {
        return;
    }
    for (int c = 0; c < C; c++)
    {
        double in_block = in[c] * i_points;
        double full_block = total_avg[c] * full_points;
        out[c] = (full_block - in_block) / (full_points - i_points);
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
    Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};
    double min_err = infty;
    double total_avg[C] = {0.0};
    calculate_total_avg(ny, nx, data, total_avg);
    vector<double> avg_from_zero(C * nx * ny, 0.0);
    calculate_avg_from_zero(ny, nx, data, avg_from_zero);

    // TODO delete
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
            printf("%f %f %f\n", avg_from_zero[id(x, y, 0, nx)], avg_from_zero[id(x, y, 1, nx)], avg_from_zero[id(x, y, 2, nx)]);
    printf("\n");

    for (int xdim = 1; xdim <= nx; xdim++)
    {
        for (int ydim = 1; ydim <= ny; ydim++)
        {
            for (int x0 = 0; x0 <= nx - xdim; x0++)
            {
                for (int y0 = 0; y0 <= ny - ydim; y0++)
                {
                    int x1 = x0 + xdim;
                    int y1 = y0 + ydim;
                    double outer[C] = {0.0};
                    double inner[C] = {0.0};

                    // printf("x0: %d, x1: %d, y0: %d, y1: %d\n", x0, x1, y0, y1);
                    avg_colors(x0, x1, y0, y1, nx, ny, data, outer, inner);
                    double sq_err = squared_error(x0, x1, y0, y1, nx, ny, data, outer, inner);
                    // printf("inner: %f,%f,%f\n", inner[0], inner[1], inner[2]);
                    // printf("outer: %f,%f,%f\n", outer[0], outer[1], outer[2]);
                    // printf("error = %f\n", sq_err);
                    // printf("\n");
                    if (sq_err < min_err)
                    {
                        min_err = sq_err;
                        result.x0 = x0;
                        result.x1 = x1;
                        result.y0 = y0;
                        result.y1 = y1;
                        for (int c = 0; c < C; c++)
                        {
                            result.inner[c] = inner[c];
                            result.outer[c] = outer[c];
                        }
                    }

                    // TODO delete
                    printf("x0: %d, x1: %d, y0: %d, y1: %d\n", x0, x1, y0, y1);

                    double outer2[C] = {0.0};
                    double inner2[C] = {0.0};
                    calculate_avg_in_color(x0, x1 - 1, y0, y1 - 1, nx, avg_from_zero, inner2);
                    calculate_avg_out_color(x0, x1 - 1, y0, y1 - 1, nx, ny, inner2, total_avg, outer2);

                    for (int c = 0; c < C; c++)
                        printf("%d: %f-%f = %f-%f\n", c, inner[c], inner2[c], outer[c], outer2[c]);
                    for (int c = 0; c < C; c++)
                    {
                        if (inner2[c] - inner[c] > epsilon || outer2[c] - outer[c] > epsilon)
                        {
                            printf("ERROR in %d\n", c);
                            Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};
                            return result;
                        }
                    }
                    printf("ok\n\n");
                }
            }
        }
    }

    return result;
}