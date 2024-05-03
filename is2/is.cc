// #include <limits>
// #include <cstdio>
using namespace std;

// constexpr double infty = numeric_limits<double>::infinity();
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

double color(int x, int y, int c, int nx, const float *data)
{
    return data[c + C * x + C * nx * y];
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
                    double outer[C];
                    double inner[C];

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
                }
            }
        }
    }

    return result;
}