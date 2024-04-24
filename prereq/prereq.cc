struct Result
{
    float avg[3];
};

/*
This is the function you need to implement. Quick reference:
- x coordinates: 0 <= x < nx
- y coordinates: 0 <= y < ny
- horizontal position: 0 <= x0 < x1 <= nx
- vertical position: 0 <= y0 < y1 <= ny
- color components: 0 <= c < 3
- input: data[c + 3 * x + 3 * nx * y]
- output: avg[c]
*/
Result calculate(int ny, int nx, const float *data, int y0, int x0, int y1, int x1)
{
    Result result{{0.0f, 0.0f, 0.0f}};
    for (int c = 0; c < 3; c++)
    {
        int count = (x1 - x0) * (y1 - y0);
        double avg = 0.0;
        for (int x = x0; x < x1; x++)
        {
            for (int y = y0; y < y1; y++)
            {
                double value = (double)data[c + 3 * x + 3 * nx * y];
                avg += value / count;
            }
        }
        result.avg[c] = (float)avg;
    }
    return result;
}
