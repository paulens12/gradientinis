#include <iostream>
#include <mgl2/mgl.h>
#include <tuple>

#define ZMIN (-0.005)
#define ZMAX 0.125

struct triangle
{
    mglPoint points[3];
    double values[3];
};

double dfndx(double x, double y)
{
    if(x == 0 && y == 0)
        return -0.01;
    return (2*x*y+y*y-y)/8;
}

double dfndy(double x, double y)
{
    if(x == 0 && y == 0)
        return -0.01;
    return (2*x*y+x*x-x)/8;
}

double fn(double x, double y)
{
    //return sin(5*x)*cos(5*y)/5;
    return (x+y-1)*x*y/8;
}

double fnGreic(double x, double y, double g)
{
    return fn(x - g * dfndx(x, y), y - g * dfndy(x, y));
}

const double fib = (sqrt(5) - 1) / 2; // NOLINT(cert-err58-cpp)

double auks_pjuvis(double x0, double y0, double min, double max, double e, std::vector<double> &points)
{
    double l = max - min;
    double g1 = max - fib * l;
    double g2 = min + fib * l;

    if(l < e)
    {
        if(fnGreic(x0, y0, g1) < fnGreic(x0, y0, g2))
            return g1;
        return g2;
    }

    if(fnGreic(x0, y0, g2) < fnGreic(x0, y0, g1))
    {
        points.push_back(g1);
        return auks_pjuvis(x0, y0, g1, max, e, points);
    }

    points.push_back(g2);
    return auks_pjuvis(x0, y0, min, g2, e, points);
}

mglPoint gradientinis(double x0, double y0, double e, double g, std::vector<std::pair<mglPoint, mglPoint> > &points)
{
    double stp;
    double x=x0, y=y0;
    do
    {
        double dx = dfndx(x, y);
        double dy = dfndy(x, y);
        points.emplace_back(std::pair<mglPoint, mglPoint>(mglPoint(x, y, ZMAX), mglPoint(dx, dy)));
        x -= g * dx;
        y -= g * dy;
        stp = sqrt(dx*dx+dy*dy);
    } while(stp > e);

    return mglPoint(x, y, ZMAX);
}

//mglPoint greiciausias(double x0, double y0, double e, std::vector<mglPoint> &points)
//{
//    double dx = dfndx(x0, y0);
//    double dy = dfndy(x0, y0);
//    double maxX = dx >= 0 ? x0 / dx : (x0 - 1) / dx;
//    double maxY = dy >= 0 ? y0 / dy : (y0 - 1) / dy;
//
//
//    std::vector<double> gamos;
//    double g = auks_pjuvis(x0, y0, 0, std::min(maxX, maxY), e, gamos);
//
//    for(auto & _g : gamos)
//    {
//        points.emplace_back(mglPoint(x0 - _g * dfndx(x0, y0), y0 - _g * dfndy(x0, y0)));
//    }
//
//    return mglPoint(x0 - g * dfndx(x0, y0), y0 - g * dfndy(x0, y0));
//}

mglPoint greiciausias(double x0, double y0, double e, std::vector<std::pair<mglPoint, mglPoint> > &points)
{
    double stp;
    double x=x0, y=y0;
    do
    {
        double dx = dfndx(x, y);
        double dy = dfndy(x, y);
        if(dx == 0 && dy == 0)
            break;

        double maxX = dx >= 0 ? x / dx : (x - 1) / dx;
        double maxY = dy >= 0 ? y / dy : (y - 1) / dy;

        std::vector<double> gamos;
        double g = auks_pjuvis(x, y, 0, std::min(maxX, maxY), e, gamos);
        points.emplace_back(std::pair<mglPoint, mglPoint>(mglPoint(x, y, ZMAX), mglPoint(dx, dy)));

        x -= g * dx;
        y -= g * dy;
        stp = sqrt(dx*dx+dy*dy);
    }while(stp > e);

    return mglPoint(x, y, ZMAX);
}

void sort(int &i0, int &i1, int &i2, const double f[])
{
    int tmp;
    if(f[i1] > f[i0])
    {
        tmp = i0;
        i0 = i1;
        i1 = tmp;
    }
    if(f[i2] > f[i0])
    {
        tmp = i2;
        i2 = i0;
        i0 = tmp;
    }
    if(f[i2] > f[i1])
    {
        tmp = i2;
        i2 = i1;
        i1 = tmp;
    }
}

void resize(const mglPoint& pntStays, mglPoint &pnt1, mglPoint &pnt2, double k, double k1)
{
    pnt1.x = k * pntStays.x + k1 * pnt1.x;
    pnt1.y = k * pntStays.y + k1 * pnt1.y;

    pnt2.x = k * pntStays.x + k1 * pnt2.x;
    pnt2.y = k * pntStays.y + k1 * pnt2.y;
}

mglPoint simpleksas(double x0, double y0, double a, double e, double fe, int m, double k, int maxStp, std::vector<triangle> &points)
{
    int steps = 0;
    double k1 = 1 - k;
    // n = 2
    double d1 = a*(sqrt(3) + 1)/(2*sqrt(2));
    double d2 = a*(sqrt(3) - 1)/(2*sqrt(2));

    mglPoint pnts[3];
    int tmp;
    pnts[0] = mglPoint(x0, y0, ZMAX);
    pnts[1] = mglPoint(x0 + d2, y0 + d1, ZMAX);
    pnts[2] = mglPoint(x0 + d1, y0 + d2, ZMAX);
    int i0 = 0, i1 = 1, i2 = 2;
    double f[3] = //funkciju reiksmes
    {
        fn(pnts[0].x, pnts[0].y),
        fn(pnts[1].x, pnts[1].y),
        fn(pnts[2].x, pnts[2].y)
    };
    int h[3] = {0, 0, 0}; //kiek laiko taskas neismestas

    points.emplace_back(triangle{.points = {pnts[0], pnts[1], pnts[2]}, .values = {f[0], f[1], f[2]}});

    sort(i0, i1, i2, f); // i0 - didziausias, i2 - maziausias

    while(a > e) {
        steps++;
        h[0]++;
        h[1]++;
        h[2]++;

        if(h[0] > m || h[1] > m || h[2] > m)
        {
            a *= k;
            resize(pnts[i2], pnts[i0], pnts[i1], k, k1);

            f[i0] = fn(pnts[i0].x, pnts[i0].y);
            f[i1] = fn(pnts[i1].x, pnts[i1].y);
            h[i0] = 0;
            h[i1] = 0;
            h[i2] = 0;
            sort(i0, i1, i2, f);
            points.emplace_back(triangle{.points = {pnts[0], pnts[1], pnts[2]}, .values = {f[0], f[1], f[2]}});
        }

        mglPoint pnew = mglPoint(pnts[i1].x + pnts[i2].x - pnts[i0].x, pnts[i1].y + pnts[i2].y - pnts[i0].y, ZMAX);
        double fnew = fn(pnew.x, pnew.y);
        if (fnew > f[i1]) {
            pnew.x = pnts[i0].x + pnts[i2].x - pnts[i1].x;
            pnew.y = pnts[i0].y + pnts[i2].y - pnts[i1].y;
            fnew = fn(pnew.x, pnew.y);
            tmp = i0;
            i0 = i1;
            i1 = tmp;
        }
        if (fnew > f[i1]) {
            // abu variantai netiko, mazinam trikampi
            a *= k;
            resize(pnts[i2], pnts[i1], pnts[i0], k, k1);

            f[i0] = fn(pnts[i0].x, pnts[i0].y);
            f[i1] = fn(pnts[i1].x, pnts[i1].y);
            h[i0] = 0;
            h[i1] = 0;
            sort(i0, i1, i2, f);
            points.emplace_back(triangle{.points = {pnts[0], pnts[1], pnts[2]}, .values = {f[0], f[1], f[2]}});
        } else {
            //kazkuris tiko, verciam trikampi...
            pnts[i0] = pnew;
            h[i0] = 0;
            f[i0] = fnew;
            sort(i0, i1, i2, f);
            points.emplace_back(triangle{.points = {pnts[0], pnts[1], pnts[2]}, .values = {f[0], f[1], f[2]}});
        }

        if((f[i0]-f[i2]) < fe || steps >= maxStp)
        {
            return pnts[i2];
        }
    }
    return pnts[i2];
}


int sample(mglGraph *gr, double step)
{
    gr->Aspect(NAN, 1);
    gr->SetTicks('x', 1);
    gr->SetTicks('y', 1);
    gr->Axis();
    //gr->Rotate(50,60);
    int steps = ceil(1.0 / step);

    mglData xd(steps), yd(steps), zd(steps, steps);

    double min = 0;

    for(int i=0;i<steps;i++)
    {
        double x = i * step;
        xd.a[i] = x;
        for(int j=0;j<steps;j++)
        {
            double y = j * step;
            yd.a[j] = y;

            double xx = fn(x, y);
            if(xx < min)
                min = xx;
            zd.SetVal(xx, i, j);
        }
    }

    gr->Surf(xd, yd, zd);
    std::cout << min << std::endl;
    return 0;
}

void printpoints(const std::vector<std::pair<mglPoint, mglPoint> >& points, mglGraph &gr)
{
    for(auto & point : points)
    {
        gr.Mark(point.first, "bk");
        std::cout << point.first.x << '\t' << point.first.y << '\t'
                  << point.second.x << '\t' << point.second.y << std::endl;
    }
}

void printpoints(const std::vector<mglPoint> &points, mglGraph &gr)
{
    for(auto & point : points)
    {
        gr.Mark(point, "bk");
        std::cout << point.x << '\t' << point.y << std::endl;
    }
}

void printpoints(const std::vector<triangle> &triangles, mglGraph &gr)
{
    for(auto & triangle : triangles)
    {
        gr.Line(triangle.points[0], triangle.points[1], "bk");
        gr.Line(triangle.points[1], triangle.points[2], "bk");
        gr.Line(triangle.points[0], triangle.points[2], "bk");
        std::cout << "(" << triangle.points[0].x << ", " << triangle.points[0].y << ")\t" << triangle.values[0] << "\t" <<
                     "(" << triangle.points[1].x << ", " << triangle.points[1].y << ")\t" << triangle.values[1] << "\t" <<
                     "(" << triangle.points[2].x << ", " << triangle.points[2].y << ")\t" << triangle.values[2] << "\t" << std::endl;
    }
}

int main() {
    std::cout << "taskas (0,0): f(X) = " <<  fn(0,0) << "; gradientas = (" << dfndx(0,0) << "," << dfndy(0,0) << ")." << std::endl;
    std::cout << "taskas (1,1): f(X) = " <<  fn(1,1) << "; gradientas = (" << dfndx(1,1) << "," << dfndy(1,1) << ")." << std::endl;
    std::cout << "taskas (0,0.8): f(X) = " <<  fn(0,0.8) << "; gradientas = (" << dfndx(0,0.8) << "," << dfndy(0,0.8) << ")." << std::endl;
    std::vector<std::pair<mglPoint, mglPoint> > pointsWithDs;
    mglGraph gr;
    gr.SetRange('x', 0, 1);
    gr.SetRange('y', 0, 1);
    gr.SetRange('z', ZMIN, ZMAX);
    gr.SetRange('c', ZMIN - 0.001, 0.001);
    gr.SetSize(2400, 1000);
    std::cout << "GRADIENTINIS NUSILEIDIMAS" << std::endl << std::endl;

    gr.SubPlot(3, 1, 0, "");
    gr.Title("(0, 0)");
    sample(&gr, 0.001);
    mglPoint ats = gradientinis(0, 0, 0.00001, 2, pointsWithDs);
    std::cout << std::endl << "taskas (0, 0): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(pointsWithDs, gr);
    pointsWithDs.clear();

    gr.SubPlot(3, 1, 1, "");
    gr.Title("(1, 1)");
    sample(&gr, 0.001);
    ats = gradientinis(1, 1, 0.00001, 2, pointsWithDs);
    std::cout << std::endl << "taskas (1, 1): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(pointsWithDs, gr);
    pointsWithDs.clear();

    gr.SubPlot(3, 1, 2, "");
    gr.Title("(0, 0.8)");
    sample(&gr, 0.001);
    ats = gradientinis(0, 0.8, 0.00001, 2, pointsWithDs);
    std::cout << std::endl << "taskas (0, 0.8): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(pointsWithDs, gr);
    pointsWithDs.clear();

    gr.WritePNG("gradientinis.png");

    std::vector<triangle> triangles;
    gr.SetSize(2400, 1000);
    std::cout << "GREICIAUSIAS NUSILEIDIMAS" << std::endl << std::endl;

    gr.SubPlot(3, 1, 0, "");
    gr.Title("(0, 0)");
    sample(&gr, 0.001);
    ats = greiciausias(0, 0, 0.00001, pointsWithDs);
    std::cout << std::endl << "taskas (0, 0): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(pointsWithDs, gr);
    pointsWithDs.clear();

    gr.SubPlot(3, 1, 1, "");
    gr.Title("(1, 1)");
    sample(&gr, 0.001);
    ats = greiciausias(1, 1, 0.00001, pointsWithDs);
    std::cout << std::endl << "taskas (1, 1): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(pointsWithDs, gr);
    pointsWithDs.clear();

    gr.SubPlot(3, 1, 2, "");
    gr.Title("(0, 0.8)");
    sample(&gr, 0.001);
    ats = greiciausias(0, 0.8, 0.00001, pointsWithDs);
    std::cout << std::endl << "taskas (0, 0.8): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(pointsWithDs, gr);
    pointsWithDs.clear();

    gr.WritePNG("greiciausias.png");

    std::cout << "DEFORMUOJAMASIS SIMPLEKSAS" << std::endl << std::endl;
    gr.SetSize(2400, 1000);
    gr.SubPlot(3, 1, 0, "");
    gr.Title("(0, 0)");
    sample(&gr, 0.001);
    ats = simpleksas(0, 0, 0.3, 0.001, 0.00001, 6, 0.5, 1000, triangles);
    std::cout << std::endl << "taskas (0, 0): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(triangles, gr);
    triangles.clear();

    gr.SubPlot(3, 1, 1, "");
    gr.Title("(1, 1)");
    sample(&gr, 0.001);
    ats = simpleksas(1, 1, 0.3, 0.001, 0.00001, 6, 0.5, 1000, triangles);
    std::cout << std::endl << "taskas (1, 1): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(triangles, gr);
    triangles.clear();

    gr.SubPlot(3, 1, 2, "");
    gr.Title("(0, 0.8)");
    sample(&gr, 0.001);
    ats = simpleksas(0, 0.8, 0.3, 0.001, 0.00001, 6, 0.5, 1000, triangles);
    std::cout << std::endl << "taskas (0, 0.8): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(triangles, gr);
    triangles.clear();
    gr.WritePNG("simpleksas.png");

    std::cout << "done!" << std::endl;
    return 0;
}