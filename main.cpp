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

double fnb(double x, double y, double z)
{
    return x+y+z-1;
}

double fnba(double x, double y, double z)
{
    return abs(fnb(x, y, z));
}

double fnb2(double x, double y, double z)
{
    double a = fnb(x, y, z);
    return a*a;
}

double dfndx(double x, double y, double z, double r)
{
    return 2*(x+y+z-1)*r - y*z/8;
}

double dfndy(double x, double y, double z, double r)
{
    return 2*(x+y+z-1)*r - x*z/8;
}

double dfndz(double x, double y, double z, double r)
{
    return 2*(x+y+z-1)*r - x*y/8;
}

double fn(double x, double y, double z)
{
    return -x*y*z/8;
}

double fnB(double x, double y, double z, double r)
{
    return fn(x, y, z) + fnb2(x, y, z) * r;
}

mglPoint gradientinis(double x0, double y0, double z0, double r, double e, double g, std::vector<std::pair<mglPoint, mglPoint> > &points)
{
    double stp;
    double x=x0, y=y0, z=z0;
    do
    {
        double dx = dfndx(x, y, z, r);
        double dy = dfndy(x, y, z, r);
        double dz = dfndz(x, y, z, r);
        points.emplace_back(std::pair<mglPoint, mglPoint>(mglPoint(x, y, ZMAX), mglPoint(dx, dy)));
        x -= g * dx;
        y -= g * dy;
        z -= g * dz;
        stp = sqrt(dx*dx+dy*dy+dz*dz);
    } while(stp > e);

    return mglPoint(x, y, ZMAX);
}

int sample(mglGraph *gr, double step)
{
    gr->Aspect(NAN, 1, 1);
    gr->SetTicks('x', 1);
    gr->SetTicks('y', 1);
    gr->SetTicks('z', 1);
    gr->Rotate(50,60);
    gr->Axis();
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

            double xx = 1 - x - y;
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

int main() {
    std::cout << "taskas (0,0): f(X) = " <<  fn(0,0,0) << "; gradientas = (" << dfndx(0,0,0,0) << "," << dfndy(0,0,0,0) << "," << dfndz(0,0,0,0) << ")." << std::endl;
    std::cout << "taskas (1,1): f(X) = " <<  fn(1,1,1) << "; gradientas = (" << dfndx(1,1,1,0) << "," << dfndy(1,1,1,0) << "," << dfndz(1,1,1,0) << ")." << std::endl;
    std::cout << "taskas (0,0.8): f(X) = " <<  fn(0.6,0,0.8) << "; gradientas = (" << dfndx(0.6,0,0.8,0) << "," << dfndy(0.6,0,0.8,0) << "," << dfndz(0.6,0,0.8,0) << ")." << std::endl;
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
    mglPoint ats = gradientinis(0, 0, 0, 100, 0.00001, 2, pointsWithDs);
    std::cout << std::endl << "taskas (0, 0): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(pointsWithDs, gr);
    pointsWithDs.clear();

    gr.SubPlot(3, 1, 1, "");
    gr.Title("(1, 1)");
    sample(&gr, 0.001);
    ats = gradientinis(1, 1, 1, 100, 0.00001, 2, pointsWithDs);
    std::cout << std::endl << "taskas (1, 1): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(pointsWithDs, gr);
    pointsWithDs.clear();

    gr.SubPlot(3, 1, 2, "");
    gr.Title("(0, 0.8)");
    sample(&gr, 0.001);
    ats = gradientinis(0.6, 0, 0.8, 100, 0.00001, 2, pointsWithDs);
    std::cout << std::endl << "taskas (0, 0.8): ats (" << ats.x << ", " << ats.y << ")" << std::endl;
    printpoints(pointsWithDs, gr);
    pointsWithDs.clear();

    gr.WritePNG("gradientinis.png");

    std::cout << "done!" << std::endl;
    return 0;
}