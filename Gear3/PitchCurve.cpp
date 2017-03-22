#include "PitchCurve.h"
#include <qmath.h>


QVector2D EllipseFunction::operator()(double t) const 
{
  double phi = t*M_PI*2;
  double r = m_r/(1.0 + m_e*cos(phi));
  return r*QVector2D(cos(phi), sin(phi));
}

double EllipseFunction::computePerimeter(double r, double e)
{
  // see https://www.mathsisfun.com/geometry/ellipse-perimeter.html
  // verificato con http://www.had2know.com/academics/ellipse-area-perimeter-calculator.html
  double a = r / (1.0 - e*e); // semi-major axis length
  double f = a - r/(1+e); // 
  double b = sqrt(a*a-f*f); // semi-minor axis length

  double h = pow((a-b)/(a+b),2.0);
  double s = 1 + 0.25*h;
  double q = 1.0/4.0;
  double he = h;
  for(int n=2;n<10;n++)
  {
    q *= pow((2.0*n-3.0)/(2.0*n),2);
    he*=h;
    s += q*he;
  }
  double length = s* M_PI*(a+b);
  return length;
}

//=============================================================================

SquareFunction::SquareFunction(double length, double r) : m_r(r) 
{
  m_l = m_r + (length - 2*M_PI*m_r)/8;
  m_a = m_l - m_r;
  m_b = sqrt(m_a*m_a+m_l*m_l);  
}

QVector2D SquareFunction::operator()(double t) const 
{
    double phi = t*M_PI*2;

    double cs = cos(phi), sn = sin(phi);
    double absCs = fabs(cs), absSn = fabs(sn); 
    double r = m_l/qMax(absCs, absSn);
    if(r > m_b)
    {
      QVector2D q(cs>0 ? m_a : -m_a, sn>0 ? m_a : -m_a);
      double cu = cs*q.x() + sn*q.y(), cc = q.lengthSquared();
      double dsc = cu*cu-cc+m_r*m_r;
      Q_ASSERT(dsc>0.0);
      r = cu+sqrt(dsc);
    }
    return QVector2D(cos(phi), sin(phi))*r;
  }



//=============================================================================

double spiralArcLength(double phi0, double r0, double phi1, double r1)
{
  /*
  // r=p+theta*q, theta=[0,T]
  Q_ASSERT(phi1>phi0);
  double T = phi1-phi0, p = r0, q = (r1-r0)/T;

  double ds = 0.5*(q*T*T + (sqrt(p*p+q*q)+p)*T);

  QVector2D p0 = r0*QVector2D(cos(phi0), sin(phi0));
  QVector2D p1 = r1*QVector2D(cos(phi1), sin(phi1));
  double ds1 = (p0-p1).length();
  double err = fabs(ds1-ds);
  Q_ASSERT(err<3.0);
  return ds1;
  */

  double s = 0.0;
  int m = 100;
  QVector2D oldp = r0 * QVector2D(cos(phi0),sin(phi0));
  for(int i = 1; i<=m; i++) 
  {
    double t = (double)i/(double)m;
    double phi = (1-t)*phi0+t*phi1;
    QVector2D p = (r0*(1-t)+r1*t) * QVector2D(cos(phi),sin(phi));
    s += (p-oldp).length();
    oldp=p;
  }
  return s;
}

//=============================================================================

PitchCurve::PitchCurve(const CurveFunction &f, int n)
  : m_pts(n)
  , m_length(0)
{ 
  for(int i=0;i<n;i++)
  {
    Point &pt = m_pts[i];
    pt.pos = f((double)i/(double)n);
    pt.phi = atan2(pt.pos.y(),pt.pos.x()); 
    if(pt.phi<0)pt.phi += 2*M_PI;
    Q_ASSERT(0<=pt.phi && pt.phi<2*M_PI);
    Q_ASSERT(i==0 || m_pts[i-1].phi < pt.phi);    
    pt.r = pt.pos.length();
  }
  for(int i=0;i<n;i++)
  {
    QVector2D v = m_pts[(i+1)%n].pos - m_pts[(i+n-1)%n].pos;
    m_pts[i].right = QVector2D(v.y(),-v.x()).normalized();
  }
  m_pts[0].s = 0.0;
  for(int i=1;i<n;i++)
  {
    m_pts[i].s = m_pts[i-1].s + spiralArcLength(m_pts[i-1].phi, m_pts[i-1].r, m_pts[i].phi, m_pts[i].r);
  }
  m_length = m_pts.back().s + spiralArcLength(m_pts.back().phi, m_pts.back().r, 2*M_PI, m_pts[0].r);



}


PitchCurve::~PitchCurve()
{
}


