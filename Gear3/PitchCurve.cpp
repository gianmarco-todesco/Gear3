#include "PitchCurve.h"
#include <qmath.h>



double normalizePeriodicValue(double value, double period, int *q)
{
  double originalValue = value;
  *q=0;
  if(value>=period) { *q = (int)floor(value/period); value -= *q * period; }
  else if(value<0) { *q = -(int)ceil(-value/period); value -= *q * period; }
  if(value>=period) { while(value>=period) { value-=period; *q += 1; } }
  else if(value<0) { while(value<0) { value += period; *q -= 1; } }
  Q_ASSERT(0<=value && value<period);
  Q_ASSERT(fabs(value + (*q)*period - originalValue) < 1.0e-8);
  return value;
}


//=============================================================================


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

double SquareFunction::computePerimeter(double e, double r)
{
  return 2*M_PI*r + (e-r)*8;
}


QVector2D SpiralFunction::operator()(double t) const
{ 
  double phi = 2*M_PI*t; 
  return (m_r0 * (1-t) + m_r1 * t) * QVector2D(cos(phi), sin(phi));
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
  m_isOpen = (f(0)-f(1)).length()>2.0;
  double dt = m_isOpen ? 1.0/(double)(n-1) : 1.0/(double)n;
  for(int i=0;i<n;i++)
  {
    Point &pt = m_pts[i];
    pt.pos = f(i*dt);
    pt.phi = atan2(pt.pos.y(),pt.pos.x()); 
    if(pt.phi<0)pt.phi += 2*M_PI;
    Q_ASSERT(0<=pt.phi && pt.phi<=2*M_PI);
    Q_ASSERT(i==0 || m_pts[i-1].phi < pt.phi);    
    pt.r = pt.pos.length();
  }
  for(int i=0;i<n;i++)
  {
    QVector2D v = m_isOpen ?  (m_pts[qMin(n-1,i+1)].pos - m_pts[qMax(0,i-1)].pos) : (m_pts[(i+1)%n].pos - m_pts[(i+n-1)%n].pos);
    m_pts[i].right = QVector2D(v.y(),-v.x()).normalized();
  }
  m_pts[0].s = 0.0;
  for(int i=1;i<n;i++)
  {
    m_pts[i].s = m_pts[i-1].s + spiralArcLength(m_pts[i-1].phi, m_pts[i-1].r, m_pts[i].phi, m_pts[i].r);
  }
  m_length = m_pts.back().s;
  if(!m_isOpen)
    m_length += spiralArcLength(m_pts.back().phi, m_pts.back().r, 2*M_PI, m_pts[0].r);
}


PitchCurve::~PitchCurve()
{
}


double PitchCurve::getPhiFromS(double s) const
{
  int a,b,q; double t;
  if(m_isOpen) {
    getIndexFromS(s, a,b, t);
    if(a<0) return m_pts[b].phi;
    else if(b<0) return m_pts[a].phi;
    else return m_pts[a].phi*(1-t) + m_pts[b].phi*t;
  }
  else
  {
    getIndexFromS(s, a,b, q, t);
    return m_pts[a].phi*(1-t) + (b==0 ? 2*M_PI : m_pts[b].phi)*t + q*2*M_PI;
  }
}

double PitchCurve::getSFromPhi(double phi) const
{
  int a,b,q; double t;
  if(m_isOpen) {
    getIndexFromS(phi, a,b, t);
    if(a<0) return m_pts[b].s;
    else if(b<0) return m_pts[a].s;
    else return m_pts[a].s*(1-t) + m_pts[b].s*t;
  }
  else
  {
    getIndexFromS(phi, a,b, q, t);
    return m_pts[a].s*(1-t) + (b==0 ? m_length : m_pts[b].s)*t + q*m_length;
  }
}


double PitchCurve::getRFromS(double s) const
{
  int a,b,q; double t;
  if(m_isOpen) {
    getIndexFromS(s, a,b, t);
    if(a<0) return m_pts[b].r;
    else if(b<0) return m_pts[a].r;
    else return m_pts[a].r*(1-t) + m_pts[b].r*t;
  }
  else
  {
    getIndexFromS(s, a,b, q, t);
    return m_pts[a].r*(1-t) + m_pts[b].r*t;
  }
}

double PitchCurve::getRFromPhi(double phi) const
{
  int a,b,q; double t;
  if(m_isOpen) {
    getIndexFromS(phi, a,b, t);
    if(a<0) return m_pts[b].r;
    else if(b<0) return m_pts[a].r;
    else return m_pts[a].r*(1-t) + m_pts[b].r*t;
  }
  else
  {
    getIndexFromS(phi, a,b, q, t);
    return m_pts[a].r*(1-t) + m_pts[b].r*t;
  }
}


QVector2D PitchCurve::getPosFromS(double s, double y) const
{
  int a,b,q; double t;
  QVector2D pos,rt;
  if(m_isOpen) {
    getIndexFromS(s, a,b, t);
    if(a<0) { pos = m_pts[b].pos; rt = m_pts[b].right; }
    else if(b<0) { pos = m_pts[a].pos; rt = m_pts[a].right; }
    else { pos = m_pts[a].pos*(1-t) + m_pts[b].pos*t; rt =  m_pts[a].right*(1-t) + m_pts[b].right*t; }
  }
  else
  {
    getIndexFromS(s, a,b, q, t);
    pos = m_pts[a].pos*(1-t) + m_pts[b].pos*t;
    rt = m_pts[a].right*(1-t) + m_pts[b].right*t;
  }
  if(y != 0.0) pos += y * rt.normalized();
  return pos;
}

QVector2D PitchCurve::getPosFromPhi(double phi) const
{
  int a,b,q; double t;
  if(m_isOpen) {
    getIndexFromPhi(phi, a,b, t);
    if(a<0) return m_pts[b].pos;
    else if(b<0) return m_pts[a].pos;
    else return m_pts[a].pos*(1-t) + m_pts[b].pos*t;
  }
  else
  {
    getIndexFromPhi(phi, a,b, q, t);
    return m_pts[a].pos*(1-t) + m_pts[b].pos*t;
  }
}

// for close curves
void PitchCurve::getIndexFromS(double s, int &a, int &b, int &q, double &t) const
{
  s = normalizePeriodicValue(s, m_length, &q);
  if(s>=m_pts.back().s) 
  {
    a = m_pts.count()-1;
    b = 0;
    t = (s - m_pts.back().s)/(m_length - m_pts.back().s);
  }
  else
  {
    a=0;
    b=m_pts.count()-1;
    Q_ASSERT(m_pts[a].s<=s && s<m_pts[b].s);
    while(b-a>1) 
    {
      int c = (a+b)/2;
      if(m_pts[c].s<=s)a=c; else b=c;
    }
    Q_ASSERT(m_pts[a].s<=s && s<m_pts[b].s);
    t = (s-m_pts[a].s)/(m_pts[b].s-m_pts[a].s);
  }
}

void PitchCurve::getIndexFromPhi(double phi, int &a, int &b, int &q, double &t) const
{
  phi = normalizePeriodicValue(phi, 2*M_PI, &q);
  if(phi>=m_pts.back().phi) 
  {
    a = m_pts.count()-1;
    b = 0;
    t = (phi - m_pts.back().phi)/(2*M_PI - m_pts.back().phi);
  }
  else
  {
    a=0;
    b=m_pts.count()-1;
    Q_ASSERT(m_pts[a].phi<=phi && phi<m_pts[b].phi);
    while(b-a>1) 
    {
      int c = (a+b)/2;
      if(m_pts[c].phi<=phi)a=c; else b=c;
    }
    Q_ASSERT(m_pts[a].phi<=phi && phi<m_pts[b].phi);
    t = (phi-m_pts[a].phi)/(m_pts[b].phi-m_pts[a].phi);
  }
}


// for open curves
void PitchCurve::getIndexFromS(double s, int &a, int &b, double &t) const
{
  if(s<=0.0) { a=-1; b=0; t = 0; }
  else if(s>=m_length) { a=m_pts.count()-1; b=-1; t=0; }
  a=0;
  b=m_pts.count()-1;
  Q_ASSERT(m_pts[a].s<=s && s<m_pts[b].s);
  while(b-a>1) 
  {
    int c = (a+b)/2;
    if(m_pts[c].s<=s)a=c; else b=c;
  }
  Q_ASSERT(m_pts[a].s<=s && s<m_pts[b].s);
  t = (s-m_pts[a].s)/(m_pts[b].s-m_pts[a].s);  
}

void PitchCurve::getIndexFromPhi(double phi, int &a, int &b, double &t) const
{
  if(phi<0.0) { a=-1; b=0; t = 0; }
  else if(phi>=m_pts.back().phi) { a=m_pts.count()-1; b=-1; t=0; }
  a=0;
  b=m_pts.count()-1;
  Q_ASSERT(m_pts[a].phi<=phi && phi<m_pts[b].phi);
  while(b-a>1) 
  {
    int c = (a+b)/2;
    if(m_pts[c].phi<=phi)a=c; else b=c;
  }
  Q_ASSERT(m_pts[a].phi<=phi && phi<m_pts[b].phi);
  t = (phi-m_pts[a].phi)/(m_pts[b].phi-m_pts[a].phi);
}

//=============================================================================

void testPitchCurve()
{
  PitchCurve *crv;
  double err;
  double radius = 100.0;
  crv = new PitchCurve(EllipseFunction(radius,0.0), 100);
  err = fabs(crv->getLength() - 2*M_PI*radius);
  Q_ASSERT(err<1.0e-2);

  int n = 10000;
  for(int i=0;i<n;i++)
  {
    for(int q= - 3; q<=3; q++)
    {
      double s = 2*M_PI*radius*(q + (double)i/(double)n);
      double phi = 2*M_PI*(q + (double)i/(double)n);
      double s1,phi1;
      phi1 = crv->getPhiFromS(s);
      err = fabs(phi1 - phi);
      Q_ASSERT(err<1.0e-2);

      s1 = crv->getSFromPhi(s);
      err = fabs(s1 - s);
      Q_ASSERT(err<1.0e-2);

      QVector2D p = crv->getPosFromPhi(phi);
      QVector2D p1 = radius * QVector2D(cos(phi), sin(phi));
      err = (p-p1).length();
      Q_ASSERT(err<1.0);


    }
  }

}
