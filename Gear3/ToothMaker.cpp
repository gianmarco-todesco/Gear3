#include "ToothMaker.h"
#include "PitchCurve.h"
#include "ClipperWrapper.h"
#include "PolygonSweeper.h"


#include <qmath.h>


void CircularToothMaker::makeToothSide(QList<QPair<double, double> > &qs, double r0, double r1, double r2)
{
  double rr[2] = {r1,r2}, tt[2];
  int i=0;
  double oldr=r0, oldt=0.0;
  for(;;)
  {
    double t1 = oldt + 0.01;
    double r1 = evolute(t1,r0).length();
    if(r1>rr[i])
    {
      Q_ASSERT(oldr<=rr[i] && rr[i]<r1);
      double t = oldt + (t1-oldt) * (rr[i]-oldr) / (r1-oldr);
      tt[i] = t;
      if(i==1) break;
      i++;
    }
    oldr=r1;
    oldt=t1;
  }

  int m = 10;
  int ma = 2 + (int)(m*(r1-r0)/(r2-r0));
  int mb = 2 + (int)(m*(r2-r1)/(r2-r0));
  m = ma+mb-1;

  int j = 0;
  for(int i=0;i<m;i++)
  {
    double t;
    if(i<=ma)
      t = tt[0] * (double)i/(double)ma;
    else
      t = tt[0] + (tt[1]-tt[0]) * (double)(i-ma)/(double)(m-1-ma);
    QVector2D pos = evolute(t, r0);
    double r = pos.length();
    double phi = atan2(pos.y(),pos.x());
    if(phi<0.0)phi += 2*M_PI;
    qs.append(qMakePair(phi, r));
  }
  Q_ASSERT(fabs(qs[0].second-r0)<1.0e-2);
  Q_ASSERT(fabs(qs[ma].second-r1)<1.0e-2);
  Q_ASSERT(fabs(qs[m-1].second-r2)<1.0e-2);
  double dphi = qs[ma].first;
  for(int i=0;i<m;i++) qs[i].first -= dphi;
}

void CircularToothMaker::completeTooth(QList<QPair<double, double> > &qs, int toothCount)
{
  int m = qs.count();
  double dphi = M_PI/(toothCount*2);
  for(int i=0;i<m;i++) qs[i].first -= dphi;
  for(int i=m-1;i>=0;i--)
    qs.append(qMakePair(-qs[i].first, qs[i].second));
  
}

void CircularToothMaker::makeTeeth(QVector<QVector2D> &pts, const Params &params)
{
  QList<QPair<double, double> > qs;
  double r1 = params.pitchRadius;
  double r0 = r1 - params.toothDedendum, r2 = r1 + params.toothAddendum;
  makeToothSide(qs, r0,r1,r2);
  completeTooth(qs, params.toothCount);
  int n = params.toothCount;
  for(int i=0;i<n;i++)
  {
    double phiOff = 2*M_PI*i/n;
    for(int j=0;j<qs.count();j++)
    {
      double phi = phiOff + qs[j].first;
      double r = qs[j].second;
      pts.append(r*QVector2D(cos(phi), sin(phi)));
    }
  }
}


//=========================================================================================================

SimpleToothMaker::SimpleToothMaker()
{
}

void SimpleToothMaker::makeTeeth(QVector<QVector2D> &pts, const PitchCurve *curve, const Params &params)
{
  int m = 5, a = 3;
  int m4 = m*4;
  int n = params.toothCount * m4;
  double ds = curve->getLength() / n;
  double h = params.toothHeight * 0.5;
  QVector<double> ys(m4);
  int i1 = a, i2 = 2*m-a, i3 = 2*m+a, i4 = 4*m-a, i5 = 4*m-1;
  
  for(int i=0;i<=i1;i++) ys[i] = h;
  for(int i=i1+1;i<=i2;i++) { double t = (double)(i-i1)/(double)(i2-i1); ys[i] = h * (1-2*t); }
  for(int i=i2+1;i<=i3;i++) ys[i] = -h;
  for(int i=i3+1;i<=i4;i++) { double t = (double)(i-i3)/(double)(i4-i3); ys[i] = h * (-1+2*t); }
  for(int i=i4+1;i<=i5;i++) ys[i] = h;  
  int joff = ((int)(0.5+params.toothOffset * m4))%m4;

  for(int i=0;i<n;i++)
  {
    int j = (i+joff)%m4;
    pts.append(curve->getPosFromS(ds*i,ys[j]));
  }
}


//=========================================================================================================

void SquareToothMaker::makeTeeth(QVector<QVector2D> &pts, const PitchCurve *curve, const Params &params)
{
  double ds = curve->getLength() / params.toothCount;

  double y0 = -params.toothHeight * 0.5;
  double y1 = y0 + params.toothHeight * 0.9;

  int ma = 5, mb = 10;
  double a = 0.2;

  for(int i=0;i<params.toothCount;i++) 
  {
    double s0 = ds*i;
    double s1 = s0 + a * ds;
    double s2 = s0 + ds - a * ds;
    double s3 = s0 + ds;
    for(int i=0;i<ma;i++) {
      double t = (double)i/(double)(ma-1);
      pts.append(curve->getPosFromS(s0*(1-t)+s1*t, y1));
    }
    for(int i=0;i<mb;i++) {
      double t = (double)i/(double)(mb-1);
      pts.append(curve->getPosFromS(s1*(1-t)+s2*t, y0));
    }
    for(int i=0;i<ma-1;i++) {
      double t = (double)i/(double)(ma-1);
      pts.append(curve->getPosFromS(s2*(1-t)+s3*t, y1));
    }
  }    
}


//=========================================================================================================


void MagicToothMaker::makeTeeth(QVector<QVector2D> &pts, const PitchCurve *crv)
{
  QList<QPointF> shape;
  //shape << QPointF(20,-20) << QPointF(5,20) << QPointF(-5,20) << QPointF(-20,-20);

  int toothCount = 40;
  double slope = 0.3;
  double y0 = 10, y1 = 10;
  double addendum = 7.5;

  double ds = crv->getLength()/(2*toothCount);
  double x0 = -ds*0.5, x1 = ds*0.5;

  shape 
    << QPointF(x0-slope*y0,-y0) 
    << QPointF(x0+slope*y1,y1) 
    << QPointF(x1-slope*y1,y1) 
    << QPointF(x1+slope*y0,-y0); 

  ClipperWrapper cw;
 
  int n = toothCount;
  for(int k=0;k<n;k++)
  {
    double s0 = crv->getLength()*k/n;

    QList<QPolygonF> polygons;
    PolygonSweeper ps;
    ps.setShape(shape);

    int m = 40;
    for(int i=0;i<m;i++)
    {
      double t = (double)i/(double)(m-1);
      double s = s0 + 200*(-1+2*t);
      PitchCurve::Point pt = crv->getPointFromS(s);
      QMatrix matrix;
      matrix.translate(pt.pos.x(), pt.pos.y());
      matrix.rotate(180.0*atan2(pt.right.y(), pt.right.x())/M_PI + 90);
      matrix.translate(-s + s0 + 0.5*ds, 0);
      ps.add(matrix);
      QPolygonF polygon;
      for(int i=0;i<shape.count();i++) polygon.append(matrix.map(shape[i]));
      polygons.append(polygon);
    }

    QVector<QVector<QPointF > > strip;
    ps.getOutline(strip);
    cw.add(strip);
  }

  QVector<QVector<QPointF > > outline;
  cw.getResult(outline);

  QVector<QPointF> bound;
  for(int i=0;i<crv->getPointCount(); i++)
  {
    PitchCurve::Point pt = crv->getPoint(i);
    bound.append((pt.pos + pt.right*addendum).toPointF());
  }

  cw.antiSubtract(bound);    
  outline.clear();

  cw.getResult(outline);
  for(int i=0;i<outline[0].count();i++) pts.append(QVector2D(outline[0][i]));
}

