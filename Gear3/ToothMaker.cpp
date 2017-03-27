#include "ToothMaker.h"
#include "PitchCurve.h"
#include "ClipperWrapper.h"
#include "PolygonSweeper.h"


#include <qmath.h>


void CircularToothMaker::makeToothSide(QList<QPair<double, double> > &qs, double pitchRadius, double toothHeight)
{
  double r0 = pitchRadius - toothHeight*0.5;
  double r1 = pitchRadius;
  double r2 = r0 + toothHeight;
  struct Node { double t, r; QVector2D pos; };
  QList<Node> nodes;
  double t = 0.0;
  for(;;)
  {
    Node node;
    node.t = t;
    node.pos = evolute(t, r0);
    node.r = node.pos.length();
    nodes.append(node);
    if(node.r>r2) break;
    t += 0.02;
  }

  int j = 0;
  int m = 7;
  for(int i=0;i<m;i++)
  {
    double r = r0 + (r2-r0)*i/(double)(m-1);
    if(r>r2)r=r2;
    while(nodes[j+1].r<=r) j++;
    double w = (r-nodes[j].r)/(nodes[j+1].r-nodes[j].r);
    double t = nodes[j].t *(1-w) + nodes[j+1].t * w;
    QVector2D pos = evolute(t, r0);
    double err = fabs(r-pos.length());
    Q_ASSERT(err<0.05);
    double phi = atan2(pos.y(),pos.x());
    if(phi<0.0)phi += 2*M_PI;
    qs.append(qMakePair(phi, r));
  }
  double dphi = qs[m/2].first;
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
  makeToothSide(qs, params.pitchRadius, params.toothHeight);
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

void SquareToothMaker::makeTeeth(QVector<QVector2D> &pts, const PitchCurve *curve, const Params &params)
{
  double ds = curve->getLength() / params.toothCount;

  double y0 = -params.toothHeight * 0.5;
  double y1 = y0 + params.toothHeight;

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



void MagicToothMaker::makeTeeth(QVector<QVector2D> &pts, const PitchCurve *crv)
{
  QList<QPointF> shape;
  shape << QPointF(20,-20) << QPointF(5,20) << QPointF(-5,20) << QPointF(-20,-20);


  ClipperWrapper cw;
 
  int n = 10;
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
      matrix.translate(-s + s0,0);
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
    bound.append((pt.pos + pt.right*20).toPointF());
  }

  cw.antiSubtract(bound);    
  outline.clear();

  cw.getResult(outline);
  for(int i=0;i<outline[0].count();i++) pts.append(QVector2D(outline[0][i]));
}

