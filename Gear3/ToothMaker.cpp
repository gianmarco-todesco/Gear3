#include "ToothMaker.h"
#include "PitchCurve.h"


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
  int m = 4;
  double ds = params.toothLength/m;
  double h = params.toothHeight * 0.5;
  double s = 0;
  int i = 0;
  while(s<curve->getLength())
  {
    double y = 0.0;
    if(i!=0 && i!=m) y = i<m ? h : -h; 
    pts.append(curve->getPosFromS(s,y));
    s += ds;
    i = (i+1)%(2*m);
  }
}

