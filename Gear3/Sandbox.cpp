#include "Sandbox.h"
#include "PitchCurve.h"
#include <QPainter>
#include <QPainterPath>
#include <qmath.h>
#include <QDebug>


class CircularToothMaker {

public:

  struct Params {
    double pitchRadius;
    int toothCount;
    double toothHeight;
  };


  CircularToothMaker() {}

  void makeTeeth(QVector<QVector2D> &pts, const Params &params);

private:
  QVector2D evolute(double phi, double radius) const { 
      double cs = cos(phi), sn = sin(phi);
      return radius * QVector2D(cs + phi*sn, sn - phi*cs);
  }

  void makeToothSide(QList<QPair<double, double> > &qs, double pitchRadius, double toothHeight);
  void completeTooth(QList<QPair<double, double> > &qs, int toothCount);
};

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
    Q_ASSERT(err<0.02);
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


struct Sandbox::Imp {
  double parameter;
  PitchCurve *curve;

};

Sandbox::Sandbox()
  : m_imp(new Imp)
{
  m_imp->parameter = 0.0;
  double r=200.0, e = 0.6;

  PitchCurve *crv = m_imp->curve = new PitchCurve(EllipseFunction(r,e),200);
  double err = crv->getLength() - EllipseFunction::computePerimeter(r,e);
  qDebug() << "Err = " << err; // 0.137783; -0.0300675
}

Sandbox::~Sandbox()
{
  delete m_imp;
}

void Sandbox::setParameter(double v)
{
  m_imp->parameter = v;
}

double Sandbox::getParameter() const
{
  return m_imp->parameter;
}

void Sandbox::paint(QPainter &pa)
{
  pa.setPen(Qt::gray);
  double r = 500;
  pa.drawLine(-r,0,r,0);
  pa.drawLine(0,-r,0,r);


  /*

  PitchCurve *crv = m_imp->curve;
  QPainterPath pp;
  pp.moveTo(crv->getPoint(0).pos.toPointF());
  for(int i=1;i<crv->getPointCount();i++)
    pp.lineTo(crv->getPoint(i).pos.toPointF());
  pp.closeSubpath();
  pa.setPen(Qt::black);
  pa.drawPath(pp);

  pp = QPainterPath();
  for(int i=0;i<crv->getPointCount();i++)
  {
    pp.moveTo(crv->getPoint(i).pos.toPointF());
    pp.lineTo((crv->getPoint(i).pos+30*crv->getPoint(i).right).toPointF());
    //pp.addEllipse(crv->getPoint(i).pos.toPointF(), 30,30);
  }
  pa.setPen(Qt::cyan);
  pa.drawPath(pp);
  */
  int n1 = 23, n2 = 23;
  double r1 = 300.0, r2 = r1*n2/n1;

  CircularToothMaker ctm;
  CircularToothMaker::Params params;

  params.toothHeight = 40;
  params.pitchRadius = r1;
  params.toothCount = n1;

  QVector<QVector2D> pts1, pts2;  
  ctm.makeTeeth(pts1, params);

  params.pitchRadius = r2;
  params.toothCount = n2;
  ctm.makeTeeth(pts2, params);



  QPainterPath pps[2];
  for(int i=0;i<2;i++)
  {
    QVector<QVector2D> &pts = i==0 ? pts1 : pts2;
    pps[i].moveTo(pts[0].toPointF());
    for(int j=0;j<pts.count();j++) pps[i].lineTo(pts[j].toPointF());
    pps[i].closeSubpath();
  }


  double rotation = getParameter() * 1.0; 
  for(int i=0;i<2;i++)
  {
    pa.save();
    pa.translate(-(r1+r2)*i,0);
    pa.rotate(i==0 ? rotation : -rotation*n1/n2);
    pa.setPen(Qt::black);
    pa.setBrush(Qt::cyan);
    pa.drawPath(pps[i]);

    pa.setBrush(Qt::NoBrush);
    pa.setPen(QPen(Qt::black, 0, Qt::DotLine));
    double r = i==0 ? r1 : r2;
    pa.drawEllipse(QPointF(0,0), r,r);
    pa.restore();

  }


  pa.setPen(Qt::black);
  pa.drawLine(-(r1+r2),0,0,0);
  /*

  double r0 = 300;
  double r1 = 330;
  double r2 = 360;

  r = 200;
  pa.setPen(Qt::black);
  pa.drawEllipse(QPointF(0,0), r0,r0);
  pa.setPen(Qt::gray);
  pa.drawEllipse(QPointF(0,0), r1,r1);
  pa.drawEllipse(QPointF(0,0), r2,r2);
  
  QVector<QVector2D> model;
  double phi = 0.0;
  double dphi = 0.1;
  double oldr = 0.0;
  double psi = 0.0;
  for(;;)
  {
    double cs = cos(phi), sn = sin(phi);
    QVector2D p = r0*QVector2D(cs,sn) + r0*phi*QVector2D(sn,-cs);
    r = p.length();
    if(oldr < r1 && r1<=r)
    {
      double t = (r1-oldr)/(r-oldr);
      phi = phi - dphi + t*dphi;
      cs = cos(phi), sn = sin(phi);
      p = r0*QVector2D(cs,sn) + r0*phi*QVector2D(sn,-cs);
      r = p.length();
      double err = fabs(r-r1);
      Q_ASSERT(err<0.5);      
      r = r1;
      psi = atan2(p.y(),p.x());
    }
    else if(oldr < r2 && r2<=r)
    {
      double t = (r2-oldr)/(r-oldr);
      phi = phi - dphi + t*dphi;
      cs = cos(phi), sn = sin(phi);
      p = r0*QVector2D(cs,sn) + r0*phi*QVector2D(sn,-cs);
      r = p.length();
      double err = fabs(r-r2);
      Q_ASSERT(err<0.5);
      model.append(p);
      break;
    }
    model.append(p);
    oldr = r;
    phi += dphi;
  }

  int toothCount = 17;

  double rotation = getParameter() * 1.0; 
  QPainterPath pp;
  pp.moveTo(model[0].toPointF());
  for(int i=1;i<model.count();i++)
    pp.lineTo(model[i].toPointF());

  pa.setPen(Qt::red);
  for(int i=0;i<toothCount;i++)
  {
    pa.save();
    pa.rotate(-90.0/toothCount + rotation);
    pa.rotate(360.0*i/toothCount - psi*180/M_PI);
    pa.drawPath(pp);
    pa.restore();

    pa.save();
    pa.rotate( 90.0/toothCount + rotation);
    pa.rotate(360.0*i/toothCount);
    pa.scale(1,-1);
    pa.rotate( - psi*180/M_PI - 0/toothCount);
    pa.drawPath(pp);
    pa.restore();

  }
  
  
 

  pa.save();
  pa.translate(-2*r1,0);

  pa.setPen(Qt::black);
  pa.drawEllipse(QPointF(0,0), r0,r0);
  pa.setPen(Qt::gray);
  pa.drawEllipse(QPointF(0,0), r1,r1);
  pa.drawEllipse(QPointF(0,0), r2,r2);
  

  pa.setPen(Qt::red);
  for(int i=0;i<toothCount;i++)
  {
    pa.save();
    pa.rotate(-90.0/toothCount - rotation);
    pa.rotate(360.0*i/toothCount - psi*180/M_PI);
    pa.drawPath(pp);
    pa.restore();

    pa.save();
    pa.rotate( 90.0/toothCount - rotation);
    pa.rotate(360.0*i/toothCount);
    pa.scale(1,-1);
    pa.rotate( - psi*180/M_PI - 0/toothCount);
    pa.drawPath(pp);
    pa.restore();

  }

  pa.restore();
  */



}

