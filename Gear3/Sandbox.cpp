#include "Sandbox.h"
#include "PitchCurve.h"
#include "ToothMaker.h"
#include "Gear.h"

#include <QPainter>
#include <QPainterPath>
#include <qmath.h>
#include <QDebug>

void draw(QPainter &pa, PitchCurve *crv)
{
  QPainterPath pp;

  for(int i=0;i<crv->getPointCount();i++)
  {
    const PitchCurve::Point &pt = crv->getPoint(i);
    QPointF p = pt.pos.toPointF();
    pp.moveTo(p);
    pp.lineTo(p+(pt.right*20).toPointF());
  }
  pa.setPen(Qt::magenta);
  pa.drawPath(pp);

  pp = QPainterPath();
  pp.moveTo(crv->getPoint(0).pos.toPointF());
  for(int i=1;i<crv->getPointCount();i++)
    pp.lineTo(crv->getPoint(i).pos.toPointF());
  pp.closeSubpath();
  pa.setPen(Qt::blue);
  pa.setBrush(Qt::NoBrush);
  pa.drawPath(pp);

}


//=========================================================================================================


struct Sandbox::Imp {
  double parameter;
  // PitchCurve *curve;
  Gear *m_gear1, *m_gear2;

};

Sandbox::Sandbox()
  : m_imp(new Imp)
{
  m_imp->parameter = 0.0;
  // testPitchCurve();
  /*
  double r=200.0, e = 0.6;
  PitchCurve *crv = m_imp->curve = new PitchCurve(EllipseFunction(r,e),200);
  double err = crv->getLength() - EllipseFunction::computePerimeter(r,e);
  qDebug() << "Err = " << err; // 0.137783; -0.0300675
  */
  /*
  double r=600.0, corner = 20;
  PitchCurve *crv = m_imp->curve = new PitchCurve(SquareFunction(r,corner),200);
  double err = crv->getLength() - SquareFunction::computePerimeter(r,corner);
  qDebug() << "Err = " << err; // 0.137783; -0.0300675
  */
  // PitchCurve *crv = m_imp->curve = new PitchCurve(SpiralFunction(100,200),200);

  Gear *gear = new Gear(new PitchCurve(SpiralFunction(200,500)));
  QVector<QVector2D> pts;
  SimpleToothMaker::Params params;
  params.toothLength = 0.5*gear->getCurve()->getLength()/50;
  params.toothHeight = 10;
  SimpleToothMaker().makeTeeth(pts, gear->getCurve(), params);
  gear->setBodyPath(pts);

  m_imp->m_gear1 = gear;
  
  
  gear = new Gear(new PitchCurve(SpiralFunction(200,500)));
  pts.clear();
  SimpleToothMaker().makeTeeth(pts, gear->getCurve(), params);
  gear->setBodyPath(pts);
  m_imp->m_gear2 = gear;
  gear->setPosition(-(200+500),0);
  
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
  m_imp->m_gear1->setAngle(M_PI - getParameter()*0.01);

  GearLink(m_imp->m_gear1, m_imp->m_gear2).update();

  
  m_imp->m_gear1->draw(pa);
  m_imp->m_gear2->draw(pa);
  */






  double r0 = 60, r1 = 200, d = 20;
  double q = r1-3*d;

  double r2 = r1 + 50;

  double angle1 = getParameter();

  pa.save();
  pa.translate(r1,r1);
  pa.rotate(angle1);

  pa.setPen(QPen(Qt::black,2));
  pa.setBrush(Qt::cyan);
  pa.drawEllipse(QPointF(0,0),r2,r2);

  pa.setBrush(Qt::blue);  
  {
  double r = r1-3*d-2;
  QRectF box(-r,-r,2*r,2*r);
  QRectF box2(-r,-r-r*sqrt(2.0),2*r,2*r);

  QPainterPath pp;
  pp.arcMoveTo(box, 135);
  pp.arcTo(box, 135, 270); 
  pp.arcTo(box2, -45, -90); 
  pp.closeSubpath();
  pa.drawPath(pp);
  }

  pa.setBrush(Qt::red);
  pa.drawEllipse(QPointF(0,-r1), 15,15);

  pa.setBrush(Qt::white);
  pa.drawEllipse(QPointF(0,0), d,d);

  pa.restore();

  double angle2 = 0;
  if(-90<=angle1 && angle1<0) 
  {
    double phi = M_PI*angle1/180.0;
    QPointF p(r1+r1*sin(phi), r1-r1*cos(phi));
    double psi = atan2(p.y(),p.x());
    angle2 = 180.0*psi/M_PI;
  }


  pa.setPen(QPen(Qt::black,2));
  QPainterPath pp;
  pp.moveTo(d, r0);

  QPointF e0(1,0), e1(0,1);
  double phi = 0.0;
  QPointF p;

  for(int i=0;i<4;i++)
  {
    pp.lineTo(d*e0+r1*e1);
    p = e0*(2*d)+r1*e1;
    pp.arcTo(p.x()-d,p.y()-d,2*d,2*d,180+phi,180);
    p = r1*(e0+e1);
    pp.arcTo(p.x()-q,p.y()-q,2*q,2*q,180+phi,-90);
    p = e1*(2*d)+r1*e0;
    pp.arcTo(p.x()-d,p.y()-d,2*d,2*d,-90+phi,180);
    p = r0*e0;
    pp.arcTo(p.x()-d,p.y()-d,2*d,2*d,-90+phi,-180);
    phi += 90;
    QPointF tmp = e1; e1 = e0; e0 = -tmp; 
  }

  pp.closeSubpath();
  pp.addEllipse(QPointF(0,0), d,d);

  pa.setPen(QPen(Qt::blue, 2));
  pa.setBrush(Qt::yellow);
  pa.save();
  pa.rotate(angle2);
  pa.drawPath(pp);
  pa.restore();

  
  /*

  draw(pa, m_imp->curve);

  QVector<QVector2D> pts;
  SimpleToothMaker::Params params;
  params.toothLength = 50;
  params.toothHeight = 30;

  SimpleToothMaker().makeTeeth(pts, m_imp->curve, params);

  QPainterPath pp;
  pp.moveTo(pts[0].toPointF());
  for(int i=1;i<pts.count();i++) pp.lineTo(pts[i].toPointF());
  pa.setPen(Qt::black);
  pa.setBrush(Qt::yellow);
  pa.drawPath(pp);
  */


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

  /*
  int n1 = 23, n2 = 43;
  double r1 = 300.0, r2 = r1*n2/n1;
  double alpha = 20 * M_PI/180.0;
  double csAlpha = cos(alpha);

  CircularToothMaker ctm;
  CircularToothMaker::Params params;

  params.toothHeight = 2 * r1 * (1-csAlpha);
  params.pitchRadius = r1;
  params.toothCount = n1;

  QVector<QVector2D> pts1, pts2;  
  ctm.makeTeeth(pts1, params);

  params.toothHeight = 2 * r2 * (1-csAlpha);
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
    pa.translate(i==0 ? r1 : -r2, 0);
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
  pa.drawLine(-r2, 0, r1, 0);
  pa.drawLine(-sin(alpha)*300, -300, sin(alpha)*300,300);
  */


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

