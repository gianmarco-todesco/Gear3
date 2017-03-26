#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"



class MalteseCrossPage : public Page, public Pannable {
  Gear *m_gear1, *m_gear2;
public:
  MalteseCrossPage() : Page("malteseCross") { 
  }

  ~MalteseCrossPage() {  }

  void draw(QPainter &pa);
} titlePage;


void MalteseCrossPage::draw(QPainter &pa)
{
  

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

}


