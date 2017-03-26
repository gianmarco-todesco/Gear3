#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"

#include <QLineF>

namespace {
Gear *makeSpiralGear(int toothCount)
{
  double radius = 200;
  Gear *gear = new Gear(new PitchCurve(SpiralFunction(300,500)));
  SimpleToothMaker ctm;
  SimpleToothMaker::Params params;
  params.toothHeight = 40;
  params.toothCount = toothCount;

  QVector<QVector2D> pts;
  ctm.makeTeeth(pts, gear->getCurve(), params);

  gear->setBodyPath(pts);
  return gear;
}

}


Gear *makeCircularGear(int toothLength, int toothCount, int flag = 0) {
  double radius = (toothLength * toothCount) / (2*M_PI);
  Gear *gear = new Gear(new PitchCurve(EllipseFunction(radius,0.0),200));

  
  QVector<QVector2D> pts;
  if(flag == 0)
  {
    SimpleToothMaker::Params params;
    params.toothHeight = toothLength*0.75*0.5;
    params.toothCount = toothCount;
    SimpleToothMaker().makeTeeth(pts, gear->getCurve(), params);
  }
  else if(flag==1)
  {
    SquareToothMaker::Params params;
    params.toothHeight = toothLength*0.75*0.5;
    params.toothCount = toothCount;
    SquareToothMaker().makeTeeth(pts, gear->getCurve(), params);
  }

  gear->setBodyPath(pts);
  
  return gear;

}


Gear *makeEllipticGear(int toothLength, int toothCount, double e) {
  double radius = (toothLength * toothCount) / (2*M_PI);
  Gear *gear = new Gear(new PitchCurve(EllipseFunction(radius,e),200));

  SimpleToothMaker::Params params;
  params.toothHeight = toothLength*0.75*0.5;
  params.toothCount = toothCount;
  
  QVector<QVector2D> pts;
  SimpleToothMaker().makeTeeth(pts, gear->getCurve(), params);

  gear->setBodyPath(pts);
  
  return gear;

}


Gear *makeSpiralGear(int toothCount, double r0, double r1, double toothOffset) {
  Gear *gear = new Gear(new PitchCurve(SpiralFunction(r0,r1),200));

  SimpleToothMaker::Params params;
  params.toothHeight = 20.0;
  params.toothCount = toothCount;
  params.toothOffset = toothOffset;
  
  QVector<QVector2D> pts;
  SimpleToothMaker().makeTeeth(pts, gear->getCurve(), params);

  gear->setBodyPath(pts);
  
  return gear;

}


class ArcFunction : public CurveFunction {
public:
  QVector2D operator()(double t) const { double phi = M_PI*0.5*t; return 200*QVector2D(cos(phi), sin(phi)); }
};


class MyGearBox : public GearBox {
public:
  virtual void setParams(double param1, double param2) {}
};


class MyGearBox1 : public MyGearBox {
public:
  MyGearBox1() {
    double toothLength = 30.0;
    addGear(makeCircularGear(toothLength, 53));
    addGear(makeCircularGear(toothLength, 21));
    addGear(makeCircularGear(toothLength, 31));
    GearLink *link1 = addLink(0,1);
    GearLink *link2 = addLink(1,2);
    link1->moveDriven(M_PI);
    link2->moveDriven(M_PI);
  }
  void setParams(double param1, double param2) {
    double psi = M_PI + param1;
    getLink(0)->moveDriven(psi);
    getLink(1)->moveDriven(psi);

    getGear(0)->setAngle(param2);
  }
};


class MyGearBox2 : public MyGearBox {
public:
  MyGearBox2() {
    double toothLength = 30.0;
    int toothCount = 23;
    double eccentricity = 0.6;
    for(int i=0;i<4;i++)
      addGear(makeEllipticGear(toothLength, toothCount, eccentricity));
    getGear(0)->setPosition(600,0);
    for(int i=1;i<4;i++)
      addLink(i-1,i)->moveDriven(M_PI);
  }
  void setParams(double param1, double param2) {
    double psi = M_PI + param1;
    
    getGear(0)->setAngle(param2);
  }
};


class MyGearBox3 : public MyGearBox {
public:
  MyGearBox3() {
    double r0 = 100, r1 = 350;
    int toothCount = 16;

    Gear *gear1 = addGear(makeSpiralGear(toothCount, r0, r1, 0.0));
    Gear *gear2 = addGear(makeSpiralGear(toothCount, r0, r1, 0.5));

    double dist = r0+r1;
    gear1->setPosition( dist*0.5,0);
    gear2->setPosition(-dist*0.5,0);

    double L1 = gear1->getCurve()->getLength();
    double L2 = gear2->getCurve()->getLength();
    double s1 = L1/2 ;
    double s2 = L2/2 ;


    gear1->setAngle(M_PI - gear1->getCurve()->getPhiFromS(s1));
    gear2->setAngle(-gear2->getCurve()->getPhiFromS(s2));

    addLink(0,1);
  }
  void setParams(double param1, double param2) {
    double psi = M_PI + param1;
    
    getGear(0)->setAngle(param2);

  }
};


class MyGearBox4 : public MyGearBox {
  double m_r0, m_r1;
  int m_toothCount;
  double m_theta;

public:
  MyGearBox4() : m_toothCount(11) {
    double toothLength = 80.0;
    Gear *gear1 = addGear(makeCircularGear(toothLength, m_toothCount, 1));
    Gear *gear2 = addGear(makeCircularGear(toothLength, m_toothCount, 1));
    double dist = gear1->getCurve()->getRFromS(0) + gear2->getCurve()->getRFromS(0);
    double radius = dist*0.5;
    double h = toothLength*0.5*0.75;
    m_r0 = radius - h/2;
    m_r1 = m_r0 + h;
    m_theta = 0.2*(2*M_PI/m_toothCount);
    gear2->setPosition(-dist,0);
  }


  void setParams(double param1, double param2) {
    
    getGear(0)->setAngle(param2 * 0.01);
  }

  void getSegments(Gear *gear, QList<QLineF> &pts) {
    for(int i=0;i<m_toothCount;i++)
    {
      double phi_i = 2*M_PI*i/m_toothCount + gear->getAngle();
      for(int j=0;j<2;j++)
      {
        double phi = phi_i + (-1+2*j)*m_theta;
        QPointF u(cos(phi), sin(phi));
        pts.append( QLineF(gear->getPosition() + m_r0 * u,  gear->getPosition() + m_r1 * u ) );
      }
    }
  }

  void drawSegments(QPainter &pa, const QList<QLineF> &pts) {
    QPainterPath pp;
    for(int i=0;i<pts.count();i++)
    {
      pp.moveTo(pts[i].p1()); pp.lineTo(pts[i].p2());
    }
    pa.drawPath(pp);
  }

  bool intesectCircle(const QLineF &line, Gear *gear, QPointF *p) {
    // p = line.p1() + (line.p2()-line.p1())*lambda 
    // |A + B*lambda| = m_r1
    // lambda^2*<B,B> + 2*lambda*<A,B> +<A,A>-m_r1 = 0

    QPointF A = line.p1() - gear->getPosition();
    QPointF B = line.p2()-line.p1();
    double a = B.x()*B.x()+B.y()*B.y();
    double b = A.x()*B.x()+A.y()*B.y();
    double c = A.x()*A.x()+A.y()*A.y() - m_r1*m_r1;
    double dsc = b*b-a*c;
    if(dsc<0) return false;
    double lambda = (-b-sqrt(dsc))/a;
    if(0<=lambda && lambda<1.0) 
    {
      *p = line.p1() + B*lambda;
      return true;
    }
    else return false;

  }
  

  void draw(QPainter &pa) {

    QList<QLineF> pts1, pts2;
    getSegments(getGear(0), pts1);
    getSegments(getGear(1), pts2);
    
    double phi1,phi2;
    bool flag = false;
    QPointF q1,q2;
    QLineF line1,line2;
    QPointF intersection;
    bool intersectionFlag = false;

    for(int i=0;i<pts1.count(); i++)
    {
      for(int j=0;j<pts2.count(); j++) 
      {
        if(pts1[i].intersect(pts2[j], &intersection) == QLineF::BoundedIntersection)
        {
          intersectionFlag=true;
          line1 = pts1[i];
          line2 = pts2[j];
          break;
        }
      }
      if(intersectionFlag) break;
    }

    if(intersectionFlag)
    {
      QPointF gearPos1 = getGear(0)->getPosition();
      QPointF gearPos2 = getGear(1)->getPosition();
      q2 = line1.p2() - gearPos2;
      q1 = line2.p2() - gearPos2;
      QPointF dq = q1-q2;
      if(dq.y()*q2.x()-dq.x()*q2.y()>0)
      {
        phi1 = atan2(q1.y(),q1.x());
        phi2 = atan2(q2.y(),q2.x());
        flag = true;

        getGear(1)->setAngle(getGear(1)->getAngle()+phi2-phi1);

      }
      else
      {
        QPointF q;
        if(intesectCircle(line1, getGear(1), &q))
        {
          q1 = line2.p2() - gearPos2;
          phi1 = atan2(q1.y(),q1.x());
          q2 = q - gearPos2;
          phi2 = atan2(q2.y(),q2.x());
          flag = true;

          getGear(1)->setAngle(getGear(1)->getAngle()+phi2-phi1);


        }
      }

    }



    MyGearBox::draw(pa);

    

    
    pa.setPen(Qt::red);
    pa.setBrush(Qt::NoBrush);
    drawSegments(pa, pts1);
    drawSegments(pa, pts2);

    pa.setPen(Qt::black);
    
    if(flag)
    {
      QPointF p = getGear(1)->getPosition();
      pa.setPen(Qt::blue);
      pa.drawLine(p, p+300*QPointF(cos(phi1),sin(phi1)));
      pa.setPen(Qt::green);
      pa.drawLine(p, p+300*QPointF(cos(phi2),sin(phi2)));

      pa.setPen(Qt::magenta);
      pa.drawEllipse(p+q1,2,2);
      pa.drawEllipse(p+q2,2,2);

    }
    
  }

};


class SandboxPage : public Page, public Pannable {
  int m_mode;
  double m_param1, m_param2;
  MyGearBox* m_gearBox;
public:
  SandboxPage() : Page("sandbox"), m_mode(0), m_param1(0), m_param2(0) { 
    m_gearBox = new MyGearBox4();
  }
  ~SandboxPage() {  }

  void draw(QPainter &pa);
  bool onKey(int key) {
    if(Qt::Key_1 == key) { delete m_gearBox; m_gearBox = new MyGearBox1(); }
    else if(Qt::Key_2 == key) { delete m_gearBox; m_gearBox = new MyGearBox2(); }
    else if(Qt::Key_3 == key) { delete m_gearBox; m_gearBox = new MyGearBox3(); }
    else if(Qt::Key_4 == key) { delete m_gearBox; m_gearBox = new MyGearBox4(); }
    else return false;
    return true;
  }
  void drawCurve(QPainter &pa, PitchCurve *crv);
  void drag(int dx, int dy, int modifier) {
    m_param1 += dx*0.01;
    m_param2 += dy*0.01;
  }
} sandboxPage;

/*
void SandboxPage::draw(QPainter &pa)
{
  pa.setPen(Qt::gray);
  QPainterPath pp;
  pp.moveTo(-500,0);pp.lineTo(500,0);
  pp.moveTo(0,-500);pp.lineTo(0,500);
  pa.drawPath(pp);


  PitchCurve *crv = 0;
  if(m_mode==0)      crv = new PitchCurve(EllipseFunction(200.0,0.0),100);  
  else if(m_mode==1) crv = new PitchCurve(EllipseFunction(200.0,0.6),100); 
  else if(m_mode==2) crv = new PitchCurve(SpiralFunction(100.0,250.0),100); 
  else if(m_mode==3) crv = new PitchCurve(ArcFunction(),25); 

  if(crv)
  {
   
     drawCurve(pa,crv);
    delete crv;
    

    
  }
}

*/

void SandboxPage::draw(QPainter &pa)
{
  pa.setPen(Qt::gray);
  QPainterPath pp;
  pp.moveTo(-500,0);pp.lineTo(500,0);
  pp.moveTo(0,-500);pp.lineTo(0,500);
  pa.drawPath(pp);

  m_gearBox->setParams(m_param1, m_param2);
  m_gearBox->draw(pa);

}



void SandboxPage::drawCurve(QPainter &pa, PitchCurve *crv)
{
  QPainterPath pp;
  pp.moveTo(crv->getPoint(0).pos.toPointF());
  for(int i=1;i<crv->getPointCount();i++) pp.lineTo(crv->getPoint(i).pos.toPointF());
  pa.setPen(Qt::red);
  pa.setBrush(Qt::NoBrush);
  pa.drawPath(pp);

  pp = QPainterPath();
  for(int i=0;i<crv->getPointCount();i++) 
  {
    PitchCurve::Point pt = crv->getPoint(i);
    pp.moveTo(pt.pos.toPointF());
    pp.lineTo((pt.pos + pt.right*10).toPointF());
  }
  pa.setPen(Qt::black);
  pa.drawPath(pp);


  pa.setPen(Qt::blue);
  int m = 20;
  for(int i=0;i<m;i++)
  {
    double s = crv->getLength() * i / (m-1);
    QPointF p = crv->getPosFromS(s).toPointF();
    pa.drawEllipse(p,3,3);

  }

}

