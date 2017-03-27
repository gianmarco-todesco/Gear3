#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"
#include "ClipperWrapper.h"
#include "PolygonSweeper.h"

#include <QLineF>

#ifdef CICCIO
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

  double qqq[360];
  QPointF m_otherPoint;

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
    for(int i=0;i<360;i++) qqq[i]=0.0;

  }


  void setParams(double param1, double param2) {
    
    getGear(0)->setAngle(param1 * 0.1);
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
  
  bool push1(const QPointF &p1, const QPointF &p2, double direction, double &angle, QPointF &contactPoint) {
    QPointF c = getGear(1)->getPosition();
    QPointF v1 = p1-c, v2 = p2-c;
    if(QVector2D(v1).length()>m_r1) return false;
    QPointF w1(-v2.y(), v2.x());
    double pd = w1.x()*v1.x()+w1.y()*v1.y();
    if(pd*direction>=0) return false;

    double dphi = atan2(v1.y(),v1.x()) - atan2(v2.y(),v2.x());
    
    Q_ASSERT(dphi * direction < 0.0);
    if(fabs(dphi) > fabs(angle)) 
    {
      angle = dphi;
      contactPoint = p1;
      m_otherPoint = p2;
    }
    return true;
  }

  bool push2(const QPointF &p1, const QPointF &p2, double direction, double &angle, QPointF &contactPoint) {
 
    QPointF p;
    QPointF c = getGear(1)->getPosition();
    
    if(!intesectCircle(QLineF(getGear(0)->getPosition(), p1), getGear(1), &p)) return false;
    QVector2D e0 = QVector2D(p2-c).normalized();
    QVector2D e1(-e0.y(),e0.x());
    QVector2D q(p-c);
    double dphi = atan2(QVector2D::dotProduct(e1,q), QVector2D::dotProduct(e0,q));
   //  Q_ASSERT(dphi * direction > 0.0);
    if(fabs(dphi) > fabs(angle)) 
    {
      angle = dphi;
      contactPoint = p;
    }
    return true;
  }


  void draw(QPainter &pa) {

    Gear *gear1 = getGear(0);
    Gear *gear2 = getGear(1);

    QPointF pos1 = gear1->getPosition();
    QPointF pos2 = gear2->getPosition();

    QList<QLineF> segm1, segm2;
    getSegments(gear1, segm1);
    getSegments(gear2, segm2);

    QPointF contactPoint;
    double angle = 0.0;
    bool contact = false;

    for(int i=0;i+1<segm1.count();i+=2)
    {
      QPointF p1 = (segm1[i].p2()+segm1[i+1].p2())*0.5;
      if((p1-pos1).x()>0) continue;
      for(int j=0;j+1<segm2.count();j+=2)
      {
        QPointF p2 = (segm2[j].p2()+segm2[j+1].p2())*0.5;
        if((p2-pos2).x()<0) continue;
        if(QVector2D(p1-p2).length() > (m_r1-m_r0)) continue;
        bool ret;
        ret = push1(segm1[i+1].p2(), segm2[j+1].p2(), 1.0, angle, contactPoint);
        contact = contact || ret;
       // ret = push2(segm1[i+1].p2(), segm2[j+1].p2(), 1.0, angle, contactPoint);
       // contact = contact || ret;
      }
    }
    // Q_ASSERT(fabs(angle)<0.1);
    //gear2->setAngle(gear2->getAngle() + angle);
    //qqq[(int)(360.0*gear1->getAngle())%360] = 360*gear2->getAngle();

    MyGearBox::draw(pa);

    if(contact)
    {
      pa.setPen(Qt::magenta);
      pa.drawEllipse(contactPoint, 3,3);
      pa.setPen(Qt::blue);
      pa.drawEllipse(m_otherPoint, 3,3);
    }

    /*
    QPainterPath pp2;
    for(int i=0;i<360;i++)
    {
      if(i==0) pp2.moveTo(i*0.2,qqq[i]*0.2); else pp2.lineTo(i*0.2,qqq[i]*0.2);
    }
    pa.setPen(Qt::red);
    pa.setBrush(Qt::NoBrush);
    pa.drawPath(pp2);
    */
    /*

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
    */

    
  }

};


class SandboxPage2 : public Page, public Pannable {
  int m_mode;
  double m_param1, m_param2;
  MyGearBox* m_gearBox;
public:
  SandboxPage2() : Page("sandbox2"), m_mode(0), m_param1(0), m_param2(0) { 
    m_gearBox = new MyGearBox4();
  }
  ~SandboxPage2() {  }

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
} sandboxPage2;

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

void SandboxPage2::draw(QPainter &pa)
{
  pa.setPen(Qt::gray);
  QPainterPath pp;
  pp.moveTo(-500,0);pp.lineTo(500,0);
  pp.moveTo(0,-500);pp.lineTo(0,500);
  pa.drawPath(pp);

  m_gearBox->setParams(m_param1, m_param2);
  m_gearBox->draw(pa);

}



void SandboxPage2::drawCurve(QPainter &pa, PitchCurve *crv)
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


class Sweeper {
  QVector<QVector2D> m_shape;
  QPainterPath m_mainPath, m_borderPath;
  QList<QVector2D> m_lastPolygon;
  QVector2D m_lastCenter;
  QVector<QLineF> m_lines[2];
  QVector<QPointF> m_outline;

public:
  Sweeper() {}

  void setShape(const QVector<QVector2D> &shape) { m_shape = shape; }

  void add(const QMatrix &matrix);

  void buildOutline() {
    
    QVector<QPointF> borders[2];      
    for(int i=0;i<2;i++)
    {
      QVector<QPointF> &border = borders[i];
      border.append(m_lines[i][0].p1());
      for(int j=0;j+1<m_lines[i].count();j++)
      {
        QLineF line1 = m_lines[i][j];
        QLineF line2 = m_lines[i][j+1];
        if(QVector2D(line1.p2()-line2.p1()).lengthSquared()<4)
        {
          border.append(line1.p2());
        }
        else
        {
          QPointF p;
          if(line1.intersect(line2, &p)!=QLineF::BoundedIntersection)
          {
            border.append(line1.p2());
            border.append(line2.p1());
          }
          else
          {
            border.append(p);
          }
        }
      }
      border.append(m_lines[i].back().p2());
    }
    for(int i=0;i<borders[0].count();i++) m_outline.append(borders[0][i]);
    for(int i=borders[1].count()-1;i>=0;i--) m_outline.append(borders[1][i]);
  }

  void getOutline(QVector<QPointF> &outline) { outline = m_outline; }

  void draw(QPainter &pa) {

    pa.setPen(Qt::gray);
    pa.drawPath(m_mainPath);
    
    for(int i=0;i<2;i++)
    {
      QVector<QPointF> border;
      border.append(m_lines[i][0].p1());
      for(int j=0;j+1<m_lines[i].count();j++)
      {
        QLineF line1 = m_lines[i][j];
        QLineF line2 = m_lines[i][j+1];
        if(QVector2D(line1.p2()-line2.p1()).lengthSquared()<4)
        {
          border.append(line1.p2());
        }
        else
        {
          QPointF p;
          if(line1.intersect(line2, &p)!=QLineF::BoundedIntersection)
          {
            border.append(line1.p2());
            border.append(line2.p1());
          }
          else
          {
            border.append(p);
          }
        }
      }
      border.append(m_lines[i].back().p2());
      
      QPainterPath pp;
      pp.moveTo(border[0]);
      for(int j=1;j<border.count();j++) 
        pp.lineTo(border[j]);
      pa.setPen(i==0 ? Qt::red : Qt::green);
      pa.drawPath(pp);
    }


    pa.setPen(Qt::red);
    for(int i=0;i<2;i++)
    {
      for(int j=0;j<m_lines[i].count();j++)
        pa.drawLine(m_lines[i][j].p1(), m_lines[i][j].p2());

    }
    
  }
};


void Sweeper::add(const QMatrix &matrix)
{
  QList<QVector2D> polygon;
  QVector2D center;
  int n = m_shape.count();
  for(int i=0;i<n;i++) { QVector2D p(matrix.map(m_shape[i].toPointF())); polygon.append(p); center += p;} 
  center *= 1.0/n;

  m_mainPath.moveTo(polygon[0].toPointF());
  for(int i=1;i<m_shape.count();i++) m_mainPath.lineTo(polygon[i].toPointF());

  
  QVector2D movement = center - m_lastCenter;
  if(!m_lastPolygon.empty() && movement.length()>1)
  {
    for(int side=0; side<2; side++) {
      double sgn = 1-2*side;

      int k = -1;
      double maxDist = 0;
      QVector2D perp(-movement.y(), movement.x());

      for(int i=0;i<n;i++) 
      {
        double d = sgn*QVector2D::dotProduct(polygon[i]-m_lastCenter, perp);
        if(k<0 || d>maxDist) { maxDist = d; k = i; }
      }
      int k1 = k;
      int k2 = 0;
      int count=0;
    
      for(;;)
      {
        QVector2D v = m_lastPolygon[k2] - polygon[k1];
        for(int i=0;i<n;i++)
        {
          QVector2D v2 = m_lastPolygon[i] - polygon[k1];
          if(sgn*(-v2.x()*v.y()+v2.y()*v.x())<0) { v=v2; k2 = i; }
        }
        bool done = true;
        v = polygon[k1] - m_lastPolygon[k2];
        for(int i=0;i<n;i++)
        {
          QVector2D v2 = polygon[i] - m_lastPolygon[k2]; 
          if(sgn*(-v2.x()*v.y()+v2.y()*v.x())>0) { v=v2; k1 = i; done = false;}
        }
        if(done|| ++count>100) break;
      }
      QVector<QLineF> &lines = m_lines[side];
      lines.append(QLineF(m_lastPolygon[k2].toPointF(), polygon[k1].toPointF()));
    }


    /*
    k1 = k;
    k2 = 0;
    int count=0;
    for(;;)
    {
      QVector2D v = m_lastPolygon[k2] - polygon[k1];
      for(int i=0;i<n;i++)
      {
        QVector2D v2 = m_lastPolygon[i] - polygon[k1];
        if(v2.y()*v.x()-v2.x()*v.y()>0) { v=v2; k2 = i; }
      }
      bool done = true;
      v = polygon[k1] - m_lastPolygon[k2];
      for(int i=0;i<n;i++)
      {
        QVector2D v2 = polygon[i] - m_lastPolygon[k2]; 
        if(v2.y()*v.x()-v2.x()*v.y()<0) { v=v2; k1 = i; done = false;}
      }
      if(done || ++count>10) break;
    }

    m_borderPath.moveTo(polygon[k1].toPointF());
    m_borderPath.lineTo(m_lastPolygon[k2].toPointF());
    */
    

  }
  m_lastPolygon.swap(polygon);
  m_lastCenter = center;


}




class SandboxPage : public Page, public Pannable {
public:
  SandboxPage() : Page("sandbox") { 
  }
  ~SandboxPage() {  }

  void draw(QPainter &pa, const QVector<QPointF> &pts) {
     QPainterPath pp;
     pp.moveTo(pts.back());
     for(int i=0;i<pts.count();i++) pp.lineTo(pts[i]);
     pa.drawPath(pp);
  }

  void draw(QPainter &pa) {
    PitchCurve *crv = new PitchCurve(EllipseFunction(100,0.6),200);
    QPainterPath pp;
    pp.moveTo(crv->getPoint(0).pos.toPointF());
    for(int i=1;i<crv->getPointCount();i++) pp.lineTo(crv->getPoint(i).pos.toPointF());
    pa.setPen(Qt::black);
    pa.drawPath(pp);

    QVector<QVector2D> pts;
    pts << QVector2D(20,-20) << QVector2D(5,20) << QVector2D(-5,20) << QVector2D(-20,-20);

    ClipperWrapper cw;
 
    int n = 10;
    for(int k=0;k<1;k++)
    {
      double s0 = crv->getLength()*k/n;
      for(int d=0; d<1;d++)
      {
        //Sweeper sweeper;
        //sweeper.setShape(pts);
        QVector<QPointF> pts2;

        int m = 20;
        for(int i=0;i<m;i++)
        {
          double t = (double)i/(double)(m-1);
          double s = s0 + 200*t*(1-2*d);
          PitchCurve::Point pt = crv->getPointFromS(s);
          QMatrix matrix;
          matrix.translate(pt.pos.x(), pt.pos.y());
          matrix.rotate(180.0*atan2(pt.right.y(), pt.right.x())/M_PI + 90);
          matrix.translate(-s + s0,0);
          for(int j=0;j<pts.count();j++) pts2.append(matrix.map(pts[j].toPointF()));
          //sweeper.add(matrix);
          pa.setPen(Qt::gray);
          QVector<QPointF> trpts;
          for(int j=0;j<pts.count();j++) trpts.append(matrix.map(pts[j].toPointF()));
          draw(pa,trpts);
          
        }

        pa.setPen(Qt::black);
        QVector<QPointF> strip;
        for(int j=0;j<m;j++) strip.append(pts2[pts.count()*j+0]);
        for(int j=m-1;j>=0;j--) strip.append(pts2[pts.count()*j+1]);

        draw(pa,strip);
        /*
        //sweeper.buildOutline();
        QVector<QPointF> outline;
        sweeper.getOutline(outline);
        sweeper.draw(pa);
        QPainterPath pp;
        pp.moveTo(outline.back());
        for(int i=0;i<outline.count();i++) pp.lineTo(outline[i]);
        pa.setPen(Qt::magenta);
        pa.drawPath(pp);
        cw.add(outline);
        */

      }
    }

    /*
    QVector<QVector<QPointF> > outlines;
    cw.getOutline(outlines);
    for(int i=0;i<outlines.count();i++)
    {
      QVector<QPointF> &outline = outlines[i];
      QPainterPath pp;
      pp.moveTo(outline.back());
      for(int i=0;i<outline.count();i++) pp.lineTo(outline[i]);
      pa.setPen(Qt::cyan);
      pa.drawPath(pp);

    }


    QVector<QPointF> urka;
    for(int i=0;i<crv->getPointCount();i++)
    {
      PitchCurve::Point pt = crv->getPoint(i);
      urka.append((pt.pos + pt.right*30).toPointF());
    }

    pp = QPainterPath() ;
    pp.moveTo(urka.back());
    for(int i=0;i<urka.count();i++) pp.lineTo(urka[i]);
    pa.setPen(Qt::blue);
    pa.drawPath(pp);

    cw.intersect(urka);
    outlines.clear();
    cw.getOutline(outlines);
    pp = QPainterPath ();      
    for(int i=0;i<outlines.count();i++)
    {
      QVector<QPointF> &outline = outlines[i];
      pp.moveTo(outline.back());
      for(int i=0;i<outline.count();i++) pp.lineTo(outline[i]);
      pp.closeSubpath();
    }
    pa.setPen(Qt::blue);
    pa.setBrush(Qt::yellow);
    //pa.drawPath(pp);


    double s = getParameter();
    PitchCurve::Point pt = crv->getPointFromS(s);
    pp = QPainterPath();
    for(int i= -10; i<=10; i++)
    {
      QMatrix matrix;
      matrix.translate(pt.pos.x(), pt.pos.y());
      matrix.rotate(180.0*atan2(pt.right.y(), pt.right.x())/M_PI + 90);
      matrix.translate(-s + crv->getLength()*i/n,0);
      pp.moveTo(matrix.map(pts.back().toPointF()));
      for(int j=0;j<pts.count();j++)
        pp.lineTo(matrix.map(pts[j].toPointF()));
    }
    pa.setBrush(Qt::cyan);
    pa.setPen(Qt::blue);
    pa.drawPath(pp);
    */

  }



  bool onKey(int key) {
    return false;
  }
} sandboxPage;

#endif


class SandboxPage : public Page, public Pannable {
public:
  SandboxPage() : Page("sandbox") { 
  }
  ~SandboxPage() {  }

  void draw(QPainter &pa) {


    PitchCurve *crv = new PitchCurve(EllipseFunction(100,0.6),200);

    QPainterPath pp;
    pp.moveTo(crv->getPoint(0).pos.toPointF());
    for(int i=1;i<crv->getPointCount();i++) pp.lineTo(crv->getPoint(i).pos.toPointF());
    pa.setPen(Qt::black);
    pa.drawPath(pp);

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

      pa.setBrush(Qt::yellow);
      pa.setPen(QPen(Qt::black, 0.5));

      QPainterPath pp;
      for(int i=0;i<strip.count();i++)
      {
        pp.addPolygon(QPolygonF(strip[i]));
        pp.closeSubpath();
      }
      pa.drawPath(pp);

      pa.setBrush(Qt::NoBrush);
      pa.setPen(Qt::gray);
      pp = QPainterPath();
      for(int i=0;i<polygons.count();i++) 
      {
        pp.addPolygon(polygons[i]);
        pp.closeSubpath();
      }
      pa.drawPath(pp);
    }  

    QVector<QVector<QPointF > > outline;
    cw.getResult(outline);
    pp = QPainterPath ();
    for(int i=0;i<outline.count();i++) { pp.addPolygon(QPolygonF(outline[i])); pp.closeSubpath(); }
    pa.setPen(QPen(Qt::blue, 2));
    pa.setBrush(Qt::NoBrush);
    pa.drawPath(pp);
  }

  bool onKey(int key) {
    return false;
  }
} sandboxPage;
