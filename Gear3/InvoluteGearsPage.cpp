#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"
#include <QVector2D>

/*
namespace {


Gear *makeSimpleGear2(int toothCount)
{
  double radius = toothCount * 30 / (2*M_PI);
  Gear *gear = new Gear(new PitchCurve(EllipseFunction(radius,0.0)));
  SimpleToothMaker ctm;
  SimpleToothMaker::Params params;
  params.toothHeight = 20;
  params.toothCount = toothCount;

  QVector<QVector2D> pts;
  ctm.makeTeeth(pts, gear->getCurve(), params);

  gear->setBodyPath(pts);
  return gear;
}

}

*/


class InvoluteGearsPage : public Page, public Pannable {
  GearBox m_gearBox;
  double m_r0,m_r1,m_r2;
  int m_toothCount;
  double m_theta;

  double m_w, m_d;
  QPointF m_center1, m_center2, m_q1, m_q2;
  int m_flags;
  int m_mode;
  bool m_showContactPoints;
  bool m_showPitchLine;
  double m_maxParameter;

  int m_toothIndex;


public:
  InvoluteGearsPage() 
    : Page("involuteGears")
    , m_flags(0)
    , m_mode(1)
    , m_showContactPoints(false)
    , m_showPitchLine(false)
    , m_maxParameter(0)
    , m_toothIndex(0)
  { 

    m_toothCount = 17;
    double radius = m_toothCount * 90 / (2*M_PI);

    CircularToothMaker ctm;
    CircularToothMaker::Params params;
    params.pitchRadius = radius;
    params.toothCount = m_toothCount;
    params.toothDedendum = 35;
    params.toothAddendum = 20;

    QVector<QVector2D> pts;
    ctm.makeTeeth(pts, params);

    m_theta = atan2(pts[0].y(), pts[0].x());

    m_r0 = params.pitchRadius - params.toothDedendum;
    m_r1 = params.pitchRadius;
    m_r2 = params.pitchRadius + params.toothAddendum;

    Gear *gear1 = new Gear(new PitchCurve(EllipseFunction(radius,0.0)));
    gear1->setBodyPath(pts);

    Gear *gear2 = new Gear(new PitchCurve(EllipseFunction(radius,0.0)));
    gear2->setBodyPath(pts);

    m_gearBox.addGear(gear1); 
    m_gearBox.addGear(gear2); 
    m_gearBox.addLink(0,1);
    double dist = m_gearBox.getLink(0)->getDistance();
    gear1->setPosition(dist/2,0);
    gear2->setPosition(-dist/2,0);
    gear1->setBrush(Qt::yellow);
    gear2->setBrush(QColor(200,123,10));

    m_center1 = gear1->getPosition();
    m_center2 = gear2->getPosition();

    // semidistanza fra i centri
    m_w = 0.5*m_gearBox.getLink(0)->getDistance();
    // semilunghezza pitchline
    m_d = sqrt(m_w*m_w-m_r0*m_r0);
    // coordinate rispetto a center2 di un'estremità della pitchline
    QPointF q(m_r0*m_r0/m_w, -m_d*m_r0/m_w);
    m_q1 = m_center1 - q;
    m_q2 = m_center2 + q;
  }
  ~InvoluteGearsPage() { }

  void draw(QPainter &pa);
  void draw2(QPainter &pa);

  void drawGear(QPainter &pa, Gear *gear);

  void drawContactPoints(QPainter &pa);

  void drawString1(QPainter &pa);
  void drawString2(QPainter &pa, const QPointF &center, const QPointF &p1, const QPointF &p2);

  void setToothIndex();

  QPointF involute(const Gear *gear, int side, int toothIndex, double t);
  
  void drawInvoluteCurve(QPainter &pa);

  void drag(int dx, int dy, int modifiers) {
    Gear *gear = m_gearBox.getGear(0);
    gear->setAngle(gear->getAngle() + dx * 0.01);

    setParameter(getParameter() + dy * 0.01);
  }


  bool onKey(int key) {
    if(key == Qt::Key_1) m_mode = 1;
    else if(key == Qt::Key_2) m_mode = 2;
    else if(key == Qt::Key_3) { m_mode = 3; m_maxParameter = 0.0; setParameter(0.0); setToothIndex(); }
    else if(key == Qt::Key_4) { m_mode = 4; setToothIndex(); }

    else if(key == Qt::Key_C) 
      m_showContactPoints = !m_showContactPoints;
    else if(key == Qt::Key_P) 
      m_showPitchLine = !m_showPitchLine;
    else return false;
    return true;
  }


} involuteGearsPage;


void InvoluteGearsPage::draw(QPainter &pa)
{
  Gear *gear1 = m_gearBox.getGear(0);
  Gear *gear2 = m_gearBox.getGear(1);

  // gear1->setAngle(getParameter()*0.01);
  drawGear(pa, gear1);

  if(m_mode==2 || m_mode==4)
  {
    m_gearBox.getLink(0)->update();
    drawGear(pa, gear2);

    if(m_showPitchLine)
    {
      pa.setPen(Qt::black);
      pa.drawLine(m_q1,m_q2);
    }

    if(m_showContactPoints) 
      drawContactPoints(pa);
  }

  if(m_mode==3 || m_mode==4)
  {
    drawString1(pa);
  }


}

void InvoluteGearsPage::drawGear(QPainter &pa, Gear *gear)
{
  pa.save();
  gear->rotoTranslate(pa);
  pa.setPen(Qt::black);
  pa.setBrush(gear->getBrush());
  pa.drawPath(gear->getBodyPath());
  pa.setBrush(Qt::NoBrush);
  pa.setPen(QColor(120,120,120));
  double r = m_r0-10;
  pa.drawEllipse(QPointF(0,0), r,r);
  r-=4;
  pa.drawEllipse(QPointF(0,0), r,r);

  pa.setPen(Qt::black);
  pa.drawEllipse(QPointF(0,0), m_r0,m_r0);

  pa.restore();
}



void InvoluteGearsPage::draw2(QPainter &pa)
{
  Gear *gear1 = m_gearBox.getGear(0);
  Gear *gear2 = m_gearBox.getGear(1);

  // gear1->setAngle(getParameter() * 0.01);
  m_gearBox.draw(pa);

  drawInvoluteCurve(pa);

  QPointF center1 = gear1->getPosition();
  QPointF center2 = gear2->getPosition();

  pa.setPen(QPen(Qt::black,1));
  pa.setBrush(Qt::NoBrush);
  pa.drawEllipse(center2, m_r0, m_r0);
  pa.drawEllipse(center2, m_r2, m_r2);

  pa.drawEllipse(center1, m_r0, m_r0);
  pa.drawEllipse(center1, m_r0, m_r0);

  pa.drawLine(m_q1,m_q2);

  double theta = -m_theta + gear2->getAngle();
  pa.drawLine(center2, center2 + 400*QPointF(cos(theta),sin(theta)));

  drawContactPoints(pa);

  /*
  for(int i=0;i<2;i++)
  {
    for(int j=0;j<m_toothCount;j++)
    {
      for(int side = 0; side<2; side++)
      {
        QPainterPath pp;
        int m = 100;
        for(int k=0;k<m;k++)
        {
          double t = 1.0*k/(m-1);
          QPointF p = involute(m_gearBox.getGear(i),side,j,t);
          if(k==0) pp.moveTo(p); else pp.lineTo(p);
        }
        if(side==0) pa.setPen(QPen(Qt::red)); else pa.setPen(QPen(Qt::green));
        pa.drawPath(pp);
      }
    }
  }
  */
}


void InvoluteGearsPage::setToothIndex()
{
  int side = 0;
  double sgn = 1-2*side;
  int toothIndex = 0;
  Gear *gear = m_gearBox.getGear(0);

  // cerco il dente più adatto
  double alpha =  M_PI + sgn * m_theta - gear->getAngle();
  if(alpha<0) while(alpha<0) alpha += 2*M_PI;
  else if(alpha>=2*M_PI) while(alpha>=2*M_PI) alpha -= 2*M_PI;

  toothIndex = (int)floor(0.5 + m_toothCount * alpha / (2* M_PI));
  if(toothIndex<0) toothIndex = 0;
  else if(toothIndex>=m_toothCount) toothIndex %= m_toothCount;

  m_toothIndex = toothIndex;
}


void InvoluteGearsPage::drawString1(QPainter &pa)
{
  int side = 0;
  double sgn = 1-2*side;
  int toothIndex = m_toothIndex;
  Gear *gear = m_gearBox.getGear(0);

  double t;
  {
  double outPhi = M_PI - qMax(0.0, qMin(M_PI/2, getParameter()*0.6));
  if(m_mode == 4) { 
    QPointF q = m_q2 - m_gearBox.getGear(1)->getPosition(); 
    outPhi = M_PI - atan2(-q.y(),q.x()); 
  }
  
  double theta0 = -sgn * m_theta + gear->getAngle();
  double phi = theta0 + (toothIndex * 2* M_PI ) /m_toothCount;
  t = (phi - outPhi)*sgn;
  if(t<0.0) t = 0.0;
  }

  m_maxParameter = qMin(1.5, qMax(m_maxParameter,t));
  

  if(m_mode==3 && m_maxParameter>0.0)
  {
    QPainterPath curve;
    curve.moveTo(involute(gear, side, toothIndex, 0.0));
    int m = qMax(4, (int)(m_maxParameter*50.0));
    for(int i=1;i<m;i++)
    {
      curve.lineTo(involute(gear, side, toothIndex, m_maxParameter * i/(m-1)));
    }
    pa.setPen(QPen(Qt::red, 2));
    pa.drawPath(curve);
  }

  double theta0 = -sgn * m_theta + gear->getAngle();
  double phi = theta0 + (toothIndex * 2* M_PI ) /m_toothCount;
  double phi1 = phi - t * sgn;
  double cs = cos(phi1), sn = sin(phi1);
  QPointF p1 = gear->getPosition() + m_r0 * QPointF(cs,sn);
  QPointF p2 = p1 + (sgn * m_r0 * t) * QPointF(-sn,cs);
  double len = QVector2D(p1-p2).length();

  drawString2(pa, gear->getPosition(), p1,p2);
  

  if(m_mode==4)  
    drawString2(pa, m_gearBox.getGear(1)->getPosition(), m_q2, p2);

  pa.setPen(Qt::black);
  pa.setBrush(Qt::black);
  pa.drawEllipse(p2,5,5);
  // pa.drawEllipse(pts.back(),5,5);

  /*
  double phi = -0.5;
  double param = getParameter()*0.01;
  if(param<0.0)param=0.0;

  double theta = phi - param;
  QPointF p1 = m_center2 + m_r0 * QPointF(cos(theta), sin(theta));
  QPointF p2 = p1 + (m_r0 * param) * QPointF(-sin(theta), cos(theta));
  
  pa.setPen(QPen(Qt::black,2));
  pa.drawLine(p1,p2);
  pa.setPen(Qt::black);
  
  pa.setBrush(Qt::yellow);
  double s = 0;
  for(;;)
  {
    s += 6;  
    if(s>m_r0*param) break;
    double t = s/(m_r0 * param);
    QPointF p = p1*(1-t) + p2*t;
    pa.drawEllipse(p,3,3);
  }
  for(;;)
  {
    if(theta<-2*M_PI) break;
    QPointF p = m_center2 + m_r0 * QPointF(cos(theta), sin(theta));
    pa.drawEllipse(p,3,3);
    theta -= 6.0/m_r0;
  }


  pa.setBrush(Qt::red);
  pa.drawEllipse(p2,4,4);
  */

}

void InvoluteGearsPage::drawString2(QPainter &pa, const QPointF &center, const QPointF &p1, const QPointF &p2)
{
  double len = QVector2D(p2-p1).length();
  QPointF v1 = p1-center;
  double phi1 = atan2(v1.y(),v1.x());
  QVector<QPointF> pts;
  double s = 0.0;
  double ds = 2.0;
  while(s<m_r0*M_PI*1.5)
  {
    QPointF p;
    if(s<len) { p = p2 + (p1-p2)*(s/len); }
    else
    {
      double a = phi1 - (s-len)/m_r0;
      p = center + m_r0 * QPointF(cos(a),sin(a));
    }
    pts.append(p);
    s += ds;
  }
  
  pa.setBrush(Qt::NoBrush);
  
  QPainterPath stripPath;
  stripPath.moveTo(pts[0]);
  for(int i=1;i<pts.count();i++) stripPath.lineTo(pts[i]);
  QPen pen(QColor(200,0,200), 5);
  pen.setDashPattern(QVector<qreal>() << 10 << 10);
  pa.setPen(pen);
  pa.drawPath(stripPath);
  pen.setColor(QColor(0,200,100));
  pen.setDashOffset(10);
  pa.setPen(pen);
  pa.drawPath(stripPath);
}

QPointF InvoluteGearsPage::involute(const Gear *gear, int side, int toothIndex, double t)
{
  double sgn = 1.0 - 2*side;

  double theta0 = -sgn * m_theta + gear->getAngle();
  double phi = theta0 + (toothIndex * 2* M_PI ) /m_toothCount;
  double phi1 = phi - t * sgn;
  double cs = cos(phi1), sn = sin(phi1);
  QPointF p = gear->getPosition() + m_r0 * QPointF(cs,sn) + (sgn * m_r0 * t) * QPointF(-sn,cs);
  return p;
}

class SmallFoo {
  QPointF m_refPoint;
  QVector2D m_direction;
  QPointF m_oldp;
  double m_olddot;
  int m_count;
public:
  SmallFoo(const QPointF &pos, const QVector2D &direction) : m_refPoint(pos), m_direction(direction), m_olddot(0), m_count(0) {}

  bool addPoint(const QPointF &p, QPointF &rp) {
    double dot = QVector2D::dotProduct(QVector2D(p-m_refPoint), m_direction);
    if(dot>0.0)
    {
      if(m_count==0) rp = p;
      else
      {
        double s = (0.0-m_olddot)/(dot-m_olddot);
        rp = (1-s)*m_oldp + s*p;
      }
      return true;
    }
    m_olddot = dot;
    m_oldp = p;
    m_count++;
    return false;
  }

  void clear()
  {
    m_count = 0;
  }

};



void InvoluteGearsPage::drawContactPoints(QPainter &pa)
{
  Gear *gear1 = m_gearBox.getGear(0);
  Gear *gear2 = m_gearBox.getGear(1);

  double theta0 = -m_theta + gear2->getAngle();
  
  SmallFoo smallFoo(m_q2, QVector2D(m_q2 - m_center2).normalized()); 

  for(int i=0; i<m_toothCount; i++)
  {
    double phi = theta0 + i * 2*M_PI/m_toothCount;
    if(cos(phi)<0.75) continue;
    int m = 50;
    QPointF ip;
    bool found = false;
    for(int j=0; j<m;j++)
    {
      if(smallFoo.addPoint(involute(gear2, 0, i, (double)j/(double)(m-1)), ip))
      {
        found = true;
        break;
      }
    }

    if(found)
    {
      if(QVector2D(ip-m_center2).length()<m_r2 && QVector2D(ip-m_center1).length()<m_r2)
      {
        QPainterPath pp;
        pp.addEllipse(ip,6,6);
        pp.addEllipse(ip,3,3);

        pa.setPen(Qt::black);
        pa.setBrush(Qt::red);
        pa.drawPath(pp);
      }
    }
  }
  /*
  SmallFoo smallFoo2(m_q2, -QVector2D(m_q2 - m_center2).normalized()); 

  for(int i=0; i<m_toothCount; i++)
  {
    double phi = (-m_theta + gear1->getAngle()) + i * 2*M_PI/m_toothCount;
    if(cos(phi)>-0.75) continue;
    int m = 50;
    QPointF ip;
    bool found = false;
    for(int j=0; j<m;j++)
    {
      if(smallFoo2.addPoint(involute(gear1, 0, i, (double)j/(double)(m-1)), ip))
      {
        found = true;
        break;
      }
    }

    if(found)
    {
      pa.setPen(Qt::green);
      pa.drawEllipse(ip,5,5);
    }
  }
  */

}

void InvoluteGearsPage::drawInvoluteCurve(QPainter &pa)
{
  double phi = -0.5;
  double param = getParameter()*0.01;
  if(param<0.0)param=0.0;

  double theta = phi - param;
  QPointF p1 = m_center2 + m_r0 * QPointF(cos(theta), sin(theta));
  QPointF p2 = p1 + (m_r0 * param) * QPointF(-sin(theta), cos(theta));
  
  pa.setPen(QPen(Qt::black,2));
  pa.drawLine(p1,p2);
  pa.setPen(Qt::black);
  
  pa.setBrush(Qt::yellow);
  double s = 0;
  for(;;)
  {
    s += 6;  
    if(s>m_r0*param) break;
    double t = s/(m_r0 * param);
    QPointF p = p1*(1-t) + p2*t;
    pa.drawEllipse(p,3,3);
  }
  for(;;)
  {
    if(theta<-2*M_PI) break;
    QPointF p = m_center2 + m_r0 * QPointF(cos(theta), sin(theta));
    pa.drawEllipse(p,3,3);
    theta -= 6.0/m_r0;
  }


  pa.setBrush(Qt::red);
  pa.drawEllipse(p2,4,4);
}

