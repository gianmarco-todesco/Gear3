#include "Gear.h"
#include "PitchCurve.h"
#include <qmath.h>


Gear::Gear(PitchCurve *curve) 
  : m_curve(curve)
  , m_angle(0)
{
  updatePitchPath();
}

Gear::~Gear()
{
  delete m_curve;
}


void Gear::draw(QPainter &pa)
{
  pa.save();
  rotoTranslate(pa);

  pa.setPen(QPen(Qt::blue, 2));
  pa.setBrush(Qt::cyan);
  pa.drawPath(m_bodyPath);

  pa.setPen(QPen(Qt::magenta, 0, Qt::DotLine));
  pa.setBrush(Qt::NoBrush);
  pa.drawPath(m_pitchLinePath);

  QPointF p = m_curve->getPoint(0).pos.toPointF();
  pa.setPen(Qt::black);
  pa.drawLine(p,-p*10);

  pa.restore();
}

void Gear::rotoTranslate(QPainter &pa)
{
  pa.translate(getPosition());
  pa.rotate(getAngle()*180.0/M_PI);
}


void Gear::setBodyPath(const QVector<QVector2D> &pts)
{
  QVector2D oldp = pts[0];
  m_bodyPath.moveTo(oldp.toPointF());
  for(int i=1;i<pts.count();i++) 
  {
    QVector2D p = pts[i];
    if((p-oldp).lengthSquared()<1) continue;
    oldp = p;
    m_bodyPath.lineTo(p.toPointF());
  }
  m_bodyPath.closeSubpath();
  m_bodyPath.addEllipse(QPointF(0,0), 5,5);
}


void Gear::updatePitchPath()
{
  QVector2D lastPos = m_curve->getPoint(0).pos;
  m_pitchLinePath.moveTo(lastPos.toPointF());
  double ds = 5, s = ds;
  while(s<m_curve->getLength())
  {
    lastPos = m_curve->getPosFromS(s);
    m_pitchLinePath.lineTo(lastPos.toPointF());
    s += ds;
  }
}


//=============================================================================


GearLink::GearLink(Gear*driver, Gear *driven) 
  : m_gear1(driver)
  , m_gear2(driven)
  , m_theta1(driver->getAngle())
  , m_theta2(driven->getAngle())
{ 
  m_distance = computeDistance();
}

GearLink::~GearLink()
{
}

double GearLink::computeDistance() const
{
  double r0 = m_gear1->getCurve()->getRFromPhi(M_PI - m_theta1);
  double r1 = m_gear2->getCurve()->getRFromPhi(-m_theta2);
  return r0 + r1;
}


void GearLink::update()
{
  const PitchCurve *crv1 = m_gear1->getCurve();
  const PitchCurve *crv2 = m_gear2->getCurve();
  
  QPointF d = m_gear2->getPosition() - m_gear1->getPosition();
  double psi = atan2(d.y(),d.x()) - M_PI;

  double angle1 = m_gear1->getAngle() - psi;

  /*
  QPointF p = m_driven->getPosition() - m_driver->getPosition();
  double psi = atan2(p.y(),p.x()); if(psi<0)psi+=2*M_PI;
  double s = driverCrv->getSFromPhi(M_PI*0) - driverCrv->getSFromPhi(psi-m_driver->getAngle());
  m_driven->setAngle(drivenCrv->getPhiFromS(-s) + psi - M_PI); 
  */

  double s1a = crv1->getSFromPhi(M_PI-m_theta1);
  double s1b = crv1->getSFromPhi(M_PI-angle1);

  double s1 = s1a-s1b;

  double s2a = crv2->getSFromPhi(-m_theta2) ;
  double s = s1 - s2a;

  double angle2 = - crv2->getPhiFromS( 
       crv1->getSFromPhi(M_PI-m_theta1) 
     - crv1->getSFromPhi(M_PI-angle1) 
     + crv2->getSFromPhi(-m_theta2) );
  m_gear2->setAngle(angle2 + psi);

}


void GearLink::moveDriven(double psi)
{
  m_gear2->setPosition(m_gear1->getPosition() + m_distance * QPointF(cos(psi), sin(psi)));
}

//=============================================================================

GearBox::GearBox()
{
}

GearBox::~GearBox()
{
  clear();
}

void GearBox::clear()
{
  for(int i=0;i<m_gears.count();i++) delete m_gears[i];
  for(int i=0;i<m_links.count();i++) delete m_links[i];
}

Gear * GearBox::addGear(Gear*gear)
{
  m_gears.append(gear);
  return gear;
}

GearLink*  GearBox::addLink(GearLink*link)
{
  m_links.append(link);
  return link;
}

GearLink *GearBox::addLink(int a, int b)
{
  return GearBox::addLink(new GearLink(m_gears[a], m_gears[b]));
}


void GearBox::draw(QPainter &pa)
{
  for(int i=0;i<m_links.count();i++) m_links[i]->update();
  pa.setPen(Qt::black);
  pa.setBrush(Qt::yellow);
  for(int i=0;i<m_gears.count();i++)
  {
    Gear *gear = m_gears[i];
    pa.save();
    gear->rotoTranslate(pa);
    pa.drawPath(gear->getBodyPath());
    pa.restore();
  }
  pa.setPen(QPen(Qt::red,0,Qt::DotLine));
  pa.setBrush(Qt::NoBrush);
  for(int i=0;i<m_gears.count();i++)
  {
    Gear *gear = m_gears[i];
    pa.save();
    gear->rotoTranslate(pa);
    pa.drawPath(gear->getPitchLinePath());
    pa.restore();
  }

}
