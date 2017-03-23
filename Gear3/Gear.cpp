#include "Gear.h"
#include "PitchCurve.h"
#include <qmath.h>


Gear::Gear(PitchCurve *curve) 
  : m_curve(curve)
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
  pa.translate(getPosition());
  pa.rotate(getAngle()*180.0/M_PI);

  pa.setPen(QPen(Qt::blue, 2));
  pa.setBrush(Qt::cyan);
  pa.drawPath(m_bodyPath);

  pa.setPen(QPen(Qt::magenta, 0, Qt::DotLine));
  pa.setBrush(Qt::NoBrush);
  pa.drawPath(m_pitchLinePath);
  pa.restore();
}


void Gear::setBodyPath(const QVector<QVector2D> &pts)
{
  m_bodyPath.moveTo(pts[0].toPointF());
  for(int i=1;i<pts.count();i++)
    m_bodyPath.lineTo(pts[i].toPointF());
  m_bodyPath.closeSubpath();
  m_bodyPath.addEllipse(QPointF(0,0), 3,3);
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


GearLink::GearLink(Gear*driver, Gear *driven) : m_driver(driver), m_driven(driven) 
{
}

GearLink::~GearLink()
{
}

void GearLink::update()
{
  const PitchCurve *driverCrv = m_driver->getCurve();
  const PitchCurve *drivenCrv = m_driven->getCurve();
  
  /*
  QPointF p = m_driven->getPosition() - m_driver->getPosition();
  double psi = atan2(p.y(),p.x()); if(psi<0)psi+=2*M_PI;
  double s = driverCrv->getSFromPhi(M_PI*0) - driverCrv->getSFromPhi(psi-m_driver->getAngle());
  m_driven->setAngle(drivenCrv->getPhiFromS(-s) + psi - M_PI); 
  */

  double driverPhi0 = M_PI;
  double drivenPhi0 = 0;

  double driverAngle = m_driver->getAngle();
  double s = driverCrv->getSFromPhi(M_PI - driverPhi0) - driverCrv->getSFromPhi(M_PI - driverAngle);

  double q = drivenCrv->getSFromPhi(- drivenPhi0);
  double drivenAngle = drivenCrv->getPhiFromS(drivenCrv->getLength() + s + q);
  m_driven->setAngle(-drivenAngle);

}


