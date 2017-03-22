#include "Gear.h"
#include "PitchCurve.h"


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
  pa.setPen(QPen(Qt::blue, 2));
  pa.setBrush(Qt::cyan);
  pa.drawPath(m_bodyPath);

  pa.setPen(QPen(Qt::magenta, 0, Qt::DotLine));
  pa.setBrush(Qt::NoBrush);
  pa.drawPath(m_pitchLinePath);
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

