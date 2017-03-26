#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"

namespace {
Gear *makeSimpleGear(int toothCount)
{
  double radius = toothCount * 30 / (2*M_PI);
  Gear *gear = new Gear(new PitchCurve(EllipseFunction(radius,0.0)));
  CircularToothMaker ctm;
  CircularToothMaker::Params params;
  params.pitchRadius = radius;
  params.toothCount = toothCount;
  params.toothHeight = 20;

  QVector<QVector2D> pts;
  ctm.makeTeeth(pts, params);

  gear->setBodyPath(pts);
  return gear;
}

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



class InvoluteGearsPage : public Page, public Pannable {
  Gear *m_gear1, *m_gear2;
public:
  InvoluteGearsPage() : Page("involuteGears") { 
      m_gear1 = makeSimpleGear2(21); 
      m_gear2 = makeSimpleGear2(21); 
      m_gear1->setPosition(600,400);
      m_gear2->setPosition(m_gear1->getPosition() + QPointF(-GearLink(m_gear1,m_gear2).getDistance(),0));

  }
  ~InvoluteGearsPage() { delete m_gear1; delete m_gear2; }

  void draw(QPainter &pa);
} titlePage;


void InvoluteGearsPage::draw(QPainter &pa)
{

  m_gear1->setAngle(getParameter() * 0.01);
  GearLink(m_gear1, m_gear2).update();
  // m_gear1->getAngle() + 0.01);

  pa.setBrush(Qt::yellow);
  pa.setPen(Qt::black);
  
  pa.save();
  pa.translate(m_gear1->getPosition());
  pa.rotate(m_gear1->getAngle()*180.0/M_PI);
  pa.drawPath(m_gear1->getBodyPath());
  pa.restore();

  pa.save();
  pa.translate(m_gear2->getPosition());
  pa.rotate(m_gear2->getAngle()*180.0/M_PI);
  pa.drawPath(m_gear2->getBodyPath());
  pa.restore();

  pa.setPen(QPen(Qt::red, 0, Qt::DotLine));
  pa.setBrush(Qt::NoBrush);

  pa.save();
  pa.translate(m_gear1->getPosition());
  pa.rotate(m_gear1->getAngle()*180.0/M_PI);
  pa.drawPath(m_gear1->getPitchLinePath());
  pa.restore();

  pa.save();
  pa.translate(m_gear2->getPosition());
  pa.rotate(m_gear2->getAngle()*180.0/M_PI);
  pa.drawPath(m_gear2->getPitchLinePath());
  pa.restore();
}


