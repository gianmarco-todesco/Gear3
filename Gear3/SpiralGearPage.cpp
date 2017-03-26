#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"

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



class SpiralGearsPage : public Page, public Pannable {
  Gear *m_gear1, *m_gear2;
  GearLink *m_link;
public:
  SpiralGearsPage() : Page("spiralGear") { 
      m_gear1 = makeSpiralGear(20); 
      m_gear2 = makeSpiralGear(20); 
      m_gear1->setAngle(M_PI - m_gear1->getCurve()->getPhiFromS(m_gear1->getCurve()->getLength()/2));
      m_gear2->setAngle(-m_gear2->getCurve()->getPhiFromS(m_gear2->getCurve()->getLength()/2));
      m_gear1->setPosition(600,400);
      m_link = new GearLink(m_gear1, m_gear2);
      //m_gear2->setPosition(m_gear1->getPosition() + QPointF(-GearLink(m_gear1,m_gear2).getDistance(),0));
      m_gear2->setPosition(m_gear1->getPosition() + QPointF(-800,0));

  }
  ~SpiralGearsPage() { delete m_gear1; delete m_gear2; }

  void draw(QPainter &pa);
} titlePage;


void SpiralGearsPage::draw(QPainter &pa)
{

  m_gear1->setAngle(getParameter() * 0.01 );
  m_link->update();
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


