#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"

Gear *makeEllipticGear() 
{
  int toothCount = 19;
  double radius = 100;
  Gear *gear = new Gear(new PitchCurve(EllipseFunction(radius,0.6)));
  SimpleToothMaker ctm;
  SimpleToothMaker::Params params;
  params.toothHeight = 20;
  params.toothCount = toothCount;
  
  QVector<QVector2D> pts;
  ctm.makeTeeth(pts, gear->getCurve(), params);

  gear->setBodyPath(pts);
  return gear;
}

class EllipticGearPage : public Page, public Pannable {
  Gear *m_gear1, *m_gear2;
public:
  EllipticGearPage() : Page("ellipticGear") { 
    m_gear1 = makeEllipticGear();
    m_gear2 = makeEllipticGear();
    m_gear2->setPosition(m_gear1->getPosition() - QPointF(GearLink(m_gear1,m_gear2).getDistance(),0));
  }

  ~EllipticGearPage() {  }

  void draw(QPainter &pa);
} titlePage;


void EllipticGearPage::draw(QPainter &pa)
{
  m_gear1->setAngle(getParameter()*0.01);
  GearLink(m_gear1,m_gear2).update();
  m_gear1->draw(pa);
  m_gear2->draw(pa);
}


