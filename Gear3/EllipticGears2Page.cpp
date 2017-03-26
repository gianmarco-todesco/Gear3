#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"

Gear *makeEllipticGear2() 
{
  int toothCount = 19;
  double radius = 100;
  Gear *gear = new Gear(new PitchCurve(EllipseFunction(radius,0.6)));
  SimpleToothMaker ctm;
  SimpleToothMaker::Params params;
  params.toothHeight = 20;
  params.toothCount = toothCount;
  params.toothOffset = 0.0;
  
  QVector<QVector2D> pts;
  ctm.makeTeeth(pts, gear->getCurve(), params);

  gear->setBodyPath(pts);
  return gear;
}

class EllipticGears2Page : public Page, public Pannable {
  GearBox m_gearBox;
public:
  EllipticGears2Page() : Page("ellipticGears2") { 
    m_gearBox.addGear(makeEllipticGear2());
    m_gearBox.addGear(makeEllipticGear2());
    m_gearBox.addGear(makeEllipticGear2());
    m_gearBox.addGear(makeEllipticGear2());
    for(int i=1;i<m_gearBox.getGearCount();i++)
    {
      Gear *g1 = m_gearBox.getGear(i-1);
      Gear *g2 = m_gearBox.getGear(i);
      GearLink *link = new GearLink(g1,g2);
      g2->setPosition(g1->getPosition() -  QPointF(link->getDistance(),0));
      m_gearBox.addLink(link);
    }
  }

  ~EllipticGears2Page() {  }

  void draw(QPainter &pa);
} titlePage;


void EllipticGears2Page::draw(QPainter &pa)
{
  m_gearBox.getGear(0)->setAngle(getParameter()*0.01);

  m_gearBox.draw(pa);
  /*

  GearLink(m_gearBox.getGear(0), m_gearBox.getGear(1), M_PI).update();
  m_gearBox.getGear(0)->draw(pa);
  m_gearBox.getGear(1)->draw(pa);
  */

}


