#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"



class EllipticGearPage : public Page, public Pannable {
  GearBox m_gearBox;
public:
  EllipticGearPage() : Page("ellipticGear") { 
    Gear*gear1 = m_gearBox.addGear(makeEllipticGear());
    Gear*gear2 = m_gearBox.addGear(makeEllipticGear());
    GearLink *link = m_gearBox.addLink(0,1);
    double x = link->getDistance()*0.5;
    gear1->setPosition(x,0);
    gear2->setPosition(-x,0);
    gear1->setBrush(Qt::cyan);
    gear2->setBrush(Qt::green);

  }

  ~EllipticGearPage() {  }

  void draw(QPainter &pa);
} titlePage;


void EllipticGearPage::draw(QPainter &pa)
{
  m_gearBox.getGear(0)->setAngle(getParameter()*0.01);
  m_gearBox.draw(pa);
}


