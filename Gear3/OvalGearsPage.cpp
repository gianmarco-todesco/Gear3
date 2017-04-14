#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"

#include <QElapsedTimer>
#include <QDebug>

class OvalGearsPage : public Page, public Pannable {
  GearBox m_gearBox;
public:
  OvalGearsPage() : Page("ovalGears") { 
    setDefaultScale(0.8);
  }

  void build() {
    QElapsedTimer timer;
    timer.start();
    Gear*gear1 = m_gearBox.addGear(makeOvalSelfMatchingGear());
    Gear*gear2 = m_gearBox.addGear(makeOvalSelfMatchingGear());
    gear1->setBrush(Qt::yellow);
    gear2->setBrush(Qt::cyan);
    
    GearLink *link = m_gearBox.addLink(0,1);
    double dist = link->getDistance();
    gear1->setPosition( dist/2,0);
    gear2->setPosition(-dist/2,0);
    qDebug() << "oval gears done: " << timer.elapsed() << "ms";
  }

  ~OvalGearsPage() {  }

  void draw(QPainter &pa);
} ovalGears;


void OvalGearsPage::draw(QPainter &pa)
{
  if(m_gearBox.getGearCount()==0) build();
  m_gearBox.getGear(0)->setAngle(getParameter()*0.01);
  m_gearBox.draw(pa);
}


