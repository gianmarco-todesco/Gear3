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

class SquareGearsPage : public Page, public Pannable {
  GearBox m_gearBox;
  bool m_spinning;
public:
  SquareGearsPage() : Page("squareGears"), m_spinning(false) { 
  }

  ~SquareGearsPage() {  }

  bool onKey(int key) {
    if(key == Qt::Key_P) { m_spinning = !m_spinning; }
    else return false;
    return true;
  }

  void build() {
    QElapsedTimer timer; timer.start();
    Gear*gear1 = m_gearBox.addGear(makeSquareSelfMatchingGear());
    Gear*gear2 = m_gearBox.addGear(makeSquareSelfMatchingGear());
    gear1->setBrush(Qt::yellow);
    gear2->setBrush(Qt::cyan);
    

    GearLink *link = m_gearBox.addLink(0,1);
    double dist = link->getDistance();
    gear1->setPosition( dist/2,0);
    gear2->setPosition(-dist/2,0);
    qDebug() << "Square gear built: " << timer.elapsed() << "ms";
  }

  void rotateFirstGear(double delta) {
    Gear *gear = m_gearBox.getGear(0);
    double angle = gear->getAngle();
    angle += delta;
    if(angle>2*M_PI) while(angle>2*M_PI) angle -= 2*M_PI;
    else if(angle<0) while(angle<0) angle += 2*M_PI;
    gear->setAngle(angle);
  }

  void drag(int dx, int dy, int modifiers) {
    rotateFirstGear(dx*0.01);
  }

  void draw(QPainter &pa);
} squareGears;


void SquareGearsPage::draw(QPainter &pa)
{
  if(m_gearBox.getGearCount()==0) build();

  if(m_spinning) rotateFirstGear(0.001*getElapsedTime());
  m_gearBox.getLink(0)->update();

  for(int i=0;i<m_gearBox.getGearCount();i++)
  {
    Gear *gear = m_gearBox.getGear(i);
    pa.save();
    gear->rotoTranslate(pa);
    const QPainterPath &pp = gear->getBodyPath();
    pa.setPen(QPen(Qt::black, 2));
    pa.setBrush(gear->getBrush());
    pa.drawPath(pp);
    pa.setBrush(Qt::NoBrush);
    pa.setPen(Qt::black);
    pa.rotate(-25);
    double r = 180;
    pa.drawRect(-r,-r,2*r,2*r);
    r-=5;
    pa.drawRect(-r,-r,2*r,2*r);
    pa.restore();
  }
}


