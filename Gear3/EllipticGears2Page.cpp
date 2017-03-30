#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"

class EllipticGears2Page : public Page, public Pannable {
  GearBox m_gearBox;
public:
  EllipticGears2Page() : Page("ellipticGears2") { 
    Gear *gear1 = m_gearBox.addGear(makeEllipticGear());
    Gear *gear2 = m_gearBox.addGear(makeEllipticGear());
    Gear *gear3 = m_gearBox.addGear(makeEllipticGear());
    Gear *gear4 = m_gearBox.addGear(makeEllipticGear());
    Gear *gear5 = m_gearBox.addGear(makeEllipticGear());
    //Gear *gear6 = m_gearBox.addGear(makeEllipticGear());
    gear1->setPosition(400,0);
    double theta = M_PI*0.8;
    m_gearBox.addLink(0,1)->moveDriven(theta);
    m_gearBox.addLink(1,2)->moveDriven(-theta);
    m_gearBox.addLink(2,3)->moveDriven(theta);
    m_gearBox.addLink(3,4)->moveDriven(-theta);
    //m_gearBox.addLink(4,5)->moveDriven(5*M_PI/3);
/*
    double dist = m_gearBox.getLink(0)->getDistance();
    gear1->setPosition( 1.5*dist,0);
    gear2->setPosition( 0.5*dist,0);
    gear3->setPosition(-0.5*dist,0);
    gear4->setPosition(-1.5*dist,0);
    m_gearBox.getLink(3)->moveDriven(M_PI/3);
  */  
    for(int i=0;i<m_gearBox.getGearCount();i++)
    {
      m_gearBox.getGear(i)->setBrush(Qt::cyan);
    }
    setDefaultScale(0.8);
  }

  ~EllipticGears2Page() {  }

  void draw(QPainter &pa);
} titlePage;


void EllipticGears2Page::draw(QPainter &pa)
{
  setParameter(getParameter() + 0.1*getElapsedTime());
  m_gearBox.getGear(0)->setAngle(getParameter()*0.005);

  m_gearBox.draw(pa);
 

}


