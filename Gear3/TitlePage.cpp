#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"

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



class TitlePage : public Page {
  Gear *m_gear1, *m_gear2;
public:
  TitlePage() : Page("title") { 
      m_gear1 = makeSimpleGear(21); 
      m_gear2 = makeSimpleGear(31); 
      m_gear1->setPosition(600,400);
      m_gear2->setPosition(m_gear1->getPosition() + QPointF(-GearLink(m_gear1,m_gear2).getDistance(),0));

  }
  ~TitlePage() { delete m_gear1; delete m_gear2; }

  void draw(QPainter &pa);
} titlePage;


void TitlePage::draw(QPainter &pa)
{
  QLinearGradient myGradient;
  
  QString texts[2] = { "Weird", "Gears" };
  
  
  
  QFont font("AR HERMANN",200);
  QFontMetrics fm(font);
  

  QPointF p(100,200);
  QPainterPath pp;
  pp.addText(p, font, texts[0]);
  p += QPointF(150,fm.height()-20);
  pp.addText(p, font, texts[1]);
  
  pa.drawPath(pp);
  //pa.fillPath(pp, Qt::yellow);
  //pa.strokePath(pp, QPen(Qt::black, 2));
  

  GearLink(m_gear1, m_gear2).update();
  m_gear1->setAngle(m_gear1->getAngle() + 0.01);

  m_gear1->draw(pa);
  m_gear2->draw(pa);

}


