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
  params.toothAddendum = 10;
  params.toothDedendum = 10;

  QVector<QVector2D> pts;
  ctm.makeTeeth(pts, params);

  gear->setBodyPath(pts);
  return gear;
}



class TitlePage : public Page {
  GearBox m_gearBox;
public:
  TitlePage() : Page("title") { 
      Gear *gear1 = m_gearBox.addGear(makeSimpleGear(21)); 
      gear1->setPosition(120,350);
      gear1->setBrush(QColor(100,230,150));
      m_gearBox.addGear(makeSimpleGear(31))->setBrush(QColor(200,100,50)); 
      m_gearBox.addLink(0, 1)->moveDriven(M_PI*0.45);
      m_gearBox.addGear(makeSimpleGear(17))->setBrush(QColor(200,230,50)); 
      m_gearBox.addLink(1, 2)->moveDriven(0.3);
      m_gearBox.addGear(makeSimpleGear(17))->setBrush(QColor(200,230,50)); 
      m_gearBox.addLink(2, 3)->moveDriven(-0.4);
      m_gearBox.addGear(makeSimpleGear(13))->setBrush(QColor(100,130,250)); 
      m_gearBox.addLink(3, 4)->moveDriven(-0.1);

  }
  ~TitlePage() { }

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
  p += QPointF(150,fm.height()-40);
  pp.addText(p, font, texts[1]);
  
  pa.setPen(QPen(Qt::black,1));
  pa.setBrush(Qt::yellow);
  pa.drawPath(pp);
  //pa.fillPath(pp, Qt::yellow);
  //pa.strokePath(pp, QPen(Qt::black, 2));
  
  m_gearBox.getGear(0)->setAngle(getTime() * 0.001);
  m_gearBox.draw(pa);
}


