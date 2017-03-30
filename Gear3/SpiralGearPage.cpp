#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"
#include <QDebug>

#include <QFile>
#include <QDataStream>



namespace {
Gear *makeSpiralGear(int toothCount)
{
  double radius = 200;
  Gear *gear = new Gear(new PitchCurve(SpiralFunction(100,300)));

  /*
  SimpleToothMaker ctm;
  SimpleToothMaker::Params params;
  params.toothHeight = 20;
  params.toothCount = toothCount;
  params.toothOffset = 0.25;

  QVector<QVector2D> pts;
  ctm.makeTeeth(pts, gear->getCurve(), params);
  */
  
  QVector<QVector2D> pts;
  if(false)
  {
    MagicToothMaker ctm;
    QElapsedTimer clock;
    clock.start();
    ctm.makeTeeth(pts, gear->getCurve());
    int time = clock.elapsed();
    qDebug() << "time = " << time;
  
    QFile file("spiralGearTeeth.dat");
    file.open(QIODevice::WriteOnly);
    QDataStream out(&file);
    out << pts;
    file.close();
  }
  else
  {
    QFile file("spiralGearTeeth.dat");
    file.open(QIODevice::ReadOnly);
    QDataStream out(&file);
    out >> pts;
    file.close();
  }

  gear->setBodyPath(pts);

  QPainterPath pp = gear->getBodyPath();
  const PitchCurve *crv = gear->getCurve();
  double s = 10.0;
  while(s + 150.0 < crv->getLength())
  {
    double ss[4] = {s,s+5, s+60-5, s+60};
    QVector2D vv1[4],vv2[4];
    for(int i=0;i<4;i++) { vv1[i] = crv->getPosFromS(ss[i],-20.0); vv2[i] = 50.0 * vv1[i].normalized(); }
    pp.moveTo(vv1[0].toPointF());
    pp.lineTo(vv1[3].toPointF());
    pp.cubicTo(
      (vv1[3]+(vv1[3]-vv1[2])*10).toPointF(), 
      (0.5*(vv1[3]+vv2[3])).toPointF(), 
      vv2[3].toPointF() );

    pp.lineTo(vv2[0].toPointF());
    pp.closeSubpath();
    s += 100.0;
  }
  gear->setBodyPath(pp);
  
  return gear;
}

}



class SpiralGearsPage : public Page, public Pannable {
  Gear *m_gear1, *m_gear2;
  GearLink *m_link;
public:
  SpiralGearsPage() : Page("spiralGear") { 
      m_gear1 = makeSpiralGear(21); 
      m_gear2 = makeSpiralGear(21); 
      m_gear1->setAngle(M_PI - m_gear1->getCurve()->getPhiFromS(m_gear1->getCurve()->getLength()/2));
      m_gear2->setAngle(-m_gear2->getCurve()->getPhiFromS(m_gear2->getCurve()->getLength()/2));
      m_gear1->setPosition(300,100);
      m_link = new GearLink(m_gear1, m_gear2);
      //m_gear2->setPosition(m_gear1->getPosition() + QPointF(-GearLink(m_gear1,m_gear2).getDistance(),0));
      m_gear2->setPosition(m_gear1->getPosition() + QPointF(-400,0));

  }
  ~SpiralGearsPage() { delete m_gear1; delete m_gear2; }

  void draw(QPainter &pa);
} titlePage;


void SpiralGearsPage::draw(QPainter &pa)
{
  double angle = getParameter() * 0.01;
  if(angle>=M_PI) while(angle>=M_PI) angle -= 2*M_PI;
  else if(angle<-M_PI) while(angle<-M_PI) angle += 2*M_PI;
  m_gear1->setAngle( angle );
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


