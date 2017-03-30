#include "Gear.h"
#include "PitchCurve.h"
#include <qmath.h>
#include "ToothMaker.h"
#include <QDebug>
#include <QElapsedTimer>


Gear::Gear(PitchCurve *curve) 
  : m_curve(curve)
  , m_angle(0)
  , m_brush(Qt::white)
{
  updatePitchPath();
}

Gear::~Gear()
{
  delete m_curve;
}


void Gear::draw(QPainter &pa)
{
  pa.save();
  rotoTranslate(pa);

  pa.setPen(QPen(Qt::blue, 2));
  pa.setBrush(m_brush);
  pa.drawPath(m_bodyPath);

  pa.setPen(QPen(Qt::magenta, 0, Qt::DotLine));
  pa.setBrush(Qt::NoBrush);
  pa.drawPath(m_pitchLinePath);

  QPointF p = m_curve->getPoint(0).pos.toPointF();
  pa.setPen(Qt::black);
  pa.drawLine(p*0.9, p*0.5);

  pa.restore();
}

void Gear::rotoTranslate(QPainter &pa)
{
  pa.translate(getPosition());
  pa.rotate(getAngle()*180.0/M_PI);
}


void Gear::setBodyPath(const QVector<QVector2D> &pts)
{
  QVector2D oldp = pts[0];
  m_bodyPath.moveTo(oldp.toPointF());
  for(int i=1;i<pts.count();i++) 
  {
    QVector2D p = pts[i];
    if((p-oldp).lengthSquared()<1) continue;
    oldp = p;
    m_bodyPath.lineTo(p.toPointF());
  }
  m_bodyPath.closeSubpath();
  m_bodyPath.addEllipse(QPointF(0,0), 5,5);
}


void Gear::updatePitchPath()
{
  QVector2D lastPos = m_curve->getPoint(0).pos;
  m_pitchLinePath.moveTo(lastPos.toPointF());
  double ds = 5, s = ds;
  while(s<m_curve->getLength())
  {
    lastPos = m_curve->getPosFromS(s);
    m_pitchLinePath.lineTo(lastPos.toPointF());
    s += ds;
  }
}


//=============================================================================


GearLink::GearLink(Gear*driver, Gear *driven) 
  : m_gear1(driver)
  , m_gear2(driven)
  , m_theta1(driver->getAngle())
  , m_theta2(driven->getAngle())
{ 
  m_distance = computeDistance();
}

GearLink::~GearLink()
{
}

double GearLink::computeDistance() const
{
  double r0 = m_gear1->getCurve()->getRFromPhi(M_PI - m_theta1);
  double r1 = m_gear2->getCurve()->getRFromPhi(-m_theta2);
  return r0 + r1;
}


void GearLink::update()
{
  const PitchCurve *crv1 = m_gear1->getCurve();
  const PitchCurve *crv2 = m_gear2->getCurve();
  
  QPointF d = m_gear2->getPosition() - m_gear1->getPosition();
  double psi = atan2(d.y(),d.x()) - M_PI;

  double angle1 = m_gear1->getAngle() - psi;

  /*
  QPointF p = m_driven->getPosition() - m_driver->getPosition();
  double psi = atan2(p.y(),p.x()); if(psi<0)psi+=2*M_PI;
  double s = driverCrv->getSFromPhi(M_PI*0) - driverCrv->getSFromPhi(psi-m_driver->getAngle());
  m_driven->setAngle(drivenCrv->getPhiFromS(-s) + psi - M_PI); 
  */

  double s1a = crv1->getSFromPhi(M_PI-m_theta1);
  double s1b = crv1->getSFromPhi(M_PI-angle1);

  double s1 = s1a-s1b;

  double s2a = crv2->getSFromPhi(-m_theta2) ;
  double s = s1 - s2a;

  double angle2 = - crv2->getPhiFromS( 
       crv1->getSFromPhi(M_PI-m_theta1) 
     - crv1->getSFromPhi(M_PI-angle1) 
     + crv2->getSFromPhi(-m_theta2) );
  m_gear2->setAngle(angle2 + psi);

}


void GearLink::moveDriven(double psi)
{
  m_gear2->setPosition(m_gear1->getPosition() + m_distance * QPointF(cos(psi), sin(psi)));
}

//=============================================================================

GearBox::GearBox()
{
}

GearBox::~GearBox()
{
  clear();
}

void GearBox::clear()
{
  for(int i=0;i<m_gears.count();i++) delete m_gears[i];
  for(int i=0;i<m_links.count();i++) delete m_links[i];
}

Gear * GearBox::addGear(Gear*gear)
{
  m_gears.append(gear);
  return gear;
}

GearLink*  GearBox::addLink(GearLink*link)
{
  m_links.append(link);
  return link;
}

void GearBox::draw(QPainter &pa)
{
  for(int i=0;i<m_links.count();i++) m_links[i]->update();
  pa.setPen(Qt::black);
  for(int i=0;i<m_gears.count();i++)
  {
    Gear *gear = m_gears[i];
    pa.save();
    gear->rotoTranslate(pa);
    pa.setBrush(gear->getBrush());
    pa.drawPath(gear->getBodyPath());
    pa.restore();
  }
  pa.setPen(QPen(Qt::red,0,Qt::DotLine));
  pa.setBrush(Qt::NoBrush);
  for(int i=0;i<m_gears.count();i++)
  {
    Gear *gear = m_gears[i];
    pa.save();
    gear->rotoTranslate(pa);
    pa.drawPath(gear->getPitchLinePath());
    pa.restore();
  }

}



//=============================================================================


Gear *makeCircularGear(int toothLength, int toothCount, int flag) 
{
  double radius = (toothLength * toothCount) / (2*M_PI);
  Gear *gear = new Gear(new PitchCurve(EllipseFunction(radius,0.0),200));

  
  QVector<QVector2D> pts;
  if(flag == 0)
  {
    SimpleToothMaker::Params params;
    params.toothHeight = toothLength*0.75*0.5;
    params.toothCount = toothCount;
    SimpleToothMaker().makeTeeth(pts, gear->getCurve(), params);
  }
  else if(flag==1)
  {
    SquareToothMaker::Params params;
    params.toothHeight = toothLength*0.75*0.5 * 1.3;
    params.toothCount = toothCount;
    SquareToothMaker().makeTeeth(pts, gear->getCurve(), params);
  }

  gear->setBodyPath(pts);
  
  return gear;

}


//=============================================================================


class SelfMatchingGearBuilder {
public:
  struct PolarPoint {
    double phi, radius;
    PolarPoint() : phi(0), radius(0) {}
    PolarPoint(double _phi, double _radius) : phi(_phi), radius(_radius) {}
  };

  double getDeltaPhi(double r0, double r1, double ds) const { return acos((r0*r0+r1*r1-ds*ds)/(2.0*r0*r1)); }

  double foobar(QList<PolarPoint> &polarPoints, double x, double y, int m);

  void makeSeed1(QList<PolarPoint> &polarPoints, int &off);
  void makeSeed2(QList<PolarPoint> &polarPoints, int &off);
  void duplicate(QList<PolarPoint> &polarPoints);
  void adjust(QList<PolarPoint> &polarPoints, int off);

  void createPitchCurvePoints(QVector<PitchCurve::Point> &pts, const QList<PolarPoint> &polarPoints);
  PitchCurve *createPitchCurve(const QList<PolarPoint> &polarPoints);

  PitchCurve *buildPitchCurve1();
  PitchCurve *buildPitchCurve2();

  Gear *makeGear();
};


void SelfMatchingGearBuilder::makeSeed1(QList<PolarPoint> &polarPoints, int &off)
{
  double dist = 300.0;
  double r = dist*0.5;
  
  polarPoints.clear();
  polarPoints.append(PolarPoint(0.0,r));
  double phi0 = 0.0, phi1 = 0.0;
  double r0 = r, r1 = r;
  double ds = 1;
  off = 0;
  for(;;)
  {
    double ra0 = r0 - 0.2;
    double ra1 = dist - ra0;
    phi0 -= acos((ra0*ra0+r0*r0-ds*ds)/(2.0*ra0*r0));
    phi1 += acos((ra1*ra1+r1*r1-ds*ds)/(2.0*ra1*r1));
    polarPoints.push_back(PolarPoint(phi1,ra1));
    polarPoints.push_front(PolarPoint(phi0,ra0));
    off++;
    if(phi1-phi0>=M_PI/2) break;
    r0=ra0;
    r1=ra1;
  }
  for(int i=1;i<polarPoints.count();i++) 
  { 
    Q_ASSERT(polarPoints[i-1].phi < polarPoints[i].phi); 
  }
}



PitchCurve *SelfMatchingGearBuilder::buildPitchCurve1()
{
  int off;
  QList<PolarPoint> polarPoints;
  makeSeed1(polarPoints, off);
  duplicate(polarPoints);
  duplicate(polarPoints);
  adjust(polarPoints, off);

  return createPitchCurve(polarPoints);
}

double SelfMatchingGearBuilder::foobar(QList<PolarPoint> &polarPoints, double x, double ymax, int m)
{
  double ds = ymax/(m-1);
  polarPoints.clear();
  for(int i=0;i<m;i++)
  {
    double y = ds * i;
    PolarPoint pt(atan2(y,x), sqrt(x*x+y*y));
    polarPoints.append(pt);
  }
  double dist = 2.0*polarPoints.back().radius;
  PolarPoint pt;
  pt.phi = polarPoints.back().phi;
  for(int i=m-2;i>=0;i--)
  {
    pt.radius = dist - polarPoints[i].radius;
    pt.phi += getDeltaPhi(polarPoints.back().radius, pt.radius, ds);
    polarPoints.append(pt);
  }
  return polarPoints.back().phi;
}

void SelfMatchingGearBuilder::makeSeed2(QList<PolarPoint> &polarPoints, int &off)
{
  QElapsedTimer clock;
  clock.start();
  double theta = M_PI/4;
  double x = 200;
  double y = 10;
  QList<PolarPoint> pts;
  double phi = foobar(pts, x,y,(int)y);
  Q_ASSERT(phi <= theta);
  double y1, phi1;
  for(;;)
  {
    y1 = y + 5; 
    phi1 = foobar(pts, x,y1,(int)y1);
    if(phi1>theta) break;
    y=y1; phi=phi1;
  }
  Q_ASSERT(phi<=theta && theta<phi1);
  double t = (theta-phi)/(phi1-phi);
  y = y*(1-t)+y1*t;
  phi = foobar(pts, x,y,(int)y);
  off=pts.count()/2;

  double err = fabs(phi - theta);
  qDebug() << err;
  qDebug() << "time " << clock.elapsed();
  polarPoints.swap(pts);

  for(int i=1;i<polarPoints.count();i++) 
  { 
    Q_ASSERT(polarPoints[i-1].phi < polarPoints[i].phi); 
  }
}


PitchCurve *SelfMatchingGearBuilder::buildPitchCurve2()
{
  int off;
  QList<PolarPoint> polarPoints;
  makeSeed2(polarPoints, off);
  duplicate(polarPoints);
  duplicate(polarPoints);
  duplicate(polarPoints);
  adjust(polarPoints, off);

  return createPitchCurve(polarPoints);
}


void SelfMatchingGearBuilder::duplicate(QList<PolarPoint> &polarPoints)
{
  int n = polarPoints.count();
  double phi0 = polarPoints.back().phi;
  for(int i=n-2;i>=0;i--) 
  {
    double dphi = polarPoints[n-1].phi-polarPoints[i].phi;
    Q_ASSERT(dphi>0.0);
    double r = polarPoints[i].radius;
    polarPoints.append( PolarPoint(phi0+dphi, r));
    Q_ASSERT(polarPoints.back().phi > polarPoints[polarPoints.count()-2].phi);    
  }
  for(int i=1;i<polarPoints.count();i++) 
  { 
    Q_ASSERT(polarPoints[i-1].phi < polarPoints[i].phi); 
  }
}


void SelfMatchingGearBuilder::adjust(QList<PolarPoint> &polarPoints, int off)
{
  // remove extra points
  while(polarPoints.back().phi > polarPoints[0].phi + 2*M_PI) 
    polarPoints.pop_back();

  // move the first off point to the end
  for(int i=0;i<off;i++) 
  {
    polarPoints[0].phi += 2*M_PI; 
    Q_ASSERT(polarPoints[0].phi>polarPoints.back().phi);
    polarPoints.push_back(polarPoints[0]); 
    polarPoints.pop_front(); 
  }

  // rotating points ( polarPoints[0].phi must be 0 )
  double dphi = -polarPoints[0].phi;
  for(int i=0;i<polarPoints.count();i++) polarPoints[i].phi += dphi;

  // shouldn't need this
  while(polarPoints.back().phi>=2*M_PI) 
    polarPoints.pop_back();
}

void SelfMatchingGearBuilder::createPitchCurvePoints(QVector<PitchCurve::Point> &pts, const QList<PolarPoint> &polarPoints)
{
  pts.resize(polarPoints.count());
  for(int i=0;i<polarPoints.count();i++)
  {
    PitchCurve::Point pt;
    pt.phi = polarPoints[i].phi;
    pt.r = polarPoints[i].radius;
    pt.pos = pt.r * QVector2D(cos(pt.phi), sin(pt.phi));
    pt.s = i==0 ? 0.0 : pts[i-1].s + (pts[i-1].pos-pt.pos).length();
    if(i>0)
    {
      Q_ASSERT(pt.phi>pts.back().phi);
    }
    pts[i] = pt;
  }
  // build 'right' vectors
  for(int i=0;i<pts.count();i++)
  {
    QVector2D p1 = pts[(i+1)%pts.count()].pos;
    QVector2D p0 = pts[(i+pts.count()-1)%pts.count()].pos;
    QVector2D v = (p1-p0).normalized();
    pts[i].right = QVector2D(-v.y(),v.x());
  }
}

PitchCurve *SelfMatchingGearBuilder::createPitchCurve(const QList<PolarPoint> &polarPoints)
{
  QVector<PitchCurve::Point> pts;
  createPitchCurvePoints(pts, polarPoints);
  return new PitchCurve(pts);
}



Gear *SelfMatchingGearBuilder::makeGear()
{
  Gear *gear = new Gear(buildPitchCurve2());
  SimpleToothMaker::Params params;
  params.toothCount = 70;
  params.toothHeight = 12;
  params.toothOffset = 0.25;
  QVector<QVector2D> teethPts;
  SimpleToothMaker().makeTeeth(teethPts, gear->getCurve(), params);
  gear->setBodyPath(teethPts);
  return gear;
}


Gear *makeSquareSelfMatchingGear()
{
  return SelfMatchingGearBuilder().makeGear();
}




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