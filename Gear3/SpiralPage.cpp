#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>
#include <QVector2D>
#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"
#include <QMatrix>
#include <QDebug>

double atanp(const QPointF &p)
{
  return atan2(p.y(),p.x());
}

QPointF normalize(const QPointF &p) { return p * (1.0/sqrt(p.x()*p.x()+p.y()*p.y())); }

void drawAngle(QPainter &pa, QPointF p0, QPointF p1, QPointF p2)
{
  double r = QVector2D(p1-p0).length();
  double phi = -180*atanp(p1-p0)/M_PI;
  QVector2D v1(p1-p0); v1 *= (1.0/r);
  QVector2D v2(p2-p0); v2.normalize();

  double dphi = acos(QVector2D::dotProduct(v1,v2)) * 180.0/M_PI;
  QRectF rect(p0.x()-r,p0.y()-r,2*r,2*r);
  QPainterPath pp;
  pp.moveTo(p0);
  pp.arcTo(rect, phi, dphi);
  pp.closeSubpath();
  QPen oldPen = pa.pen();
  pa.setPen(Qt::NoPen);
  pa.drawPath(pp);
  pa.setPen(oldPen);
  pp=QPainterPath();
  pp.arcMoveTo(rect,phi);
  pp.arcTo(rect, phi, dphi);
  QBrush oldBrush = pa.brush();
  pa.setBrush(Qt::NoBrush);
  pa.drawPath(pp);
  pa.setBrush(oldBrush);
}



class SpiralPage : public Page, public Pannable {
  Gear* m_gear1, *m_gear2;
  double m_par1, m_par2;
  int m_mode;
  double m_dist;

public:
  SpiralPage() : Page("spiral"), m_par1(0), m_par2(0), m_mode(5), m_dist(400) {
    buildGears();
  }

  ~SpiralPage() {  
    delete m_gear1;
    delete m_gear2;
  }

  bool onKey(int key) {
    if(key == Qt::Key_1) { m_mode = 1; }
    else if(key == Qt::Key_2) m_mode = 2;
    else if(key == Qt::Key_3) m_mode = 3;
    else if(key == Qt::Key_4) m_mode = 4;
    else if(key == Qt::Key_5) m_mode = 5;
    else if(key == Qt::Key_6) m_mode = 6;
    else return false;
    return true;
  }

  void drag(int dx, int dy, int modifiers) {
    m_par1 += dx*0.01;
    m_par2 += dy*0.01;
  }

  void buildGears() {
    m_gear1 = makeSpiral1();
    m_gear2 = makeSpiral2();
  }

  void makeOutline(Gear*gear);
  Gear *makeSpiral1();
  Gear *makeSpiral2();

  void draw(QPainter &pa);
  void draw1(QPainter &pa);
  void draw2(QPainter &pa);
  void draw3(QPainter &pa);

  void drawUnfolding(QPainter &pa);

  void drawRuler(QPainter &pa, const PitchCurve *curve);

} spiralPage;

void SpiralPage::makeOutline(Gear*gear)
{
  const PitchCurve *crv = gear->getCurve();
  QPainterPath pp;
  pp.moveTo(crv->getPoint(0).pos.toPointF());
  for(int i=1;i<crv->getPointCount();i++)
    pp.lineTo(crv->getPoint(i).pos.toPointF());
  pp.closeSubpath();
  double r = 15;
  pp.addEllipse(QPointF(0,0), r,r);
  gear->setBodyPath(pp);
}


Gear *SpiralPage::makeSpiral1()
{
  Gear *gear = new Gear(new PitchCurve(ArchimedeanSpiralFunction(100,300), 200));
  makeOutline(gear);
  return gear;
}

Gear *SpiralPage::makeSpiral2()
{
  Gear *gear = new Gear(new PitchCurve(SpiralFunction(100,300), 200));
  makeOutline(gear);
  return gear;
}

void SpiralPage::drawRuler(QPainter &pa, const PitchCurve *curve)
{
  int n = curve->getPointCount();
  QPainterPath pp;
  pp.moveTo(curve->getPoint(0).pos.toPointF());
  for(int i=1;i<n;i++) 
    pp.lineTo(curve->getPoint(i).pos.toPointF());
  for(int i=n-1;i>=0;i--) 
    pp.lineTo((curve->getPoint(i).pos - curve->getPoint(i).right * 30).toPointF());
  pp.closeSubpath();
  pa.setBrush(QColor(200,200,200,100));
  pa.setPen(Qt::black);
  pa.drawPath(pp);

  pp = QPainterPath();  
  double sunit = curve->getLength()/100;
  int k = 1;
  double s = sunit;
  while(s+sunit<=curve->getLength())
  {
    double h = k%10 == 0 ? 14 : 7;
    pp.moveTo(curve->getPosFromS(s,0).toPointF());
    pp.lineTo(curve->getPosFromS(s,-h).toPointF());
    if(k==50)
    {
      pp.addEllipse(curve->getPosFromS(s,-15).toPointF(),5,5);
    }
    s += sunit;
    k++;
  }
  pa.drawPath(pp);
}

double normalizeAngle(double angle)
{
  if(angle>=2*M_PI) while(angle>=2*M_PI) angle -= 2*M_PI;
  else if(angle<0) while(angle<0.0) angle += 2*M_PI;
  return angle;
}

void SpiralPage::draw(QPainter &pa)
{
  if(m_mode<=3) draw1(pa);
  else if(m_mode==4) draw2(pa);
  else if(m_mode==5) draw3(pa);
}


void SpiralPage::draw1(QPainter &pa)
{
  Gear *gear = (m_mode<=2 ? m_gear1 : m_gear2);
  const PitchCurve *crv = gear->getCurve(); 
  double halfLength = crv->getLength()*0.5;
  double dist = crv->getRFromS(halfLength);
  double phi1, phi2;

  phi1 = m_par1-crv->getPhiFromS(halfLength) + M_PI;
  if(m_mode==1)
  {
    phi2 = -m_par1-crv->getPhiFromS(halfLength);
  }
  else
  {
    double s = crv->getSFromPhi(normalizeAngle(M_PI - phi1) );
    qDebug() << s;
    phi2 = -crv->getPhiFromS(crv->getLength()-s);
  }

  pa.save();
  pa.translate(dist, 0);
  pa.rotate(phi1*180/M_PI);
  pa.setPen(Qt::black);
  pa.setBrush(Qt::yellow);
  pa.drawPath(gear->getBodyPath());
  drawRuler(pa, gear->getCurve());
  pa.restore();

  pa.save();
  pa.translate(-dist, 0);
  pa.rotate(180*phi2/M_PI);
    
  pa.setPen(Qt::black);
  pa.setBrush(Qt::cyan);
  pa.drawPath(gear->getBodyPath());
  drawRuler(pa, gear->getCurve());
  pa.restore();
  pa.setPen(Qt::black);
  pa.setBrush(Qt::NoBrush);
  pa.drawLine(-dist,0, dist,0);
}

void SpiralPage::draw2(QPainter &pa)
{
  Gear *gear = m_gear2;
  const PitchCurve *crv = gear->getCurve(); 
  double halfLength = crv->getLength()*0.5;
  double dist = crv->getRFromS(halfLength);
  double phi1, phi2;

  phi1 = -crv->getPhiFromS(halfLength) + M_PI;

  pa.save();
  pa.translate(dist, 0);
  pa.rotate(phi1*180/M_PI);
  pa.setPen(Qt::black);
  pa.setBrush(Qt::yellow);
  pa.drawPath(gear->getBodyPath());

  QPointF center(dist,0);

  phi2 = m_par1 * M_PI;
  if(phi2<0)phi2=0; else if(phi2>2*M_PI)phi2=2*M_PI;

  PitchCurve::Point pt = crv->getPointFromS(crv->getSFromPhi(phi2));

  pa.setPen(Qt::black);
  pa.setBrush(Qt::NoBrush);
  QPointF p = pt.pos.toPointF();
  pa.drawLine(QPointF(0,0), p*2);

  QPointF v(-pt.right.y(),pt.right.x());

  pa.drawLine(p - v*200, p + v*200);

  pa.setBrush(QColor(200,150,20,127));
  pa.setPen(QPen(Qt::black,1.5));

  drawAngle(pa, p, p + v*100, p*2);

  pa.setPen(Qt::black);
  pa.setBrush(Qt::black);
  pa.drawEllipse(p,3,3);
  pa.setBrush(Qt::NoBrush);
  pa.restore();

}

void SpiralPage::draw3(QPainter &pa)
{
  Gear *gear = m_gear2;
  const PitchCurve *crv = gear->getCurve(); 
  double halfLength = crv->getLength()*0.5;
  double dist = crv->getRFromS(halfLength);
  
  pa.save();
  pa.translate(dist,0);
  drawUnfolding(pa);
  pa.restore();


}

void SpiralPage::drawUnfolding(QPainter &pa)
{
  double parameter = m_par2;
  if(parameter<0)parameter=0; else if(parameter>1)parameter=1;
  const PitchCurve *crv = m_gear2->getCurve();
  double s0 = crv->getLength()*0.5;
  double ds = 5;
  QList<QPointF> pts;
  double s = s0;
  while(s<=crv->getLength())
  {
    pts.push_back(crv->getPosFromS(s).toPointF());
    s+=ds;
  }
  int i0 = 0;
  s = s0;
  while(s>=0.0)
  {
    pts.push_front(crv->getPosFromS(s).toPointF());
    s-=ds;
    i0++;
  }

  QVector<QPointF> pts2(pts.count());
  QVector<QPointF> pts3(pts.count());


  
  double pegRadius = 10.0;
  double direction = atanp(pts[i0+1]-pts[i0-1]);


  QMatrix matrix;
  double oldTheta = direction;
  pts2[i0] = pts[i0];
  pts3[i0] = normalize(pts[i0])*pegRadius;

  for(int i=i0+1;i<pts.count();i++)
  {
    double theta = atanp(pts[i]-pts[i-1]);
    double rot = (oldTheta-theta);
    if(rot>M_PI) rot-=2*M_PI;
    else if(rot<-M_PI) rot += 2*M_PI;
    QMatrix rotMatrix;
    rotMatrix.translate(pts[i-1].x(),pts[i-1].y());
    rotMatrix.rotate(180.0*rot*parameter/M_PI);
    rotMatrix.translate(-pts[i-1].x(),-pts[i-1].y());
    matrix = rotMatrix * matrix;
    pts2[i] = matrix.map(pts[i]);
    oldTheta = theta;
    pts3[i] = matrix.map(normalize(pts2[i])*pegRadius);
  }
  matrix = QMatrix();
  oldTheta = direction;
  for(int i=i0-1;i>=0;i--)
  {
    QPointF oldp = pts[i+1];
    double theta = atanp(oldp-pts[i]);
    double rot = (oldTheta-theta);
    if(rot>M_PI) rot-=2*M_PI;
    else if(rot<-M_PI) rot += 2*M_PI;
    QMatrix rotMatrix;
    rotMatrix.translate(oldp.x(),oldp.y());
    rotMatrix.rotate(180.0*rot*parameter/M_PI);
    rotMatrix.translate(-oldp.x(),-oldp.y());
    matrix = rotMatrix * matrix;
    pts2[i] = matrix.map(pts[i]);
    oldTheta = theta;
    pts3[i] = matrix.map(normalize(pts2[i])*pegRadius);
  }


  QPainterPath pp;

  pp.moveTo(pts[0]);
  for(int i=1;i<pts.count();i++) pp.lineTo(pts[i]);
  pa.setPen(QPen(Qt::red,3));
  pa.drawPath(pp);

  pp = QPainterPath();
  pp.moveTo(pts2[0]);
  for(int i=1;i+1<pts.count();i++) pp.lineTo(pts2[i]);
  pa.setPen(QPen(Qt::black,3));
  pa.drawPath(pp);

  pp = QPainterPath();
  int m = 10;
  int k = i0%m;
  if(k>0) k-=m;
  for(;k+1<pts.count();k+=m)
  {
    int a = qMax(0,k), b = qMin(pts.count()-1,k+m);
    QPointF pa = pts3[a], pb = pts3[b];
    QPointF pab = 0.5*(pa+pb);

    double t = parameter * (1.0 - 10.0/QVector2D(pa-pb).length());
    pp.moveTo(pa*(1-t) + pab*t);
    for(int i=a;i<=b;i++) pp.lineTo(pts2[i]);
    pp.lineTo(pb*(1-t) + pab*t);
    pp.closeSubpath();
  }
  if(parameter==0.0) 
  {
    pa.setPen(Qt::NoPen);
    pa.setBrush(Qt::yellow);
    pa.drawPath(pp);
  }
  else
  {
    pa.setPen(Qt::black);
    pa.setBrush(Qt::yellow);
    pa.drawPath(pp);
  }  

}



#ifdef DOPO

class SpiralPage : public Page, public Pannable {
  QVector<QPointF> m_pts;
  double m_par1, m_par2;
  int m_mode;

public:
  SpiralPage() : Page("spiral"), m_par1(0.5), m_par2(0), m_mode(1) {
    buildPoints();
  }

  ~SpiralPage() {  }

  void buildPoints();

  void drawFoldingSpiral(QPainter &pa, int k, double param);
  QMatrix getMatrixAt(int k);

  void draw(QPainter &pa);
  void draw1(QPainter &pa);
  void drawSpiral(QPainter &pa);

  void drawArchimedeanSpirals(QPainter &pa);
  void drawArchimedeanSpiral(QPainter &pa, const QBrush &brush);

  void drag(int dx, int dy, int modifiers) {
    m_par1 += 0.003*dx; if(m_par1<0)m_par1=0; else if(m_par1>1)m_par1=1.0;
    m_par2 += 0.01*dy; if(m_par2<0)m_par2=0; else if(m_par2>1)m_par2=1.0;

  }

  bool onKey(int key) {
    if(key == Qt::Key_1) { m_par1 = 0.5; m_mode = 1; }
    else if(key == Qt::Key_2) m_mode = 2;
    else if(key == Qt::Key_3) m_mode = 3;
    else if(key == Qt::Key_4) m_mode = 4;
    else return false;
    return true;
  }

} spiralPage;

void SpiralPage::buildPoints()
{
  PitchCurve crv(SpiralFunction(100,300));
  int n = 250;
  for(int i=0;i<n;i++)
  {
    double s = crv.getLength()*i/(n-1);
    m_pts.append(crv.getPosFromS(s).toPointF());
  }
}

void SpiralPage::draw(QPainter &pa)
{
  if(m_mode==1 || m_mode == 2)
  {
    drawArchimedeanSpirals(pa);
    return;
  }
  else
  {
    /*
    drawSpiral(pa);
    pa.save();

    pa.restore();
    */
  }
 

  int k = (int)(m_pts.count()*m_par1);
  drawSpiral(pa);

  pa.save();
  pa.setMatrix(getMatrixAt(k), true);

  pa.save();
  pa.rotate(180);
  pa.setMatrix(getMatrixAt(m_pts.count()-2-k).inverted(), true);
  drawSpiral(pa);
  //if(m_par2>0) 
  //  drawFoldingSpiral(pa, m_pts.count()-2-k, m_par2);
  pa.restore();


  pa.restore();

  /*
  draw1(pa);
  int k = (int)(1 * (1-m_par1) + (m_pts.count()-2) * m_par1);
  QPointF p = m_pts[k];
  QPointF dir = m_pts[k+1]-p;
  double phi = atan2(dir.y(),dir.x());
  pa.save();
  pa.translate(p.x(),p.y());
  pa.rotate(180);
  pa.translate(-p.x(),-p.y());

  draw1(pa);
  pa.restore();
  */

}

void SpiralPage::drawArchimedeanSpirals(QPainter &pa)
{
  double angle = -180 + m_par1*360;
  pa.save();
  pa.translate(200,0);
  pa.rotate(angle);
  drawArchimedeanSpiral(pa, Qt::yellow);
  pa.restore();

  pa.save();
  pa.translate(-200,0);
  pa.rotate(180-angle);
  drawArchimedeanSpiral(pa, Qt::cyan);
  pa.restore();

}

void SpiralPage::drawArchimedeanSpiral(QPainter &pa, const QBrush &brush)
{
  PitchCurve crv(ArchimedeanSpiralFunction(100,300), 200);
  QPainterPath pp;
  pp.moveTo(crv.getPoint(0).pos.toPointF());
  for(int i=1;i<crv.getPointCount();i++)
    pp.lineTo(crv.getPoint(i).pos.toPointF());
  pp.closeSubpath();
  pp.addEllipse(QPointF(0,0), 10,10);
  pa.setPen(QPen(Qt::black, 1.5));
  pa.setBrush(brush);
  pa.drawPath(pp);

  if(m_mode==1) return;
  int n = crv.getPointCount();
  pp = QPainterPath();
  pp.moveTo(crv.getPoint(0).pos.toPointF());
  for(int i=1;i<n;i++) 
    pp.lineTo(crv.getPoint(i).pos.toPointF());
  for(int i=n-1;i>=0;i--) 
    pp.lineTo((crv.getPoint(i).pos - crv.getPoint(i).right * 30).toPointF());
  pp.closeSubpath();
  pa.setBrush(QColor(200,200,200,100));
  pa.setPen(Qt::black);
  pa.drawPath(pp);

  pp = QPainterPath();  
  double s1 = crv.getSFromPhi(M_PI);
  double sunit = 10.0;
  int k = floor(s1/sunit);
  double s0 = s1 - sunit * k;
  double s = s0;
  while(s+10<crv.getLength())
  {
    double h = ((k+1000)%10) == 0 ? 20 : 10;
    pp.moveTo(crv.getPosFromS(s,0).toPointF());
    pp.lineTo(crv.getPosFromS(s,-h).toPointF());
    s += sunit;
    k--;
  }
  pa.drawPath(pp);
}


void SpiralPage::drawSpiral(QPainter &pa)
{
  pa.setPen(QPen(Qt::black, 3));
  pa.setBrush(Qt::yellow);
  QPainterPath pp;
  pp.addPolygon(QPolygonF(m_pts));
  pp.closeSubpath();
  pp.addEllipse(QPointF(0,0), 5,5);
  pa.drawPath(pp);
  pp = QPainterPath();
  for(int i=0;i<m_pts.count();i++)
  {
    pp.moveTo(m_pts[i]*0.95);
    pp.lineTo(m_pts[i]*0.80);
  }
  pa.setPen(Qt::black);
  pa.setBrush(Qt::NoBrush);
  pa.drawPath(pp);
}


void SpiralPage::drawFoldingSpiral(QPainter &pa, int k0, double parameter)
{
  if(parameter>1.0) parameter = 1.0;
  else if(parameter<0.0) parameter = 0.0;

  if(parameter>0) pa.setPen(Qt::black);
  else pa.setPen(Qt::NoPen);
  pa.setBrush(Qt::cyan);

  QMatrix matrix;

  int k = k0;

  for(;;)
  {
    QPointF p0 = m_pts[k];
    QPointF p1 = m_pts[k+1];

    QPolygonF poly; poly << matrix.map(p0) << matrix.map(p1) << matrix.map(QPointF(0,0));
    pa.drawConvexPolygon(poly);

    if(k+2>=m_pts.count()) break;
    QPointF p2 = m_pts[k+2];

    double cs = QVector2D::dotProduct(
    QVector2D(p2-p1).normalized(), 
    QVector2D(p1-p0).normalized());
    double theta = acos(cs) * parameter;

    QMatrix rot;
    rot.translate(p1.x(),p1.y());
    rot.rotate(-180.0*theta/M_PI);
    rot.translate(-p1.x(),-p1.y());

    matrix = rot * matrix;
    k++;
  }

  k = k0+1;
  matrix = QMatrix();
  for(;;)
  {
    if(k-2<0) break;
    QPointF p0 = m_pts[k];
    QPointF p1 = m_pts[k-1];
    QPointF p2 = m_pts[k-2];

    double cs = QVector2D::dotProduct(
    QVector2D(p2-p1).normalized(), 
    QVector2D(p1-p0).normalized());
    double theta = -acos(cs) * parameter;

    QMatrix rot;
    rot.translate(p1.x(),p1.y());
    rot.rotate(-180.0*theta/M_PI);
    rot.translate(-p1.x(),-p1.y());
    matrix = rot * matrix;

    QPolygonF poly; poly << matrix.map(p1) << matrix.map(p2) << matrix.map(QPointF(0,0));
    pa.drawConvexPolygon(poly);


    k--;
  }
}

QMatrix SpiralPage::getMatrixAt(int k)
{
  if(k<0)k=0; else if(k>=m_pts.count()-1) k = m_pts.count()-2;
  QPointF p = 0.5*(m_pts[k]+m_pts[k+1]);
  QPointF d = m_pts[k+1]-m_pts[k];
  double theta = 180.0 * atan2(d.y(),d.x())/M_PI;
  QMatrix matrix;
  matrix.translate(p.x(),p.y());
  matrix.rotate(180.0 + theta);
  return matrix;
}

void SpiralPage::draw1(QPainter &pa)
{
  QPainterPath pp;
  pp.moveTo(m_pts[0]);
  for(int i=0; i<m_pts.count(); i++)
  {
    pp.lineTo(m_pts[i]);
  }
  pp.closeSubpath();
  pa.setPen(Qt::black);
  pa.setBrush(Qt::yellow);
  pa.drawPath(pp);

  double parameter = getParameter()*0.01;
  if(parameter>1.0) parameter = 1.0;
  else if(parameter<0.0) parameter = 0.0;
  
  int k = (int)(1 * (1-m_par1) + (m_pts.count()-2) * m_par1);
  QMatrix matrix;

  int k0 = k;

  if(m_par2>0) pa.setPen(Qt::black);
  else pa.setPen(Qt::NoPen);
  pa.setBrush(Qt::cyan);
  

  for(;;)
  {
    QPointF p0 = m_pts[k];
    QPointF p1 = m_pts[k+1];

    QPolygonF poly; poly << matrix.map(p0) << matrix.map(p1) << matrix.map(QPointF(0,0));
    pa.drawConvexPolygon(poly);

    if(k+2>=m_pts.count()) break;
    QPointF p2 = m_pts[k+2];

    double cs = QVector2D::dotProduct(
    QVector2D(p2-p1).normalized(), 
    QVector2D(p1-p0).normalized());
    double theta = acos(cs) * m_par2;

    QMatrix rot;
    rot.translate(p1.x(),p1.y());
    rot.rotate(-180.0*theta/M_PI);
    rot.translate(-p1.x(),-p1.y());

    matrix = rot * matrix;
    k++;
  }

  k = k0+1;
  matrix = QMatrix();
  for(;;)
  {
    if(k-2<0) break;
    QPointF p0 = m_pts[k];
    QPointF p1 = m_pts[k-1];
    QPointF p2 = m_pts[k-2];

    double cs = QVector2D::dotProduct(
    QVector2D(p2-p1).normalized(), 
    QVector2D(p1-p0).normalized());
    double theta = -acos(cs) * m_par2;

    QMatrix rot;
    rot.translate(p1.x(),p1.y());
    rot.rotate(-180.0*theta/M_PI);
    rot.translate(-p1.x(),-p1.y());
    matrix = rot * matrix;

    QPolygonF poly; poly << matrix.map(p1) << matrix.map(p2) << matrix.map(QPointF(0,0));
    pa.drawConvexPolygon(poly);


    k--;
  }

  /*
  pa.setPen(Qt::black);
  pa.setBrush(Qt::cyan);
  pa.drawPath(pp);
  */


}

#endif

