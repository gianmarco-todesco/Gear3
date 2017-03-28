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


class SpiralPage : public Page, public Pannable {
  QVector<QPointF> m_pts;
  double m_par1, m_par2;

public:
  SpiralPage() : Page("spiral"), m_par1(0), m_par2(0) {
    buildPoints();
  }

  ~SpiralPage() {  }

  void buildPoints();

  void drawFoldingSpiral(QPainter &pa, int k, double param);
  QMatrix getMatrixAt(int k);

  void draw(QPainter &pa);
  void draw1(QPainter &pa);
  void drawSpiral(QPainter &pa);

  void drag(int dx, int dy, int modifiers) {
    m_par1 += 0.003*dx; if(m_par1<0)m_par1=0; else if(m_par1>1)m_par1=1.0;
    m_par2 += 0.01*dy; if(m_par2<0)m_par2=0; else if(m_par2>1)m_par2=1.0;

  }

  bool onKey(int key) {
    return false;
  }

} spiralPage;

void SpiralPage::buildPoints()
{
  PitchCurve crv(SpiralFunction(200,500));
  int n = 250;
  for(int i=0;i<n;i++)
  {
    double s = crv.getLength()*i/(n-1);
    m_pts.append(crv.getPosFromS(s).toPointF());
  }
}

void SpiralPage::draw(QPainter &pa)
{
  int k = (int)(m_pts.count()*m_par1);
  drawSpiral(pa);

  pa.save();
  pa.setMatrix(getMatrixAt(k), true);

  pa.save();
  pa.rotate(180);
  pa.setMatrix(getMatrixAt(m_pts.count()-2-k).inverted(), true);
  drawSpiral(pa);
  if(m_par2>0) 
    drawFoldingSpiral(pa, m_pts.count()-2-k, m_par2);
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

