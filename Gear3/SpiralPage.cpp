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

  void draw(QPainter &pa);
  void draw1(QPainter &pa);

  void drag(int dx, int dy, int modifiers) {
    m_par1 += 0.01*dx; if(m_par1<0)m_par1=0; else if(m_par1>1)m_par1=1.0;
    m_par2 += 0.01*dy; if(m_par2<0)m_par2=0; else if(m_par2>1)m_par2=1.0;

  }

  bool onKey(int key) {
    return false;
  }

} spiralPage;

void SpiralPage::buildPoints()
{
  PitchCurve crv(SpiralFunction(200,500));
  int n = 50;
  for(int i=0;i<n;i++)
  {
    double s = crv.getLength()*i/(n-1);
    m_pts.append(crv.getPosFromS(s).toPointF());
  }
}

void SpiralPage::draw(QPainter &pa)
{
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

