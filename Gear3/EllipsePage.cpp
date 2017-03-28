#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>
#include <QVector2D>
#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"

#include <qmath.h>
namespace {
  double atan(const QPointF &p) { return atan2(p.y(),p.x()); }
}

class EllipsePage : public Page, public Pannable {
  double m_a, m_b;
  double m_f;
  int m_mode;
  QVector<QPointF> m_ticks; 
public:
  EllipsePage() : Page("ellipse"), m_a(300), m_b(200), m_mode(1) {
    m_f = sqrt(m_a*m_a-m_b*m_b);
    computeTicks();
  }

  ~EllipsePage() {  }

  void computeTicks();

  void draw(QPainter &pa);
  void drawEllipse(QPainter &pa);

  void drawSegments(QPainter &pa, const QPointF &pc);
  void drawTangent(QPainter &pa, double phi);

  bool onKey(int key) {
    if(key == Qt::Key_1) {m_mode=0;}
    else if(key == Qt::Key_2) {m_mode=1;}
    else if(key == Qt::Key_3) {m_mode=2;}
    else return false;
    return true;

  }

} ellipsePage;


void EllipsePage::computeTicks()
{
  PitchCurve crv(CenteredEllipseFunction(m_a,m_b),1000);
  double length = crv.getLength();
  int m = 100;
  for(int i=0;i<m;i++)
  {
    m_ticks.append(crv.getPosFromS(i*length/m, -20).toPointF());
  }

}

void EllipsePage::draw(QPainter &pa)
{
  double phi = getParameter() * 0.01;
  QPointF pc(m_a*cos(phi), m_b*sin(phi));  
  drawEllipse(pa);
  
  if(m_mode >= 1)
  {
    drawTangent(pa,phi);

    if(m_mode >= 2)
    {
      pa.save();
      pa.translate(pc);

      {
        double cs = cos(phi), sn = sin(phi);
        QPointF pc = QPointF(m_a*cs, m_b*sn);
        QPointF dir = QPointF(-m_a*sn, m_b*cs);
        double theta = atan(dir);

        pa.rotate(180.0*theta/M_PI);
        pa.scale(1,-1);
        pa.rotate(-180.0*theta/M_PI);

      }
      pa.translate(-pc);
      drawEllipse(pa);
      pa.restore();
    }
  }

}

void EllipsePage::drawEllipse(QPainter &pa)
{
  pa.setPen(Qt::black);
  pa.setBrush(QColor(200,200,200));
  
  
  QPainterPath pp;
  pp.addEllipse(QPointF(0,0), m_a, m_b);
  for(int i=0;i<m_ticks.count();i++) 
    pp.addEllipse(m_ticks[i], 4,4);
  pa.drawPath(pp);

  pa.setBrush(Qt::NoBrush);
  pa.setPen(QPen(Qt::black, 3));
  pa.drawEllipse(QPointF(0,0), m_a, m_b);

  double phi = getParameter() * 0.01;
  QPointF pc(m_a*cos(phi), m_b*sin(phi));
  pa.setBrush(Qt::black);


  QPointF f0(-m_f,0);
  QPointF f1( m_f,0);

  if(m_mode==0) drawSegments(pa, pc);
  else
  {
    pa.setPen(QPen(Qt::red,2));
    pa.drawLine(f0,pc);
    pa.setPen(QPen(Qt::green,2));
    pa.drawLine(f1,pc);
    
  }

  pa.setPen(Qt::black);
  pa.setBrush(Qt::black);
  pa.drawEllipse(pc,5,5);

  pa.drawEllipse(QPointF(-m_f,0),5,5);
  pa.drawEllipse(QPointF( m_f,0),5,5);


}



void EllipsePage::drawSegments(QPainter &pa, const QPointF &pc)
{
  QPointF f0(-m_f,0);
  QPointF f1( m_f,0);


  double LL[2];
  LL[0] = QVector2D(pc-f0).length();
  LL[1] = QVector2D(pc-f1).length();
  double L = LL[0]+LL[1];

  int m = 11;
  QVector<QPointF> pts[2];
  for(int i=0;i<2;i++)
  {
    QPointF pf = i==0 ? f0 : f1;
    QVector<QPointF> &pp = pts[i];
    for(int j=0;;j++)
    {
      double s = L*j/(m-1);
      if(s>=LL[i])
      {
        if(QVector2D(pc-pp.back()).length()>2) pp.append(pc);
        break;
      }
      double t = s/LL[i]; 
      pp.append((1-t)*pf+t*pc); 
    }
  }

  for(int colorIndex=0;colorIndex<2;colorIndex++)
  {
    QPainterPath pp;
    for(int i=0;i<2;i++)
    {
      for(int j=(colorIndex+i)%2;j+1<pts[i].count();j+=2)
      {
        pp.moveTo(pts[i][j]); 
        pp.lineTo(pts[i][j+1]);
      }
    }
    pa.setPen(QPen(colorIndex==0 ? Qt::magenta : Qt::cyan, 4));
    pa.drawPath(pp);

  }


}

void EllipsePage::drawTangent(QPainter &pa, double phi)
{
  double cs = cos(phi), sn = sin(phi);
  QPointF pc = QPointF(m_a*cs, m_b*sn);
  QPointF dir = QPointF(-m_a*sn, m_b*cs);
  pa.setPen(Qt::gray);

  QPointF f1(-m_f,0), f2(m_f,0);
  QPointF p1 = -500*dir+pc, p2 = 500*dir+pc;
  pa.drawLine(-500*dir+pc, 500*dir+pc);

  if(QVector2D::dotProduct(QVector2D(p1-pc),QVector2D(f1-pc))<0) qSwap(p1,p2);


  QPainterPath pp;
  double r = 60;
  QRectF box(pc.x()-r,pc.y()-r,2*r,2*r);
  
  QColor color(200,100,10,127);

  double theta = atan(f1-pc);
  double dtheta = atan(p1-pc) - theta;
  if(dtheta<-M_PI) dtheta += 2*M_PI;
  else if(dtheta>M_PI) dtheta -= 2*M_PI;

  pp.moveTo(pc);
  pp.arcTo(box, -(int)(180.0*theta/M_PI), -(int)(180.0*dtheta/M_PI));
  pa.setBrush(color);
  pa.setPen(Qt::NoPen);
  pa.drawPath(pp);

  pp = QPainterPath();
  pp.arcMoveTo(box, -(int)(180.0*theta/M_PI));
  pp.arcTo(box, -(int)(180.0*theta/M_PI), -(int)(180.0*dtheta/M_PI));
  pa.setBrush(Qt::NoBrush);
  pa.setPen(Qt::black);
  pa.drawPath(pp);

  theta = atan(f2-pc);
  dtheta = atan(p2-pc) - theta;
  if(dtheta<-M_PI) dtheta += 2*M_PI;
  else if(dtheta>M_PI) dtheta -= 2*M_PI;

 
  pp = QPainterPath();
  pp.moveTo(pc);
  pp.arcTo(box, -(int)(180.0*theta/M_PI), -(int)(180.0*dtheta/M_PI));
  pa.setBrush(color);
  pa.setPen(Qt::NoPen);
  pa.drawPath(pp);


  pp = QPainterPath();
  pp.arcMoveTo(box, -(int)(180.0*theta/M_PI));
  pp.arcTo(box, -(int)(180.0*theta/M_PI), -(int)(180.0*dtheta/M_PI));
  pa.setBrush(Qt::NoBrush);
  pa.setPen(Qt::black);
  pa.drawPath(pp);
  

}
