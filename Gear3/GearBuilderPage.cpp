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
#include <QDebug>


class CurveBuilder {
  QVector<QPointF> m_pts;
public:
  CurveBuilder(const QPointF &p) { m_pts.append(p); }

  void addLine(const QPointF &p) {
    QPointF oldp = m_pts.last();
    int m = 2 + (int)(QVector2D(oldp-p).length()*0.1);
    for(int i=1;i<=m;i++)
    {
      double t = (double)i/(double)m;
      m_pts.append(oldp*(1-t) + p*t);
    }
  }
  void addArc(const QPointF &p1, double f) {
    QPointF p0 = m_pts.last();
    QPointF u = (p1-p0)*0.5;
    QPointF v(-u.y()*f,u.x()*f);
    QPointF pc = (p0+p1)*0.5 - v;
    double r = QVector2D(pc-p0).length();
    u = QVector2D(u).normalized().toPointF();
    v = QVector2D(v).normalized().toPointF();
    int m = 2 + (int)(QVector2D(p1-p0).length()*0.1);
    for(int i=1;i<=m;i++)
    {
      double t = (double)i/(double)m;
      double phi = -M_PI/4 + M_PI/2*t;
      m_pts.append(pc + r*(u*sin(phi) + v*cos(phi)));
    }
  }


  const QVector<QPointF> getPoints() const { return m_pts; }
  
};

//----------------------------------------------------

class GearBuilderPage : public Page, public Pannable {
  int m_mode;
  PitchCurve *m_source, *m_target;
  double m_distance;
  double m_angle;

public:
  GearBuilderPage() : Page("gearBuilder"), m_mode(1), m_distance(400.0), m_angle(0) {
    m_source = m_target = 0;
  }

  ~GearBuilderPage() {  }


  PitchCurve *buildSource();

  


  PitchCurve *makeConjugate(const PitchCurve *curve, double distance);

  void drag(int dx, int dy, int modifiers) {
    if((modifiers & Qt::ShiftModifier)==0 || m_mode==5)
    {
      m_angle += 0.01*dx;
      if(m_angle<0) m_angle=0;
    }
    else
    {
      m_distance += dx * 0.1;
      if(m_source)
      {
        delete m_target; 
        m_target = makeConjugate(m_source, m_distance);
      }
    }
  }


  void draw(QPainter &pa);

  void drawTile(QPainter &pa);
  void drawMosaic(QPainter &pa);

  bool onKey(int key) {
    if(key == Qt::Key_1) {m_mode=1;}
    else if(key == Qt::Key_2) {m_mode=2;}
    else if(key == Qt::Key_3) {m_mode=3;}
    else if(key == Qt::Key_4) {m_mode=4;}
    else if(key == Qt::Key_5) {m_mode=5;}
    else return false;
    return true;

  }

} gearBuilder;


class PointsFunction : public CurveFunction {
  QVector<QPointF> m_pts; 
public:
  PointsFunction(const QVector<QPointF> &pts) : m_pts(pts) {}
  QVector2D operator()(double t) const {
    int n = m_pts.count();
    t *= (n-1);
    int k = (int)t;
    if(k>=n-1) k = n-2;
    double s = t - k;
    return QVector2D(m_pts[k]*(1-s)+m_pts[k+1]*s);
  }
};

PitchCurve *GearBuilderPage::buildSource()
{
/*  
  QVector<QPointF> pts;
  int n = 100;
  pts.resize(n);
  for(int i=0;i<n;i++)
  {
    double phi = 2*M_PI*i/(n-1);
    double r = 100 + 30*sin(phi*3);
    pts[i] = r*QPointF(cos(phi), sin(phi));
  }
  */

  double lx = 100, ly = lx*tan(M_PI/6);
  CurveBuilder builder(QPointF(2*lx,0));
  builder.addLine(QPointF(lx,ly));
  builder.addArc(QPointF(0,2*ly),1);
  builder.addArc(QPointF(-lx,ly),-1);
  builder.addLine(QPointF(-2*lx,0));
  builder.addLine(QPointF(-lx,-ly));
  builder.addArc(QPointF(0,-2*ly),-1);
  builder.addArc(QPointF(lx,-ly),1);
  builder.addLine(QPointF(2*lx,0));
  return new PitchCurve(PointsFunction(builder.getPoints()));
 
}

PitchCurve *GearBuilderPage::makeConjugate(const PitchCurve *curve, double distance)
{
  int n = curve->getPointCount();
  for(int i=0;i<n;i++)
  {
    double r = curve->getPoint(i).r + 50;
    if(r>distance)distance = r;
  }

  QVector<PitchCurve::Point> pts;
  pts.reserve(curve->getPointCount());

  int i = 0;
  PitchCurve::Point pt;
  pt.r = distance - curve->getPoint(0).r;
  pt.s = 0.0;
  pt.phi = 0;

  pts.append(pt);

  for(;;)
  {
    double ds;
    if(i==0) {ds = curve->getLength() - curve->getPoint(n-1).s; i=n-1; }
    else { ds = curve->getPoint(i).s - curve->getPoint(i-1).s; i--; }
    double r0 = pt.r;
    double r1 = distance - curve->getPoint(i).r;
    double dphi = acos((r0*r0+r1*r1-ds*ds)/(2*r0*r1));
    pt.phi += dphi;
    pt.s += ds;
    pt.r = r1;
    if(pt.phi<2*M_PI) 
    {
      pts.append(pt);
    }
    else 
    {
      if(pt.phi>2*M_PI)
      {
        double f = (pt.phi - 2*M_PI)/dphi;
        pt.r = pt.r * (1-f) + r0 * f;
        pt.phi = 2*M_PI;
        pt.s -= ds*f;
      }
      pts.append(pt);
      break;
    }
  }

  for(int i=0;i<pts.count();i++)
  {
    double phi = pts[i].phi;
    pts[i].pos = pts[i].r * QVector2D(cos(phi), sin(phi));
    pts[i].right = QVector2D(pts[(i+1)%pts.count()].pos - pts[(i+pts.count()-1)%pts.count()].pos).normalized();
  }
  return new PitchCurve(pts);
}


void GearBuilderPage::drawTile(QPainter &pa)
{
  QPainterPath pp;
  pp.moveTo(m_source->getPoint(0).pos.toPointF());
  for(int i=4;i<m_source->getPointCount();i+=4)
    pp.lineTo(m_source->getPoint(i).pos.toPointF());
  pp.closeSubpath();
  pa.drawPath(pp);
}

void GearBuilderPage::drawMosaic(QPainter &pa)
{
  double lx = 100, ly = lx*tan(M_PI/6.0);

  pa.setPen(QPen(Qt::black, 2));
  QColor colors[3] = {
    QColor(250,250,250),
    QColor(90,90,90),
    QColor(200,180,90)};
  double d = sqrt(lx*lx+ly*ly);

  for(int iy= -4;iy<=4;iy++)
  for(int ix= -4;ix<=4;ix++)
  {
    for(int i=0;i<3;i++)
    {
      pa.save();
      pa.translate(0,ly*2);
      pa.translate(4*lx*ix + 2*lx*iy, (2*d+2*ly)*iy);
      pa.rotate(120*i);
      pa.translate(0,-ly*2);
      pa.setBrush(colors[(i+ix-iy+1000)%3]);
      drawTile(pa);
      pa.restore();
    }
  }


  /*
  pa.setBrush(QColor(200,160,20));
  for(int i=0;i<6;i++)
  {
    pa.save();
    pa.rotate(120*(i/2));
    pa.translate(-400,0);
    pa.rotate(120*(1-2*(i%2)));
    pa.translate(-200,0);
    drawTile(pa);
    pa.restore();
  }
  */
}

void GearBuilderPage::draw(QPainter &pa)
{
  const double distance = m_distance;
  if(!m_source)
  {
    QElapsedTimer timer;
    timer.start();
    m_source = buildSource();
    qDebug() << "source " << timer.elapsed();
    timer.start();
    m_target = makeConjugate(m_source, distance);
    qDebug() << "target " << timer.elapsed();
  }
  if(m_mode==1) { drawMosaic(pa); return; }


  double angle = -2*M_PI + m_angle;
  if(angle<-2*M_PI) angle = -2*M_PI;

  double maxAngle = -m_source->getPhiFromS(m_source->getLength() - m_target->getLength());
  if(angle>maxAngle) angle = maxAngle;

  // draw first gear
  pa.save();
  pa.rotate(angle*180/M_PI);

  QPainterPath pp;
  pp.moveTo(m_source->getPoint(0).pos.toPointF());
  for(int i=1;i<m_source->getPointCount();i++)
    pp.lineTo(m_source->getPoint(i).pos.toPointF());

  pp.addEllipse(QPointF(0,0), 5,5);

  pa.setPen(QPen(Qt::black, 3));
  pa.setBrush(QColor(90,90,90));
  pa.drawPath(pp);

  pa.restore();

  pa.setPen(Qt::black);
  pa.setBrush(Qt::NoBrush);
  pa.drawLine( 6, 0, 100, 0);
  pa.drawLine( -100,0,-6,0);
  pa.drawLine( 0, -100, 0,-6);
  pa.drawLine( 0, 6, 0, 100);


  // draw second gear

  if(m_mode > 2)
  {



    pa.save();
    pa.translate(distance,0);


    double s = m_source->getLength() - m_source->getSFromPhi(-angle);
    double angle2 = M_PI - m_target->getPhiFromS(s);
    pa.rotate(angle2 * 180/M_PI);

    if(m_mode==3)
    {
      // disegno progressivo 


      double angle0 = -2*M_PI;
      double da = 0.01;
      if(angle-angle0>da)
      {
        pp = QPainterPath();
        pp.moveTo(QPointF(0,0));
        double a = angle0;
        for(;;)
        {
          if(a>angle)a=angle;
          double si = m_source->getLength() - m_source->getSFromPhi(-a);
          pp.lineTo(m_target->getPosFromS(si).toPointF());
          if(a==angle) break;
          a+=da;
        }
        pa.setPen(QPen(Qt::black, 3));
        pa.setBrush(Qt::cyan);
        pa.drawPath(pp);

        pa.setPen(QPen(Qt::red, 3));
        pa.drawLine(
          QPointF(0,0), 
          m_target->getPosFromS(m_source->getLength() - m_source->getSFromPhi(-angle0)).toPointF());
        pa.drawLine(
          QPointF(0,0), 
          m_target->getPosFromS(m_source->getLength() - m_source->getSFromPhi(-angle)).toPointF());

                  
      }
      pa.setPen(Qt::black);
      pa.setBrush(Qt::white);
      pa.drawEllipse(QPointF(0,0), 10,10);
    }
    else if(m_mode==4 || m_mode==5)
    {
      QPainterPath pp;
      pp.moveTo(m_target->getPoint(0).pos.toPointF());
      for(int i=1;i<m_target->getPointCount();i++)
        pp.lineTo(m_target->getPoint(i).pos.toPointF());
      pp.addEllipse(QPointF(0,0), 5,5);
      pa.setPen(QPen(Qt::black, 3));
      pa.setBrush(Qt::cyan);
      pa.drawPath(pp);
    }
    pa.restore();

    pa.drawLine(distance,-50, distance,50);
    pa.drawLine(distance-50,0,distance+50,0);
  }

  




  
  

}


