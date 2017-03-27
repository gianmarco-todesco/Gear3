#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>

#include <qmath.h>
#include "Gear.h"
#include "ToothMaker.h"
#include "PitchCurve.h"
#include "ClipperWrapper.h"
#include "PolygonSweeper.h"

#include "ToothMaker.h"

#include <QLineF>
#include <QDebug>


namespace {
  double getNorm2(const QPointF &p) { return p.x()*p.x()+p.y()*p.y(); }
  double atan2p(const QPointF &p) { return atan2(p.y(),p.x()); }
}

class SquareTeethGearPage : public Page, public Pannable {
  double m_r0, m_r1;
  int m_toothCount;
  double m_theta;
  Gear *m_gear1, *m_gear2;

  QPointF m_contactPoint;
  bool m_contact;
  double m_contactAngle;

  QList<QPointF> m_graph;

public:
  SquareTeethGearPage() : Page("squareGears"), m_toothCount(11), m_contact(false), m_contactAngle(0.0) { 
    double toothLength = 80.0;
    Gear *gear1 = m_gear1 = makeCircularGear(toothLength, m_toothCount, 1);
    Gear *gear2 = m_gear2 = makeCircularGear(toothLength, m_toothCount, 1);
    double dist = gear1->getCurve()->getRFromS(0) + gear2->getCurve()->getRFromS(0);
    double radius = dist*0.5;
    double h = toothLength*0.5*0.75;
    m_r0 = radius - h/2;
    m_r1 = m_r0 + h;
    m_theta = 0.2*(2*M_PI/m_toothCount);
    gear2->setPosition(-dist,0);
  }


  ~SquareTeethGearPage() {  }

   void draw(QPainter &pa) {
     m_gear1->setAngle(getParameter()*0.01);

     QPainterPath pp1, pp2;
     QPointF contactPoint;
     bool contact = foo(contactPoint);
     if(contact) {m_contactPoint = contactPoint; m_contactAngle = m_gear1->getAngle(); }
     else {
       double d = m_contactAngle - m_gear1->getAngle();
       if(d<-M_PI) d+=2*M_PI; else if(d>M_PI) d-=2*M_PI;
       if(fabs(d)>0.05) {} else contact = true;
     }

     m_gear1->draw(pa);
     m_gear2->draw(pa);

     if(contact)
     {
       QPainterPath pp;
       pp.addEllipse(m_contactPoint, 5,5);
       pp.closeSubpath();
       pp.addEllipse(m_contactPoint, 3,3);
       pp.closeSubpath();

       pa.setPen(Qt::black);
       pa.setBrush(Qt::red);
       pa.drawPath(pp);
     }
     qDebug() << m_gear1->getAngle() + m_gear2->getAngle();

     double angle1 = m_gear1->getAngle(), angle2 = m_gear2->getAngle();
     if(m_graph.empty() || angle1>m_graph.back().x()) m_graph.append(QPointF(angle1, angle1+angle2));
     else {
       while(!m_graph.empty() && m_graph.back().x()>angle1) m_graph.pop_back();
     }
   }

   void drawOverlay(QPainter &pa)
   {
     int w = getWidth(), h = getHeight();
     pa.setPen(Qt::black);
     pa.setBrush(Qt::NoBrush);
     pa.drawRect(10,h-110,100,100);

     QPainterPath pp;
     for(int i=0;i<m_graph.count();i++)
     {
       double angle1 = m_graph[i].x();
       double da = m_graph[i].y();
       double x = 10 + 100 * angle1 / (2*M_PI);
       double y = h-10-50 - da * 250;
       if(i==0) pp.moveTo(x,y); else pp.lineTo(x,y);
     }
     pa.drawPath(pp);


   }

   bool foo(QPointF &contactPoint) {
     QList<QLineF> l1;
     QList<QLineF> l2;
     getSegments(m_gear1, l1);
     getSegments(m_gear2, l2);
     QPointF center1 = m_gear1->getPosition();
     QPointF center2 = m_gear2->getPosition();
     double r1_2 = m_r1*m_r1;

     double maxPositiveAngle = 0.0;
     double maxNegativeAngle = 0.0;

     for(int i=0;i+1<l1.count();i+=2)
     {
       QPointF p1a = l1[i].p2();
       QPointF p1b = l1[i+1].p2();

       QPointF p1 = 0.5*(p1a+p1b);
       if(getNorm2(p1 - center2)>r1_2) continue;
       
       double phi1a = atan2p(p1a - center2);
       double phi1b = atan2p(p1b - center2);

       double theta1a = atan2p(center1 - p1a);
       double theta1b = atan2p(center1 - p1b);
       
       
       for(int j=0; j+1<l2.count();j+=2)
       {
         QPointF p2a = l2[j].p2();
         QPointF p2b = l2[j+1].p2();

         QPointF p2 = 0.5*(p2a+p2b);
         if(getNorm2(p2 - center1)>r1_2) continue;
         double phi2a = atan2p(p2a - center2);
         double phi2b = atan2p(p2b - center2);
         double theta2a = atan2p(center1 - p2a);
         double theta2b = atan2p(center1 - p2b);

         if((phi1a-phi2a)*(phi1a-phi2b)<0.0)           
         {
           double theta = phi1a-phi2a;
           if(theta>0.0 && theta>maxPositiveAngle) { maxPositiveAngle = theta; contactPoint = p1a; }
         }
         else if((phi1b-phi2a)*(phi1b-phi2b)<0.0 )
         {
           double theta = phi1b-phi2b;
           if(theta<0.0 && -theta>maxNegativeAngle) { maxNegativeAngle = -theta;contactPoint = p1b; }
         }
         QPointF p;           
         if((theta2a-theta1a)*(theta2a-theta1b)<0.0)
         {
           intesectCircle(l1[i], center2, &p);
           double theta = atan2p(p-center2) - phi2a;
           if(theta>0.0 && theta>maxPositiveAngle) { maxPositiveAngle = theta;   contactPoint = p; }   
         }
         else if((theta2b-theta1a)*(theta2b-theta1b)<0.0)
         {
           intesectCircle(l1[i+1], center2, &p);
           
           double theta = atan2p(p-center2) - phi2b;
           if(theta<0.0 && -theta>maxNegativeAngle) { maxNegativeAngle = -theta;    contactPoint = p; }       
         }
       }
     }

     bool flag = true;
     if(maxPositiveAngle>0.0) 
       m_gear2->setAngle(m_gear2->getAngle() + maxPositiveAngle);
     else if(maxNegativeAngle>0.0)
       m_gear2->setAngle(m_gear2->getAngle() - maxNegativeAngle);
     else
       flag = false;
     return flag;

    
   }

   bool intesectCircle(const QLineF &line, const QPointF &center, QPointF *p) {
    // p = line.p1() + (line.p2()-line.p1())*lambda 
    // |A + B*lambda| = m_r1
    // lambda^2*<B,B> + 2*lambda*<A,B> +<A,A>-m_r1 = 0

    QPointF A = line.p1() - center;
    QPointF B = line.p2()-line.p1();
    double a = B.x()*B.x()+B.y()*B.y();
    double b = A.x()*B.x()+A.y()*B.y();
    double c = A.x()*A.x()+A.y()*A.y() - m_r1*m_r1;
    double dsc = b*b-a*c;
    if(dsc<0) return false;
    double lambda = (-b-sqrt(dsc))/a;
    if(0<=lambda && lambda<1.0) 
    {
      *p = line.p1() + B*lambda;
      return true;
    }
    else return false;
  }

   void getSegments(Gear *gear, QList<QLineF> &pts) {
    for(int i=0;i<m_toothCount;i++)
    {
      double phi_i = 2*M_PI*i/m_toothCount + gear->getAngle();
      for(int j=0;j<2;j++)
      {
        double phi = phi_i + (-1+2*j)*m_theta;
        QPointF u(cos(phi), sin(phi));
        pts.append( QLineF(gear->getPosition() + m_r0 * u,  gear->getPosition() + m_r1 * u ) );
      }
    }
  }

} squareTeethGearPage;
