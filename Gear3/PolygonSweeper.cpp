#include "PolygonSweeper.h"
#include "ClipperWrapper.h"
#include <QPainterPath>
#include <QPainter>
#include <qmath.h>
#include <QtAlgorithms>
#include <QPolygonF>


static double getNorm2(const QPointF &p) { return p.x()*p.x()+p.y()*p.y(); }

PolygonSweeper::PolygonSweeper()
{
}


PolygonSweeper::~PolygonSweeper()
{
}


void PolygonSweeper::add(const QMatrix &matrix)
{
  QVector<QPointF> polygon; polygon.reserve(m_shape.count());
  for(int i=0;i<m_shape.count();i++) polygon.append(matrix.map(m_shape[i]));
  m_cw.add(polygon);
  if(!m_lastPolygon.empty())
  {
    QList<QPointF> hull;
    makeConvexHull((polygon + m_lastPolygon).toList(), hull);
    m_cw.add(hull.toVector());
  }
  m_lastPolygon = polygon;
}

/*
void PolygonSweeper::draw(QPainter &pa)
{
  ClipperWrapper cw;
  for(int i=0;i+1<m_polygons.count();i++)
  {
    QList<QPointF> hull;
    makeConvexHull(m_polygons[i] + m_polygons[i+1], hull);
    cw.add(hull.toVector());
    / *
    QPainterPath pp;
    pp.addPolygon(QPolygonF(hull.toVector()));
    pp.closeSubpath();
    pa.setBrush(Qt::yellow);
    pa.setPen(QPen(Qt::black, 1));
    pa.drawPath(pp);
    * /
  }
  QVector<QVector<QPointF> > lines;
  cw.getResult(lines);
  QPainterPath pp;
  for(int i=0;i<lines.count();i++) { pp.addPolygon(QPolygonF(lines[i])); pp.closeSubpath(); }
  pa.setBrush(Qt::yellow);
  pa.setPen(QPen(Qt::black, 1));
  pa.drawPath(pp);

  for(int i=0;i<m_polygons.count();i++)
  {
    const QList<QPointF> &polygon = m_polygons[i];
    QPainterPath pp;
    pp.moveTo(polygon.last());
    for(int j=0;j<polygon.count();j++) pp.lineTo(polygon[j]);
    pa.setBrush(Qt::NoBrush);
    pa.setPen(Qt::red);
    pa.drawPath(pp);
  }
}
*/


void PolygonSweeper::makeConvexHull(const QList<QPointF> &pts, QList<QPointF> &hull)
{
  int n = pts.count();
  int k = 0;
  for(int i=1;i<n;i++)
  {
    if(pts[i].y()<pts[k].y() || pts[i].y()==pts[k].y() && pts[i].x()<pts[k].x()) k=i;
  }
  QVector<QPair<double, int> > lst;
  lst.reserve(n);
  for(int i=0;i<n;i++)
  {
    if(i==k) continue;
    QPointF p = pts[i]-pts[k];
    double phi = atan2(p.y(),p.x());
    Q_ASSERT(0<=phi && phi<M_PI);
    lst.append(qMakePair(phi, i));
  }
  qSort(lst);
  QList<QPointF> pts1; pts1.reserve(n);
  for(int i=0;i<lst.count();i++)
  {
    QPointF p = pts[lst[i].second];
    while(i+1<lst.count() && lst[i].first == lst[i+1].first)
    {
      QPointF p1 = pts[lst[i+1].second];
      if(getNorm2(p1)>getNorm2(p)) p=p1;
      i++;
    }
    pts1.append(p);
  }
  pts1.append(pts[k]);

  hull.reserve(n);
  hull.append(pts[k]);
  hull.append(pts1[0]);
  for(int i=1;i<pts1.count();i++)
  {
    QPointF p2 = pts1[i];
    for(;;)
    {
      Q_ASSERT(hull.count()>=2);
      QPointF p1 = hull.back();
      QPointF p0 = hull[hull.count()-2];
      QPointF p01 = p1-p0, p12 = p2-p1;
      double q = -p12.x()*p01.y()+p12.y()*p01.x();
      if(q>0.0) break;
      hull.pop_back();
    }
    hull.append(p2);
  }
}



