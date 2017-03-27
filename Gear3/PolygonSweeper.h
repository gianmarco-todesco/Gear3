#pragma once

#include <QList>
#include <QVector>
#include <QLineF>
#include <QMatrix>
#include "ClipperWrapper.h"

class QPainter;

class PolygonSweeper
{
  QList<QPointF> m_shape;
  QVector<QPointF> m_lastPolygon;
  ClipperWrapper m_cw;
  

public:
  PolygonSweeper();
  ~PolygonSweeper();

  void setShape(const QList<QPointF> &shape) { m_shape = shape; }

  void add(const QMatrix &matrix);

  // void draw(QPainter &pa);

  void getOutline(QVector<QVector<QPointF> > &outline) { m_cw.getResult(outline); }


  void makeConvexHull(const QList<QPointF> &pts, QList<QPointF> &hull);


};

