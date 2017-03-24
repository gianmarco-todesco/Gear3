#pragma once

#include <QPointF>
#include <QVector>


class ClipperWrapper
{
  class Imp;
  Imp *m_imp;
public:
  ClipperWrapper();
  ~ClipperWrapper();

  void add(QPointF &p0, QPointF &p1, QPointF &p2);

  void getOutline(QVector<QVector<QPointF> > &lines);
};

