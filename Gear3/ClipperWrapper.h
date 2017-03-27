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
  void add(const QVector<QPointF> &outline);
  void sub(const QVector<QPointF> &outline);
  void intersect(const QVector<QPointF> &outline);

  void getOutline(QVector<QVector<QPointF> > &lines);
};

