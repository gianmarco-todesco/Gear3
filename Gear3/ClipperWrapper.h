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

  void add(const QVector<QPointF> &outline);
  void add(const QVector<QVector<QPointF> > &outline);
  void add(QPointF &p0, QPointF &p1, QPointF &p2);

  void subtract(const QVector<QPointF> &outline);
  void subtract(const QVector<QVector<QPointF> > &outline);

  void getResult(QVector<QVector<QPointF> > &lines);
};

