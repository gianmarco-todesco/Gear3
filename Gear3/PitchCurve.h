#pragma once

#include <QVector>
#include <QVector2D>

class CurveFunction {
public:
  virtual QVector2D operator()(double t) const = 0;
};

class EllipseFunction : public CurveFunction {
  double m_r,m_e;
public:
  EllipseFunction(double r, double e) : m_r(r), m_e(e) {}
  QVector2D operator()(double t) const;
  static double computePerimeter(double r, double e);
};

class SquareFunction : public CurveFunction {
  double m_l, m_r; // half square edge & corner radius
  double m_a, m_b; // coord of corner circles and length of diagonal from center to corner first point
public:
  SquareFunction(double halfSquareEdge, double cornerRadius);
  QVector2D operator()(double t) const;
  static double computePerimeter(double e, double r);
};

class PitchCurve
{
public:
  struct Point {
    double s,r,phi;
    QVector2D pos,right;
  };

private:
  QVector<Point> m_pts;
  double m_length;

public:
  PitchCurve(const CurveFunction &f, int n = 500);
  ~PitchCurve();

  double getLength() const { return m_length; }
  
  int getPointCount() const { return m_pts.count(); }
  const Point &getPoint(int index) const { return m_pts[index]; }
};


