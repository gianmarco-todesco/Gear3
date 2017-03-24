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

class SpiralFunction : public CurveFunction {
  double m_r0, m_r1;
public:
  SpiralFunction(double r0, double r1) : m_r0(r0), m_r1(r1) {}
  QVector2D operator()(double t) const;
};


//=============================================================================

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
  bool m_isOpen;

public:
  PitchCurve(const CurveFunction &f, int n = 500);
  ~PitchCurve();

  bool isOpen() const { return m_isOpen; }
  double getLength() const { return m_length; }
  
  int getPointCount() const { return m_pts.count(); }
  const Point &getPoint(int index) const { return m_pts[index]; }

  double getPhiFromS(double s) const;
  double getSFromPhi(double phi) const;

  double getRFromS(double s) const;
  double getRFromPhi(double phi) const;

  QVector2D getPosFromS(double s, double y = 0.0) const;
  QVector2D getPosFromPhi(double phi) const;

  Point getPointFromS(double s) const;

  // for close curves
  void getIndexFromS(double s, int &a, int &b, int &q, double &t) const;
  void getIndexFromPhi(double phi, int &a, int &b, int &q, double &t) const;

  // for open curves
  void getIndexFromS(double s, int &a, int &b, double &t) const;
  void getIndexFromPhi(double phi, int &a, int &b, double &t) const;

};


void testPitchCurve();


