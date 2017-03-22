#pragma once

#include <QVector>
#include <QVector2D>
#include <QPair>

class PitchCurve;


class CircularToothMaker {
public:

  struct Params {
    double pitchRadius;
    int toothCount;
    double toothHeight;
  };

  CircularToothMaker() {}

  void makeTeeth(QVector<QVector2D> &pts, const Params &params);

private:
  QVector2D evolute(double phi, double radius) const { 
      double cs = cos(phi), sn = sin(phi);
      return radius * QVector2D(cs + phi*sn, sn - phi*cs);
  }

  void makeToothSide(QList<QPair<double, double> > &qs, double pitchRadius, double toothHeight);
  void completeTooth(QList<QPair<double, double> > &qs, int toothCount);
};

class SimpleToothMaker {
public:

  struct Params {
    int toothLength;
    double toothHeight;
  };

  SimpleToothMaker();
  void makeTeeth(QVector<QVector2D> &pts, const PitchCurve *curve, const Params &params);

};
