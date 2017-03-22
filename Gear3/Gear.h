#pragma once

#include <QPointF>
#include <QPainterPath>
#include <QPainter>

class PitchCurve;


class Gear {
  PitchCurve *m_curve;
  QPointF m_position;
  double m_angle;
  QPainterPath m_pitchLinePath;
  QPainterPath m_bodyPath;

public:
  Gear(PitchCurve *curve);
  ~Gear();

  const PitchCurve *getCurve() const { return m_curve; }

  QPointF getPosition() const { return m_position; }
  void setPosition(const QPointF &p) { m_position = p; }

  double getAngle() const { return m_angle; }
  void setAngle(double angle) { m_angle = angle; }

  const QPainterPath &getBodyPath() const { return m_bodyPath; }
  void setBodyPath(const QPainterPath &path) { m_bodyPath = path; }
  void setBodyPath(const QVector<QVector2D> &pts);

  void draw(QPainter &pa);

private:
  void updatePitchPath();
};
