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
  QBrush m_brush;

public:
  Gear(PitchCurve *curve);
  ~Gear();

  void setBrush(const QBrush &brush) { m_brush = brush; }
  const QBrush &getBrush() const { return m_brush; }

  const PitchCurve *getCurve() const { return m_curve; }

  QPointF getPosition() const { return m_position; }
  void setPosition(const QPointF &p) { m_position = p; }
  void setPosition(double x, double y) { m_position = QPointF(x,y); }

  double getAngle() const { return m_angle; }
  void setAngle(double angle) { m_angle = angle; }

  const QPainterPath &getBodyPath() const { return m_bodyPath; }
  void setBodyPath(const QPainterPath &path) { m_bodyPath = path; }
  void setBodyPath(const QVector<QVector2D> &pts);

  const QPainterPath &getPitchLinePath() const { return m_pitchLinePath; }
  
  void draw(QPainter &pa);
  void rotoTranslate(QPainter &pa);

private:
  void updatePitchPath();
};



class GearLink {
  Gear *m_gear1, *m_gear2; // 1 = driver, 2 = driven
  double m_theta1, m_theta2; // the gears match when: m_gear1->getAngle() == m_theta1 && m_gear2->getAngle() == m_theta2;
  double m_distance;
public:
  GearLink(Gear*driver, Gear *driven);
  ~GearLink();
  double getDistance() const { return m_distance; }
  double computeDistance() const;
  void update();

  void moveDriven(double psi);

};


class GearBox {
  QList<Gear*> m_gears;
  QList<GearLink*> m_links;

public:
  GearBox();
  virtual ~GearBox();

  void clear();

  Gear *addGear(Gear*gear);
  GearLink *addLink(GearLink*link);
  GearLink *addLink(Gear *a, Gear *b) { return addLink(new GearLink(a,b)); }
  GearLink *addLink(int a, int b) {return addLink(m_gears[a], m_gears[b]); }

  int getGearCount() const { return m_gears.count(); }
  Gear *getGear(int index) const { return m_gears.at(index); }

  int getLinkCount() const { return m_links.count(); }
  GearLink *getLink(int index) const { return m_links.at(index); }

  virtual void draw(QPainter &pa);
};


Gear *makeCircularGear(int toothLength, int toothCount, int flag = 0);
