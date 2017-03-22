#pragma once

class QPainter;


class Sandbox {
  struct Imp;
  Imp*m_imp;

public:
  Sandbox();
  ~Sandbox();

  void setParameter(double v);
  double getParameter() const;
  void changeParameter(double dv) { setParameter(getParameter()+dv); }

  void paint(QPainter &pa);
};
