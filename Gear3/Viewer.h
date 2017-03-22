#ifndef VIEWER_H
#define VIEWER_H

#include <QWidget>

class Sandbox;

class Viewer : public QWidget
{
  Q_OBJECT
  QPoint m_lastMousePos, m_firstMousePos;
  int m_button;
  QPointF m_pan;
  double m_scale;

  Sandbox *m_sandbox;

public:
  Viewer(QWidget *parent = 0);
  ~Viewer();

protected:
  QSize sizeHint() const;
  void paintEvent(QPaintEvent*);

  void mousePressEvent(QMouseEvent*);
  void mouseMoveEvent(QMouseEvent*);
  void mouseReleaseEvent(QMouseEvent*);

  void wheelEvent(QWheelEvent*);

  void keyPressEvent(QKeyEvent*);

};

#endif // VIEWER_H
