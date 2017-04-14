#ifndef VIEWER_H
#define VIEWER_H

#include <QWidget>
#include "Page.h"
#include <QElapsedTimer>
#include <QVector>
#include <QPair>
#include <QGLWidget>

class Sandbox;

template<class T>
class MyQueue {
  QVector<T> m_queue;
  int m_index, m_count;

public:
  MyQueue(int n) : m_queue(n), m_index(0), m_count(0) {}
  ~MyQueue() {}

  void add(const T &v);
  int getCount() const { return m_count; }
  const T &getValue(int index) const;
};



class Viewer : public QWidget
{
  Q_OBJECT
  QPoint m_lastMousePos, m_firstMousePos;
  int m_button;
  QPointF m_pan;
  double m_scale;

  Sandbox *m_sandbox;
  PageManager m_pageMngr;

  QElapsedTimer m_timer;
  MyQueue<QPair<int,int> > m_timesRecord;
  bool m_firstDraw;
  bool m_showFps;

public:
  Viewer(QWidget *parent = 0);
  ~Viewer();

  Page *getCurrentPage() { return m_pageMngr.getCurrentPage(); }

protected:
  QSize sizeHint() const;
  void paintEvent(QPaintEvent*);
  void resizeEvent(QResizeEvent*);

  void mousePressEvent(QMouseEvent*);
  void mouseMoveEvent(QMouseEvent*);
  void mouseReleaseEvent(QMouseEvent*);

  void wheelEvent(QWheelEvent*);

  void keyPressEvent(QKeyEvent*);

  void timerEvent(QTimerEvent*);


  void drawFps(QPainter &);
};

#endif // VIEWER_H
