#pragma once

#include <QPainter>
#include <QMap>
#include <QList>
#include <QString>
#include <QElapsedTimer>

class QMouseEvent;

class Pannable {
  QPointF m_pan;
  double m_scale;
  double m_defaultScale;

public:
  Pannable() : m_scale(1.0), m_pan(0,0), m_defaultScale(1.0) {}
  virtual ~Pannable() {}

  QPointF getPanOffset() const { return m_pan; }
  void setPanOffset(const QPointF &pan) { m_pan = pan; }
  void pan(const QPointF &p) { m_pan += p; }
  void pan(double dx, double dy) { m_pan += QPointF(dx,dy); }

  double getScale() const { return m_scale; }
  void setScale(double scale) { m_scale = scale; }

  void zoom(double scale) { m_scale *= scale; }
  void zoom(const QPoint &center, double scale);
  void zoom(int cx, int cy, double scale) { zoom(QPoint(cx,cy), scale); }

  double getDefaultScale() const { return m_defaultScale; }
  void setDefaultScale(double scale) { m_defaultScale = scale; }
};

class Page
{
  QString m_name;
  static QMap<QString, Page*> *m_pages;
  int m_width, m_height;
  double m_parameter;
  QElapsedTimer m_clock;
  int m_elapsedTime;

protected:
  QPoint m_firstMousePos, m_lastMousePos;

public:
  Page(const QString &name);
  virtual ~Page();

  const QString &getName() const { return m_name; }

  int getWidth() const { return m_width; }
  void setWidth(int width) { m_width=width; }

  int getHeight() const { return m_height; }
  void setHeight(int height) { m_height=height; }

  QRect getBounds() const { return QRect(0,0,m_width,m_height); }

  virtual void draw(QPainter &pa) = 0;
  virtual void drawOverlay(QPainter &pa) {}
  virtual void tick() {}
  virtual bool onKey(int key) { return false; }

  static Page *getPage(const QString &name);

  double getParameter() const { return m_parameter; }
  void setParameter(double parameter) { m_parameter = parameter; }
  
  virtual void mousePressEvent(QMouseEvent*);
  virtual void mouseMoveEvent(QMouseEvent*);
  virtual void mouseReleaseEvent(QMouseEvent*);

  void resetTimer() { m_clock.start(); }

  virtual void drag(int dx, int dy, int modifiers);

  int getTime() const { return m_clock.elapsed(); }
  int getElapsedTime() const { return m_elapsedTime; }
  void setElapsedTime(int time) { m_elapsedTime = time; }
};

class PageManager {
  QList<Page*> m_pages;
  int m_currentIndex;
  QList<Page*> m_createdPages;
  QSize m_viewerSize;
  QElapsedTimer m_frameClock;
  
public:
  PageManager(const QString &fileName);
  ~PageManager();

  void setSize(int w, int h) { m_viewerSize = QSize(w,h); }

  Page *getCurrentPage() const { return m_pages[m_currentIndex]; }

  void goToPrevPage();
  void goToNextPage();
  void goToPage(int index);
  void tick();

  void draw(QPainter &pa, int w, int h);


private:
  void load(const QString &fileName);
  Page *createErrorPage(const QString &name);

  
};
