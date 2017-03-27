#include "Viewer.h"
#include <QPainter>
#include <QPainterPath>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QWheelEvent>
#include <qmath.h>

#include "PitchCurve.h"
#include "Sandbox.h"


Viewer::Viewer(QWidget *parent)
: QWidget(parent)
, m_button(Qt::NoButton)
, m_pan(0,0)
, m_scale(1.0)
, m_pageMngr("pageList.txt")
, m_timesRecord(100)
, m_firstDraw(true)
{
  m_sandbox = new Sandbox();
  startTimer(40);
  m_timer.start();
}

Viewer::~Viewer()
{
  delete m_sandbox;
}


QSize Viewer::sizeHint() const
{
  return QSize(1024,768);
}



void Viewer::paintEvent(QPaintEvent*)
{
  if(m_firstDraw) { m_firstDraw = false; m_pageMngr.goToPage(0); }
  int ms2 = m_timer.restart();
  QElapsedTimer timer1; 
  timer1.start();

  QPainter pa(this);
  pa.setRenderHints(QPainter::Antialiasing);

  Pannable *pannable = dynamic_cast<Pannable*>(getCurrentPage());
  if(pannable) 
  {
    pa.save();
    QPointF p = pannable->getPanOffset();
    pa.translate(p.x(),p.y());
    double s = pannable->getScale();
    pa.scale(s,s);
  }
  //pa.translate(m_pan + QPointF(width(),height())*0.5);
  //pa.scale(m_scale,m_scale);

  // m_sandbox->paint(pa);  
  m_pageMngr.draw(pa, width(), height());

  if(pannable) { pa.restore(); getCurrentPage()->drawOverlay(pa); }

  drawFps(pa);

  int ms1 = timer1.elapsed();
  m_timesRecord.add(qMakePair(ms1,ms2));
}

void Viewer::resizeEvent(QResizeEvent*)
{
  m_pageMngr.setSize(width(), height());
}


void Viewer::drawFps(QPainter &pa)
{
  int x0 = 110, y0 = 110;

  int n = m_timesRecord.getCount();
  pa.setPen(Qt::cyan);
  for(int i=0;i<n;i++)
  {
    int x = x0-i;
    pa.drawLine(x,y0,x,y0-m_timesRecord.getValue(i).second);
  }
  pa.setPen(Qt::blue);
  for(int i=0;i<n;i++)
  {
    int x = x0-i;
    pa.drawLine(x,y0,x,y0-m_timesRecord.getValue(i).first);
  }
  pa.setPen(Qt::gray);
  for(int i=0;i<100;i+=10) pa.drawLine(x0-100,y0-i,x0,y0-i);
  
}



void Viewer::mousePressEvent(QMouseEvent*e)
{
  m_lastMousePos = m_firstMousePos = e->pos();
  m_button = e->button();
  if(m_button == Qt::LeftButton) getCurrentPage()->mousePressEvent(e);
}


void Viewer::mouseMoveEvent(QMouseEvent*e)
{
  if(m_button == Qt::LeftButton) getCurrentPage()->mouseMoveEvent(e);
  else
  {
    QPoint delta = e->pos() - m_lastMousePos;
    m_lastMousePos = e->pos();
    Pannable *pannable = dynamic_cast<Pannable*>(getCurrentPage());
    if(pannable) pannable->pan(delta);
  }
}


void Viewer::mouseReleaseEvent(QMouseEvent*e)
{
  if(m_button == Qt::LeftButton) getCurrentPage()->mouseReleaseEvent(e);
  m_button = Qt::NoButton;
}


void Viewer::wheelEvent(QWheelEvent*e)
{
  Pannable *pannable = dynamic_cast<Pannable*>(getCurrentPage());
  if(pannable) pannable->zoom(e->pos(),  exp(e->delta()*0.001));

  /*
  {
    pan
    
  QPointF center = QPointF(width(), height())*0.5;
  
  */

  // update();
}


void Viewer::keyPressEvent(QKeyEvent*e)
{
  if(e->key() == Qt::Key_Down) m_pageMngr.goToNextPage();
  else if(e->key() == Qt::Key_Up) m_pageMngr.goToPrevPage();
  else 
  {
    if(!getCurrentPage()->onKey(e->key())) e->ignore();
  }
}

void Viewer::timerEvent(QTimerEvent*)
{
  m_pageMngr.tick();
  update();
}



//=============================================================================

template<class T>
void MyQueue<T>::add(const T &v)
{
  const int n = m_queue.count();
  if(m_count<n) m_count++;
  m_queue[m_index] = v;
  m_index = (m_index + 1) % n; 
}

template<class T>
const T &MyQueue<T>::getValue(int index) const
{
  const int n = m_queue.count();
  if(index>=m_count) index = m_count-1;
  int j = (m_index + n - index-1) % n;
  return m_queue[j];
}
