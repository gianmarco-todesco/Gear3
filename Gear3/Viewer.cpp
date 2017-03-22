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
{
  m_sandbox = new Sandbox();
}

Viewer::~Viewer()
{
  delete m_sandbox;
}


QSize Viewer::sizeHint() const
{
  return QSize(1000,800);
}



void Viewer::paintEvent(QPaintEvent*)
{
  QPainter pa(this);
  pa.setRenderHints(QPainter::Antialiasing);

  pa.translate(m_pan + QPointF(width(),height())*0.5);
  pa.scale(m_scale,m_scale);

  m_sandbox->paint(pa);  
}


void Viewer::mousePressEvent(QMouseEvent*e)
{
  m_lastMousePos = m_firstMousePos = e->pos();
  m_button = e->button();
}


void Viewer::mouseMoveEvent(QMouseEvent*e)
{
  QPoint delta = e->pos() - m_lastMousePos;
  m_lastMousePos = e->pos();
  if(m_button == Qt::LeftButton)
  {
    m_sandbox->changeParameter(delta.x());
  }
  else
  {
    m_pan += delta;
  }
  update();
}


void Viewer::mouseReleaseEvent(QMouseEvent*e)
{
  m_button = Qt::NoButton;
}


void Viewer::wheelEvent(QWheelEvent*e)
{
  QPointF center = QPointF(width(), height())*0.5;
  QPointF p = e->pos();
  QPointF wp = (p-(m_pan+center))/m_scale;
  m_scale *= exp(e->delta()*0.001);
  m_pan = p - wp*m_scale - center;
  update();
}


void Viewer::keyPressEvent(QKeyEvent*e)
{
  e->ignore();
}
