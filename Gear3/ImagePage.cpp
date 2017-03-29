#include "Page.h"
#include <QLinearGradient>
#include <QPainterPath>
#include <QFontMetrics>
#include <QVector2D>
#include <QImage>
#include <qpixmap.h>



class SimpleImagePage : public Page {
  QString m_filepath;
  
public:
  SimpleImagePage(const QString &name, const QString &path) : Page(name) {
    m_filepath = path;
  }

  ~SimpleImagePage() {  }

  void draw(QPainter &pa) {
    QPixmap m_image(m_filepath);
    pa.setRenderHints(QPainter::SmoothPixmapTransform | QPainter::Antialiasing);
    double w = getWidth(), h = getHeight();
    double sc = qMin(w/m_image.width(), h/m_image.height());
    QPointF p(0.5*(w-m_image.width()*sc),0.5*(h-m_image.height()*sc));
    
    pa.drawPixmap(QRectF(p,QSizeF(m_image.width()*sc, m_image.height()*sc)), m_image, QRectF());
  }

};


SimpleImagePage fig1("fig1", "images\\geneve_wheel_patent.png");
