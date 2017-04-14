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


SimpleImagePage fig1("fig_geneve_wheel_patent", "images\\geneve_wheel_patent.png");
SimpleImagePage fig2("fig_oval_patent", "images\\oval_patent.png");
SimpleImagePage fig3("fig_fluxometer", "images\\fluxometer.png");
SimpleImagePage fig4("fig_maltese_cross", "images\\maltese_cross2.png");
SimpleImagePage fig6("fig_laser_cutter", "images\\laser_cutter.png");
SimpleImagePage fig7("fig_upfixes", "images\\upfixes.png");
SimpleImagePage fig8("fig_italy", "images\\italy.png");
SimpleImagePage fig9("fig_gearify", "images\\gearify.png");

SimpleImagePage fig10("fig_cat1", "images\\cat_1.png");
SimpleImagePage fig11("fig_cat2", "images\\cat_2.png");
SimpleImagePage fig12("fig_cat3", "images\\cat_3.png");
SimpleImagePage fig13("fig_cat4", "images\\cat_4.png");
SimpleImagePage fig14("fig_cat5", "images\\cat_5.png");
SimpleImagePage fig15("fig_cat6", "images\\cat_6.png");
SimpleImagePage fig16("fig_cat7", "images\\cat_7.png");
