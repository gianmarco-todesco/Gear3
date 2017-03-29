#include "Page.h"
#include <QFile>
#include <QTextStream>
#include <QPainter>
#include <QDebug>
#include <QMouseEvent>

Page::Page(const QString &name)
  : m_name(name)
  , m_parameter(0)
{
  if(!m_pages) m_pages = new QMap<QString, Page*>();
  (*m_pages)[name] = this;
}


Page::~Page()
{
}

Page *Page::getPage(const QString &name)
{
  if(!m_pages) return 0;
  return (*m_pages)[name];
}

void Page::mousePressEvent(QMouseEvent*e)
{
  m_lastMousePos = m_firstMousePos = e->pos();
}

void Page::mouseMoveEvent(QMouseEvent*e)
{
  QPoint delta = e->pos() - m_lastMousePos; 
  m_lastMousePos = e->pos();
  drag(delta.x(), delta.y(), e->modifiers());
}

void Page::mouseReleaseEvent(QMouseEvent*e)
{
}

void Page::drag(int dx, int dy, int modifiers)
{
  m_parameter += dx;
}



QMap<QString, Page*> *Page::m_pages = 0;

//=============================================================================

class DummyPage : public Page {
public:
  DummyPage(const QString &name) : Page(name) {}
  void draw(QPainter &pa) {
    pa.setFont(QFont("Arial", 50));
    pa.setPen(Qt::black);
    int w = getWidth(), h = getHeight();
    pa.drawText(QRect(-w/2,-h/2,w,h), Qt::AlignCenter, getName());
  }
};

//=============================================================================


void Pannable::zoom(const QPoint &center, double scale)
{
  QPointF wp = QPointF(center - m_pan) * (1.0/ m_scale);
  m_scale *= scale;
  m_pan = center - wp*m_scale;
  
}


//=============================================================================

PageManager::PageManager(const QString &fileName)
  : m_currentIndex(0)
{
  load(fileName);
  m_frameClock.start();
}

PageManager::~PageManager()
{
  for(int i=0;i<m_createdPages.count(); i++) delete m_createdPages.at(i);
  m_createdPages.clear();
}

Page *PageManager::createErrorPage(const QString &name)
{
  Page *page = Page::getPage(name);
  if(!page) { page = new DummyPage(name); m_createdPages.append(page); }
  return page;
}

void PageManager::load(const QString &fileName)
{
  QFile inputFile(fileName);
  if (inputFile.open(QIODevice::ReadOnly))
  {
    QTextStream in(&inputFile);
    while (!in.atEnd())
    {
      QString name = in.readLine();
      if(name == "" || name.startsWith("!")) continue;
      Page *page = Page::getPage(name);
      if(!page) page = createErrorPage(name + " not found");
      m_pages.append(page);
    }
    inputFile.close();
  }
  else
  {
    qDebug() << "Could not open " << fileName;
    m_pages.append(createErrorPage("Page list not found"));
  }
}


void PageManager::goToPrevPage()
{
  goToPage(m_currentIndex-1);
}

void PageManager::goToNextPage()
{
  goToPage(m_currentIndex+1);
}

void PageManager::tick()
{
  if(m_currentIndex>=0)
    m_pages[m_currentIndex]->tick();
}

void PageManager::draw(QPainter &pa, int w, int h)
{
  Page *page = getCurrentPage();
  page->setElapsedTime(m_frameClock.restart());
  page->setWidth(w);
  page->setHeight(h);
  page->draw(pa);
}

void PageManager::goToPage(int index)
{
  if(index<0 || index>=m_pages.count()) return;
  Page *page = m_pages[index];

  Pannable *pannable = dynamic_cast<Pannable *>(page);
  if(pannable)
  {
    pannable->setScale(pannable->getDefaultScale());
    pannable->setPanOffset(QPointF(m_viewerSize.width()*0.5, m_viewerSize.height()*0.5));
  }
  page->resetTimer();
  page->setElapsedTime(0);
  
  m_currentIndex = index;
}

