#include "ClipperWrapper.h"
#include "clipper.hpp"

class ClipperWrapper::Imp {
public:

  ClipperLib::Paths result;
  double scale, iscale;

  Imp() : scale(1000.0) { iscale = 1.0/scale; }

  ClipperLib::IntPoint toIntPoint(const QPointF &p) const { return ClipperLib::IntPoint((int)(0.5+p.x()*scale), (int)(0.5+p.y()*scale)); }
  QPointF toQPointF(const ClipperLib::IntPoint &p) const { return iscale * QPointF(p.X,p.Y); }

  void toPath(ClipperLib::Path &path, const QVector<QPointF> &pts) { for(int i=0;i<pts.count();i++) path.push_back(toIntPoint(pts[i])); } 
  void toPaths(ClipperLib::Paths &paths, const QVector<QVector<QPointF> > &lines) { 
    for(int i=0;i<lines.count();i++) {
      ClipperLib::Path path; toPath(path, lines[i]);
      paths.push_back(path);
    }
  } 

  void execute(ClipperLib::ClipType op, const QVector<QPointF> &outline);
  void execute(ClipperLib::ClipType op, const  QVector<QVector<QPointF> > &outline);


};

void ClipperWrapper::Imp::execute(ClipperLib::ClipType op, const QVector<QPointF> &outline)
{
  using namespace ClipperLib;
  Path path; toPath(path, outline);
  if(result.empty()) 
  {
    Q_ASSERT(op == ctUnion);
    result.push_back(path);
  }
  else
  {
    Clipper clpr;
    clpr.AddPaths(result, ptSubject, true);
    clpr.AddPath(path, ptClip, true);
    Paths tmp;
    clpr.Execute(op, tmp, pftEvenOdd, pftEvenOdd);   
    result.swap(tmp);
  }
}


void ClipperWrapper::Imp::execute(ClipperLib::ClipType op, const QVector<QVector<QPointF> > &outline)
{
  using namespace ClipperLib;
  Paths paths; toPaths(paths, outline);
  if(result.empty()) 
  {
    Q_ASSERT(op == ctUnion);
    result = paths;
  }
  else
  {
    Clipper clpr;
    clpr.AddPaths(result, ptSubject, true);
    clpr.AddPaths(paths, ptClip, true);
    Paths tmp;
    clpr.Execute(op, tmp, pftEvenOdd, pftEvenOdd);   
    result.swap(tmp);
  }
}



ClipperWrapper::ClipperWrapper(void)
  : m_imp(new Imp)
{
}


ClipperWrapper::~ClipperWrapper(void)
{
  delete m_imp;
}



void ClipperWrapper::add(const QVector<QPointF> &outline)
{
  m_imp->execute(ClipperLib::ctUnion, outline);
}

void ClipperWrapper::add(const QVector<QVector<QPointF> > &outline)
{
  m_imp->execute(ClipperLib::ctUnion, outline);
}

void ClipperWrapper::subtract(const QVector<QPointF> &outline)
{
  m_imp->execute(ClipperLib::ctDifference, outline);
}

void ClipperWrapper::subtract(const QVector<QVector<QPointF> > &outline)
{
  m_imp->execute(ClipperLib::ctDifference, outline);
}

void ClipperWrapper::add(QPointF &p0, QPointF &p1, QPointF &p2)
{
  QVector<QPointF> pts; pts << p0 << p1 << p2;
  add(pts);
}


void ClipperWrapper::getResult(QVector<QVector<QPointF> > &lines)
{
  using namespace ClipperLib;
  for(int i=0; i<(int)m_imp->result.size();i++)
  {
    Path &path = m_imp->result[i];
    QVector<QPointF> line(path.size());
    for(int j=0;j<path.size();j++) line[j] = m_imp->toQPointF(path[j]);
    lines.push_back(line);
  }
}

