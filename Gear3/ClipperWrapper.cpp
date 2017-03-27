#include "ClipperWrapper.h"
#include "clipper.hpp"

class ClipperWrapper::Imp {
public:

  ClipperLib::Paths result;
  double scale, iscale;

  Imp() : scale(1000.0) { iscale = 1.0/scale; }

  ClipperLib::IntPoint toIntPoint(const QPointF &p) const { return ClipperLib::IntPoint((int)(0.5+p.x()*scale), (int)(0.5+p.y()*scale)); }
  QPointF toQPointF(const ClipperLib::IntPoint &p) const { return iscale * QPointF(p.X,p.Y); }

};

ClipperWrapper::ClipperWrapper(void)
  : m_imp(new Imp)
{
}


ClipperWrapper::~ClipperWrapper(void)
{
  delete m_imp;
}


void ClipperWrapper::add(QPointF &p0, QPointF &p1, QPointF &p2)
{
  using namespace ClipperLib;
  Path triangle; 
  triangle.push_back(m_imp->toIntPoint(p0));
  triangle.push_back(m_imp->toIntPoint(p1));
  triangle.push_back(m_imp->toIntPoint(p2));
  triangle.push_back(m_imp->toIntPoint(p0));
  if(m_imp->result.empty()) 
  {
    m_imp->result.push_back(triangle);
  }
  else
  {
    Clipper clpr;
    clpr.AddPaths(m_imp->result, ptSubject, true);
    clpr.AddPath(triangle, ptClip, true);
    Paths result;
    clpr.Execute(ctUnion, result, pftEvenOdd, pftEvenOdd);   
    result.swap(m_imp->result);
  }
}

void ClipperWrapper::add(const QVector<QPointF> &outline)
{
  using namespace ClipperLib;
  Path path; 
  for(int i=0;i<outline.count();i++) path.push_back(m_imp->toIntPoint(outline[i]));
  path.push_back(m_imp->toIntPoint(outline[0]));
  if(m_imp->result.empty()) 
  {
    m_imp->result.push_back(path);
  }
  else
  {
    Clipper clpr;
    clpr.AddPaths(m_imp->result, ptSubject, true);
    clpr.AddPath(path, ptClip, true);
    Paths result;
    clpr.Execute(ctUnion, result, pftEvenOdd, pftEvenOdd);   
    result.swap(m_imp->result);
  }
}

void ClipperWrapper::sub(const QVector<QPointF> &outline)
{
 using namespace ClipperLib;
  Path path; 
  for(int i=0;i<outline.count();i++) path.push_back(m_imp->toIntPoint(outline[i]));
  path.push_back(m_imp->toIntPoint(outline[0]));

  Clipper clpr;
  clpr.AddPaths(m_imp->result, ptSubject, true);
  clpr.AddPath(path, ptClip, true);
  Paths result;
  clpr.Execute(ctDifference, result, pftEvenOdd, pftEvenOdd);   
  result.swap(m_imp->result);
}

void ClipperWrapper::intersect(const QVector<QPointF> &outline)
{
  using namespace ClipperLib;
  Path path; 
  for(int i=0;i<outline.count();i++) path.push_back(m_imp->toIntPoint(outline[i]));
  path.push_back(m_imp->toIntPoint(outline[0]));

  Clipper clpr;
  clpr.AddPaths(m_imp->result, ptClip, true);
  clpr.AddPath(path, ptSubject, true);
  Paths result;
  clpr.Execute(ctDifference, result, pftEvenOdd, pftEvenOdd);   
  result.swap(m_imp->result);
}


void ClipperWrapper::getOutline(QVector<QVector<QPointF> > &lines)
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


