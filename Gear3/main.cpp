#include "Viewer.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Viewer w;
    w.show();
    return a.exec();
}
