#include "Viewer.h"
#include <QtGui/QApplication>
#include <Phonon/VideoPlayer>
#include <QUrl>
#include <QFontDatabase>


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    using namespace Phonon;

    int ret = QFontDatabase::addApplicationFont("ARHERMANN.ttf");

    Viewer w;
    w.show();
    w.resize(1024,768);
    /*
    // see http://doc.qt.io/qt-4.8/qt-multimedia-videographicsitem-videoplayer-cpp.html
VideoPlayer *player = new VideoPlayer(Phonon::VideoCategory, 0);
QObject::connect(player, SIGNAL(finished()), player, SLOT(deleteLater()));
player->play(QUrl::fromLocalFile("C:\\Users\\fw552fw131\\Videos\\ennagon1.avi"));

player->show();
*/
    return a.exec();
}
