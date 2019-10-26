TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CFLAGS += -std=c99
QMAKE_CXXFLAGS += -std=c++98

SOURCES += \
        main.c \
        yins_core/ins.c \
        yins_core/inscmn.c \
        yins_core/insio.c \
        yins_core/inskf.c \
        yinsapp.c

config_nmea2ycsv{
    TARGET = nmea2ycsv
    SOURCES +=
}

config_main{
    TARGET = main
    SOURCES +=
}

config_UT{
    TARGET = unittest
    SOURCES += yins_core/unittest.c
}


HEADERS += \
    yins_core/ins.h \
    yins_core/insmacro.h \
    yinsapp.h

unix:!macx: LIBS += -lcriterion

DISTFILES += \
    Doxyfile
