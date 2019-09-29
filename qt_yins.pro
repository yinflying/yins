TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CFLAGS += -std=c99
QMAKE_CXXFLAGS += -std=c++98

SOURCES += \
        yins_core/ins.c \
        yins_core/inscmn.c \
        yins_core/insio.c \
        yins_core/inskf.c

config_nmea2ycsv{
    TARGET = nmea2ycsv
    SOURCES +=
}

config_UT{
    TARGET = unittest
    SOURCES += yins/unittest.c
}

config_main{
    TARGET = main
    SOURCES += main.cpp
}

HEADERS += \
    yins/ins.h

unix:!macx: LIBS += -lcriterion

DISTFILES += \
    Doxyfile
