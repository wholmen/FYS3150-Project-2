TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -larmadillo -lblas -llapack
SOURCES += main.cpp \
    lib.cpp

HEADERS += \
    lib.h

