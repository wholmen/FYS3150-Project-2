TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -larmadillo -lblas -llapack
SOURCES += main.cpp \
    lib.cpp \
    Jacobi_function.cpp

HEADERS += \
    lib.h \
    Jacobi_function.h

