TGT = 1.6
DESTDIR = .
INSTALLDIR = ../lib
CLASSPATH = .:../lib/TauP-2.1.1.jar:../lib/seisFile-1.5.1.jar

DISTOBJS =	MatCurve.class \
	MatTauP_Curve.class \
	MatPath.class \
	MatArrival.class \
	LatLon.class \
	MatTauP_Path.class

OBJS = ${DISTOBJS}

.SUFFIXES : .class .html .java

.java.class :
	javac -target ${TGT} -classpath ${CLASSPATH} \
		-d ${DESTDIR} $<

dist:	${OBJS}

install: jar
	cp -f $(DESTDIR)/MatTauP-2.1.1.jar ${INSTALLDIR}/MatTauP-2.1.1.jar

jar:	dist
	cd ${DESTDIR} ; jar -cf MatTauP-2.1.1.jar edu



