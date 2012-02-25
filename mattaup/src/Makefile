TGT = 1.5
DESTDIR = .
INSTALLDIR = ../lib
CLASSPATH = .:/opt/taup/1.2beta4/lib/TauP-1.2beta4.jar:/opt/taup/1.2beta4/lib/seisFile-1.0.8.jar

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
	cp -f $(DESTDIR)/MatTauP-1.2beta4.jar ${INSTALLDIR}/MatTauP-1.2beta4.jar

jar:	dist
	cd ${DESTDIR} ; jar -cf MatTauP-1.2beta4.jar edu



