TGT = 1.5
DESTDIR = .
INSTALLDIR = ../lib
#CLASSPATH = .:/opt/taup/1.1.7/lib/TauP-1.1.7.jar:/opt/taup/1.1.7/lib/seisFile-1.0.1.jar
CLASSPATH = .:/opt/taup/2.0/lib/TauP-2.0.jar:/opt/taup/2.0/lib/seisFile-1.2.jar

DISTOBJS =	LatLon.class \
	MatPath.class \
	MatArrival.class \
	TT_Curve.class \
	MatTauP_Curve.class \
	MatTauP_Path.class \
	MatTauP_Pierce.class \
	MatTauP_Time.class

OBJS = ${DISTOBJS}

.SUFFIXES : .class .html .java

.java.class :
	javac -target ${TGT} -classpath ${CLASSPATH} \
		-d ${DESTDIR} $<

dist:	${OBJS}

install: jar
#	cp -f $(DESTDIR)/MatTauP-1.1.7.jar ${INSTALLDIR}/MatTauP-1.1.7.jar
	cp -f $(DESTDIR)/MatTauP-2.0.jar ${INSTALLDIR}/MatTauP-2.0.jar

jar:	dist
#	cd ${DESTDIR} ; jar -cf MatTauP-1.1.7.jar edu
	cd ${DESTDIR} ; jar -cf MatTauP-2.0.jar edu


