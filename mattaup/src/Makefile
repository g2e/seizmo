TGT = 1.5
DESTDIR = .
INSTALLDIR = ../lib
#CLASSPATH = .:/opt/taup/1.1.7/lib/TauP-1.1.7.jar:/opt/taup/1.1.7/lib/seisFile-1.0.1.jar
#CLASSPATH = .:/opt/taup/1.2beta2/lib/TauP-1.2beta2.jar:/opt/taup/1.2beta2/lib/seisFile-1.0.6.jar
CLASSPATH = .:/opt/taup/1.2beta3/lib/TauP-1.2beta3.jar:/opt/taup/1.2beta3/lib/seisFile-1.0.7.jar

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
#	cp -f $(DESTDIR)/MatTauP-1.2beta2.jar ${INSTALLDIR}/MatTauP-1.2beta2.jar
	cp -f $(DESTDIR)/MatTauP-1.2beta3.jar ${INSTALLDIR}/MatTauP-1.2beta3.jar

jar:	dist
#	cd ${DESTDIR} ; jar -cf MatTauP-1.1.7.jar edu
#	cd ${DESTDIR} ; jar -cf MatTauP-1.2beta2.jar edu
	cd ${DESTDIR} ; jar -cf MatTauP-1.2beta3.jar edu


