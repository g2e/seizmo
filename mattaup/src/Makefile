TGT = 1.5
DESTDIR = .
INSTALLDIR = ../lib
CLASSPATH = .:/opt/taup/1.2beta2/lib/TauP-1.2beta2.jar:/opt/taup/1.2beta2/lib/seisFile-1.0.6.jar

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
	cp -f $(DESTDIR)/MatTauP-1.2beta2.jar ${INSTALLDIR}/MatTauP-1.2beta2.jar

jar:	dist
	cd ${DESTDIR} ; jar -cf MatTauP-1.2beta2.jar edu



