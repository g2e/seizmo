/*
  The TauP Toolkit: Flexible Seismic Travel-Time and Raypath Utilities.
  Copyright (C) 1998-2000 University of South Carolina

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  The current version can be found at 
  <A HREF="www.seis.sc.edu">http://www.seis.sc.edu</A>

  Bug reports and comments should be directed to 
  H. Philip Crotwell, crotwell@seis.sc.edu or
  Tom Owens, owens@seis.sc.edu

*/
 
 /**
 * MATLAB TauP_Path interface
 * Writen by Qin Li, qinli@u.washington.edu
 * University of Washington
 * Nov, 2002
 */

package edu.sc.seis.TauP;
 
import java.io.*;
import java.util.Vector;
import edu.sc.seis.TauP.*;

/**
  * Calculate travel paths for different phases using a linear interpolated
  * ray parameter between known slowness samples.
  *
  * @version 1.1 Wed Feb  2 20:40:49 GMT 2000
  * @author H. Philip Crotwell
  *
  */
public class MatTauP_Path extends TauP_Path 
{
    protected static double maxPathInc = 1.0;
    protected static double eventLat = 0.0;
    protected static double eventLon = 0.0;
    protected static double azimuth = 0.0;

    protected MatArrival matArrivals[];
    
    protected MatTauP_Path() 
    {
        super();
        outFile = null;
    }

    public MatTauP_Path(TauModel tMod) throws TauModelException {
        super(tMod);
        outFile = null;
    }

    public MatTauP_Path(String modelName) throws TauModelException {
        super(modelName);
        outFile = null;
    }

    public MatArrival getMatArrival(int i)
    {
        return matArrivals[i];
    }

    public void init() throws IOException 
    {
        if (phaseNames.size() == 0) {
            if ( toolProps.containsKey("taup.phase.file")) {
                if ( toolProps.containsKey("taup.phase.list")) {
                    parsePhaseList( toolProps.getProperty("taup.phase.list"));
                }
                try {
                    readPhaseFile(toolProps.getProperty("taup.phase.file"));
                } catch (IOException e) {
                    Alert.warning(
                                  "Caught IOException while attempting to reading phase file "+
                                  toolProps.getProperty("taup.phase.file"), e.getMessage());
                    if (phaseNames.size() <= 0) {
                        parsePhaseList(toolProps.getProperty("taup.phase.list",
                                                             "p,s,P,S,Pn,Sn,PcP,ScS,Pdiff,Sdiff,PKP,SKS,PKiKP,SKiKS,PKIKP,SKIKS"));
                    }
                }
            } else {
                parsePhaseList( toolProps.getProperty("taup.phase.list",
                                                      "p,s,P,S,Pn,Sn,PcP,ScS,Pdiff,Sdiff,PKP,SKS,PKiKP,SKiKS,PKIKP,SKIKS"));
            }
        }

        depth = Double.valueOf(toolProps.getProperty("taup.source.depth", "0.0")
                               ).doubleValue();

        if (tMod == null) {
            modelName = toolProps.getProperty("taup.model.name", "iasp91");
            try {
                readTauModel();
            } catch (TauModelException ee) {
                Alert.error("Caught TauModelException", ee.getMessage());
            }
        }
    }

    public void pathInterpolate() 
    {
        double calcTime, calcDist, calcDepth;
        LatLon coord;
        Arrival currArrival;
        int interpNum, maxInterpNum, count, k;
        double interpInc;
        double radiusOfEarth = tModDepth.getRadiusOfEarth();
        boolean longWayRound;

        matArrivals=new MatArrival[arrivals.size()];
        for (int i=0;i<arrivals.size();i++) 
        {
            currArrival = (Arrival)arrivals.get(i);
            matArrivals[i] = new MatArrival(currArrival);

            // we need to figure out how many interpolated points will be added (gge)
            count=currArrival.path.length;
            for (int j=0; j< currArrival.path.length; j++) 
            {
                if(j < currArrival.path.length - 1
                        && (currArrival.rayParam != 0.0 && 
                                (currArrival.path[j + 1].getDistDeg() - currArrival.path[j].getDistDeg()) > maxPathInc))
                {
                    count=count-1+(int)Math.ceil((currArrival.path[j + 1].getDistDeg() - currArrival.path[j].getDistDeg())
                            / maxPathInc);
                }
            }

            matArrivals[i].matPath.init_path(count);
            longWayRound = false;
            if ((currArrival.dist*180/Math.PI) % 360 > 180) 
            {
                longWayRound = true;
            }

            k=0; // count number of extra points for interpolation (gge)
            calcTime = 0.0;
            calcDist = 0.0;
            calcDepth = currArrival.sourceDepth;
            for (int j=0; j< currArrival.path.length; j++) 
            {
                calcTime = currArrival.path[j].time;
                calcDepth = currArrival.path[j].depth;
                calcDist = currArrival.path[j].getDistDeg();
                if (longWayRound && calcDist != 0.0) 
                {
                   calcDist =  -1.0*calcDist;
                }
                matArrivals[i].matPath.time[j+k]=calcTime;
                matArrivals[i].matPath.dist[j+k]=calcDist;
                matArrivals[i].matPath.depth[j+k]=calcDepth;
                matArrivals[i].matPath.lat[j+k]=(calcLatLon(calcDist)).lat;
                matArrivals[i].matPath.lon[j+k]=(calcLatLon(calcDist)).lon;

                // interpolate to steps of at most maxPathInc degrees for path (gge)
                if(j < currArrival.path.length - 1
                        && (currArrival.rayParam != 0.0 && 
                                (currArrival.path[j + 1].getDistDeg() - currArrival.path[j].getDistDeg()) > maxPathInc)) {
                    maxInterpNum = (int)Math.ceil((currArrival.path[j + 1].getDistDeg() - currArrival.path[j].getDistDeg())
                            / maxPathInc);
                    for(interpNum = 1; interpNum < maxInterpNum; interpNum++) {
                        k=k+1;
                        calcTime += (currArrival.path[j + 1].time - currArrival.path[j].time)
                                / maxInterpNum;
                        if(longWayRound) {
                            calcDist -= (currArrival.path[j + 1].getDistDeg() - currArrival.path[j].getDistDeg())
                                    / maxInterpNum;
                        } else {
                            calcDist += (currArrival.path[j + 1].getDistDeg() - currArrival.path[j].getDistDeg())
                                    / maxInterpNum;
                        }
                        calcDepth = currArrival.path[j].depth + interpNum
                                * (currArrival.path[j + 1].depth - currArrival.path[j].depth)
                                / maxInterpNum;

                        matArrivals[i].matPath.time[j+k]=calcTime;
                        matArrivals[i].matPath.dist[j+k]=calcDist;
                        matArrivals[i].matPath.depth[j+k]=calcDepth;
                        matArrivals[i].matPath.lat[j+k]=(calcLatLon(calcDist)).lat;
                        matArrivals[i].matPath.lon[j+k]=(calcLatLon(calcDist)).lon;
                    }
                }
            } //for
        }
    }

    protected LatLon calcLatLon(double calcDist) 
    {
	double lat=0, lon=0;
    LatLon pathCoord=new LatLon();
    		
	if ( eventLat != Double.MAX_VALUE &&
	     eventLon != Double.MAX_VALUE &&
	     azimuth != Double.MAX_VALUE) 
    {
	    lat = SphericalCoords.latFor(eventLat, eventLon, 
					 calcDist, azimuth);
	    lon = SphericalCoords.lonFor(eventLat, eventLon,
					 calcDist, azimuth);
	} else if ( stationLat != Double.MAX_VALUE &&
		    stationLon != Double.MAX_VALUE &&
		    backAzimuth != Double.MAX_VALUE) {
	    lat = SphericalCoords.latFor(stationLat, stationLon,
					 degrees-calcDist, backAzimuth);
	    lon = SphericalCoords.lonFor(stationLat, stationLon,
					 degrees-calcDist, backAzimuth);
        } else if (stationLat != Double.MAX_VALUE &&  stationLon != Double.MAX_VALUE &&
                    eventLat != Double.MAX_VALUE && eventLon != Double.MAX_VALUE) 
        {
            azimuth = SphericalCoords.azimuth(eventLat, eventLon, stationLat, stationLon);
            backAzimuth = SphericalCoords.azimuth(stationLat, stationLon, eventLat, eventLon);
            lat = SphericalCoords.latFor(eventLat, eventLon,calcDist, azimuth);
            lon = SphericalCoords.lonFor(eventLat, eventLon, calcDist, azimuth);
        }
        pathCoord.lat=lat;
        pathCoord.lon=lon;
        return pathCoord;
    }

    public void setEv(double eventLat, double eventLon) {
        this.eventLat=eventLat;
        this.eventLon=eventLon;
    }
    public void setAz(double azimuth) {
        this.azimuth=azimuth;
    }
}
