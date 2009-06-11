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
    protected MatArrival matArrivals[];
    
    protected MatTauP_Path() 
    {
        super();
        outFile = null;
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
            } catch(FileNotFoundException ee) {
                Alert.error("Can't find saved model file for model "+
                            modelName+".", "");
                return;
            } catch (InvalidClassException ee) {
                Alert.error("Model file "+
                            modelName+" is not compatible with the current version.",
                            "Recreate using taup_create.");
                return;
            }
        }
        dos=null;
    }

    public void printResult(DataOutputStream dos) throws IOException 
    {
        Writer s = null;
        printResult(s);
    }

    public void printResult(Writer out) throws IOException 
    {
        double calcTime, calcDist, calcDepth;
        LatLon coord;
        Arrival currArrival;
        int interpNum, maxInterpNum;
        double interpInc;
        double radiusOfEarth = tModDepth.getRadiusOfEarth();
        boolean longWayRound;

        matArrivals=new MatArrival[arrivals.size()];
        for (int i=0;i<arrivals.size();i++) 
        {
            currArrival = (Arrival)arrivals.elementAt(i);
            matArrivals[i] = new MatArrival(currArrival);
            matArrivals[i].matPath.init_path(currArrival.path.length);
            longWayRound = false;
            if ((currArrival.dist*180/Math.PI) % 360 > 180) 
            {
                longWayRound = true;
            }

            calcTime = 0.0;
            calcDist = 0.0;
            calcDepth = currArrival.sourceDepth;
            for (int j=0; j< currArrival.path.length; j++) 
            {
                calcTime = currArrival.path[j].time;
                calcDepth = currArrival.path[j].depth;
                calcDist = currArrival.path[j].dist *180.0/Math.PI;
                if (longWayRound && calcDist != 0.0) 
                {
                   calcDist =  -1.0*calcDist;
                }
                matArrivals[i].matPath.p[j]=currArrival.path[j].p;
                matArrivals[i].matPath.time[j]=calcTime;
                matArrivals[i].matPath.dist[j]=calcDist;
                matArrivals[i].matPath.depth[j]=calcDepth;
                matArrivals[i].matPath.lat[j]=(calcLatLon(calcDist)).lat;
                matArrivals[i].matPath.lon[j]=(calcLatLon(calcDist)).lon;
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
   
    /** Allows TauP_Path to run as an application. Creates an instance of 
     * TauP_Path and calls TauP_Path.init() and TauP_Path.start(). */
    public static MatArrival[] run_path(String[] args)
        throws FileNotFoundException,
        IOException,
        StreamCorruptedException,
        ClassNotFoundException,
        OptionalDataException
    {
        MatArrival[] matArrivals=null;
	try {
	    MatTauP_Path tauPPath = new MatTauP_Path();
	    String[] noComprendoArgs = tauPPath.parseCmdLineArgs(args);
	    if (noComprendoArgs.length > 0) {
		for (int i=0;i<noComprendoArgs.length;i++) {
		    if (noComprendoArgs[i].equals("-help") ||
			noComprendoArgs[i].equals("-version")) {
			return matArrivals;
		    }
		}
		System.out.println("I don't understand the following arguments, continuing:");
		for (int i=0;i<noComprendoArgs.length;i++) {
		    System.out.print(noComprendoArgs[i]+" ");
		    if (noComprendoArgs[i].equals("-help") ||
			noComprendoArgs[i].equals("-version")) {
			System.out.println();
			return matArrivals;
		    }
		}
		System.out.println();
		noComprendoArgs = null;
	    }

	    tauPPath.init();
	    if (tauPPath.DEBUG) {
		System.out.println("Done reading "+tauPPath.modelName);
	    }
 
	    tauPPath.start();
        matArrivals=tauPPath.matArrivals;
	    tauPPath.destroy();
        return matArrivals;
 
	} catch (TauModelException e) {
	    System.out.println("Caught TauModelException: "+e.getMessage());
	    e.printStackTrace();
	} catch (TauPException e) {
	    System.out.println("Caught TauPException: "+e.getMessage());
	    e.printStackTrace();
	}
    return matArrivals;
    }

    
    public static void main(String[] args)
        throws FileNotFoundException,
        IOException,
        StreamCorruptedException,
        ClassNotFoundException,
        OptionalDataException
    {
        MatArrival[] matArrivals=null;
	try {
	    MatTauP_Path tauPPath = new MatTauP_Path();
	    String[] noComprendoArgs = tauPPath.parseCmdLineArgs(args);
	    if (noComprendoArgs.length > 0) {
		for (int i=0;i<noComprendoArgs.length;i++) {
		    if (noComprendoArgs[i].equals("-help") ||
			noComprendoArgs[i].equals("-version")) {
			return ;
		    }
		}
		System.out.println("I don't understand the following arguments, continuing:");
		for (int i=0;i<noComprendoArgs.length;i++) {
		    System.out.print(noComprendoArgs[i]+" ");
		    if (noComprendoArgs[i].equals("-help") ||
			noComprendoArgs[i].equals("-version")) {
			System.out.println();
			return ;
		    }
		}
		System.out.println();
		noComprendoArgs = null;
	    }

	    tauPPath.init();
	    if (tauPPath.DEBUG) {
		System.out.println("Done reading "+tauPPath.modelName);
	    }
 
	    tauPPath.start();
        matArrivals=tauPPath.matArrivals;
	    tauPPath.destroy();
        return ;
 
	} catch (TauModelException e) {
	    System.out.println("Caught TauModelException: "+e.getMessage());
	    e.printStackTrace();
	} catch (TauPException e) {
	    System.out.println("Caught TauPException: "+e.getMessage());
	    e.printStackTrace();
	}
    return ;
    }

}
