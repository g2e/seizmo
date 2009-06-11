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
 * MATLAB TauP_Curve interface
 * Writen by Qin Li, qinli@u.washington.edu
 * University of Washington
 * Nov, 2002
 */
 
package edu.sc.seis.TauP;
 
import java.io.*;
import java.util.Vector;

/**
  * Calculates travel time curves 
  * at known slowness samples.
  *
  * @version 1.1 Wed Feb  2 20:40:49 GMT 2000
  * @author H. Philip Crotwell
  *
  */
public class MatTauP_Curve extends TauP_Curve 
{
    protected TT_Curve tt_curve[];

    protected MatTauP_Curve() {
        super();
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
        SeismicPhase phase;
        Arrival currArrival;
        double[] dist, time, rayParams;
        double arcDistance, timeReduced;
        Arrival[] phaseArrivals;
        double maxTime = -1*Double.MAX_VALUE, minTime = Double.MAX_VALUE;

        tt_curve=new TT_Curve[phases.size()];
        for (int phaseNum=0;phaseNum<phases.size(); phaseNum++) {
            phase = (SeismicPhase)phases.elementAt(phaseNum);
            dist = phase.getDist();
            time = phase.getTime();
            rayParams = phase.getRayParams();

            tt_curve[phaseNum]=new TT_Curve(dist.length);
            tt_curve[phaseNum].phaseName=phase.getName();
            tt_curve[phaseNum].sourceDepth=depth;
            tt_curve[phaseNum].rayParam=rayParams;
            for (int i=0; i<dist.length; i++) {

                /* Here we use a trig trick to make sure the dist is 0 to PI. */
                arcDistance = Math.acos(Math.cos(dist[i]));
                if (reduceTime) {
                    timeReduced = time[i]-arcDistance/reduceVel;
                } else {
                    timeReduced = time[i];
                }
                tt_curve[phaseNum].time[i]=timeReduced;
                tt_curve[phaseNum].dist[i]=180.0/Math.PI*arcDistance;

                if (i<dist.length-1 && (rayParams[i] == rayParams[i+1]) && 
                    rayParams.length > 2) {
                    /* Here we have a shadow zone, so put a break in the curve. */
                    System.out.println("  Shadow Zone detected for phase "+phase.getName()+"\n");
                    continue;
                }

                /* Here we check to see if we cross a 180 degree mark, in which
                 * case sin changes sign. */
                if (i<dist.length-1 && Math.sin(dist[i]) > 0 && 
                    Math.sin(dist[i+1]) <0) {
                    phase.calcTime(180.0);
                    phaseArrivals = phase.getArrivals();
                    int j=0;
                    while (j<phaseArrivals.length) {
                        if ((phase.rayParams[i]-phaseArrivals[j].rayParam)*
                            (phaseArrivals[j].rayParam-phase.rayParams[i+1])>0) {
                            if (reduceTime) {
                                System.out.println("180.0  "+
                                          (float)(phaseArrivals[j].time- Math.PI/reduceVel)+"\n");
                            } else {
                                System.out.println("180.0  "+(float)phaseArrivals[j].time+"\n");
                            }
                            break;
                        }
                        j++;
                    }
                }

                /* Here we check to see if we cross a 0 degree mark, in which
                 * case sin changes sign. No need for reduce vel at 0 distance.*/
                if (i<dist.length-1 && Math.sin(dist[i]) < 0 && 
                    Math.sin(dist[i+1]) > 0) {
                    phase.calcTime(0.0);
                    phaseArrivals = phase.getArrivals();
                    int j=0;
                    while (j<phaseArrivals.length) {
                        if ((phase.rayParams[i]-phaseArrivals[j].rayParam)*
                            (phaseArrivals[j].rayParam-phase.rayParams[i+1])>0) {
                            System.out.println("0.0  "+(float)phaseArrivals[j].time+"\n");
                            break;
                        }
                        j++;
                    }
                }
            }

        }
    }

    public static TT_Curve[] run_curve(String[] args)
        throws FileNotFoundException,
        IOException,
        StreamCorruptedException,
        ClassNotFoundException,
        OptionalDataException
    {
        TT_Curve[] tt_curve=null;
        boolean doInteractive = true;
        try {
            MatTauP_Curve tauPCurve = new MatTauP_Curve();
            tauPCurve.outFile = "taup_curve.gmt";

            String[] noComprendoArgs = tauPCurve.parseCmdLineArgs(args);
            if (noComprendoArgs.length > 0) {
                for (int i=0;i<noComprendoArgs.length;i++) {
                    if (noComprendoArgs[i].equals("-help") ||
                        noComprendoArgs[i].equals("-version")) {
                        return tt_curve;
                    }
                }
                System.out.println("I don't understand the following arguments, continuing:");
                for (int i=0;i<noComprendoArgs.length;i++) {
                    System.out.print(noComprendoArgs[i]+" ");
                    if (noComprendoArgs[i].equals("-help")) {
                        System.out.println();
                        return tt_curve;
                    }
                }
                System.out.println();
                noComprendoArgs = null;
            }

            for (int i=0; i<args.length; i++) {
                if (args[i] == "-h") {
                    doInteractive = false;
                }
            }

            if (tauPCurve.DEBUG) {
                System.out.println("Done reading "+tauPCurve.modelName);
            }

            tauPCurve.init();
            if (doInteractive) {
                tauPCurve.start();
            } else {
                /* enough info given on cmd line, so just do one calc. */
                tauPCurve.depthCorrect( Double.valueOf(
                                                       tauPCurve.toolProps.getProperty("taup.source.depth", "0.0")
                                                       ).doubleValue());
                tauPCurve.calculate(tauPCurve.degrees);
                tauPCurve.printResult(tauPCurve.dos);
            }
            tt_curve=tauPCurve.tt_curve;
            tauPCurve.destroy();
            return tt_curve;
        } catch (TauModelException e) {
            System.out.println("Caught TauModelException: "+e.getMessage());
            e.printStackTrace();
        } 
        return tt_curve;
    }
    
    /** Allows TauP_Curve to run as an application. Creates an instance 
     * of TauP_Curve.
     * . */
    public static void main(String[] args)
        throws FileNotFoundException,
        IOException,
        StreamCorruptedException,
        ClassNotFoundException,
        OptionalDataException
    {
        TT_Curve[] tt_curve;
        boolean doInteractive = true;
        try {
            MatTauP_Curve tauPCurve = new MatTauP_Curve();
            tauPCurve.outFile = "taup_curve.gmt";

            String[] noComprendoArgs = tauPCurve.parseCmdLineArgs(args);
            if (noComprendoArgs.length > 0) {
                for (int i=0;i<noComprendoArgs.length;i++) {
                    if (noComprendoArgs[i].equals("-help") ||
                        noComprendoArgs[i].equals("-version")) {
                        System.exit(0);
                    }
                }
                System.out.println("I don't understand the following arguments, continuing:");
                for (int i=0;i<noComprendoArgs.length;i++) {
                    System.out.print(noComprendoArgs[i]+" ");
                    if (noComprendoArgs[i].equals("-help")) {
                        System.out.println();
                        System.exit(0);
                    }
                }
                System.out.println();
                noComprendoArgs = null;
            }

            for (int i=0; i<args.length; i++) {
                if (args[i] == "-h") {
                    doInteractive = false;
                }
            }

            if (tauPCurve.DEBUG) {
                System.out.println("Done reading "+tauPCurve.modelName);
            }

            tauPCurve.init();
            if (doInteractive) {
                tauPCurve.start();
            } else {
                /* enough info given on cmd line, so just do one calc. */
                tauPCurve.depthCorrect( Double.valueOf(
                                                       tauPCurve.toolProps.getProperty("taup.source.depth", "0.0")
                                                       ).doubleValue());
                tauPCurve.calculate(tauPCurve.degrees);
                tauPCurve.printResult(tauPCurve.dos);
            }
            tt_curve=tauPCurve.tt_curve;
            tauPCurve.destroy();
        } catch (TauModelException e) {
            System.out.println("Caught TauModelException: "+e.getMessage());
            e.printStackTrace();
        } 

    }
}
