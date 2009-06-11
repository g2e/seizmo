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
 * MATLAB TauP_Time interface
 * Writen by Qin Li, qinli@u.washington.edu
 * University of Washington
 * Nov, 2002
 */

package edu.sc.seis.TauP;

import java.io.*;
import java.util.*;
import java.util.zip.*;
import edu.sc.seis.TauP.*;

/**
  *  Calculate travel times for different branches using linear interpolation
  * between known slowness samples.
  *
  * @version 1.1 Wed Feb  2 20:40:49 GMT 2000
  * @author H. Philip Crotwell
  *
  */
public class MatTauP_Time extends TauP_Time{

    /* Constructors  */
    public MatTauP_Time() 
    {
        super();
    }

    /** Reads the velocity model, slowness model, and tau model from
     * a file saved using Java's Serializable interface.
     * Performs a depth correction if the current depth is not 0.0
     */
    protected void readTauModel()
        throws FileNotFoundException,
        InvalidClassException,
        IOException,
        StreamCorruptedException,
        OptionalDataException,
        TauModelException
    {
        try {
            TauModel tModLoad = TauModelLoader.load(modelName,
                                                    toolProps.getProperty("taup.model.path"));
            if (tModLoad != null) {
                tMod = tModLoad;
                tModDepth = tMod;
                      this.modelName = tMod.sMod.vMod.getModelName();
            }

        } catch (ClassNotFoundException e) {
            Alert.error("Caught ClassNotFoundException", e.getMessage()+
                        "\nThere must be something wrong with your installation of TauP.\n"+
                        "Exiting.");
            return;
        } catch (InvalidClassException e) {
            Alert.error("Model file "+
                        modelName+ " is not compatible with the current version.",
                        "Recreate using taup_create.");
        }
    }


    /** preforms intialization of the tool. Properties are queried for
     *  the the default model to load, source depth to use, phases to use,
     *  etc. Note that because of the IO inherent in these operations, this
     *  method is not appropriate for Applets. Applets should load
     *  TauModels themselves and use the setTauModel(TauModel) method.
     */
    public void init() throws IOException {

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
        if (dos==null) return;
        Writer s = new BufferedWriter(new OutputStreamWriter(dos));
        printResult(s);
        s.flush();
    }

    public static Arrival[] run_time(String[] args)
        throws FileNotFoundException,
        IOException,
        StreamCorruptedException,
        ClassNotFoundException,
        OptionalDataException
    {
        Arrival[] arrivals=null;
        try {
            long prevTime = 0;
            long currTime;

            prevTime = System.currentTimeMillis();
            MatTauP_Time tauPTime = new MatTauP_Time();
            String[] noComprendoArgs = tauPTime.parseCmdLineArgs(args);
            if (noComprendoArgs.length > 0) {
                for (int i=0;i<noComprendoArgs.length;i++) {
                    if (noComprendoArgs[i].equals("-help") ||
                        noComprendoArgs[i].equals("-version")) {
                        return arrivals;
                    }
                }
                String outStringA =
                    "I don't understand the following arguments, continuing:";
                String outStringB = "";
                for (int i=0;i<noComprendoArgs.length;i++) {
                    outStringB += noComprendoArgs[i]+" ";
                }
                Alert.warning(outStringA, outStringB);
                noComprendoArgs = null;
            }
            currTime = System.currentTimeMillis();

            prevTime = System.currentTimeMillis();
            tauPTime.init();
            currTime = System.currentTimeMillis();
            if (tauPTime.DEBUG) {
                Alert.info("taup model read time="+(currTime-prevTime));
            }

            tauPTime.start();
            arrivals=tauPTime.getArrivals();
            tauPTime.destroy();
            return arrivals;

        } catch (TauModelException e) {
            Alert.error("Caught TauModelException", e.getMessage());
            e.printStackTrace();
        } catch (TauPException e) {
            Alert.error("Caught TauPException", e.getMessage());
            e.printStackTrace();
        }
        return arrivals;
    }
}
