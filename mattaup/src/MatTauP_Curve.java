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
    public MatCurve matCurves[];

    protected MatTauP_Curve() {
        super();
    }

    public MatTauP_Curve(TauModel tMod) throws TauModelException {
        super(tMod);
    }

    public MatTauP_Curve(String modelName) throws TauModelException {
        super(modelName);
    }

    public MatCurve getMatCurve(int i)
    {
        return matCurves[i];
    }

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
            }
        }
    }

    public void curvecalculate(double degrees) {
        /* no need to do any calculations, just check the phases since
         * they have already been
         * done within the seismic phase. So, this just overrides
         * TauP_Time.calculate. printResult handles everything else. */
        recalcPhases();
        SeismicPhase phase;
        double[] dist, time, rayParams;
        double arcDistance, timeReduced, rayParamsReduced;
        double maxTime = -1*Double.MAX_VALUE, minTime = Double.MAX_VALUE;

        matCurves=new MatCurve[phases.size()];
        
        for (int phaseNum=0;phaseNum<phases.size(); phaseNum++) {
            phase = phases.get(phaseNum);
		dist = phase.getDist();
		time = phase.getTime();
		rayParams = phase.getRayParams();

		matCurves[phaseNum]=new MatCurve(dist.length);
		matCurves[phaseNum].phaseName=phase.getName();
		matCurves[phaseNum].puristPhaseName=phase.getPuristName();
		matCurves[phaseNum].sourceDepth=depth;
		
		for (int i=0; i<dist.length; i++) {

		    /* trig trick to make sure the dist is 0 to PI. */
		    arcDistance = Math.acos(Math.cos(dist[i]));
		    if (reduceTime) {
			timeReduced = time[i]-dist[i]/reduceVel;
			rayParamsReduced = rayParams[i]-1/reduceVel;
		    } else {
			timeReduced = time[i];
			rayParamsReduced = rayParams[i];
		    }
		    matCurves[phaseNum].time[i]=timeReduced;
		    matCurves[phaseNum].distance[i]=180.0/Math.PI*dist[i];
		    matCurves[phaseNum].mindistance[i]=180.0/Math.PI*arcDistance;
		    matCurves[phaseNum].rayParam[i]=rayParamsReduced;
		}
        }
    }
}
