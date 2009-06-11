/**
 * MATLAB TauP interface
 * Writen by Qin Li, qinli@u.washington.edu
 * University of Washington
 * Nov, 2002
 */

package edu.sc.seis.TauP;
import edu.sc.seis.TauP.*;

public class MatArrival extends Arrival 
{
    protected MatPath matPath;

    public MatArrival(SeismicPhase phase, double time, double dist, double rayParam, 
		   int rayParamIndex, String name, String puristName, double sourceDepth) 
    {
        super(phase,time,dist,rayParam,rayParamIndex,name,puristName,sourceDepth);
        this.matPath=new MatPath();
    }

    public MatArrival(Arrival arrival)
    {
        super(arrival.phase,arrival.time,arrival.dist,arrival.rayParam,
            arrival.rayParamIndex,arrival.name,arrival.puristName,arrival.sourceDepth);
        this.path=arrival.path;
        this.pierce=arrival.pierce;
        this.matPath=new MatPath();
    }

    public MatPath getMatPath()
    {
        return matPath;
    }
}
