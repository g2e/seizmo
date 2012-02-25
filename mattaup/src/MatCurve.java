/**
 * MATLAB TauP interface
 * Writen by Qin Li, qinli@u.washington.edu
 * University of Washington
 * Nov, 2002
 */

package edu.sc.seis.TauP;

public class MatCurve 
{
    public int n_points;
    public String phaseName;
    public String puristPhaseName;
    public double sourceDepth;
    public double time[];
    public double distance[];
    public double mindistance[];
    public double rayParam[];
    
    public MatCurve()
    {
        n_points=0;
        phaseName=null;
        puristPhaseName=null;
        sourceDepth=0;
        time=null;
        distance=null;
        mindistance=null;
        rayParam=null;
    }

    public MatCurve(int n)
    {
        n_points=n;
        phaseName=null;
        puristPhaseName=null;
        sourceDepth=0;
        time=new double[n];
        distance=new double[n];
        mindistance=new double[n];
        rayParam=new double[n];
    }

    public void init_curve(int n)
    {
        n_points=n;
        time=new double[n];
        distance=new double[n];
        mindistance=new double[n];
        rayParam=new double[n];
    }

    public int getNumPoints()
    {
        return n_points;
    }

    public String getPhaseName()
    {
        return phaseName;
    }

    public String getPuristPhaseName()
    {
        return puristPhaseName;
    }

    public double getSourceDepth()
    {
        return sourceDepth;
    }

    public double[] getDistances()
    {
        return distance;
    }

    public double[] getTimes()
    {
        return time;
    }

    public double[] getMinDistances()
    {
        return mindistance;
    }

    public double[] getRayParams()
    {
        return rayParam;
    }

    public double getDistance(int i)
    {
        return distance[i];
    }

    public double getTime(int i)
    {
        return time[i];
    }

    public double getMinDistance(int i)
    {
        return mindistance[i];
    }

    public double getRayParam(int i)
    {
        return rayParam[i];
    }
}