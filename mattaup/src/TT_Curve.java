/**
 * MATLAB TauP interface
 * Writen by Qin Li, qinli@u.washington.edu
 * University of Washington
 * Nov, 2002
 */

package edu.sc.seis.TauP;

public class TT_Curve 
{
    public int n_points;
    public String phaseName;
    public double sourceDepth;
    public double time[];
    public double dist[];
    public double rayParam[];
    
    public TT_Curve()
    {
        n_points=0;
        phaseName=null;
        sourceDepth=0;
        time=null;
        dist=null;
        rayParam=null;
    }

    public TT_Curve(int n)
    {
        n_points=n;
        phaseName=null;
        sourceDepth=0;
        time=new double[n];
        dist=new double[n];
        rayParam=new double[n];
    }

    public void init_Curve(int n)
    {
        n_points=n;
        time=new double[n];
        dist=new double[n];
        rayParam=new double[n];
    }
}