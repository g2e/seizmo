/**
 * MATLAB TauP interface
 * Writen by Qin Li, qinli@u.washington.edu
 * University of Washington
 * Nov, 2002
 */

package edu.sc.seis.TauP;

public class MatPath  
{
    public int n_point;
    public double[] p;
    public double[] depth;
    public double[] time;
    public double[] dist;
    public double[] lat;
    public double[] lon;
    
    public MatPath()
    {
        this.n_point=0;
        this.p=null;
        this.time=null;
        this.dist=null;
        this.depth=null;
        this.lat=null;
        this.lon=null;
    }

    public MatPath(int n)
    {
        this.n_point=n;
        this.p=new double[n];
        this.time=new double[n];
        this.dist=new double[n];
        this.depth=new double[n];
        this.lat=new double[n];
        this.lon=new double[n];
    }

    public void init_path(int n)
    {
        n_point=n;
        p=new double[n];
        time=new double[n];
        dist=new double[n];
        depth=new double[n];
        lat=new double[n];
        lon=new double[n];
    }
}