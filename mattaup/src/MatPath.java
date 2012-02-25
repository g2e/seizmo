/**
 * MATLAB TauP interface
 * Writen by Qin Li, qinli@u.washington.edu
 * University of Washington
 * Nov, 2002
 */

package edu.sc.seis.TauP;

public class MatPath  
{
    public int n_points;
    public double[] depth;
    public double[] time;
    public double[] dist;
    public double[] lat;
    public double[] lon;
    
    public MatPath()
    {
        this.n_points=0;
        this.time=null;
        this.dist=null;
        this.depth=null;
        this.lat=null;
        this.lon=null;
    }

    public MatPath(int n)
    {
        this.n_points=n;
        this.time=new double[n];
        this.dist=new double[n];
        this.depth=new double[n];
        this.lat=new double[n];
        this.lon=new double[n];
    }

    public void init_path(int n)
    {
        n_points=n;
        time=new double[n];
        dist=new double[n];
        depth=new double[n];
        lat=new double[n];
        lon=new double[n];
    }

    public int getNumPoints()
    {
        return n_points;
    }

    public double[] getDepths()
    {
        return depth;
    }

    public double[] getDistances()
    {
        return dist;
    }

    public double[] getTimes()
    {
        return time;
    }

    public double[] getLats()
    {
        return lat;
    }

    public double[] getLons()
    {
        return lon;
    }

    public double getDepth(int i)
    {
        return depth[i];
    }

    public double getDistance(int i)
    {
        return dist[i];
    }

    public double getTime(int i)
    {
        return time[i];
    }

    public double getLat(int i)
    {
        return lat[i];
    }

    public double getLon(int i)
    {
        return lon[i];
    }
}