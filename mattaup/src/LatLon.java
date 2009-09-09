/**
 * MATLAB TauP interface
 * Writen by Qin Li, qinli@u.washington.edu
 * University of Washington
 * Nov, 2002
 */

package edu.sc.seis.TauP;

public class LatLon
{
    public double lat;
    public double lon;

    public LatLon()
    {
        lat=0;
        lon=0;
    }
    public LatLon(double lat, double lon)
    {
        this.lat=lat;
        this.lon=lon;
    }
}
