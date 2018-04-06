/*
                Interactive Moc Program

     Program to perform two dimensional analysis of supersonic flow
         through nozzles using the method of characteristics

                     Version 1.6g   - 18 Apr 14

                         Written by Tom Benson
                       NASA Glenn Research Center

>                              NOTICE
>This software is in the Public Domain.  It may be freely copied and used in
>non-commercial products, assuming proper credit to the author is given.  IT
>MAY NOT BE RESOLD.  If you want to use the software for commercial
>products, contact the author.
>No copyright is claimed in the United States under Title 17, U. S. Code.
>This software is provided "as is" without any warranty of any kind, either
>express, implied, or statutory, including, but not limited to, any warranty
>that the software will conform to specifications, any implied warranties of
>merchantability, fitness for a particular purpose, and freedom from
>infringement, and any warranty that the documentation will conform to the
>program, or any warranty that the software will be error free.
>In no event shall NASA be liable for any damages, including, but not
>limited to direct, indirect, special or consequential damages, arising out
>of, resulting from, or in any way connected with this software, whether or
>not based on warranty, contract, tort or otherwise, whether or not injury
>was sustained by persons or property or otherwise, and whether or not loss
>was sustained from, or arose out of the results of, or use of, the software
>or services provided hereunder.

      New Test : 
                 * add cone flow - Taylor-Maccoll
                 *   use integer prob = 0  single wedge
                 *               prob = 1  single cone
                 *               prob = 2   double wedge
                 *               prob = 3  ext compression single wedge inlt
                 *               prob = 4 prandtl-meyer expansion
                 *               prob = 5  2d - isen ext comp inlt
                 *                       design - off-design option
                 *               prob = 6  2d jet exhaust
                 *               prob = 7  MOC nozzle design - zones
                 *               prob = 8  MOC nozzle - points
                 *               prob = 9  Axi MOC nozzle - points
                 *   put in limiter for detached shock
                 *   put in 10 rays of output - new output panel for cone
                 *   add double wedge problem
                 *   include PM expansion
                 *   add zoom widget and translation to graphics
                   add double cone
                 *   add cowl - 
                 *      single wedge inlet (ext compression)
                       single cone inlet (ext compression)
                       mixed compression
                 * add isentropic compression section
                   add internal compression - moc
                 * add jet problem 
                 * add nozzle problem
                 *     get flow
                 *     add geometry
                 *     add pressure - temp - other flow variables
                 *     cleanup low PM angle integration
                 *     check/cleanup geometry to match a/a*
                 * add axisymmetric nozzle
                 *     add additional terms to equations
                 *     change x-y(r) to be non-dimensionalized by throat height
                 *     add geometry output to probs 7-8-9
                 * take out earth/mars gamma button
                   correct error in axi nozzle MOC 
                   put in logic for delhi -- limit for nozzle
                   cleanup
                   
                                                     TJB 18 Apr 14
*/
package Sfs;

import java.awt.*;
import java.lang.Math ;

public class Moc extends java.applet.Applet {

   final double convdr = 3.14515926/180.;
   double gama, mach0, ang1, ang2, sepval, machpm ;
   double angmn, angmx, machlo, machhi, ang2mn, ang2mx, mlaslo, mlashi ;
   double xlong, xr1, xr1mn, xr1mx, yr1 ;
   double cwly, cwymn, cwymx, cwlx, cwxmn, cwxmx ;
   double vpold,vpnew,vpoldr,vpnewr,vpoldt,vpnewt,thetold,thetnew,machnew ;
   double macl,nzht,nzhtlo,nzhthi,nzlg,nzlghi,nzlglo;
   double nprat,prathi,pratlo ;

   double angr, delmax, gamma, thmx ;
   double thetmx, delthet, dely, machend ;
   double delx, delxlo, delxhi ;

   int prob, nramps, nshocks, nslip, planet, numray, numd, desmod, ncycle ;
   int outzn,drawray ;
              //  plot data 
   static double fact,fac1,fac2 ;
   int xt,yt,sldloc ;

              //  flow parameters
   static double[] turning = new double[30] ;
   static double[] mach1   = new double[30] ;
   static double[] mach2   = new double[30] ;
   static double[] prat    = new double[30] ;
   static double[] trat    = new double[30] ;
   static double[] ptrat   = new double[30] ;
   static double[] rhorat  = new double[30] ;
   static double[] defl    = new double[30] ;
   static double[] shkang  = new double[30] ;
   static double[] pspo    = new double[30] ;
   static double[] tsto    = new double[30] ;
   static double[] ptpto   = new double[30] ;
   static double[] rsro    = new double[30] ;
   static double[] pm      = new double[30] ;
   static double[] mang    = new double[30] ;
   static double[] ppt     = new double[30] ;
   static double[] ttt     = new double[30] ;
   static double[] rrt     = new double[30] ;
   static double[] xflow   = new double[30] ;
   static double[] yflow   = new double[30] ;
   static boolean[] detach = new boolean[30] ;
              //  wedge geometry
   static double[] ang = new double[110] ;
   static int[] wfamily = new int[110] ;
   static double[] winter = new double[110] ;
   static double[] wxbgn = new double[110] ;
   static double[] wxnd = new double[110] ;
   static double[] wybgn = new double[110] ;
   static double[] wynd = new double[110] ;
   static double[] wslope = new double[110] ;
              // shock geometry
   static double[] sang = new double[35] ;
   static int[] sfamily = new int[35] ;
   static double[] sinter = new double[35] ;
   static double[] sxbgn = new double[35] ;
   static double[] sxnd = new double[35] ;
   static double[] sybgn = new double[35] ;
   static double[] synd = new double[35] ;
   static double[] sslope = new double[35] ;
               // conical ray flowfield
   static double[] vp = new double[200] ;
   static double[] vpr = new double[200] ;
   static double[] vpt = new double[200] ;
   static double[] rmach = new double[200] ;
   static double[] rthet = new double[200] ;
   static double[] rpt = new double[20] ;
   static double[] rdel = new double[20] ;
   static double[] rps = new double[20] ;
   static double[] rtemp = new double[20] ;
   static double[] rrho = new double[20] ;
         // expansion or compression fan  geometry
   static double[] expang = new double[35] ;
   static int[] efamily = new int[35] ;
   static double[] einter = new double[35];
   static double[] exbgn = new double[35] ;
   static double[] exnd = new double[35] ;
   static double[] eybgn = new double[35] ;
   static double[] eynd = new double[35] ;
   static double[] eslope = new double[35] ;
              // slip line geometry
   static double[] slinter = new double[10] ;
   static double[] slxbgn = new double[10] ;
   static double[] slxnd = new double[10] ;
   static double[] slybgn = new double[10] ;
   static double[] slynd = new double[10] ;
   static double[] slslope = new double[10] ;
                // isentropic relations
   double poverpt,tovertt,roverrt,arat,mu,nu ;
   double nuexit,pref,tref,rref ;
              // streamlines 
   static double[][] strx = new double[26][10] ;
   static double[][] stry = new double[26][10] ;
              // moc flow
   static double[][] mcmach  = new double[55][55] ;
   static double[][] mcturn  = new double[55][55] ;
   static double[][] mcdefl  = new double[55][55] ;
   static double[][] mcpm    = new double[55][55] ;
   static double[][] mcmang  = new double[55][55] ;
   static double[][] mcQ  = new double[55][55] ;
   static double[][] mcR  = new double[55][55] ;
   static double[][] mcx  = new double[55][55] ;
   static double[][] mcy  = new double[55][55] ;
   static double[][] mcal  = new double[55][55] ;
   static double[][] mcbe  = new double[55][55] ;
   static double[][] mcaxi1  = new double[55][55] ;
   static double[][] mcaxi2  = new double[55][55] ;

   static double[][] mcxul   = new double[55][55] ;
   static double[][] mcyul   = new double[55][55] ;
   static double[][] mcxur   = new double[55][55] ;
   static double[][] mcyur   = new double[55][55] ;
   static double[][] mcxll   = new double[55][55] ;
   static double[][] mcyll   = new double[55][55] ;
   static double[][] mcxlr   = new double[55][55] ;
   static double[][] mcylr   = new double[55][55] ;

   static double[][] mcprat  = new double[55][55] ;
   static double[][] mctrat  = new double[55][55] ;
   static double[][] mcrrat  = new double[55][55] ;
   static double[][] mcpp0   = new double[55][55] ;
   static double[][] mctt0   = new double[55][55] ;
   static double[][] mcrr0   = new double[55][55] ;
   static double[][] mcarat  = new double[55][55] ;

   int row, rowlo, rowhi, col, collo, colhi ;

   CardLayout layout,layin,layin2 ;

   Viewer view ;
   Num num ;
   Image offscreenImg ;
   Graphics offsGg ;

   public void init() {
     boolean which = false ;
 
     offscreenImg = createImage(this.size().width,
                      this.size().height) ;
     offsGg = offscreenImg.getGraphics() ;
 
     setLayout(new GridLayout(2,1,0,0)) ;

     setDefaults () ;

     view = new Viewer(this) ;

     num = new Num(this) ;

     add(view) ;
     add(num) ;

     comPute( ) ;
     view.start() ;
  }
 
   public void setDefaults() {

     ncycle = 1 ;
     numd = 1;
     drawray = 0;
     prob = 0 ;  // wedge flow
     mach0 = 2.0 ;
     ang1 = 10.0 ;
     ang2 = 0.0 ;
     gamma  = 1.4 ;
     gama = 1.4 ;
     numray = 10 ;
     desmod = 1 ;

     nramps = 1 ;
     nshocks = 1 ;
     xlong = 1000. ;

     planet = 0 ;
     angmn = 4.0 ;
     angmx = 35.0 ;
     ang2mn = 0.0 ;
     ang2mx = 35.0 ;
     machlo = 1.0 ;
     machhi = 4.0 ;
     macl = 1.3 ;
     mlaslo = 1.05 ;
     mlashi = 1.6 ;

     xr1 = 500. ;
     xr1mn = 10. ;
     xr1mx = 1000. ; 

     cwly = 150. ;
     cwymn = 5. ;
     cwymx = 500. ; 
     cwlx = 300. ;
     cwxmn = 100. ;
     cwxmx = 950. ; 

     ang[0]  = 0.0 ;
     prat[0] = 1.0 ;
     trat[0] = 1.0 ;
     ptrat[0] = 1.0 ;
     rhorat[0] = 1.0 ;
     pspo[0] = 1.0 ;
     tsto[0] = 1.0 ;
     rsro[0] = 1.0 ;
     ptpto[0]= 1.0 ;
     defl[0] = 0.0 ;
     turning[0] = 0.0 ;
     shkang[0] = 0.0 ;

     nzht = 50.0 ; nzhtlo = 1.0 ; nzhthi = 200. ;
     nzlg = .2 ; nzlglo = .05 ; nzlghi = 1. ;
     nprat = 1.0 ; pratlo = .2 ;  prathi = 2.0 ;

     delx = .01 ; delxlo = .0001; delxhi = .25 ;
     row =1; rowlo=1; rowhi= 100;
     col =1; collo=1; colhi= 100;

     xt = 80; yt = 20; sldloc = 140;
     fac1 = .625; fac2 = .015;
     fact = 1.9 ;
   }

   public void comPute() {

       if (prob == 0) {   // single wedge
           mach2[0] = mach1[0] = mach0;
           ang[1] = ang1 ;
           anlSing() ;
           getStream() ;
       }
       if (prob == 1) {   // single cone
           mach2[0] = mach1[0] = mach0;
           ang[1] = ang1 ;
           tayMac() ;
           getStream() ;
       }
       if (prob == 2) {   // double wedge
           mach2[0] = mach1[0] = mach0;
           ang[1] = ang1 ;
           ang[2] = ang2 ;
           anlDub() ;
       }
       if (prob == 3) {   // external compression single wedge inlet
           mach2[0] = mach1[0] = mach0;
           ang[1] = ang1 ;
           anlSingExt() ;
       }
       if (prob == 4) {   // centered prandtl-meyer expansion
           mach2[0] = mach1[0] = mach0;
           ang[1] = -ang1 ;
           anlExpan() ;
       }
       if (prob == 5) {   // 2d - isentropic ext comp inlet
           mach2[0] = mach1[0] = mach0;
           ang[1] = ang1 ;
           anlIsnExt() ;
       }
       if (prob == 6) {   // 2d - jet exhaust - over or under expanded
           mach2[0] = mach1[0] = mach0;
           anlJet() ;
       }
       if (prob == 7) {   // 2d - nozzle design - MOC zone
           machend = mach0;
           anlNoz() ;
       }
       if (prob == 8) {   // 2d - nozzle design - MOC points
           machend = mach0;
           anlNozII() ;
       }
       if (prob == 9) {   // Axi - nozzle design - MOC points
           machend = mach0;
           anlNozIII() ;
       }

       loadOut() ;
 
       view.repaint();
   }

   public void loadZero() {
       int i ;

       for(i=1; i<=27; ++i) {
            pspo[i] = 0.0 ;
            tsto[i] = 0.0 ;
            rsro[i] = 0.0 ;
            ptpto[i] = 0.0 ;
            prat[i] = 0.0 ;
            trat[i] = 0.0 ;
            rhorat[i] = 0.0 ;
            ptrat[i] = 0.0 ;
            turning[i] = 0.0 ;
            defl[i] = 0.0 ;
            shkang[i] = 0.0 ;
            mach2[i] = 0.0 ;
            detach[i] = false ;
       }
       for(i=1; i<=14; ++i) {
            getGeomOblq(300.0,300.0,1,i) ;
       }
   }

   public void  tayMac() {     //  Analysis for Single Cone.
       double thets,delc,dtht,gm1,angerr,signa ;
       double delts,machs,vps,vpsr,vpst ;
       double pspts,tstts,rsrts ; 
       int i,iter,itrk,itermx ;

       loadZero ();

       gm1 = gamma - 1.0 ;
       getThetmax (mach0,gamma) ;
       thets = .99 * thmx ;
       delc = convdr * ang1 ;
       iter = 0 ;
       itermx = 200 ;
       dtht = -.005 ;  
       nramps = 1;
       nshocks = 1;
       getGeom() ;

       getAnglim (mach0,gamma) ;
       if (ang1 > delmax ) {   // check for detached shock
           shkang[1] = 90.0 ;
           getNorm(mach1[1],gamma,0,1) ;
           getGeomNorm(wxbgn[1],wybgn[1],wxbgn[1],xlong,wfamily[1],1) ;
           detach[1] = true ;
           turning[1] = ang[1] ;  
           rmach[numray] = mach2[1];
           rpt[numray] = ptpto[1] ;
           rps[numray] = pspo[1] ;
           rthet[numray] = shkang[1] * convdr ;
           rtemp[numray] = tsto[1] ;
           rrho[numray] = rsro[1] ;
           rdel[numray] = delc ;
           for (i=1; i<=numray-1; ++i) {
              rmach[i] = rmach[numray];
              rpt[i] = rpt[numray] ;
              rps[i] = rps[numray] ;
              rtemp[i] = rtemp[numray] ;
              rrho[i] = rrho[numray] ;
              rdel[i] = rdel[numray] ;
              rthet[i] = rthet[numray] ;
           }
           return ;         
       }

       while(true) {   //  find the shock angle - for the given cone half angle
           iter = iter + 1;
           if (iter > itermx) break;
           getOblq(mach0,thets/convdr,gamma,0,1) ;

           delts = defl[1] * convdr ;
           machs = mach2[1] ;
           vps = 1.0 / Math.sqrt(2.0/(gm1 * machs*machs) + 1.0) ;
           vpsr = vps * Math.cos(thets - delts) ;
           vpst = -vps * Math.sin(thets - delts) ;

           thetold = thets ;
           vpold = vps ;
           vpoldr = vpsr ;
           vpoldt = vpst ;

           itrk = 0 ;
           numd = iter ;
           while(vpoldt < 0.0) {  // integrate until positive look at every .01 radian (.5 degree)
               itrk = itrk + 1;

               rkCalc (gamma,dtht) ;
               thetnew = thetold + dtht ;
               vpnew = Math.sqrt(vpnewr*vpnewr + vpnewt*vpnewt) ;
               machnew = Math.sqrt(2.0*vpnew*vpnew/(gm1*(1.0-vpnew*vpnew))) ;
               vpoldt = vpnewt ;
               vpoldr = vpnewr ;
               thetold = thetnew ;
               
               if (itrk < 190) {  // diagnostics
                 rthet[itrk] = thetold / convdr ;
                 vp[itrk] = vpnew ;
                 vpr[itrk] = vpoldr ;
                 vpt[itrk] = vpoldt ;
                 rmach[itrk] = machnew ;
               }
               numray = itrk ;
           }

           angerr = Math.abs(thetnew - delc) ;
           signa = angerr / Math.abs(angerr) ;

           if(angerr <= .005) { // we found it
               break ;
           }
           else { // change thets and try again
               thets = thets - .05 * signa * angerr ;
               if (Math.abs(thets) > thmx) thets = thmx ; // protection
           }
       }
       if(iter >= itermx) {  // detached - no Taylor Maccoll possible
           shkang[1] = 90.0 ;
           getNorm(mach1[1],gamma,0,1) ;
           getGeomNorm(wxbgn[1],wybgn[1],wxbgn[1],xlong,wfamily[1],1) ;
           detach[1] = true ;
           turning[1] = ang[1] ; 
           rmach[numray] = mach2[1];
           rpt[numray] = ptpto[1] ;
           rps[numray] = pspo[1] ;
           rtemp[numray] = tsto[1] ;
           rrho[numray] = rsro[1] ;
           rdel[numray] = delc ; 
           rthet[numray] = shkang[1] * convdr ;
           for (i=1; i<=numray-1; ++i) {
              rmach[i] = rmach[numray];
              rpt[i] = rpt[numray] ;
              rps[i] = rps[numray] ;
              rtemp[i] = rtemp[numray] ;
              rrho[i] = rrho[numray] ;
              rdel[i] = rdel[numray] ;
              rthet[i] = rthet[numray] ;
           }
           return ;         
       }

       numray = 10 ;
       dtht = -(thets - delc) / (numray-1) ;
       shkang[1] = thets / convdr ;
          // get conditions across the oblique shock
       getOblq(mach0,shkang[1],gamma,0,1) ;
       turning[1] = turning[0] + defl[1] ;
       sang[1] = turning[0] + shkang[1] ;
       getGeomOblq(wxbgn[1],wybgn[1],wfamily[1],1) ;
       delts = defl[1] * convdr ;
       machs = mach2[1] ;
       vps = 1.0 / Math.sqrt(2.0/(gm1 * machs*machs) + 1.0) ;
       vpsr = vps * Math.cos(thets - delts) ;
       vpst = -vps * Math.sin(thets - delts) ;
       rmach[numray] = machs;
       thetold = thets ;
       rthet[numray] = thetold ;
       rpt[numray] = ptpto[1] ;
       rps[numray] = pspo[1] ;
       rtemp[numray] = tsto[1] ;
       rrho[numray] = rsro[1] ;
       rdel[numray] = delts ;
          // get refernece conditions
       getIsen(machs,gamma) ;
       pspts = poverpt ;
       tstts = tovertt ;
       rsrts = roverrt ;
           // initialize for itegration
       vpold = vps ;
       vpoldr = vpsr ;
       vpoldt = vpst ;
       numd = numray - 1;
       while(numd > 0) {  // integrate Taylor Maccol equations to the surface
           rkCalc (gamma,dtht) ;
           thetnew = thetold + dtht ;
           vpnew = Math.sqrt(vpnewr*vpnewr + vpnewt*vpnewt) ;
           machnew = Math.sqrt(2.0*vpnew*vpnew/(gm1*(1.0-vpnew*vpnew))) ;
           vpoldt = vpnewt ;
           vpoldr = vpnewr ;
           thetold = thetnew ;
           rmach[numd] = machnew;
           getIsen(machnew,gamma) ;
           rthet[numd] = thetold ;
           rpt[numd] = rpt[numray] ;
           rps[numd] = rps[numray] * poverpt / pspts ;
           rtemp[numd] = rtemp[numray] * tovertt / tstts ;
           rrho[numd] = rrho[numray] * roverrt / rsrts ;
           rdel[numd] = rthet[numd] - Math.acos(vpoldr / vpnew) ;
           if (numd == 1) rdel[numd] = delc;
           numd = numd -1 ;
        }
        mach2[1] = rmach[1] ;
   
        return ;
   }

   public void rkCalc (double gam, double dthet) {
//           integrated runge -kutta function for Taylor-Maccoll
       double alp1,alp2,alp3 ;
       double vpr1,vpr2,vpr3 ;
       double vpt1,vpt2,vpt3 ;
       double tht1,tht2,tht3 ;

       alp1 = .25  ;
       alp2 = .33333333 ;
       alp3 = .5 ;

       vpt1 = vpoldt + alp1 * dthet * rkfunc(gam,thetold,vpoldr,vpoldt) ;
       vpr1 = vpoldr + alp1 * dthet * vpoldt ;
       tht1 = thetold + alp1 * dthet ;

       vpt2 = vpoldt + alp2 * dthet * rkfunc(gam,tht1,vpr1,vpt1) ;
       vpr2 = vpoldr + alp2 * dthet * vpt1 ;
       tht2 = thetold + alp2 * dthet ;

       vpt3 = vpoldt + alp3 * dthet * rkfunc(gam,tht2,vpr2,vpt2) ;
       vpr3 = vpoldr + alp3 * dthet * vpt2 ;
       tht3 = thetold + alp3 * dthet ;
 
       vpnewt = vpoldt + dthet * rkfunc(gam,tht3,vpr3,vpt3) ;
       vpnewr = vpoldr + dthet * vpt3 ;

       return ;
   }

   public double rkfunc (double gam, double theta, double vpra, double vpta) {
//         runge - kutta function for Taylor Maccoll
       double a,b,c,ans;

       a = 1.0 - vpra*vpra - vpta*vpta ;
       b = 2.0 * vpra + vpta / Math.tan(theta) ;
       c = .5 * (gam - 1.0) ;

       ans = (vpra * vpta*vpta -c*a*b) / (c*a-vpta*vpta) ;

       return ans ;
   }

   public void anlSing() {     //  Analysis for Single  Wedge.
       loadZero() ;

       mach1[1] = mach2[0] ;
       nshocks = 1 ;
       nramps = 1 ;
       angr  = ang[1] * convdr ;
       detach[1] = false ;
       getGeom () ;

       getAnglim(mach1[1],gamma) ;

       if (ang[1] > delmax) {
         getNorm(mach1[1],gamma,0,1) ;
         getGeomNorm(wxbgn[1],wybgn[1],wxbgn[1],xlong,wfamily[1],1) ;
         detach[1] = true ;
         turning[1] = ang[1] ;
       }
       else {
         getShkang(mach1[1],angr,gamma,nshocks) ;
         getOblq(mach1[1],shkang[1],gamma,0,1) ;
         turning[1] = turning[0] + defl[1] ;
         sang[1] = turning[0] + shkang[1] ;
         getGeomOblq(wxbgn[1],wybgn[1],wfamily[1],1) ;
       }
       return ;
   }

   public void anlDub() {     //  Analysis for Double Wedge.
       int i ;
       loadZero() ;

       nshocks = 0 ;
       nramps = 2 ;
       getGeom () ;

       for (i=1; i<=nramps; i++) {  // ramp shocks
          mach1[i] = mach2[i-1] ;
          getAnglim(mach1[i],gamma) ;
          delmax = delmax - .01 ;
          detach[i] = false ;

          if (ang[i] > delmax) {   // shock detached
            getNorm(mach1[i],gamma,i-1,i) ;
            getGeomNorm(wxbgn[i],wybgn[i],wxbgn[1],25.,wfamily[i],i) ;
            detach[i] = true ;
            nshocks++ ;
            if (i == 1) return ;
          }
          else {
            angr = ang[i] * convdr ;
            if (angr > 0.0001) {         // oblique shock wave
              nshocks++ ;
              getShkang(mach1[i],angr,gamma,i) ;
              getOblq(mach1[i],shkang[nshocks],gamma,i-1,nshocks) ;
              turning[nshocks] = turning[nshocks-1] + defl[nshocks] ;
              sang[nshocks] = turning[nshocks-1] + shkang[nshocks] ;
              getGeomOblq(wxbgn[i],wybgn[i],wfamily[i],nshocks) ;
            }
          }
       }

       if(nshocks == 2) {   // intersecting shocks
           nslip = 1;                
           getIntersect(1,2,3,4) ;
       }

       return ;
     }

   public void anlExpan() {     //  Analysis for Prandtl-Meyer Expansion.
       loadZero() ;

       mach1[1] = mach2[0] ;
       nshocks = 1 ;
       nramps = 1 ;
       angr  = ang[1] * convdr ;
       detach[1] = false ;
       getGeom () ;

       getIsenExp(mach1[1],gamma,angr,0,1) ;
       turning[1] = turning[0] + defl[1] ;
       sang[1] = turning[0] + shkang[1] ;
       if (mach1[1] > 1.001) getGeomOblq(wxbgn[1],wybgn[1],wfamily[1],1) ;
       else getGeomNorm(wxbgn[1],wybgn[1],wxbgn[1],25.,wfamily[1],1) ;
       getGeomExpan(wxbgn[1],wybgn[1],wfamily[1],0,1) ;

       return ;
   }

   public void anlIsnExt() {     //  Analysis for 2d isentropic- ext comp inlet
       double xa,delm,bet,alp;
       double lray,lht,ry,rx;
       int i ;

       if (desmod == 1) {
          loadZero() ;
       }

       getIsenRamp(mach2[0],gamma,0,0) ;
       mach1[1] = mach2[0] ;
       nshocks = 1 ;
       nramps = 1 ;
       angr  = ang[1] * convdr ;
       detach[1] = false ;
       if (desmod == 1) {
          getGeom () ;
       }

// calculate oblique shock wave from the first wedge
       getAnglim(mach1[1],gamma) ;
       if (ang[1] > delmax) {
         getNorm(mach1[1],gamma,0,1) ;
         getGeomNorm(wxbgn[1],wybgn[1],wxbgn[1],xlong,wfamily[1],1) ;
         detach[1] = true ;
         turning[1] = ang[1] ;
         return ;
       }
       else {
         getShkang(mach1[1],angr,gamma,nshocks) ;
         getOblq(mach1[1],shkang[1],gamma,0,1) ;
         turning[1] = turning[0] + defl[1] ;
         sang[1] = turning[0] + shkang[1] ;
         getGeomOblq(wxbgn[1],wybgn[1],wfamily[1],1) ;
       }
       xflow[1] = wxbgn[1] ;
       yflow[1] = wybgn[1] ;

 // in design mode
       if (desmod == 1) {
 // find the cowl leading edge
          cwlx = cwly / Math.tan(sang[1] * convdr) ;
 // get the mach angle and PM function for flow downstream of shock
          getIsenRamp(mach2[1],gamma,0,1) ;
  // find the end of ramp 1
          xa = cwlx - cwly / Math.tan(convdr*mang[1] + angr) ;
          wxnd[1] = xa * Math.tan(convdr*mang[1] + angr) / (Math.tan(convdr*mang[1] + angr) -
             Math.tan(angr) ) ;
          wynd[1] = wxnd[1] * Math.tan(angr) ;
//  loop through isentropic compression zone
          delm = (mach2[1] - macl)/(double) (numray) ;
          for (i=2; i<=numray + 1; ++ i) {
             nramps = i ;
             mach1[i] = mach2[i-1] ;
             mach2[i] = mach1[i] - delm ;
             wxbgn[i] = wxnd[i-1] ;
             wybgn[i] = wynd[i-1] ;
             getIsenRamp(mach2[i],gamma,i-1,i) ;
             mang[i] = mu ;
             pm[i]   = nu ;
             defl[i] = pm[i-1] - pm[i] ;
             turning[i] = turning[i-1] + defl[i] ;
             xa = cwlx - (cwly - wynd[i-1])/Math.tan(convdr*(mang[i]+turning[i])) ;
             wxnd[i] = (wxnd[i-1]*Math.tan(convdr*turning[i]) - xa * Math.tan(convdr*(mang[i] + turning[i]))) /
              (Math.tan(convdr*turning[i]) - Math.tan(convdr*(mang[i] + turning[i])));
             wynd[i] = wynd[i-1] + (wxnd[i] - wxnd[i-1]) * Math.tan(convdr*turning[i]) ;
 // store results
             exbgn[i] = wxnd[i] ;
             eybgn[i] = wynd[i] ;
             exnd[i] = cwlx ;
             eynd[i] = cwly ;
             xflow[i] = wxbgn[i] ;
             yflow[i] = wybgn[i] ;
          }
       }
 // off design mode
       if (desmod == 0) {
 // get the mach angle and PM function for flow downstream of shock
          getIsenRamp(mach2[1],gamma,0,1) ;
          lray = Math.sqrt((cwly-wynd[1])*(cwly-wynd[1]) + 
                           (cwlx-wxnd[1])*(cwlx-wxnd[1])) ;
//  loop through isentropic compression zone
          for (i=2; i<=numray + 1; ++i) {
             nramps = i ;
             mach1[i] = mach2[i-1] ;
             if(mach1[i] > 1.0) {
                getIsenComp(mach1[i],gamma,convdr*defl[i],i-1,i) ;
 // store results
                exbgn[i] = wxnd[i] ;
                eybgn[i] = wynd[i] ;
                exnd[i] = exbgn[i] + lray*Math.cos(convdr*(mang[i]+turning[i])) ;
                eynd[i] = eybgn[i] + lray*Math.sin(convdr*(mang[i]+turning[i])) ;
             }
             if(mach1[i] <= 1.0) {
                exbgn[i] = wxnd[i] ;
                eybgn[i] = wynd[i] ;
                exnd[i] = exbgn[i] ;
                eynd[i] = eybgn[i] ;
             }
          }
       }
 //  add last ramp plus normal shock
       nramps = nramps + 1 ;
       wxbgn[nramps] = wxnd[nramps -1] ;
       wybgn[nramps] = wynd[nramps -1] ;
       wslope[nramps] = turning[nramps -1] * convdr ;
       turning[nramps] = turning[nramps-1] ;

       mach1[nramps] = mach2[nramps-1] ;
       if (mach1[nramps] > 1.0) {
          nshocks = 2 ;
          getNorm(mach1[nramps],gamma,nramps-1,nramps) ;
          lray = Math.sqrt((cwly - wybgn[nramps]) * (cwly - wybgn[nramps]) + 
                        (cwlx - wxbgn[nramps]) * (cwlx - wxbgn[nramps]) );
          bet = Math.atan((cwly-wybgn[nramps])/(cwlx-wxbgn[nramps]))/convdr ;
          alp = bet - turning[nramps] ;
          lht = lray * Math.sin(convdr * alp) ;
          rx = cwlx + lht * Math.sin(wslope[nramps]) ;
          ry = cwly - lht * Math.cos(wslope[nramps]) ;
          getGeomNorm(cwlx,cwly,rx,ry,-wfamily[1],nshocks) ;
          xflow[nramps] = rx ;
          yflow[nramps] = ry ;
       }

       wxnd[nramps] = xr1 ;
       wynd[nramps] = wybgn[nramps] + (wxnd[nramps]-wxbgn[nramps])*Math.tan(wslope[nramps]) ;

       return ;
   }

   public void anlSingExt() {     //  Analysis for Single Wedge External Comp Inlet
       double ly,ry,rx;

       loadZero() ;

       mach1[1] = mach2[0] ;
       nshocks = 1 ;
       nramps = 1 ;
       angr  = ang[1] * convdr ;
       detach[1] = false ;
       getGeom () ;

       getAnglim(mach1[1],gamma) ;

       if (ang[1] > delmax) {
         getNorm(mach1[1],gamma,0,1) ;
         getGeomNorm(wxbgn[1],wybgn[1],wxbgn[1],xlong,wfamily[1],1) ;
         detach[1] = true ;
         turning[1] = ang[1] ;
       }
       else {
         getShkang(mach1[1],angr,gamma,nshocks) ;
         getOblq(mach1[1],shkang[1],gamma,0,1) ;
         turning[1] = turning[0] + defl[1] ;
         sang[1] = turning[0] + shkang[1] ;
         getGeomOblq(wxbgn[1],wybgn[1],wfamily[1],1) ;

         nshocks = 2 ;
         mach1[2] = mach2[1] ;
         getNorm(mach1[2],gamma,1,2) ;
         ly = cwlx * Math.tan(convdr*ang[1]) ;
         rx = cwlx + (cwly-ly)* Math.sin(convdr*ang[1]) * Math.cos(convdr*ang[1]) ;
         ry = ly + (cwly-ly) * Math.sin(convdr*ang[1]) * Math.sin(convdr*ang[1]) ;
         getGeomNorm(cwlx,cwly,rx,ry,-wfamily[1],2) ;
         turning[2] = ang[1] ;
       }
       return ;
   }

   public void getDelmx() {     //  routine to calculate the maximum delx
       double xthrtx,xthrty,alfa,betar,betar90,xrefl,tanalf,tanbet,tanb90, x1max ;
       float fl1 ;
       int i1;

       xthrtx  = 0.0 ;
       xthrty  = .5 ;

       getIsen(machend,gamma) ;
       nuexit = nu ;
       alfa =  nuexit / numray ;

       tanalf = Math.tan(convdr*alfa) ;

       betar = mu ;
       tanbet = Math.tan(convdr*betar);
       xrefl = xthrty * tanbet ;
       
       betar90 = 90.0 - betar;
       tanb90 = Math.tan(convdr*betar90) ;

       x1max = - 2.0 * xrefl * tanb90 / (tanalf - tanb90) ;

       delxhi = x1max / (numray/2 -1) ;

       if (delx > delxhi) delx = delxhi ;
  
       fl1 = (float) delx ;
       num.inp.dwn.nzexp.inleft.f3.setText(String.valueOf(filter3(fl1))) ;
       i1 = (int) (((delx - delxlo)/(delxhi-delxlo))*1000.) ;
       num.inp.dwn.nzexp.inright.s3.setValue(i1) ;

       return ;
  }


   public void anlNozIII() {     //  Analysis for Axi Moc Nozzle design by points
       double dell1, dell2 ;
       double machold,machnew,delold, deriv ;
       int counter ;
       int i, j, k;

       loadZero() ;

       wxbgn[1] = -nzlg ;
       wybgn[1] = .5;
       wxnd[1]  = 0.0 ;
       wynd[1]  = .5 ;

       getDelmx() ;

       getIsen(machend,gamma) ;
       nuexit = nu ;
       thetmx = nuexit / 2.0 ;
 //      delthet = nuexit / numray ;
       delthet = Math.atan(Math.sqrt(arat)-1.0)/convdr / numray ;
       mcpm[0][0] = 0.0 ;
       mcmang[0][0] = 90.0 ;  
       mcmach[0][0] = 1.0 ;
       getIsen(mcmach[0][0], gamma) ;
       pref = 1.0 / poverpt ;
       tref = 1.0 / tovertt ;
       rref = 1.0 / roverrt ;
       mcpp0[0][0] = 1.0 ;
       mctt0[0][0] = 1.0 ;
       mcrr0[0][0] = 1.0 ;
       mcprat[0][0] = 1.0 ;
       mctrat[0][0] = 1.0 ;
       mcrrat[0][0] = 1.0 ;
       mcarat[0][0] = arat ;

// analysis by points 

// initialize R-K
     delold = delthet + 0.01 ;
     machold = machend + .005 ;
     machnew = machend ;
     counter = 0 ;
// R-K analysis
     while(Math.abs(machend-machold) > .001 && counter < 40) {

// 1-1  - solid boundary
       mcx[1][1] = 0.0 ;
       mcy[1][1] = wybgn[1] ;
       mcdefl[1][1] = delthet ;
       mcturn[1][1] = delthet ;
       mcpm[1][1]   = mcturn[1][1] ;
       getMOCVar(1,1,0,0) ;
// 1-2  - plane of symmetry
       mcy[1][2] = 0.0;
       mcx[1][2] = mcx[1][1] + 
          (mcy[1][1]-mcy[1][2]) / Math.tan(convdr*(mcmang[1][1]-mcturn[1][1])) ;
       dell1 = Math.sqrt(((mcx[1][2] - mcx[1][1]) * (mcx[1][2]-mcx[1][1])) 
                      +  ((mcy[1][2] - mcy[1][1]) * (mcy[1][2]-mcy[1][1]))) ;
       mcturn[1][2] = 0.0 ;
       mcdefl[1][2] = mcturn[1][2] - mcturn[1][1] ;
       mcaxi1[1][1] = (dell1 * Math.sin(convdr*mcmang[1][1]) * 
                 Math.sin(convdr*mcturn[1][1]) / mcy[1][1])/convdr ;
       mcpm[1][2]= mcpm[1][1] + mcturn[1][1] + mcaxi1[1][1] ;
       getMOCVar(1,2,1,1) ;

//  n-1 - expansion on the boundary
       for(i=2; i<=numray/2; ++i) {
          mcx[i][1] = mcx[i-1][1] + delx * wybgn[1] ;
          mcy[i][1] = mcy[i-1][1] + 
                (mcx[i][1]-mcx[i-1][1])*Math.tan(convdr*mcturn[i-1][1]);
          mcdefl[i][1] = delthet ;
          mcturn[i][1] = mcturn[i-1][1] + mcdefl[i][1] ;
          mcpm[i][1]   = mcturn[i][1] ;
          getMOCVar(i,1,i-1,1) ;
// internal points  
          for(k=2; k<=i; ++k) {
             mcal[i][k-1] = (mcmang[i][k-1] - mcturn[i][k-1]) ;
             mcbe[i-1][k] = (mcmang[i-1][k] + mcturn[i-1][k]) ;
             mcx[i][k] = (mcy[i][k-1] - mcy[i-1][k] + 
                mcx[i][k-1] * Math.tan(convdr*mcal[i][k-1]) + 
                mcx[i-1][k] * Math.tan(convdr*mcbe[i-1][k]) )/
                (Math.tan(convdr*mcal[i][k-1]) + Math.tan(convdr*mcbe[i-1][k]) );
             mcy[i][k] = mcy[i][k-1] - 
               (mcx[i][k] - mcx[i][k-1])*Math.tan(convdr*mcal[i][k-1]);
             dell1 = Math.sqrt(((mcx[i][k]-mcx[i][k-1]) * (mcx[i][k]-mcx[i][k-1])) 
                            +  ((mcy[i][k]-mcy[i][k-1]) * (mcy[i][k]-mcy[i][k-1]))) ;
             dell2 = Math.sqrt(((mcx[i][k]-mcx[i-1][k]) * (mcx[i][k]-mcx[i-1][k])) 
                            +  ((mcy[i][k]-mcy[i-1][k]) * (mcy[i][k]-mcy[i-1][k]))) ;
             mcaxi1[i][k] = (dell1 * Math.sin(convdr*mcmang[i][k-1]) * 
                       Math.sin(convdr*mcturn[i][k-1]) / mcy[i][k-1])/convdr ;
             mcaxi2[i][k] = 0.0 ;
             if(mcy[i-1][k] > 0.0) {
                 mcaxi2[i][k] = (dell2 * Math.sin(convdr*mcmang[i-1][k]) * 
                       Math.sin(convdr*mcturn[i-1][k]) / mcy[i-1][k])/convdr ;
             }  
             mcpm[i][k]   = .5*(mcpm[i][k-1] + mcpm[i-1][k])  
                          + .5*(mcturn[i][k-1] - mcturn[i-1][k])
                          + .5*(mcaxi1[i][k] + mcaxi2[i][k]) ;
             mcturn[i][k] = .5*(mcpm[i][k-1] - mcpm[i-1][k]) 
                          + .5*(mcturn[i][k-1] + mcturn[i-1][k])
                          + .5*(mcaxi1[i][k] - mcaxi2[i][k]) ;
             mcdefl[i][k] = mcturn[i][k] - mcturn[i-1][k] ;
             getMOCVar(i,k,i,k-1) ;
          }

//  plane of symmetry
          mcy[i][i+1] = 0.0;
          mcx[i][i+1] = mcx[i][i] + 
            (mcy[i][i]-mcy[i][i+1]) / Math.tan(convdr*(mcmang[i][i]-mcturn[i][i])) ;
          mcturn[i][i+1] = 0.0 ;
          mcdefl[i][i+1] = mcturn[i][i+1] - mcturn[i][i] ;
          dell1 = Math.sqrt(((mcx[i][i+1]-mcx[i][i]) * (mcx[i][i+1]-mcx[i][i])) 
                         +  ((mcy[i][i+1]-mcy[i][i]) * (mcy[i][i+1]-mcy[i][i]))) ;
          mcaxi1[i][i] = (dell1 * Math.sin(convdr*mcmang[i][i]) * 
                    Math.sin(convdr*mcturn[i][i]) / mcy[i][i])/convdr ;
          mcpm[i][i+1]= mcpm[i][i] + mcturn[i][i] + mcaxi1[i][i] ;
          getMOCVar(i,i+1,i,i) ;
          machnew = mcmach[i][i+1] ;
       }
       deriv = (machnew - machold)/(delthet - delold) ;
       machold = machnew ;
       delold = delthet ;
       delthet = delold + (machend - machold)/deriv ;
       counter = counter + 1 ;
    }  // end of R-K

 //  cancellation surface
       k = numray/2 + 1 ;
       mcdefl[k][1] = mcdefl[k-1][1] ;
       mcturn[k][1] = mcturn[k-1][1] ;
       mcpm[k][1] = mcpm[k-1][1] ;
       mcmach[k][1] = mcmach[k-1][1] ;
       getIsen(mcmach[k][1], gamma) ;
       mcmang[k][1] = mcmang[k-1][1];
       mcal[k][1] = mcal[k-1][1] ;
       mcbe[k][1] = mcbe[k-1][1] ;
       mcx[k][1] = mcx[k-1][1] ;
       mcy[k][1] = mcy[k-1][1] ; 
       mcpp0[k][1] = mcpp0[k-1][1] ;
       mctt0[k][1] = mctt0[k-1][1] ;
       mcrr0[k][1] = mcrr0[k-1][1] ;
       mcarat[k][1] = mcarat[k-1][1] ;
       mcprat[k][1] = mcprat[k-1][1] ;
       mctrat[k][1] = mctrat[k-1][1] ;
       mcrrat[k][1] = mcrrat[k-1][1] ;

       for(i=2; i<=k; ++i) {
          mcal[k][i-1] =  mcturn[k][i-1] ;
          mcbe[k-1][i] = (mcmang[k-1][i] + mcturn[k-1][i]) ;
          mcx[k][i] = (mcy[k-1][i] - mcy[k][i-1] + 
             mcx[k][i-1] * Math.tan(convdr*mcal[k][i-1]) - 
             mcx[k-1][i] * Math.tan(convdr*mcbe[k-1][i]) )/
             (Math.tan(convdr*mcal[k][i-1]) - Math.tan(convdr*mcbe[k-1][i]) );
          mcy[k][i] = mcy[k][i-1] + 
            (mcx[k][i] - mcx[k][i-1])*Math.tan(convdr*mcal[k][i-1]);
          mcdefl[k][i] = -delthet ;
          mcturn[k][i] = mcturn[k][i-1] + mcdefl[k][i] ;
          mcpm[k][i] = mcpm[k-1][i]  ;
          getMOCVar(k,i,k,i-1) ;
       }
//  wall geometry
       for(i=2; i<=numray/2; ++i) {
          wxbgn[i] = wxnd[i-1] ;
          wybgn[i] = wynd[i-1] ;
          wxnd[i]  = mcx[i][1] ;
          wynd[i]  = mcy[i][1] ;
       }
       for(i=numray/2+1; i<=numray; ++i) {
          wxbgn[i] = wxnd[i-1] ;
          wybgn[i] = wynd[i-1] ;
          j = i - numray/2 + 1 ;
          wxnd[i]  = mcx[k][j] ;
          wynd[i]  = mcy[k][j] ;
       }
       wxbgn[numray + 1] = wxnd[numray] ;
       wybgn[numray + 1] = wynd[numray] ;
       wxnd[i]  = wxbgn[numray+1] + (wxnd[numray] - wxbgn[numray]) ;
       wynd[i]  = wybgn[numray+1] ;

       return ;
   }

   public void anlNozII() {     //  Analysis for 2D Moc Nozzle design by points
       int i, j, k;

       loadZero() ;

       wxbgn[1] = -nzlg ;
       wybgn[1] = .5;
       wxnd[1]  = 0.0 ;
       wynd[1]  = .5 ;

       getDelmx() ;

       getIsen(machend,gamma) ;
       nuexit = nu ;
       thetmx = nuexit / 2.0 ;
       delthet = nuexit / numray ;
       mcpm[0][0] = 0.0 ;
       mcmang[0][0] = 90.0 ;  
       mcmach[0][0] = 1.0 ;
       getIsen(mcmach[0][0], gamma) ;
       pref = 1.0 / poverpt ;
       tref = 1.0 / tovertt ;
       rref = 1.0 / roverrt ;
       mcpp0[0][0] = 1.0 ;
       mctt0[0][0] = 1.0 ;
       mcrr0[0][0] = 1.0 ;
       mcprat[0][0] = 1.0 ;
       mctrat[0][0] = 1.0 ;
       mcrrat[0][0] = 1.0 ;
       mcarat[0][0] = arat ;

// analysis by points 
// 1-1  - solid boundary
       mcdefl[1][1] = delthet ;
       mcturn[1][1] = delthet ;
       mcpm[1][1]   = mcturn[1][1] ;
       mcQ[1][1] = mcpm[1][1] + mcturn[1][1] ;
       mcR[1][1] = mcpm[1][1] - mcturn[1][1] ;
       getMOCVar(1,1,0,0) ;
       mcx[1][1] = 0.0 ;
       mcy[1][1] = wybgn[1] ;
// 1-2  - plane of symmetry
       mcturn[1][2] = 0.0 ;
       mcdefl[1][2] = mcturn[1][2] - mcturn[1][1] ;
       mcQ[1][2] = mcQ[1][1] ;
       mcR[1][2] = mcQ[1][2] ;
       mcpm[1][2]= mcQ[1][2] ;
       getMOCVar(1,2,1,1) ;
       mcy[1][2] = 0.0;
       mcx[1][2] = mcx[1][1] + 
          (mcy[1][1]-mcy[1][2]) / Math.tan(convdr*(mcmang[1][1]-mcturn[1][1])) ;

//  n-1 - expansion on the boundary
       for(i=2; i<=numray/2; ++i) {
          mcdefl[i][1] = delthet ;
          mcturn[i][1] = mcturn[i-1][1] + mcdefl[i][1] ;
          mcpm[i][1]   = mcturn[i][1] ;
          mcQ[i][1] = mcpm[i][1] + mcturn[i][1] ;
          mcR[i][1] = mcpm[i][1] - mcturn[i][1] ;
          getMOCVar(i,1,i-1,1) ;
          mcx[i][1] = mcx[i-1][1] + delx * wybgn[1] ;
          mcy[i][1] = mcy[i-1][1] + 
                (mcx[i][1]-mcx[i-1][1])*Math.tan(convdr*mcturn[i-1][1]);
// internal points  
          for(k=2; k<=i; ++k) {
             mcQ[i][k] = mcQ[i][k-1] ;
             mcR[i][k] = mcR[i-1][k] ;
             mcpm[i][k]   = .5*(mcQ[i][k] + mcR[i][k]) ;
             mcturn[i][k] = .5*(mcQ[i][k] - mcR[i][k]) ;
             mcdefl[i][k] = mcturn[i][k] - mcturn[i-1][k] ;
             getMOCVar(i,k,i,k-1) ;
             mcal[i][k-1] = (mcmang[i][k-1] - mcturn[i][k-1]) ;
             mcbe[i-1][k] = (mcmang[i-1][k] + mcturn[i-1][k]) ;
             mcx[i][k] = (mcy[i][k-1] - mcy[i-1][k] + 
                mcx[i][k-1] * Math.tan(convdr*mcal[i][k-1]) + 
                mcx[i-1][k] * Math.tan(convdr*mcbe[i-1][k]) )/
                (Math.tan(convdr*mcal[i][k-1]) + Math.tan(convdr*mcbe[i-1][k]) );
             mcy[i][k] = mcy[i][k-1] - 
               (mcx[i][k] - mcx[i][k-1])*Math.tan(convdr*mcal[i][k-1]);
          }

//  plane of symmetry
          mcturn[i][i+1] = 0.0 ;
          mcdefl[i][i+1] = mcturn[i][i+1] - mcturn[i][i] ;
          mcQ[i][i+1] = mcQ[i][i] ;
          mcR[i][i+1] = mcQ[i][i+1] ;
          mcpm[i][i+1]= mcQ[i][i+1] ;
          getMOCVar(i,i+1,i,i) ;
          mcy[i][i+1] = 0.0;
          mcx[i][i+1] = mcx[i][i] + 
            (mcy[i][i]-mcy[i][i+1]) / Math.tan(convdr*(mcmang[i][i]-mcturn[i][i])) ;
       }
 //  cancellation surface
       k = numray/2 + 1 ;
       mcdefl[k][1] = mcdefl[k-1][1] ;
       mcturn[k][1] = mcturn[k-1][1] ;
       mcR[k][1] = mcR[k-1][1] ;
       mcpm[k][1] = mcpm[k-1][1] ;
       mcQ[k][1] = mcQ[k-1][1] ;
       mcmach[k][1] = mcmach[k-1][1] ;
       getIsen(mcmach[k][1], gamma) ;
       mcmang[k][1] = mcmang[k-1][1];
       mcal[k][1] = mcal[k-1][1] ;
       mcbe[k][1] = mcbe[k-1][1] ;
       mcx[k][1] = mcx[k-1][1] ;
       mcy[k][1] = mcy[k-1][1] ; 
       mcpp0[k][1] = mcpp0[k-1][1] ;
       mctt0[k][1] = mctt0[k-1][1] ;
       mcrr0[k][1] = mcrr0[k-1][1] ;
       mcarat[k][1] = mcarat[k-1][1] ;
       mcprat[k][1] = mcprat[k-1][1] ;
       mctrat[k][1] = mctrat[k-1][1] ;
       mcrrat[k][1] = mcrrat[k-1][1] ;

       for(i=2; i<=k; ++i) {
          mcdefl[k][i] = -delthet ;
          mcturn[k][i] = mcturn[k][i-1] + mcdefl[k][i] ;
          mcR[k][i] = mcR[k-1][i] ;
          mcpm[k][i] = mcR[k][i] + mcturn[k][i] ;
          mcQ[k][i] = mcpm[k][i] + mcturn[k][i] ;
          getMOCVar(k,i,k,i-1) ;
          mcal[k][i-1] =  mcturn[k][i-1] ;
          mcbe[k-1][i] = (mcmang[k-1][i] + mcturn[k-1][i]) ;
          mcx[k][i] = (mcy[k-1][i] - mcy[k][i-1] + 
             mcx[k][i-1] * Math.tan(convdr*mcal[k][i-1]) - 
             mcx[k-1][i] * Math.tan(convdr*mcbe[k-1][i]) )/
             (Math.tan(convdr*mcal[k][i-1]) - Math.tan(convdr*mcbe[k-1][i]) );
          mcy[k][i] = mcy[k][i-1] + 
            (mcx[k][i] - mcx[k][i-1])*Math.tan(convdr*mcal[k][i-1]);
       }
//  wall geometry
       for(i=2; i<=numray/2; ++i) {
          wxbgn[i] = wxnd[i-1] ;
          wybgn[i] = wynd[i-1] ;
          wxnd[i]  = mcx[i][1] ;
          wynd[i]  = mcy[i][1] ;
       }
       for(i=numray/2+1; i<=numray; ++i) {
          wxbgn[i] = wxnd[i-1] ;
          wybgn[i] = wynd[i-1] ;
          j = i - numray/2 + 1 ;
          wxnd[i]  = mcx[k][j] ;
          wynd[i]  = mcy[k][j] ;
       }
       wxbgn[numray + 1] = wxnd[numray] ;
       wybgn[numray + 1] = wynd[numray] ;
       wxnd[i]  = wxbgn[numray+1] + (wxnd[numray] - wxbgn[numray]) ;
       wynd[i]  = wybgn[numray+1] ;

       return ;
   }

   public void getMOCVar(int irow, int icol, int upr, int upc) {     
  // Routine to compute MOC variables
       getMachpm(mcpm[irow][icol],gamma) ;
       mcmach[irow][icol] = machpm ;
       getIsen(mcmach[irow][icol], gamma) ;
       mcmang[irow][icol] = mu ;
       mcpp0[irow][icol] = poverpt * pref ;
       mctt0[irow][icol] = tovertt * tref ;
       mcrr0[irow][icol] = roverrt * rref ;
       mcarat[irow][icol] = arat ;
       mcprat[irow][icol] = mcpp0[irow][icol] / mcpp0[upr][upc] ;
       mctrat[irow][icol] = mctt0[irow][icol] / mctt0[upr][upc] ;
       mcrrat[irow][icol] = mcrr0[irow][icol] / mcrr0[upr][upc] ;

       return ;
   }

   public void anlNoz() {     //  Analysis for 2d Moc Nozzle Design - Field Moethod
       double ang1,ang2 ;
       int i, j, k ;

       loadZero() ;

       wxbgn[1] = -nzlg ;
       wybgn[1] = .5;
       wxnd[1]  = 0.0 ;
       wynd[1]  = .5;

       getDelmx() ;

       getIsen(machend,gamma) ;
       nuexit = nu ;
       thetmx = nuexit / 2.0 ;
       delthet = nuexit / numray ;
       dely = .5 * delx * Math.tan(convdr * delthet) ;
       mcpm[0][0] = 0.0 ;
       mcmang[0][0] = 90.0 ;  
//   solve for all the flows     
// 1-1  - initial
       mcmach[1][1] = 1.0 ;
       mcdefl[1][1] = 0.0 ;
       mcturn[1][1] = 0.0 ;
       mcpm[1][1]   = 0.0 ;
       mcmang[1][1] = 90.0 ;
       getIsen(mcmach[1][1], gamma) ;
       pref = 1.0 / poverpt ;
       tref = 1.0 / tovertt ;
       rref = 1.0 / roverrt ;
       mcpp0[1][1] = 1.0 ;
       mctt0[1][1] = 1.0 ;
       mcrr0[1][1] = 1.0 ;
       mcprat[1][1] = 1.0 ;
       mctrat[1][1] = 1.0 ;
       mcrrat[1][1] = 1.0 ;
       mcarat[1][1] = arat ;
//  1-n - expansion
       for(i=2; i<=numray/2+1; ++i) {
          mcdefl[1][i] = delthet ;
          mcturn[1][i] = mcturn[1][i-1] + mcdefl[1][i] ;
          mcpm[1][i]   = mcpm[1][i-1] + mcdefl[1][i] ;
          getMOCVar(1,i,1,i-1) ;
       }
// k- columns
       for(k=2; k<=numray/2+1; ++k) {
          mcdefl[k][k] = -delthet ;
          mcturn[k][k] = 0.0 ;
          mcpm[k][k]   = mcpm[k-1][k] - mcdefl[k][k] ;
          getMOCVar(k,k,k-1,k) ;
//  k-n - expansion
          if(k <= numray/2) {
             for(i=k+1; i<=numray/2 + 1; ++i) {
                mcdefl[k][i] =  delthet ;
                mcturn[k][i] = mcturn[k][i-1] + mcdefl[k][i] ;
                mcpm[k][i]   = mcpm[k][i-1] + mcdefl[k][i] ;
                getMOCVar(k,i,k,i-1) ;
             }
          }
       }
//  solve for geometry
//  moc grid
       mcxul[1][1] = 0.0 ;
       mcyul[1][1] = wybgn[1] ;
       mcxll[1][1] = 0.0 ;
       mcyll[1][1] = 0.0 ;
       mcxur[1][1] = mcxul[1][1] ;
       mcyur[1][1] = mcyul[1][1] ;
       ang1 = convdr*(mcmang[1][2]- mcturn[1][2]) ;
       mcxlr[1][1] = (mcyul[1][1]- mcyll[1][1])/ Math.tan(ang1) ;
       mcylr[1][1] = mcyll[1][1] ;

       for (i=2; i<=numray/2; ++ i) {
         mcxul[1][i] = mcxur[1][i-1];
         mcyul[1][i] = mcyur[1][i-1];
         mcxll[1][i] = mcxlr[1][i-1];
         mcyll[1][i] = mcylr[1][i-1];
         mcxur[1][i] = mcxul[1][i] + delx ;
         mcyur[1][i] = mcyul[1][i] + delx * Math.tan(convdr*mcturn[1][i]) ;
         ang1 = convdr*(mcmang[1][i+1] - mcturn[1][i+1]) ;
         ang2 = convdr*(mcmang[2][i] + mcturn[2][i]) ;
         mcxlr[1][i] = (mcyur[1][i]-mcyll[1][i]+
                 mcxur[1][i]*Math.tan(ang1)+ mcxll[1][i]*Math.tan(ang2))/
                 (Math.tan(ang1)+ Math.tan(ang2)) ;
         mcylr[1][i] = mcyll[1][i] + (mcxlr[1][i]-mcxll[1][i])*Math.tan(ang2);
       }

       mcxul[2][2] = mcxll[1][2];
       mcyul[2][2] = mcyll[1][2];
       mcxll[2][2] = mcxll[1][2];
       mcyll[2][2] = mcyll[1][2];
       mcxur[2][2] = mcxlr[1][2] ;
       mcyur[2][2] = mcylr[1][2] ;
       mcylr[2][2] = 0.0 ;
       ang1 = convdr*(mcmang[2][3] - mcturn[2][3]) ;
       mcxlr[2][2] = mcxur[2][2] + (mcyur[2][2] - mcylr[2][2])/Math.tan(ang1) ;

       for (k=3; k<=numray/2; ++ k) {
          for (i=k; i<=numray/2; ++ i) {
            mcxul[k-1][i] = mcxll[k-2][i];
            mcyul[k-1][i] = mcyll[k-2][i];
            mcxll[k-1][i] = mcxlr[k-1][i-1];
            mcyll[k-1][i] = mcylr[k-1][i-1];
            mcxur[k-1][i] = mcxlr[k-2][i];
            mcyur[k-1][i] = mcylr[k-2][i] ;
            ang1 = convdr*(mcmang[k-1][i+1] - mcturn[k-1][i+1]) ;
            ang2 = convdr*(mcmang[k][i] + mcturn[k][i]) ;
            mcxlr[k-1][i] = (mcyur[k-1][i]-mcyll[k-1][i]+
                 mcxur[k-1][i]*Math.tan(ang1)+ mcxll[k-1][i]*Math.tan(ang2))/
                      (Math.tan(ang2)+ Math.tan(ang1)) ;
            mcylr[k-1][i] = mcyll[k-1][i] + (mcxlr[k-1][i]-mcxll[k-1][i])*Math.tan(ang2);
          }

          mcxul[k][k] = mcxll[k-1][k];
          mcyul[k][k] = mcyll[k-1][k];
          mcxll[k][k] = mcxll[k-1][k];
          mcyll[k][k] = mcyll[k-1][k];
          mcxur[k][k] = mcxlr[k-1][k] ;
          mcyur[k][k] = mcylr[k-1][k] ;
          mcylr[k][k] = 0.0 ;
          ang1 = convdr*(mcmang[k][k+1] - mcturn[k][k+1]) ;
          mcxlr[k][k] = mcxur[k][k] + (mcyur[k][k] - mcylr[k][k])/Math.tan(ang1) ;
       }
  // triangle at 1 - numray/2 + 1
       mcxul[1][numray/2 + 1] = mcxur[1][numray/2];
       mcyul[1][numray/2 + 1] = mcyur[1][numray/2];
       mcxll[1][numray/2 + 1] = mcxlr[1][numray/2];
       mcyll[1][numray/2 + 1] = mcylr[1][numray/2];
       ang1 = convdr*(mcmang[2][numray/2+1] + mcturn[2][numray/2 + 1]) ;
  //     ang1 = convdr*(mcmang[1][numray/2+1] + mcturn[1][numray/2 + 1]) ;
       ang2 = convdr*mcturn[1][numray/2+1] ;
       mcxur[1][numray/2 + 1] = (mcyll[1][numray/2 + 1] - mcyul[1][numray/2+1] + 
          mcxul[1][numray/2 +1]*Math.tan(ang2) - mcxll[1][numray/2 +1]*Math.tan(ang1)) /
          (Math.tan(ang2)-Math.tan(ang1)) ;
       mcyur[1][numray/2 + 1] = mcyul[1][numray/2 + 1] + 
          (mcxur[1][numray/2 + 1]-mcxul[1][numray/2 + 1])*Math.tan(ang2) ; 
       mcxlr[1][numray/2 + 1] = mcxur[1][numray/2 + 1] ;
       mcylr[1][numray/2 + 1] = mcyur[1][numray/2 + 1];

       for (k=2; k<=numray/2; ++k) {
         mcxul[k][numray/2 + 1] = mcxur[k][numray/2];
         mcyul[k][numray/2 + 1] = mcyur[k][numray/2];
         mcxll[k][numray/2 + 1] = mcxlr[k][numray/2];
         mcyll[k][numray/2 + 1] = mcylr[k][numray/2];
         mcxur[k][numray/2 + 1] = mcxlr[k-1][numray/2 + 1] ;
         mcyur[k][numray/2 + 1] = mcylr[k-1][numray/2 + 1];
         ang1 = convdr*(mcmang[k+1][numray/2 + 1] + mcturn[k+1][numray/2+1]) ;
  //       ang1 = convdr*(mcmang[k][numray/2 + 1] + mcturn[k][numray/2+1]) ;
         ang2 = convdr*mcturn[k][numray/2 + 1] ;
         mcxlr[k][numray/2 + 1] = (mcyll[k][numray/2 + 1] - mcyur[k][numray/2+1] + 
           mcxur[k][numray/2 +1]*Math.tan(ang2) - mcxll[k][numray/2 +1]*Math.tan(ang1)) /
           (Math.tan(ang2)-Math.tan(ang1)) ;
         mcylr[k][numray/2 + 1] =  mcyur[k][numray/2 + 1] + 
          (mcxlr[k][numray/2 + 1]-mcxur[k][numray/2 + 1])*Math.tan(ang2) ; 
       }

  // final zone at numray/2 + 1 - numray/2 + 1
       mcxul[numray/2 + 1][numray/2 + 1] = mcxlr[numray/2][numray/2 + 1];
       mcyul[numray/2 + 1][numray/2 + 1] = mcylr[numray/2][numray/2 + 1];
       mcxll[numray/2 + 1][numray/2 + 1] = mcxll[numray/2][numray/2 + 1];
       mcyll[numray/2 + 1][numray/2 + 1] = mcyll[numray/2 + 1][numray/2 + 1];
       mcxur[numray/2 + 1][numray/2 + 1] = 2.0 * mcxlr[numray/2][numray/2 + 1] - 
                            mcxlr[numray/2 - 1][numray/2 + 1]  ;
       mcyur[numray/2 + 1][numray/2 + 1] = mcylr[numray/2][numray/2 + 1] ; 
       mcxlr[numray/2 + 1][numray/2 + 1] = mcxur[numray/2 + 1][numray/2 + 1] ;
       mcylr[numray/2 + 1][numray/2 + 1] = mcyll[numray/2][numray/2 + 1];
//  wall geometry
       for(i=2; i<=numray/2; ++i) {
          wxbgn[i] = wxnd[i-1] ;
          wybgn[i] = wynd[i-1] ;
          wxnd[i]  = mcxur[1][i] ;
          wynd[i]  = mcyur[1][i] ;
       }
       for(i=numray/2+1; i<=numray; ++i) {
          wxbgn[i] = wxnd[i-1] ;
          wybgn[i] = wynd[i-1] ;
          j = i - numray/2 + 1 ;
          wxnd[i]  = mcxur[j][k] ;
          wynd[i]  = mcyur[j][k] ;
       }

       return ;
   }

   public void anlJet() {     //  Analysis for Jet exhaust - over or under-expanded
       double fac1, fac2 ;
       double xl1, xl2, xl3, yl1, yl2;

       loadZero() ;
 // geometry
       wxbgn[1] = 0.0 ;
       wybgn[1] = 0.0 ;
       wxnd[1]  = 0.0 ;
       wynd[1]  = -nzht ;

       wxbgn[2] = 0.0 ;
       wybgn[2] = 1.0 ;
       wxnd[2]  = 0.0 ;
       wynd[2]  = nzht*(1.0 + nzlg) ;

       wxbgn[3] = 0.0 ;
       wybgn[3] = 0.0 ;
       wxnd[3]  = -nzht*(1.0 + nzlg) ;
       wynd[3]  = 0.0 ;

       wxbgn[4] = 0.0 ;
       wybgn[4] = 1.0 ;
       wxnd[4]  = -nzht*(1.0 + nzlg) ;
       wynd[4]  = 1.0 ; ;

//  flow variables in zones 0 and 1
       getIsenRamp(mach2[0],gamma,0,0) ;
       mach1[1] = mach2[0] ;
       getIsenRamp(mach1[1],gamma,0,1) ;
       mach2[1] = mach1[1] ;
       getIsen(mach2[1],gamma) ;
       mang[1] = mu ;
       pm[1] = nu ;
       ppt[1] = poverpt ;
       ttt[1] = tovertt ;
       rrt[1] = roverrt ;
       ptrat[1] = 1.0 ;
       prat[1] = ppt[1] ;
       trat[1] = ttt[1] ;
       rhorat[1] = rrt[1] ;   
       pspo[1] = 1.0 ;
       tsto[1] = 1.0;
       rsro[1] = 1.0 ;
       ptpto[1] = 1.0 ;

 //  solve for flow in zone 2
       mach1[2] = mach2[1] ;
       fac1 = (gamma - 1.0) / gamma ;
       fac2 = (gamma - 1.0) * .5 ;
       mach2[2] = Math.sqrt((1.0/fac2)*(((1.0+fac2*mach1[2]*mach1[2])*
                  Math.exp(fac1*Math.log(nprat)))-1.0)) ;
       getIsen(mach2[2],gamma) ;
       mang[2] = mu ;
       pm[2] = nu ;
       defl[2] = pm[2] - pm[1] ;
       ppt[2] = poverpt ;
       ttt[2] = tovertt ;
       rrt[2] = roverrt ;
       ptrat[2] = 1.0 ;
       prat[2] = ppt[2]/ppt[1] ;
       trat[2] = ttt[2]/ttt[1] ;
       rhorat[2] = rrt[2]/rrt[1] ;   
       pspo[2] = pspo[1] * prat[2] ;
       tsto[2] = tsto[1] * trat[2] ;
       rsro[2] = rsro[1] * rhorat[2] ;
       ptpto[2] = ptpto[1] ;
// flow
       if (nprat > 1.0) {  //  underexpanded
// geometry parameters for streamlines
         yl1 = nzht * Math.tan(defl[2]*convdr) /
                    (Math.tan(mang[1]*convdr)-Math.tan(defl[2]*convdr)) ;
         xl1 = (nzht + yl1)/Math.tan(convdr* mang[1]) ;
         xl2 = (nzht + yl1)/Math.tan(convdr*(mang[2]-defl[2])) ;
 // zone 3
         mach1[3] = mach2[2] ;
         defl[3]  = 0.0 ;
         pm[3] = pm[2] + defl[2] ;
         getMachpm(pm[3],gamma) ;
         mach2[3] = machpm ;
         getIsenJet(mach2[3],gamma,2,3) ;
 // zone 4
         mach1[4] = mach2[3];
         defl[4] = - defl[2] ;
         pm[4] = pm[3] + defl[4] ;
         getMachpm(pm[4],gamma) ;
         mach2[4] = machpm ;
         getIsenJet(mach2[4],gamma,3,4) ;
 // zone 5
         mach1[5] = mach2[4];
         defl[5] = 0.0 ;
         pm[5] = pm[1] ;
         getMachpm(pm[5],gamma) ;
         mach2[5] = machpm ;
         getIsenJet(mach2[5],gamma,4,5) ;
 //  lower stream line
         strx[1][1] = 0.0 ;
         stry[1][1] = 0.0 ;
         strx[2][1] = xl1 ;
         stry[2][1] = - yl1 ;
         strx[3][1] = xl2 ;
         stry[3][1] = stry[2][1] ;
         strx[4][1] = xl2 + xl1 ;
         stry[4][1] = stry[1][1] ;
 //  middle stream line
         strx[1][2] = 0.0 ;
         stry[1][2] = .5 *(nzht)  ;
         strx[2][2] = .5 *(xl2) ;
         stry[2][2] = stry[1][2] ;
         strx[3][2] = strx[2][2] + xl2 ;
         stry[3][2] = stry[2][2] ;
         strx[4][2] = xl2 + xl1 ;
         stry[4][2] = stry[1][2] ;
 //  upper stream line
         strx[1][3] = 0.0 ;
         stry[1][3] = nzht ;
         strx[2][3] = xl1 ;
         stry[2][3] = nzht + yl1 ;
         strx[3][3] = xl2 ;
         stry[3][3] = stry[2][3] ;
         strx[4][3] = xl2 + xl1 ;
         stry[4][3] = stry[1][3] ;
// waves
         exbgn[1] = strx[1][1] ;
         eybgn[1] = stry[1][1] ;
         exnd[1]  = strx[2][3] ;
         eynd[1]  = stry[2][3] ;
         efamily[1] = 1 ;
         exbgn[2] = strx[1][3] ;
         eybgn[2] = stry[1][3] ;
         exnd[2]  = strx[2][1] ;
         eynd[2]  = stry[2][1] ;
         efamily[2] = 1 ;
         exbgn[3] = strx[1][1] ;
         eybgn[3] = stry[1][1] ;
         exnd[3]  = strx[3][3] ;
         eynd[3]  = stry[3][3] ;
         efamily[3] = 1 ;
         exbgn[4] = strx[1][3] ;
         eybgn[4] = stry[1][3] ;
         exnd[4]  = strx[3][1] ;
         eynd[4]  = stry[3][1] ;
         efamily[4] = 1;
         exbgn[5] = strx[2][1] ;
         eybgn[5] = stry[2][1] ;
         exnd[5]  = strx[4][3] ;
         eynd[5]  = stry[4][3] ;
         efamily[5] = 2 ;
         exbgn[6] = strx[2][3] ;
         eybgn[6] = stry[2][3] ;
         exnd[6]  = strx[4][1] ;
         eynd[6]  = stry[4][1] ;
         efamily[6] = 2 ;
         exbgn[7] = strx[3][1] ;
         eybgn[7] = stry[3][1] ;
         exnd[7]  = strx[4][3] ;
         eynd[7]  = stry[4][3] ;
         efamily[7] = 2 ;
         exbgn[8] = strx[3][3] ;
         eybgn[8] = stry[3][3] ;
         exnd[8]  = strx[4][1] ;
         eynd[8]  = stry[4][1] ;
         efamily[8] = 2 ;
       }
       if (nprat <= 1.0 ) { // over expanded
 // zone 3
         mach1[3] = mach2[2] ;
         defl[3]  = 0.0 ;
         pm[3] = pm[2] + defl[2] ;
         getMachpm(pm[3],gamma) ;
         mach2[3] = machpm ;
         getIsenJet(mach2[3],gamma,2,3) ;
 // zone 4
         mach1[4] = mach2[3] ;
         fac1 = (gamma - 1.0) / gamma ;
         fac2 = (gamma - 1.0) * .5 ;
         mach2[4] = Math.sqrt((1.0/fac2)*(((1.0+fac2*mach1[4]*mach1[4])*
                    Math.exp(fac1*Math.log(nprat*pspo[3])))-1.0)) ;
         getIsenJet(mach2[4],gamma,3,4) ;
         pm[4] = nu ;
         defl[4] = pm[4] - pm[3] ;
 // zone 5
         mach1[5] = mach2[4] ;
         defl[5]  = 0.0 ;
         pm[5] = pm[4] + defl[4] ;
         getMachpm(pm[5],gamma) ;
         mach2[5] = machpm ;
         getIsenJet(mach2[5],gamma,4,5) ;
 // zone 6
         mach1[6] = mach2[5];
         defl[6] = - defl[4] ;
         pm[6] = pm[5] + defl[6] ;
         getMachpm(pm[6],gamma) ;
         mach2[6] = machpm ;
         getIsenJet(mach2[6],gamma,5,6) ;
 // zone 7
         mach1[7] = mach2[6];
         defl[7] = 0.0 ;
         pm[7] = pm[3] ;
         getMachpm(pm[7],gamma) ;
         mach2[7] = machpm ;
         getIsenJet(mach2[7],gamma,6,7) ;
// geometry parameters for streamlines
         yl1 = nzht * Math.tan(defl[2]*convdr) /
                    (Math.tan(mang[1]*convdr)- Math.tan(defl[2]*convdr)) ;
         xl1 = (nzht + yl1)/Math.tan(convdr* mang[1]) ;
         yl2 = (nzht + 2.0*yl1)*Math.tan(convdr*defl[4]) /
                    (Math.tan(convdr*mang[3])-Math.tan(convdr*defl[4])) ;
         xl2 = ((nzht + 2.0*yl1) + yl2)/Math.tan(convdr*mang[3]) ;
         xl3 = ((nzht + 2.0*yl1) + yl2)/Math.tan(convdr*(mang[4]-defl[4])) ;
 //  lower stream line
         strx[1][1] = 0.0 ;
         stry[1][1] = 0.0 ;
         strx[2][1] = xl1 ;
         stry[2][1] = - yl1 ;
         strx[3][1] = xl1 + xl2 ;
         stry[3][1] = - yl1 - yl2  ;
         strx[4][1] = xl1 + xl3 ;
         stry[4][1] = stry[3][1] ;
         strx[5][1] = xl1 + xl2 + xl3 ;
         stry[5][1] = stry[2][1] ; 
 //  middle stream line
         strx[1][2] = 0.0 ;
         stry[1][2] = .5 *(nzht)  ;
         strx[2][2] = .5 *(xl1) ;
         stry[2][2] = stry[1][2] ;
         strx[3][2] = xl1 + xl2 ;
         stry[3][2] = stry[1][2] ;
         strx[4][2] = xl1 + xl3 ;
         stry[4][2] = stry[3][2] ;
         strx[5][2] = xl1 + xl2 + xl3 ;
         stry[5][2] = stry[2][2] ; 
 //  upper stream line
         strx[1][3] = 0.0 ;
         stry[1][3] = nzht ;
         strx[2][3] = xl1 ;
         stry[2][3] = nzht + yl1 ;
         strx[3][3] = xl1 + xl2 ;
         stry[3][3] = nzht + yl1 + yl2  ;
         strx[4][3] = xl1 + xl3 ;
         stry[4][3] = stry[3][3] ;
         strx[5][3] = xl1 + xl2 + xl3 ;
         stry[5][3] = stry[2][3] ; 
 // waves
         exbgn[1] = strx[1][1] ;
         eybgn[1] = stry[1][1] ;
         exnd[1]  = strx[2][3] ;
         eynd[1]  = stry[2][3] ;
         efamily[1] = 2 ;
         exbgn[2] = strx[1][3] ;
         eybgn[2] = stry[1][3] ;
         exnd[2]  = strx[2][1] ;
         eynd[2]  = stry[2][1] ;
         efamily[2] = 2 ;
         exbgn[3] = strx[2][1] ;
         eybgn[3] = stry[2][1] ;
         exnd[3]  = strx[3][3] ;
         eynd[3]  = stry[3][3] ;
         efamily[3] = 1 ;
         exbgn[4] = strx[2][1] ;
         eybgn[4] = stry[2][1] ;
         exnd[4]  = strx[4][3] ;
         eynd[4]  = stry[4][3] ;
         efamily[4] = 1 ;
         exbgn[5] = strx[2][3] ;
         eybgn[5] = stry[2][3] ;
         exnd[5]  = strx[3][1] ;
         eynd[5]  = stry[3][1] ;
         efamily[5] = 1 ;
         exbgn[6] = strx[2][3] ;
         eybgn[6] = stry[2][3] ;
         exnd[6]  = strx[4][1] ;
         eynd[6]  = stry[4][1] ;
         efamily[6] = 1 ;
         exbgn[7] = strx[3][1] ;
         eybgn[7] = stry[3][1] ;
         exnd[7]  = strx[5][3] ;
         eynd[7]  = stry[5][3] ;
         efamily[7] = 2 ;
         exbgn[8] = strx[4][1] ;
         eybgn[8] = stry[4][1] ;
         exnd[8]  = strx[5][3] ;
         eynd[8]  = stry[5][3] ;
         efamily[8] = 2 ;
         exbgn[9] = strx[3][3] ;
         eybgn[9] = stry[3][3] ;
         exnd[9]  = strx[5][1] ;
         eynd[9]  = stry[5][1] ;
         efamily[9] = 2 ;
         exbgn[10] = strx[4][3] ;
         eybgn[10] = stry[4][3] ;
         exnd[10]  = strx[5][1] ;
         eynd[10]  = stry[5][1] ;
         efamily[10] = 2 ;

       }
       return ;
   }

   public void getThetmax (double machin, double gam) {
       double a1, ab1, ac1, sints, msq, mfor ;

       msq = machin * machin ;
       mfor = msq * msq ;

       a1 = 2.0 * gam * mfor ;
       ab1 = 4.0 * msq - (gam + 1.0) * mfor ;
       ac1 = -(2.0 + (gam + 1.0) * msq) ;
       sints = (-ab1 + Math.sqrt(Math.pow(ab1,2.0)-4.0*a1*ac1))/(2.0*a1) ;
       thmx = Math.asin(Math.sqrt(sints)) ;

       return ;
   }

   public void getAnglim (double machin, double gam) {
       double a1, ab1, ac1, sints, msq, mfor, cotd ;
       double f6 = -11.5385 ;
       double f5 = 78.2051 ;
       double f4 = -190.7517 ;
       double f3 = 170.169 ;
       double f2 = 41.5865 ;
       double f1 = -79.197 ;
       double g6 = -0.1002 ;
       double g5 =  1.1018 ;
       double g4 = -3.0212 ;
       double g3 = -7.4450 ;
       double g2 = 52.7644 ;
       double g1 = -25.2683 ;

       msq = machin * machin ;
       mfor = msq * msq ;

       if(prob == 0 || prob >= 2) {  // wedge problem
          a1 = 2.0 * gam * mfor ;
          ab1 = 4.0 * msq - (gam + 1.0) * mfor ;
          ac1 = -(2.0 + (gam + 1.0) * msq) ;
          sints = (-ab1 + Math.sqrt(Math.pow(ab1,2.0)-4.0*a1*ac1))/(2.0*a1) ;
          getThetmax(machin,gam) ;

          cotd = Math.tan(thmx)*(((gam+1.0)*msq)/    
                (2.0*(msq * sints - 1.0))-1.0);
          delmax = (Math.atan(1.0/cotd))/convdr ;

          return;
       }
       if (prob == 1) { // cone problem
          if(mach0 < 2.0) {
             delmax = f6 * machin * machin * machin * machin *machin
                 + f5 * machin * machin * machin * machin
                 + f4 * machin * machin * machin
                 + f3 * machin * machin
                 + f2 * machin + f1 ;
          }
          else {
             delmax = g6 * machin * machin * machin * machin *machin
                 + g5 * machin * machin * machin * machin
                 + g4 * machin * machin * machin
                 + g3 * machin * machin
                 + g2 * machin + g1 ;
          }

           return ;
       }
   }

   public void getShkang (double machin, double delr, double gam, int index) {
        // Iterate to get Shock Angle given Wedge Angle.
      double cotd,mst,gp1,number ;    
      double theto,thetn,delo,deln,deriv ;

      gp1 = gam + 1.0 ;
      theto = Math.asin(1.0/machin) ; 
      delo = 0.0 ;
      thetn = theto + 3.0 * convdr ;
      while (Math.abs(delr - delo) > .0001) {
         mst = machin * Math.sin(thetn) ;
         cotd = Math.tan(thetn)*((gp1*machin*machin)/    
                (2.0*(mst * mst - 1.0))-1.0);
         deln = Math.atan(1.0/cotd) ;
         deriv = (deln-delo)/(thetn-theto) ;
         delo = deln ;
         theto = thetn ;
         thetn = theto + (delr-delo)/deriv ;
      }

      number = theto / convdr ;
      if (machin < 1.00001) number = 90.0 ;
      shkang[index] = number ;
      return;
   }

   public void getOblq (double machin, double shang, double gam,
                        int upstrm, int index) {
          // NACA 1135 - oblique shock relations.
       double mst, gm1, gp1, msq, m2sq, cotd ;

       mst = machin * Math.sin(shang*convdr) ;
       msq = machin * machin ;
       gm1 = gam - 1.0 ;
       gp1 = gam + 1.0 ;

       prat[index] = (2.0*gam*mst*mst - gm1)/gp1 ;
       rhorat[index] = (gp1*mst*mst)/(gm1*mst*mst + 2.0) ;
       trat[index] = (2.0*gam*mst*mst - gm1) * (gm1*mst*mst + 2.0)
                  /(mst*mst*Math.pow(gp1,2.0)) ; 
       ptrat[index] = (Math.pow(((gp1*mst*mst)/(gm1*mst*mst+2.0)),(gam/gm1)))
               * Math.pow((gp1/(2.0*gam*mst*mst - gm1)),(1.0/gm1)) ;
       m2sq = ((msq * mst * mst * Math.pow(gp1,2.0)) +
              (-4.0 * (mst*mst  - 1.0) * (gam*mst*mst + 1.0))) /
              ((2.0*gam*mst*mst - gm1) * (gm1*mst*mst + 2.0)) ;
       mach2[index] = Math.sqrt(m2sq) ;
       cotd = Math.tan(shang*convdr)*((gp1*msq)/    
                (2.0*(mst * mst - 1.0))-1.0);
       defl[index] = (Math.atan(1.0/cotd))/convdr ;
       mach1[index] = machin ;

       pspo[index] = pspo[upstrm]*prat[index] ;
       tsto[index] = tsto[upstrm]*trat[index] ;
       rsro[index] = rsro[upstrm]*rhorat[index] ;
       ptpto[index] = ptpto[upstrm]*ptrat[index] ;
  
       return;
   }

   public void getNorm (double machin, double gam, 
                         int upstrm, int index) {
          // NACA 1135 - normal shock relations.
       double gm1, gp1, msq, m2sq ;

       msq = machin * machin ;
       gm1 = gam - 1.0 ;
       gp1 = gam + 1.0 ;

       prat[index] = (2.0*gam*msq - gm1)/gp1 ;
       rhorat[index] = (gp1*msq)/(gm1*msq + 2.0) ;
       trat[index] = prat[index] / rhorat[index] ;
       ptrat[index] = (Math.pow(rhorat[index],(gam/gm1)))
               * (Math.pow((1.0/prat[index]),(1.0/gm1))) ;
       m2sq = msq / (prat[index] * rhorat[index]) ;
       mach2[index] = Math.sqrt(m2sq) ;
       defl[index] = 0.0 ;
       mach1[index] = machin ;
       shkang[index] = 90.0 ;
       sang[index] = 90.0 ;
//       detach[index] = true ;

       pspo[index] = pspo[upstrm]*prat[index] ;
       tsto[index] = tsto[upstrm]*trat[index] ;
       rsro[index] = rsro[upstrm]*rhorat[index] ;
       ptpto[index] = ptpto[upstrm]*ptrat[index] ;
 
       return;
   }

  public void getIsenExp (double machin, double gam, double delr,
                        int upstrm, int index) {
  // Centered Prandtl-Meyer Expansion
       double msq, msm1, gm1, gp1, fac1, fac2 ;
       double numax, nuo, nur, nun, m2sm1o, m2sm1n, deriv ;

       msq  = machin * machin ;
       msm1 = msq - 1.0;
       gm1 = gam - 1.0 ;
       gp1 = gam + 1.0 ;
       numax = 1.5707*(Math.sqrt(gp1/gm1) - 1.0) ;

//  limit calculations to M= 50
       numax = 124. * convdr ;

       nuo = Math.sqrt(gp1/gm1) * Math.atan( Math.sqrt(gm1*msm1/gp1)) 
             - Math.atan(Math.sqrt(msm1)) ;
       shkang[index] =  Math.asin(1.0/machin)/convdr  ;
       if (machin < 1.00001) shkang[index] = 90.0 ;

       nur = nuo - delr;
       if (nur < numax ) {

          m2sm1o = msm1;      // iterate for downstream mach 
          m2sm1n = msm1+.1 ;
          while (Math.abs(nur - nuo) > .0001) {
            nun = Math.sqrt(gp1/gm1) * Math.atan(Math.sqrt(gm1*m2sm1n/gp1)) 
                 - Math.atan(Math.sqrt(m2sm1n));
            deriv = (nun-nuo)/(m2sm1n-m2sm1o) ;
            nuo = nun ;
            m2sm1o = m2sm1n ;
            m2sm1n = m2sm1o + (nur-nuo)/deriv ;
          }

          mach2[index] = Math.sqrt(m2sm1o + 1.0) ;
          expang[index] = (Math.asin(1.0/mach2[index]) +delr)/convdr ;

          fac1 = 1.0 + .5*gm1*msq;
          fac2 = 1.0 + .5*gm1*(m2sm1o+1.0);
       }
       else {
          mach2[index] = 50.0 ;
          expang[index] = (Math.asin(1.0/mach2[index]) +delr)/convdr ;
          fac1 = 1.0 + .5*gm1*msq;
          fac2 = 1.0 + .5*gm1*(mach2[index]*mach2[index]);
       }

       prat[index] = Math.pow((fac1/fac2),(gam/gm1));
       rhorat[index] = Math.pow((fac1/fac2),(1.0/gm1));
       trat[index] = fac1/fac2 ;
       ptrat[index] = 1.0;                   /* Isentropic */
       defl[index] = delr/convdr ;
       mach1[index] = machin ;

       pspo[index] = pspo[upstrm]*prat[index] ;
       tsto[index] = tsto[upstrm]*trat[index] ;
       rsro[index] = rsro[upstrm]*rhorat[index] ;
       ptpto[index] = ptpto[upstrm]*ptrat[index] ;
  
       return;
  }

  public void getIntersect (int i, int j, int new1, int new2) {
       int k;
       double xint,yint;
       double delad,dela,delb,turno,turnn,delpo,delpn,deriv ;
                                  /* check for shock crossings */
       xint = (sinter[i] - sinter[j])/
           (sslope[j] - sslope[i]);
       if((xint > sxbgn[i] + .01) &&
          (xint < sxnd[i] -.01)) {
             yint = xint * sslope[i] + sinter[i];
             sxnd[i] = xint;
             synd[i] = yint;
             sxnd[j] = xint;
             synd[j] = yint;

             turno = Math.atan(wslope[j])/convdr ;
             delad = Math.abs(turning[i-1]-turno);
             getAnglim(mach1[i],gama);
             delmax = delmax - .001;
             if (delad > delmax)  {           /* detached - normal shock */
               getNorm (mach1[i],gama,i-1,new1) ;
               getGeomNorm (xint,yint,xint,xlong,sfamily[i],new1) ;
               detach[new1] = true ;
               nshocks++ ;
               nslip = 0 ;
               return;
             }

             dela = convdr*delad;
             getShkang(mach1[i],dela,gama,new1) ;
             getOblq(mach1[i],shkang[new1],gama,i-1,new1) ;
             nshocks++ ;
             delpo = (pspo[new1] - pspo[j])/pspo[j] ;
             delad = 0.0 ;
   
             for (k=1; k<=10; ++k) {                /* iterate for def */
                delad = delad + delpo / 2.0  * 
                    Math.sqrt(Math.pow(mach2[j],2.0)-1.0)/
                    (gama * Math.pow(mach2[j],2.0)) ;
                if (delad > 0.0) {
                   getShkang(mach2[j],delad,gama,new2) ;
                   getOblq (mach2[j],shkang[new2],gama,j,new2) ;
                }
                else {
                   getIsenExp (mach2[j],gama,delad, j, new2) ;
                }
                delb = dela - delad;
                getShkang(mach1[i],delb,gama,new1) ;
                getOblq (mach1[i],shkang[new1],gama,i-1,new1) ;
                delpo = (pspo[new1] - pspo[new2]);
                if (Math.abs(delpo) < .001) break;
                delpo = delpo/pspo[j] ;
             }
   
             nshocks++ ;
             getGeomSlip(xint, yint, turno, nslip) ;

             turning[new1] = turno - delad/convdr;
             sang[new1] = shkang[new1] ;
             getGeomOblq (xint,yint,sfamily[i],new1) ;
   
             sfamily[new2] = -sfamily[j];
             turning[new2] = turno - delad/convdr;
             sang[new2] = shkang[new2] - Math.abs(turning[j]);
             getGeomOblq (xint,yint,-sfamily[j],new2) ;
             sxnd[new2] = (sinter[new2] - winter[2])/
                            (wslope[2] - sslope[new2]);
             synd[new2] = sxnd[new2] * sslope[new2]
                                 + sinter[new2];
   
             if(prat[new2] < 1.0) {          /* expansion */
               getGeomExpan (xint,yint,-sfamily[j],j,new2) ;
               exnd[new2] = (einter[new2] - winter[2])/
                               (wslope[2] - eslope[new2]);
               eynd[new2] = exnd[new2] * eslope[new2] +
                                einter[new2] ;
             }
       }
    return;
  }
 
  public void getCrossShock (int i, int j, int new1, int new2) {

       int counter ;
       double xint,yint,y1;
       double delad,dela,turno,turnn,delpo,delpn,deriv ;
                                  /* check for shock crossings */
       xint = (sinter[i] - sinter[j])/
           (sslope[j] - sslope[i]);
       if((xint > sxbgn[i] + .01) &&
          (xint < sxnd[i] -.01)) {
             yint = xint * sslope[i] + sinter[i];
             sxnd[i] = xint;
             synd[i] = yint;
             sxnd[j] = xint;
             synd[j] = yint;
             sxbgn[new1] = xint ;
             sxbgn[new2] = xint ;
             if(nslip>1) {
               slxnd[nslip-1] = xint;
               slynd[nslip-1] = xint*slslope[nslip-1]+slinter[nslip-1];
             }
             turno = Math.atan(wslope[1])/convdr + Math.atan(wslope[2])/convdr ;
                                          /* initial guess for branch 1 */
             delad = Math.abs(turning[i]-turno);
             getAnglim(mach2[i],gama);
             delmax = delmax - .001;
             if (delad > delmax)  {           /* detached - normal shock */
               getNorm (mach2[i],gama,i,new1) ;
               yint = xint*wslope[1] + winter[1];
               y1   = xint*wslope[2] + winter[2];
               detach[new1] = true ;
               getGeomNorm (xint,yint,xint,y1,-sfamily[i],new1) ;
               return;
             }

             dela = convdr*delad;
             getShkang(mach2[i],dela,gama,new1) ;
             getOblq(mach2[i],shkang[new1],gama,i,new1) ;
                                          /* initial guess for branch 2 */
             delad = Math.abs(turning[j]-turno);
             getAnglim(mach2[j],gama) ;
             if (delad > delmax )  {
               getNorm (mach2[j],gama,j,new2) ;
               yint = xint*wslope[1] + winter[1];
               y1   = xint*wslope[2] + winter[2];
               getGeomNorm (xint,yint,xint,y1,-sfamily[j],new2) ;
               detach[new2] = true ;
               return;
             }
             dela = convdr*delad;
             getShkang(mach2[j],dela,gama,new2) ;
             getOblq (mach2[j],shkang[new2],gama,j,new2) ;
             delpo = pspo[new1] - pspo[new2];
                                      /* iterate for deflection */
             turnn = turno + .005 ;
             counter = 1;
             while (Math.abs(delpo) > .0001) {
                counter++ ;
                delad = Math.abs(turning[i]-turnn);
                dela = convdr*delad;
                getShkang (mach2[i],dela,gama,new1) ;
                getOblq (mach2[i],shkang[new1],gama,i,new1) ;
                delad = Math.abs(turning[j]-turnn);
                dela = convdr*delad;
                getShkang (mach2[j],dela,gama,new2) ;
                getOblq (mach2[j],shkang[new2],gama,j,new2) ;
                delpn = pspo[new1] - pspo[new2];
                deriv = (delpn-delpo)/(turnn-turno) ;
                delpo = delpn ;
                turno = turnn ;
                turnn = turno - delpo /deriv ;
                if (counter > 25) break;
              }

              getGeomSlip(xint, yint, turno, nslip) ;
   
              turning[new1] = turno;
              sang[new1] = shkang[new1] - Math.abs(turning[i]);
              getGeomOblq (xint,yint,-sfamily[i],new1) ;

              turning[new2] = turno;
              sang[new2] = shkang[new2] - Math.abs(turning[j]);
              getGeomOblq (xint,yint,-sfamily[j],new2) ;
       }
    return;
  }

  public void getHitWall (int j, int i, int new1) {

       double xint,yint,y1;
       double delad,dela,turno ;
                             
       xint = (sinter[i] - winter[j])/
              (wslope[j] - sslope[i]);
       if((xint > sxbgn[i] + .01) &&
          (xint < sxnd[i] -.01)) {
             sxnd[i] = xint;
             yint = xint * sslope[i] + sinter[i];
             synd[i] = yint;
             if(mach2[i] < 1.0) {
               detach[new1] = true;
               return;
             }

             getAnglim(mach2[i],gama);
             delmax = delmax - .01;

             turno = Math.atan(wslope[j])/convdr ;
             delad = Math.abs(turning[i]-turno);
             detach[new1] = false ;
             if (delad > delmax)  {  /* flow detached */
               getNorm (mach2[i],gama,i,new1) ;
               yint = xint*wslope[1] + winter[1];
               y1   = xint*wslope[2] + winter[2];
               getGeomNorm (xint,yint,xint,y1,-sfamily[i],new1) ;
               detach[new1] = true ;
               return;
             }

             dela = convdr*delad;
             getShkang(mach2[i],dela,gama,new1) ;
             getOblq(mach2[i],shkang[new1],gama,i,new1) ;
             turning[new1] = turno;
             sang[new1] = shkang[new1] - sfamily[i]*turning[i];
             getGeomOblq (xint,yint,-sfamily[i],new1) ;
       }
       return;
  }

  public void getIsen (double machin, double gam) { // isentropic relations 
     double mach1s,msm1,gm1,gp1,fac,fac1 ;   
          // poverpt and tovertt are ratiod to total conditions 
     mach1s = machin*machin ; 
     gm1 = gam - 1.0 ;
     gp1 = gam + 1.0 ;
     msm1 = mach1s - 1.0;
     fac = 1.0 + .5*gm1*mach1s;

     poverpt = Math.pow(1.0/fac,gam/gm1) ;                  /* EQ 44 */
     tovertt = 1.0 / fac ;                                   /* EQ 43 */
     roverrt = Math.pow(1.0/fac,1.0/gm1) ;                 /* EQ 45 */
     fac1 = gp1/(2.0*gm1) ;
     arat = machin * Math.pow(fac,-fac1) * Math.pow(gp1/2.0,fac1) ; /* EQ 80 */
     arat = 1.0/arat ;
     mu   = (Math.asin(1.0/machin))/convdr ;
     nu   = Math.sqrt(gp1/gm1)*Math.atan(Math.sqrt(gm1*msm1/gp1)) 
             - Math.atan(Math.sqrt(msm1)) ;
     nu   = nu / convdr;

     return;
  }

  public void getIsenJet (double machin, double gam, int upstream, int index) {
 // isentropic relations for Jet analysis
         getIsen(machin,gam) ;
         mang[index] = mu ;
         ppt[index] = poverpt ;
         ttt[index] = tovertt ;
         rrt[index] = roverrt ;
         ptrat[index] = 1.0 ;
         prat[index] = ppt[index]/ppt[upstream] ;
         trat[index] = ttt[index]/ttt[upstream] ;
         rhorat[index] = rrt[index]/rrt[upstream] ;   
         pspo[index] = pspo[upstream] * prat[index] ;
         tsto[index] = tsto[upstream] * trat[index] ;
         rsro[index] = rsro[upstream] * rhorat[index] ;
         ptpto[index] = ptpto[upstream] ;
         return ;
  }

  public void getIsenRamp (double machin, double gam, int upstream, int index) { // isentropic relations 
     double mach1s,msm1,gm1,gp1,fac,fac1 ;   
          // poverpt and tovertt are ratiod to total conditions 
     mach1s = machin*machin ; 
     gm1 = gam - 1.0 ;
     gp1 = gam + 1.0 ;
     msm1 = mach1s - 1.0;
     fac = 1.0 + .5*gm1*mach1s;

     ppt[index] = Math.pow(1.0/fac,gam/gm1) ;                  /* EQ 44 */
     ttt[index] = 1.0 / fac ;                                   /* EQ 43 */
     rrt[index] = Math.pow(1.0/fac,1.0/gm1) ;                 /* EQ 45 */
     fac1 = gp1/(2.0*gm1) ;
     arat = machin * Math.pow(fac,-fac1) * Math.pow(gp1/2.0,fac1) ; /* EQ 80 */
     arat = 1.0/arat ;
     mu = (Math.asin(1.0/machin))/convdr ;
     nu = Math.sqrt(gp1/gm1)*Math.atan(Math.sqrt(gm1*msm1/gp1)) 
             - Math.atan(Math.sqrt(msm1)) ;
     nu = nu / convdr;

     pm[index] = nu ;
     mang[index] = mu ;
     if(index >= 2) {
       ptrat[index] = 1.0 ;
       prat[index] = ppt[index]/ppt[upstream] ;
       trat[index] = ttt[index]/ttt[upstream] ;
       rhorat[index] = rrt[index]/rrt[upstream] ;
   
       pspo[index] = pspo[upstream] * prat[index] ;
       tsto[index] = tsto[upstream] * trat[index] ;
       rsro[index] = rsro[upstream] * rhorat[index] ;
       ptpto[index] = ptpto[upstream] ;
     }

     return;
  }

  public void getIsenComp (double machin, double gam, double delr,
                        int upstrm, int index) {
  // Prandtl-Meyer Compression
       double msq, msm1, gm1, gp1, fac1, fac2 ;
       double numax, nuo, nur, nun, m2sm1o, m2sm1n, deriv ;

       msq  = machin * machin ;
       msm1 = msq - 1.0;
       gm1 = gam - 1.0 ;
       gp1 = gam + 1.0 ;

       nuo = Math.sqrt(gp1/gm1) * Math.atan( Math.sqrt(gm1*msm1/gp1)) 
             - Math.atan(Math.sqrt(msm1)) ;

       nur = nuo - delr;

       m2sm1o = msm1;      // iterate for downstream mach 
       m2sm1n = msm1+.1 ;
       while (Math.abs(nur - nuo) > .000001) {
         nun = Math.sqrt(gp1/gm1) * Math.atan(Math.sqrt(gm1*m2sm1n/gp1)) 
                 - Math.atan(Math.sqrt(m2sm1n));
         deriv = (nun-nuo)/(m2sm1n-m2sm1o) ;
         nuo = nun ;
         m2sm1o = m2sm1n ;
         m2sm1n = m2sm1o + (nur-nuo)/deriv ;
       }

       mach2[index] = Math.sqrt(m2sm1o + 1.0) ;
       mang[index] =  Math.asin(1.0/mach2[index])/convdr  ;
       if (mach2[index] < 1.00001) mang[index] = 90.0 ;
       fac1 = 1.0 + .5*gm1*msq;
       fac2 = 1.0 + .5*gm1*(m2sm1o+1.0);

       prat[index] = Math.pow((fac1/fac2),(gam/gm1));
       rhorat[index] = Math.pow((fac1/fac2),(1.0/gm1));
       trat[index] = fac1/fac2 ;
       ptrat[index] = 1.0;                   /* Isentropic */
       mach1[index] = machin ;

       pspo[index] = pspo[upstrm]*prat[index] ;
       tsto[index] = tsto[upstrm]*trat[index] ;
       rsro[index] = rsro[upstrm]*rhorat[index] ;
       ptpto[index] = ptpto[upstrm]*ptrat[index] ;
  
       return;
  }

   public void getMachpm(double nuin, double gam)  {                 /* get the Mach number */
                                             /* given the Prandtl-meyer angle */
      double msm1o,msm1n,gp1,gm1;
      double nur,nuo,nun,deriv ;

      nur = nuin*convdr ;
      gm1 = gam - 1.0 ;
      gp1 = gam + 1.0 ;

      msm1o = 1.0;                                  /* iterate for mach */
      nuo = Math.sqrt(gp1/gm1)*Math.atan(Math.sqrt(gm1*msm1o/gp1)) 
             - Math.atan(Math.sqrt(msm1o));
      msm1n = msm1o+.01 ;
      while (Math.abs(nur - nuo) > .000001) {
         nun = Math.sqrt(gp1/gm1)*Math.atan(Math.sqrt(gm1*msm1n/gp1)) 
                - Math.atan(Math.sqrt(msm1n));
         deriv = (nun-nuo)/(msm1n-msm1o) ;
         nuo = nun ;
         msm1o = msm1n ;
         msm1n = msm1o + (nur-nuo)/deriv ;
       }
      machpm = Math.sqrt(msm1o + 1.0);
      return ;
   }

   public void getGeom () {
          // wedge geometry
       int i ;
       wfamily[1] = 1 ;
       wslope[1]  = Math.tan(convdr * ang[1]) ;
       wxbgn[1]   = 0.0 ;
       wybgn[1]   = 0.0 ;
       winter[1]  = wybgn[1] - wslope[1] * wxbgn[1] ;

       wfamily[2] = 1 ;
       wslope[2]  = Math.tan(convdr * (ang[1] + ang[2])) ;
       wxbgn[2]   = wxbgn[1] + xr1 ;
       wybgn[2]   = wxbgn[2] * wslope[1] + winter[1];
       winter[2]  = wybgn[2] - wslope[2] * wxbgn[2] ;
       wxnd[2]    = xlong ;
       wynd[2]    = wxnd[2] * wslope[2] + winter[2] ;

       wxnd[1]    = wxbgn[2] ;
       wynd[1]    = wybgn[2] ;
 
       return;
   }

   public void getGeomOblq (double xi, double yi, int fam, int index) {
          // oblique shock geometry
       sfamily[index] = fam ;
       sslope[index]  = fam*Math.tan(convdr * sang[index]) ;
       sxbgn[index]   = xi ;
       sybgn[index]   = yi ;
       sinter[index]  = yi - sslope[index] * xi ;
       sxnd[index]    = xlong ;
       synd[index]    = sxnd[index]*sslope[index] + sinter[index] ;
       return;
   }

   public void getGeomNorm (double xi, double yi, double xe, double ye,
                             int fam, int index) {
          // normal shock geometry
       sfamily[index] = fam ;
       sslope[index]  = 0.0 ;
       sxbgn[index]   = xi ;
       sybgn[index]   = yi ;
       sinter[index]  = yi ;
       sxnd[index]    = xe ;
       synd[index]    = ye ;
       return;
   }

  public void getGeomExpan(double xi, double yi, int fam, int prev, int index) {
          // expansion fan geometry
       efamily[index] = fam ;
       eslope[index]  = fam*Math.tan(convdr *(expang[index]
                          -Math.abs(turning[prev]))) ;
       exbgn[index]   = xi ;
       eybgn[index]   = yi ;
       einter[index]  = yi - eslope[index] * xi ;
       exnd[index]    = xlong ;
       eynd[index]    = exnd[index]*eslope[index] + einter[index] ;
       return;
  }

  public void getGeomSlip (double xi, double yi, double angle,
                             int index) {
          // slip line geometry
       slxbgn[index]   = xi ;
       slybgn[index]   = yi ;
       slslope[index]  = Math.tan(convdr*angle) ;
       slinter[index]  = yi  - slslope[index] * xi ;;
       slxnd[index]    = xlong ;
       slynd[index]    = slslope[index] * slxnd[index] + slinter[index] ;
       return;
  }

   public void loadOut() {

      if (prob == 0 || prob ==4) {  // wedge or centered expansion problem
         outzn = 1;

         num.out.wdg.o1.setText(String.valueOf(filter3(mach2[outzn]))) ;
         num.out.wdg.o3.setText(String.valueOf(filter3(shkang[outzn])));
         num.out.wdg.o2.setText(String.valueOf(filter3(prat[outzn]))) ;
         num.out.wdg.o4.setText(String.valueOf(filter3(ptrat[outzn]))) ;
         num.out.wdg.o6.setText(String.valueOf(filter3(trat[outzn]))) ;
         num.out.wdg.o7.setText(String.valueOf(filter3(rhorat[outzn])));
         return ;
      }
      if (prob ==1) {  // cone problem
         num.out.con.o3.setText(String.valueOf(filter3(shkang[1])));
         num.out.con.o1.setText(String.valueOf(filter3(mach2[1]))) ;
         num.out.con.o8.setText(String.valueOf(filter3(rmach[outzn]))) ;
         num.out.con.o2.setText(String.valueOf(filter3(rps[outzn]))) ;
         num.out.con.o4.setText(String.valueOf(filter3(rpt[outzn]))) ;
         num.out.con.o6.setText(String.valueOf(filter3(rtemp[outzn]))) ;
         num.out.con.o7.setText(String.valueOf(filter3(rrho[outzn])));
         num.out.con.o9.setText(String.valueOf(filter3(rdel[outzn]/convdr)));
         num.out.con.o10.setText(String.valueOf(filter3(rthet[outzn]/convdr)));
      }
      if (prob == 2 || prob == 3) {  // double wedge and ext comp inlet
          outzn = num.out.dub.zonch.getSelectedIndex() ;

          num.out.dub.o13.setText(String.valueOf(filter3(mach1[outzn]))) ;
          num.out.dub.o1.setText(String.valueOf(filter3(mach2[outzn]))) ;
          num.out.dub.o2.setText(String.valueOf(filter3(defl[outzn]))) ;
          num.out.dub.o3.setText(String.valueOf(filter3(shkang[outzn])));
          num.out.dub.o4.setText(String.valueOf(filter3(turning[outzn])));
          num.out.dub.o5.setText(String.valueOf(filter5(prat[outzn]))) ;
          num.out.dub.o6.setText(String.valueOf(filter5(pspo[outzn]))) ;
          num.out.dub.o7.setText(String.valueOf(filter5(ptrat[outzn]))) ;
          num.out.dub.o8.setText(String.valueOf(filter5(ptpto[outzn]))) ;
          num.out.dub.o9.setText(String.valueOf(filter5(trat[outzn]))) ;
          num.out.dub.o10.setText(String.valueOf(filter5(tsto[outzn]))) ;
          num.out.dub.o11.setText(String.valueOf(filter5(rhorat[outzn])));
          num.out.dub.o12.setText(String.valueOf(filter5(rsro[outzn])));

          num.out.dub.o18.setText(String.valueOf(filter3(turning[nshocks]))) ;
          num.out.dub.o19.setText(String.valueOf(filter3(ptpto[nshocks]))) ;
      }
      if (prob ==5) {  // 2d isentropic ext comp inlet
          outzn = num.out.isn.zonch.getSelectedIndex() ;

          num.out.isn.o1.setText(String.valueOf(filter3(mach2[outzn]))) ;
          num.out.isn.o2.setText(String.valueOf(filter3(defl[outzn]))) ;
          num.out.isn.o3.setText(String.valueOf(filter3(shkang[outzn])));
          num.out.isn.o4.setText(String.valueOf(filter3(turning[outzn])));
          num.out.isn.o5.setText(String.valueOf(filter5(prat[outzn]))) ;
          num.out.isn.o6.setText(String.valueOf(filter5(pspo[outzn]))) ;
          num.out.isn.o7.setText(String.valueOf(filter5(ptrat[outzn]))) ;
          num.out.isn.o8.setText(String.valueOf(filter5(ptpto[outzn]))) ;
          num.out.isn.o9.setText(String.valueOf(filter5(trat[outzn]))) ;
          num.out.isn.o10.setText(String.valueOf(filter5(tsto[outzn]))) ;
          num.out.isn.o11.setText(String.valueOf(filter5(rhorat[outzn])));
          num.out.isn.o12.setText(String.valueOf(filter5(rsro[outzn])));
          num.out.isn.o13.setText(String.valueOf(filter3(mach1[outzn]))) ;
          num.out.isn.o14.setText(String.valueOf(filter3(mang[outzn]))) ;
          num.out.isn.o15.setText(String.valueOf(filter3(pm[outzn]))) ;
          num.out.isn.o16.setText(String.valueOf(filter3(xflow[outzn]))) ;
          num.out.isn.o17.setText(String.valueOf(filter3(yflow[outzn]))) ;

          num.out.isn.o18.setText(String.valueOf(filter3(turning[nramps]))) ;
          num.out.isn.o19.setText(String.valueOf(filter3(ptpto[nramps]))) ;
      }
      if (prob == 6) {  // over or underexpanded jet
          outzn = num.out.jet.zonch.getSelectedIndex() ;

          num.out.jet.o1.setText(String.valueOf(filter3(mach2[outzn]))) ;
          num.out.jet.o2.setText(String.valueOf(filter3(defl[outzn]))) ;
          num.out.jet.o5.setText(String.valueOf(filter5(prat[outzn]))) ;
          num.out.jet.o6.setText(String.valueOf(filter5(pspo[outzn]))) ;
          num.out.jet.o7.setText(String.valueOf(filter5(ptrat[outzn]))) ;
          num.out.jet.o8.setText(String.valueOf(filter5(ptpto[outzn]))) ;
          num.out.jet.o9.setText(String.valueOf(filter5(trat[outzn]))) ;
          num.out.jet.o10.setText(String.valueOf(filter5(tsto[outzn]))) ;
          num.out.jet.o11.setText(String.valueOf(filter5(rhorat[outzn])));
          num.out.jet.o12.setText(String.valueOf(filter5(rsro[outzn])));
          num.out.jet.o15.setText(String.valueOf(filter3(pm[outzn]))) ;
          num.out.jet.o14.setText(String.valueOf(filter3(mang[outzn]))) ;
          num.out.jet.o13.setText(String.valueOf(filter3(mach1[outzn]))) ;
      }

      if (prob == 7) {  // moc nozzle

          num.out.noz.o1.setText(String.valueOf(filter3(mcmach[row][col]))) ;
          num.out.noz.o2.setText(String.valueOf(filter3(mcdefl[row][col]))) ;
          num.out.noz.o5.setText(String.valueOf(filter3(mcprat[row][col]))) ;
          num.out.noz.o6.setText(String.valueOf(filter3(mcpp0[row][col]))) ;
          num.out.noz.o9.setText(String.valueOf(filter3(mctrat[row][col]))) ;
          num.out.noz.o10.setText(String.valueOf(filter3(mctt0[row][col]))) ;
          num.out.noz.o11.setText(String.valueOf(filter3(mcrrat[row][col])));
          num.out.noz.o12.setText(String.valueOf(filter3(mcrr0[row][col])));
          num.out.noz.o16.setText(String.valueOf(filter3(mcarat[row][col]))) ;
          num.out.noz.o15.setText(String.valueOf(filter3(mcpm[row][col]))) ;
          num.out.noz.o13.setText(String.valueOf(filter3(mcturn[row][col]))) ;
          num.out.noz.o14.setText(String.valueOf(filter3(mcmang[row][col]))) ;

          num.out.noz.o3.setText(String.valueOf(filter3(mcxll[row][col]))) ;
          num.out.noz.o4.setText(String.valueOf(filter3(mcyll[row][col]))) ;
          num.out.noz.o7.setText(String.valueOf(filter3(mcxul[row][col]))) ;
          num.out.noz.o8.setText(String.valueOf(filter3(mcyul[row][col]))) ;
          num.out.noz.o17.setText(String.valueOf(filter3(mcxur[row][col]))) ;
          num.out.noz.o18.setText(String.valueOf(filter3(mcyur[row][col]))) ;
          num.out.noz.o19.setText(String.valueOf(filter3(mcxlr[row][col]))) ;
          num.out.noz.o20.setText(String.valueOf(filter3(mcylr[row][col]))) ;
      }

      if (prob == 8) {  // moc nozzle - points
          num.out.nozp.o1.setText(String.valueOf(filter3(mcmach[row][col]))) ;
          num.out.nozp.o2.setText(String.valueOf(filter3(mcdefl[row][col]))) ;
          num.out.nozp.o5.setText(String.valueOf(filter3(mcprat[row][col]))) ;
          num.out.nozp.o6.setText(String.valueOf(filter3(mcpp0[row][col]))) ;
          num.out.nozp.o9.setText(String.valueOf(filter3(mctrat[row][col]))) ;
          num.out.nozp.o10.setText(String.valueOf(filter3(mctt0[row][col]))) ;
          num.out.nozp.o11.setText(String.valueOf(filter3(mcrrat[row][col])));
          num.out.nozp.o12.setText(String.valueOf(filter3(mcrr0[row][col])));
          num.out.nozp.o16.setText(String.valueOf(filter3(mcarat[row][col]))) ;
          num.out.nozp.o15.setText(String.valueOf(filter3(mcpm[row][col]))) ;
          num.out.nozp.o13.setText(String.valueOf(filter3(mcturn[row][col]))) ;
          num.out.nozp.o14.setText(String.valueOf(filter3(mcmang[row][col]))) ;
          num.out.nozp.o7.setText(String.valueOf(filter3(mcQ[row][col]))) ;
          num.out.nozp.o8.setText(String.valueOf(filter3(mcR[row][col]))) ;

          num.out.nozp.o3.setText(String.valueOf(filter3(mcx[row][col]))) ;
          num.out.nozp.o4.setText(String.valueOf(filter3(mcy[row][col]))) ;
          num.out.nozp.o17.setText(String.valueOf(filter3(mcal[row][col]))) ;
          num.out.nozp.o18.setText(String.valueOf(filter3(mcbe[row][col]))) ;

      }

      if (prob == 9) {  // Axi moc nozzle - points
          num.out.nozp.o1.setText(String.valueOf(filter3(mcmach[row][col]))) ;
          num.out.nozp.o2.setText(String.valueOf(filter3(mcdefl[row][col]))) ;
          num.out.nozp.o5.setText(String.valueOf(filter3(mcprat[row][col]))) ;
          num.out.nozp.o6.setText(String.valueOf(filter3(mcpp0[row][col]))) ;
          num.out.nozp.o9.setText(String.valueOf(filter3(mctrat[row][col]))) ;
          num.out.nozp.o10.setText(String.valueOf(filter3(mctt0[row][col]))) ;
          num.out.nozp.o11.setText(String.valueOf(filter3(mcrrat[row][col])));
          num.out.nozp.o12.setText(String.valueOf(filter3(mcrr0[row][col])));
          num.out.nozp.o16.setText(String.valueOf(filter3(mcarat[row][col]))) ;
          num.out.nozp.o15.setText(String.valueOf(filter3(mcpm[row][col]))) ;
          num.out.nozp.o13.setText(String.valueOf(filter3(mcturn[row][col]))) ;
          num.out.nozp.o14.setText(String.valueOf(filter3(mcmang[row][col]))) ;
          num.out.nozp.o7.setText(String.valueOf(filter3(mcQ[row][col]))) ;
          num.out.nozp.o8.setText(String.valueOf(filter3(mcR[row][col]))) ;

          num.out.nozp.o3.setText(String.valueOf(filter3(mcx[row][col]))) ;
          num.out.nozp.o4.setText(String.valueOf(filter3(mcy[row][col]))) ;
          num.out.nozp.o17.setText(String.valueOf(filter3(mcaxi1[row][col]))) ;
          num.out.nozp.o18.setText(String.valueOf(filter3(mcaxi2[row][col]))) ;

      }

      return ;
   }

   public void getStream () {
       int i,k,p ;

      for(i=1; i<=4; ++ i) {
         strx[1][i] = wxbgn[1] - 100. ;
         stry[1][i] = 0.0 - (i-1)*100./3. ;
         strx[2][i] = wxbgn[1] ;
         stry[2][i] = stry[1][i] ;
      }
      if(detach[1]) {   // detached
         for(i=1; i<=4; ++ i) {
            stry[3][i] = stry[2][i] ;
            strx[3][i] = strx[2][i] ;
            strx[4][i] = wxnd[2] ;
            stry[4][i] = stry[3][i] - (strx[4][i]-strx[3][i]) * wslope[1]  ;
         }
      }
      else {  // attached
         if(prob == 0) { // wedge problem
             for(i=1; i<=4; ++ i) {
                stry[3][i] = stry[2][i] ;
                strx[3][i] = (stry[3][i] + sinter[1])/-sslope[1] ;
                strx[4][i] = wxnd[2] ;
                stry[4][i] = stry[3][i] - (strx[4][i]-strx[3][i]) * wslope[1]  ;
             }
          }
          else {  // cone problem
             for(i=1; i<=4; ++ i) {
                stry[3][i] = stry[2][i] ;
                strx[3][i] = (stry[3][i]) / - Math.tan(rthet[numray]) ;
             }
             for(k=4; k<=12; ++k) {
                p=k-4;
                for(i=1; i<=4; ++i) {
                    strx[k][i] = Math.abs(strx[k-1][i] * 
                        (Math.tan(rthet[numray-p])-Math.tan(rdel[numray-p]))/
                        (Math.tan(rthet[numray-p-1])-Math.tan(rdel[numray-p]))) ;
                    stry[k][i] = -strx[k][i]* Math.tan(rthet[numray-p-1]) ;
                }
             }
          }
       }

       return ;
   }

   public float filter3(double inumbr) {
     //  output only to .001
       float number ;
       int intermed ;

       intermed = (int) (inumbr * 1000.) ;
       number = (float) (intermed / 1000. );
       return number ;
  }

  public float filter5(double inumbr) {
     //  output only to .00001
       float number ;
       int intermed ;

       intermed = (int) (inumbr * 100000.) ;
       number = (float) (intermed / 100000. );
       return number ;
  }

   class Num extends Panel {
     Moc outerparent ;
     Inp inp ;
     Out out ;

     Num (Moc target) {                           
          outerparent = target ;
          setLayout(new GridLayout(1,2,10,10)) ;
 
          inp = new Inp(outerparent) ;  
          out = new Out(outerparent) ;
 
          add(inp) ;
          add(out) ;
     }

     class Inp extends Panel {
       Moc outerparent ;
       Up up ;
       Dwn dwn ;

       Inp (Moc target) {
                             
          outerparent = target ;
          setLayout(new GridLayout(2,1,2,2)) ;

          up = new Up(outerparent) ; 
          dwn = new Dwn(outerparent) ;
 
          add(up) ;
          add(dwn) ;
       }
 
       class Up extends Panel {
          Moc outerparent ;
          Con con ;
          In2 in2;

          Up (Moc target) {
                             
             outerparent = target ;
             setLayout(new GridLayout(2,1,2,2)) ;

             con = new Con(outerparent) ; 
             in2 = new In2(outerparent) ;
 
             add(con) ;
             add(in2) ;
          }

          class Con extends Panel {
             Inright inright ;
             Inleft inleft ;

             Con (Moc target) {                             
                outerparent = target ;
                setLayout(new GridLayout(1,2,2,2)) ;
   
                inleft = new Inleft(outerparent) ; 
                inright = new Inright(outerparent) ;
 
                add(inleft) ;
                add(inright) ;
             }
   
             class Inright extends Panel {
                Moc outerparent ;
                Label l1,l3,l4 ;
                Choice probch, plntch ;

                Inright (Moc target) {
                   outerparent = target ;
                   setLayout(new GridLayout(3,1,2,2)) ;

                   probch = new Choice() ;
                   probch.addItem("Single Wedge") ;
                   probch.addItem("Single Cone");
                   probch.addItem("Double Wedge");
                   probch.addItem("Centered Prandtl-Meyer Expansion");
                   probch.addItem("2D - 1 Ramp - Ext Comp Inlet");
                   probch.addItem("2D - Isentropic - Ext Comp Inlet");
                   probch.addItem("Jet - 2D Over/Under Expanded");
                   probch.addItem("Nozzle Design - MOC Field Method");
                   probch.addItem("Nozzle Design - MOC Points");
                   probch.addItem("Nozzle Design - MOC Axi Points");
//                  probch.addItem("Diagnostic");
                   probch.setBackground(Color.white) ;
                   probch.setForeground(Color.blue) ;
                   probch.select(0) ;

                   plntch = new Choice() ;
                   plntch.addItem("Earth") ;
                   plntch.addItem("Mars");
                   plntch.setBackground(Color.white) ;
                   plntch.setForeground(Color.blue) ;
                   plntch.select(0) ;

                   add(probch) ; 
                   add(new Label(" ", Label.RIGHT)) ; 
                   add(new Label(" ", Label.RIGHT)) ;  
 //                  add(plntch) ;  
                }

                public boolean handleEvent(Event evt) {
                   if(evt.target instanceof Choice) {
                      this.handleCho(evt) ;
                      return true ;
                   }
                   else return false ;
                }

                public void handleCho(Event evt) {     // problem
                   float fl2 ;
                   int probo ;

                   probo   = probch.getSelectedIndex() ;

                   if (probo == 0) {  // single wedge
                      layin2.show(in2, "first")  ;
                      layin.show(inp.dwn, "first")  ;
                      layout.show(out, "first")  ;
                      prob = 0 ;
                      desmod = 1 ;
                   }
                   if (probo == 1) {  // single cone
                      layin2.show(in2, "first")  ;
                      layin.show(inp.dwn, "first")  ;
                      layout.show(out, "second")  ;
                      prob = 1;
                      desmod = 1 ;
                   }
                   if (probo == 2) {   // double wedge
                      layin2.show(in2, "first")  ;
                      layin.show(inp.dwn, "second")  ;
                      layout.show(out, "third")  ;
                      prob = 2;
                      desmod = 1 ;
                   }
                   if (probo == 3) {   // centered Prandtl-Meyer expansion
                      layin2.show(in2, "first")  ;
                      layin.show(inp.dwn, "blank")  ;
                      layout.show(out, "first")  ;
                      prob = 4;
                      desmod = 1 ;
                   }
                   if (probo == 4) {   // 2d 1 rmp external comp inlet
                      layin2.show(in2, "first")  ;
                      layin.show(inp.dwn, "third")  ;
                      layout.show(out, "third")  ;
                      prob = 3;
                      desmod = 1 ;
                   }
                   if (probo == 5) {   //2d - isentropic ext comp inlet
                      layin2.show(in2, "first")  ;
                      layin.show(inp.dwn, "fourth")  ;
                      layout.show(out, "fourth")  ;
                      prob = 5;
                      desmod = 1 ;
                      out.isn.bt1.setBackground(Color.yellow) ;
                      out.isn.bt1.setForeground(Color.black) ;
                      out.isn.bt2.setBackground(Color.white) ;
                      out.isn.bt2.setForeground(Color.blue) ;
                   }
                   if (probo == 6) {   // 2d Jet exhaust - over or under expanded
                      layin2.show(in2, "second")  ;
                      layin.show(inp.dwn, "fifth")  ;
                      layout.show(out, "fifth")  ;
                      prob = 6;
                      desmod = 1 ;
                   }
                   if (probo == 7) {   // 2d Nozzle Design - Moc
                      layin2.show(in2, "second")  ;
                      layin.show(inp.dwn, "sixth")  ;
                      layout.show(out, "sixth")  ;
                      prob = 7;
                      desmod = 1 ;
                   }
                   if (probo == 8) {   // 2D Nozzle Design - Moc Point
                      layin2.show(in2, "second")  ;
                      layin.show(inp.dwn, "sixth")  ;
                      layout.show(out, "seventh")  ;
                      prob = 8;
                      desmod = 1 ;
                   }
                   if (probo == 9) {   // Axi Nozzle Design - Moc Point
                      layin2.show(in2, "second")  ;
                      layin.show(inp.dwn, "sixth")  ;
                      layout.show(out, "seventh")  ;
                      prob = 9;
                      desmod = 1 ;
                   }

  //                 planet   = plntch.getSelectedIndex() ;

                   if (planet == 0) gamma = 1.4;
                   if (planet == 1) gamma = 1.29;

                   fl2 = (float) gamma ;
                   inleft.f2.setText(String.valueOf(fl2)) ;

                   comPute() ;    
                }
             }  // end Right

             class Inleft extends Panel {
                Moc outerparent ;
                TextField f2;

                Inleft (Moc target) {
                   outerparent = target ;
                   setLayout(new GridLayout(3,2,2,2)) ;

                   f2 = new TextField(String.valueOf(gamma),5) ;

                   add(new Label(" ", Label.RIGHT)) ;  
                   add(new Label("Problem:", Label.RIGHT)) ;  

                   add(new Label(" ", Label.RIGHT)) ;  
                   add(new Label(" ", Label.RIGHT)) ;  

                   add(new Label("Gamma", Label.CENTER)) ;  
                   add(f2) ;  
                } 

                 public boolean handleEvent(Event evt) {
                    if(evt.id == Event.ACTION_EVENT) {
                        this.handleText(evt) ;
                        return true ;
                    }
                    else return false ;
                 }
 
                 public void handleText(Event evt) {
                      Double V1,V2,V3,V4,V5 ;
                      double v1,v2,v3,v4,v5 ;
                      float fl1,fl2 ;
                      int i1,i3,i4,i5 ;

         // Gamma - range from 1.0 to 2.0
                       V2 = Double.valueOf(f2.getText()) ;
                       v2 = V2.doubleValue() ;
                       if(v2 < 1.0) {
                          v2 = 1.02 ;
                          fl2 = (float) v2 ;
                          f2.setText(String.valueOf(fl2)) ;
                       }
                       if(v2 > 2.0) {
                          v2 = 2.0 ;
                          fl2 = (float) v2 ;
                          f2.setText(String.valueOf(fl2)) ; 
                       }       
                       gamma = v2 ;

                     comPute() ;
                 } // end handle Text
             }  // end Left
          }  // end Con

          class In2 extends Panel {
             Moc outerparent ;
             Inlt inlt ;
             Jet1 jet1 ;

             In2 (Moc target) {
                             
                outerparent = target ;
                layin2 = new CardLayout() ;
                setLayout(layin2) ;

                inlt = new Inlt(outerparent) ;
                jet1 = new Jet1(outerparent) ;

                add ("first", inlt) ;
                add ("second", jet1) ;
             }

             class Inlt extends Panel {
                Inright inright ;
                Inleft inleft ;

                Inlt (Moc target) {
                             
                   outerparent = target ;
                   setLayout(new GridLayout(1,2,2,2)) ;
   
                   inleft = new Inleft(outerparent) ; 
                   inright = new Inright(outerparent) ;
 
                   add(inleft) ;
                   add(inright) ;
                }
   
                class Inright extends Panel {
                   Moc outerparent ;
                   Scrollbar s1,s3;

                   Inright (Moc target) {
                      int i1,i3 ;
 
                      outerparent = target ;
                      setLayout(new GridLayout(3,1,2,2)) ;

                      i1 = (int) (((mach0 - machlo)/(machhi - machlo))*1000.) ;
                      i3 = (int) (((ang1 - angmn)/(angmx-angmn))*1000.) ;
 
                      s1 = new Scrollbar(Scrollbar.HORIZONTAL,i1,10,0,1000);
                      s3 = new Scrollbar(Scrollbar.HORIZONTAL,i3,10,0,1000);

                      add(new Label(" ", Label.CENTER)) ;  
                      add(s1) ;
                      add(s3) ;
                    }

                   public boolean handleEvent(Event evt) {
                      if(evt.id == Event.SCROLL_ABSOLUTE) {
                         this.handleBar(evt) ;
                         return true ;
                      }
                      if(evt.id == Event.SCROLL_LINE_DOWN) {
                          this.handleBar(evt) ;
                          return true ;
                      }
                      if(evt.id == Event.SCROLL_LINE_UP) {
                          this.handleBar(evt) ;
                          return true ;
                      }
                      if(evt.id == Event.SCROLL_PAGE_DOWN) {
                          this.handleBar(evt) ;
                          return true ;
                      }
                      if(evt.id == Event.SCROLL_PAGE_UP) {
                          this.handleBar(evt) ;
                          return true ;
                      }
                      else return false ;
                   }

                   public void handleBar(Event evt) {
                      int i1, i3, i4, i5 ;
                      Double V2 ;
                      double v1, v2, v3, v4, v5 ;
                      float fl1, fl3, fl4, fl5 ;

                      i1 = s1.getValue() ;

                      v1 = i1 * (machhi - machlo)/ 1000. + machlo ;
                      fl1 = (float) v1 ;
                      inleft.f1.setText(String.valueOf(fl1)) ;
                      mach0 = v1 ;
  
                      if(desmod == 1) {
                         i3 = s3.getValue() ;

                         v3 = i3 * (angmx - angmn)/ 1000. + angmn ;
                         fl3 = (float) v3 ;
                         inleft.f3.setText(String.valueOf(fl3)) ;
                         ang1 = v3 ;
                      }

                      comPute() ;    
                   }
                } // end of inright

                class Inleft extends Panel {
                   Moc outerparent ;
                   TextField f1, f3;
                   Label l1,l2,l3,l4 ;
   
                   Inleft (Moc target) {
            
                     outerparent = target ;
                     setLayout(new GridLayout(3,2,2,2)) ;
  
                     f1 = new TextField(String.valueOf(mach0),5) ;
                     f3 = new TextField(String.valueOf(ang1),5) ;
 
                     add(new Label(" ", Label.CENTER)) ;  
                     add(new Label(" ", Label.CENTER)) ;  
 
                     add(new Label("Mach", Label.CENTER)) ;  
                     add(f1) ;

                     add(new Label(" Angle", Label.CENTER)) ;
                     add(f3) ;
                  }
  
                  public boolean handleEvent(Event evt) {
                     if(evt.id == Event.ACTION_EVENT) {
                        this.handleText(evt) ;
                        return true ;
                     }
                     else return false ;
                  }
 
                  public void handleText(Event evt) {
                      Double V1,V2,V3,V4,V5 ;
                      double v1,v2,v3,v4,v5 ;
                      float fl1,fl2 ;
                      int i1,i3,i4,i5 ;

         // Mach number - range from machlo to machhi
                      V1 = Double.valueOf(f1.getText()) ;
                      v1 = V1.doubleValue() ;
                      if(v1 < machlo) {
                        v1 = machlo+ .02 ;
                        fl1 = (float) v1 ;
                        f1.setText(String.valueOf(fl1)) ;
                      }
                      if(v1 > machhi) {
                         v1 = machhi ;
                         fl1 = (float) v1 ;
                         f1.setText(String.valueOf(fl1)) ;
                      }
                      i1 = (int) (((v1 - machlo)/(machhi - machlo))*1000.) ;
                      inright.s1.setValue(i1) ;
                      mach0 = v1 ;

                      if(desmod == 1) {
         // ramp angle # 1  range from 0 to +30 for prob = 0
                          V3 = Double.valueOf(f3.getText()) ;
                          v3 = V3.doubleValue() ;
                          if(v3 < angmn) {
                            v3 = angmn ;
                            fl1 = (float) v3 ;
                            f3.setText(String.valueOf(fl1)) ;
                          }
                          if(v3 > angmx) {
                             v3 =  angmx ;
                             fl1 = (float) v3 ;
                             f3.setText(String.valueOf(fl1)) ;
                          }
                          i3 = (int) (((v3 - angmn)/(angmx-angmn))*1000.) ;
                          inright.s3.setValue(i3) ;
                          ang1 = v3 ;
                       }

                       comPute() ;
                  } // end handle Text
                }// end Inleft
              } // end Inlt

             class Jet1 extends Panel {
                Inright inright ;
                Inleft inleft ;

                Jet1 (Moc target) {
                             
                   outerparent = target ;
                   setLayout(new GridLayout(1,2,2,2)) ;
   
                   inleft = new Inleft(outerparent) ; 
                   inright = new Inright(outerparent) ;
 
                   add(inleft) ;
                   add(inright) ;
                }
   
                class Inright extends Panel {
                   Moc outerparent ;
                   Scrollbar s1,s3;

                   Inright (Moc target) {
                      int i1,i3 ;
 
                      outerparent = target ;
                      setLayout(new GridLayout(3,1,2,2)) ;

                      i1 = (int) (((mach0 - machlo)/(machhi - machlo))*1000.) ;
                      i3 = (int) (((nzht - nzhtlo)/(nzhthi - nzhtlo))*1000.) ;
 
                      s1 = new Scrollbar(Scrollbar.HORIZONTAL,i1,10,0,1000);
                      s3 = new Scrollbar(Scrollbar.HORIZONTAL,i3,10,0,1000);

                      add(new Label(" ", Label.CENTER)) ; 
                      add(s1) ;
                      add(s3) ;
                    }

                   public boolean handleEvent(Event evt) {
                      if(evt.id == Event.SCROLL_ABSOLUTE) {
                         this.handleBar(evt) ;
                         return true ;
                      }
                      if(evt.id == Event.SCROLL_LINE_DOWN) {
                          this.handleBar(evt) ;
                          return true ;
                      }
                      if(evt.id == Event.SCROLL_LINE_UP) {
                          this.handleBar(evt) ;
                          return true ;
                      }
                      if(evt.id == Event.SCROLL_PAGE_DOWN) {
                          this.handleBar(evt) ;
                          return true ;
                      }
                      if(evt.id == Event.SCROLL_PAGE_UP) {
                          this.handleBar(evt) ;
                          return true ;
                      }
                      else return false ;
                   }

                   public void handleBar(Event evt) {
                      int i1, i3, i4, i5 ;
                      Double V2 ;
                      double v1, v2, v3, v4, v5 ;
                      float fl1, fl3, fl4, fl5 ;

                      i1 = s1.getValue() ;
                      v1 = i1 * (machhi - machlo)/ 1000. + machlo ;
                      fl1 = (float) v1 ;
                      inleft.f1.setText(String.valueOf(fl1)) ;
                      mach0 = v1 ;
 
                      i3 = s3.getValue() ;
                      v3 = i3 * (nzhthi - nzhtlo)/ 1000. + nzhtlo ;
                      fl1 = (float) v3 ;
                      inleft.f3.setText(String.valueOf(fl1)) ;
                      nzht = v3 ;
 
                      comPute() ;    
                   }
                } // end of inright

                class Inleft extends Panel {
                   Moc outerparent ;
                   TextField f1, f2, f3;
                   Label l1,l2,l3,l4 ;
   
                   Inleft (Moc target) {
            
                     outerparent = target ;
                     setLayout(new GridLayout(3,2,2,2)) ;
  
                     f1 = new TextField(String.valueOf(mach0),5) ;
                     f3 = new TextField(String.valueOf(nzht),5) ;

                     add(new Label(" ", Label.CENTER)) ;  
                     add(new Label(" ", Label.CENTER)) ;  

                     add(new Label("Mach Exit", Label.CENTER)) ;  
                     add(f1) ;

                     add(new Label("Throat Height", Label.CENTER)) ;
                     add(f3) ;
                  }
  
                  public boolean handleEvent(Event evt) {
                     if(evt.id == Event.ACTION_EVENT) {
                        this.handleText(evt) ;
                        return true ;
                     }
                     else return false ;
                  }
 
                  public void handleText(Event evt) {
                      Double V1,V2,V3,V4,V5 ;
                      double v1,v2,v3,v4,v5 ;
                      float fl1,fl2 ;
                      int i1,i3,i4,i5 ;

         // Mach number - range from machlo to machhi
                      V1 = Double.valueOf(f1.getText()) ;
                      v1 = V1.doubleValue() ;
                      if(v1 < machlo) {
                        v1 = machlo+ .02 ;
                        fl1 = (float) v1 ;
                        f1.setText(String.valueOf(fl1)) ;
                      }
                      if(v1 > machhi) {
                         v1 = machhi ;
                         fl1 = (float) v1 ;
                         f1.setText(String.valueOf(fl1)) ;
                      }
                      i1 = (int) (((v1 - machlo)/(machhi - machlo))*1000.) ;
                      inright.s1.setValue(i1) ;
                      mach0 = v1 ;

         // Nozzle height
                       V3 = Double.valueOf(f3.getText()) ;
                       v3 = V3.doubleValue() ;
                       if(v3 < nzhtlo) {
                          v3 = nzhtlo ;
                          fl2 = (float) v3 ;
                          f3.setText(String.valueOf(fl2)) ;
                       }
                       if(v3 > nzhthi) {
                          v3 =nzhthi ;
                          fl2 = (float) v3 ;
                          f3.setText(String.valueOf(fl2)) ; 
                       }  
                       i3 = (int) (((v3 - nzhtlo)/(nzhthi - nzhtlo))*1000.) ;
                       inright.s3.setValue(i3) ;     
                       nzht = v3 ;

                       comPute() ;
                  } // end handle Text
                }// end Inleft
              } // end Jet1
           } // end In2
        } // end Up

       class Dwn extends Panel {
          Moc outerparent ;
          Oneang oneang ;
          Twoang twoang ;
          Extwdg extwdg ;
          Isnwdg isnwdg ;
          Jetinp jetinp ;
          Nzexp nzexp ;
          Blnkin blnkin;

          Dwn (Moc target) {
                             
             outerparent = target ;
             layin = new CardLayout() ;
             setLayout(layin) ;

             oneang = new Oneang(outerparent) ;
             twoang = new Twoang(outerparent) ;
             extwdg = new Extwdg(outerparent) ;
             isnwdg = new Isnwdg(outerparent) ;
             jetinp = new Jetinp(outerparent) ;
             nzexp  = new Nzexp(outerparent) ;
             blnkin = new Blnkin(outerparent) ;

             add ("first", oneang) ;
             add ("second", twoang) ;
             add ("third", extwdg) ;
             add ("fourth", isnwdg) ;
             add ("fifth", jetinp) ;
             add ("sixth", nzexp) ;
             add ("blank", blnkin);
          }
 
          class Blnkin extends Panel {
             Moc outerparent ;

             Blnkin (Moc target) {

               outerparent = target ;
               setLayout(new GridLayout(1,1,2,2)) ;

               add(new Label(" ", Label.CENTER)) ;
            }
          } // end Blnkin

          class Oneang extends Panel {
             Moc outerparent ;
             Inright inright ;
             Inleft inleft ;

             Oneang (Moc target) {

               outerparent = target ;
               setLayout(new GridLayout(1,2,2,2)) ;

               inleft = new Inleft(outerparent) ; 
               inright = new Inright(outerparent) ;
 
               add(inleft) ;
               add(inright) ;
            }

            class Inright extends Panel {
              Moc outerparent ;
              Scrollbar s3,s4 ;

              Inright (Moc target) {
                int i1, i3,i4 ;
 
                outerparent = target ;
                setLayout(new GridLayout(4,1,5,5)) ;

                i4 = (int) (((xr1 - xr1mn)/(xr1mx-xr1mn))*1000.) ;
 
                s4 = new Scrollbar(Scrollbar.HORIZONTAL,i4,10,0,1000);

                add(s4) ;
                add(new Label(" ", Label.CENTER)) ;
                add(new Label(" ", Label.CENTER)) ;
                add(new Label(" ", Label.CENTER)) ;
              }

              public boolean handleEvent(Event evt) {
                if(evt.id == Event.SCROLL_ABSOLUTE) {
                   this.handleBar(evt) ;
                   return true ;
                }
                if(evt.id == Event.SCROLL_LINE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_LINE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                else return false ;
             }

             public void handleBar(Event evt) {
                int i1, i3, i4, i5 ;
                Double V2 ;
                double v1, v2, v3, v4, v5 ;
                float fl1, fl3, fl4, fl5 ;

                i4 = s4.getValue() ;

                v4 = i4 * (xr1mx - xr1mn)/ 1000. + xr1mn ;
                fl4 = (float) v4 ;
                inleft.f4.setText(String.valueOf(fl4)) ;
                xr1 = v4 ;

                comPute() ;    
             }
           }

           class Inleft extends Panel {
             Moc outerparent ;
             TextField f3,f4;
             Label l1,l2,l3,l4 ;
   
             Inleft (Moc target) {
            
               outerparent = target ;
               setLayout(new GridLayout(4,2,5,5)) ;
  
               f4 = new TextField(String.valueOf(xr1),5) ;

               add(new Label(" X-Length ", Label.CENTER)) ;
               add(f4) ;

               add(new Label(" ", Label.CENTER)) ;
               add(new Label(" ", Label.CENTER)) ;

               add(new Label(" ", Label.CENTER)) ;
               add(new Label(" ", Label.CENTER)) ;

               add(new Label(" ", Label.CENTER)) ;
               add(new Label(" ", Label.CENTER)) ;
             }
  
             public boolean handleEvent(Event evt) {
               if(evt.id == Event.ACTION_EVENT) {
                  this.handleText(evt) ;
                  return true ;
               }
               else return false ;
             }

             public void handleText(Event evt) {
                Double V1,V2,V3,V4,V5 ;
                double v1,v2,v3,v4,v5 ;
                float fl1,fl2 ;
                int i1,i3,i4,i5 ;

 
         // ramp length # 1 for prob = 0
                V4 = Double.valueOf(f4.getText()) ;
                v4 = V4.doubleValue() ;
                if(v4 < xr1mn) {
                  v4 = xr1mn ;
                  fl1 = (float) v4 ;
                  f4.setText(String.valueOf(fl1)) ;
                }
                if(v4 > xr1mx) {
                   v4 =  xr1mx ;
                   fl1 = (float) v4 ;
                   f4.setText(String.valueOf(fl1)) ;
                } 
       
                i4 = (int) (((v4 - xr1mn)/(xr1mx-xr1mn))*1000.) ;
                inright.s4.setValue(i4) ;
                xr1 = v4 ;

                comPute() ;
             } // end handle Text
           }// end Inleft
        } // end oneang

        class Twoang extends Panel {
             Moc outerparent ;
             Inright inright ;
             Inleft inleft ;

             Twoang (Moc target) {

               outerparent = target ;
               setLayout(new GridLayout(1,2,2,2)) ;

               inleft = new Inleft(outerparent) ; 
               inright = new Inright(outerparent) ;
 
               add(inleft) ;
               add(inright) ;
            }
 
            class Inright extends Panel {
               Moc outerparent ;
               Scrollbar s3,s4,s5 ;

              Inright (Moc target) {
                int i3,i4,i5 ;
 
                outerparent = target ;
                setLayout(new GridLayout(4,1,5,5)) ;

                i4 = (int) (((ang2 - ang2mn)/(ang2mx-ang2mn))*1000.) ;
                i5 = (int) (((xr1 - xr1mn)/(xr1mx-xr1mn))*1000.) ;
 
                s4 = new Scrollbar(Scrollbar.HORIZONTAL,i4,10,0,1000);
                s5 = new Scrollbar(Scrollbar.HORIZONTAL,i5,10,0,1000);

                add(s5) ;
                add(s4) ;
                add(new Label(" ", Label.CENTER)) ;
                add(new Label(" ", Label.CENTER)) ;
              }

              public boolean handleEvent(Event evt) {
                if(evt.id == Event.SCROLL_ABSOLUTE) {
                   this.handleBar(evt) ;
                   return true ;
                }
                if(evt.id == Event.SCROLL_LINE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_LINE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                else return false ;
             }

             public void handleBar(Event evt) {
                int i1, i3, i4, i5 ;
                Double V2 ;
                double v1, v2, v3, v4, v5 ;
                float fl1, fl3, fl4, fl5 ;

                i4 = s4.getValue() ;
                i5 = s5.getValue() ;

                v4 = i4 * (ang2mx - ang2mn)/ 1000. + ang2mn ;
                fl4 = (float) v4 ;
                inleft.f4.setText(String.valueOf(fl4)) ;
                ang2 = v4 ;

                v5 = i5 * (xr1mx - xr1mn)/ 1000. + xr1mn ;
                fl5 = (float) v5 ;
                inleft.f5.setText(String.valueOf(fl5)) ;
                xr1 = v5 ;

                comPute() ;    
             }
           }

           class Inleft extends Panel {
             Moc outerparent ;
             TextField f3, f4, f5;
             Label l1,l2,l3,l4 ;
   
             Inleft (Moc target) {
            
               outerparent = target ;
               setLayout(new GridLayout(4,2,5,5)) ;
  
               f4 = new TextField(String.valueOf(ang2),5) ;
               f5 = new TextField(String.valueOf(xr1),5) ;

               add(new Label(" Length 1 ", Label.CENTER)) ;
               add(f5) ;

               add(new Label(" Angle 2", Label.CENTER)) ;
               add(f4) ;

               add(new Label(" ", Label.CENTER)) ;
               add(new Label(" ", Label.CENTER)) ;

               add(new Label(" ", Label.CENTER)) ;
               add(new Label(" ", Label.CENTER)) ;
             }
  
             public boolean handleEvent(Event evt) {
               if(evt.id == Event.ACTION_EVENT) {
                  this.handleText(evt) ;
                  return true ;
               }
               else return false ;
             }

             public void handleText(Event evt) {
                Double V1,V2,V3,V4,V5 ;
                double v1,v2,v3,v4,v5 ;
                float fl1,fl2 ;
                int i1,i3,i4,i5 ;
        
         // ramp angle # 2  range from 0 to +30 for prob = 2
                V4 = Double.valueOf(f4.getText()) ;
                v4 = V4.doubleValue() ;
                if(v4 < ang2mn) {
                  v4 = ang2mn ;
                  fl1 = (float) v4 ;
                  f4.setText(String.valueOf(fl1)) ;
                }
                if(v4 > ang2mx) {
                   v4 =  ang2mx ;
                   fl1 = (float) v4 ;
                   f4.setText(String.valueOf(fl1)) ;
                }
        
                i4 = (int) (((v4 - ang2mn)/(ang2mx-ang2mn))*1000.) ;
                inright.s4.setValue(i4) ;
                ang2 = v4 ;

         // distance between ramps for prob = 2
                V5 = Double.valueOf(f5.getText()) ;
                v5 = V5.doubleValue() ;
                if(v5 < xr1mn) {
                  v5 = xr1mn ;
                  fl1 = (float) v5 ;
                  f5.setText(String.valueOf(fl1)) ;
                }
                if(v5 > xr1mx) {
                   v5 =  xr1mx ;
                   fl1 = (float) v5 ;
                   f5.setText(String.valueOf(fl1)) ;
                }
        
                i5 = (int) (((v5 - xr1mn)/(xr1mx-xr1mn))*1000.) ;
                inright.s5.setValue(i5) ;
                xr1 = v5 ;

                comPute() ;
             } // end handle Text
           }// end Inleft
        } // end twoang

        class Extwdg extends Panel {
             Moc outerparent ;
             Inright inright ;
             Inleft inleft ;

             Extwdg (Moc target) {

               outerparent = target ;
               setLayout(new GridLayout(1,2,2,2)) ;

               inleft = new Inleft(outerparent) ; 
               inright = new Inright(outerparent) ;
 
               add(inleft) ;
               add(inright) ;
            }
 
            class Inright extends Panel {
               Moc outerparent ;
               Scrollbar s3,s4,s5,s6 ;

              Inright (Moc target) {
                int i3,i4,i5,i6 ;
 
                outerparent = target ;
                setLayout(new GridLayout(5,1,5,5)) ;

                i4 = (int) (((cwly - cwymn)/(cwymx-cwymn))*1000.) ;
                i5 = (int) (((cwlx - cwxmn)/(cwxmx-cwxmn))*1000.) ;
                i6 = (int) (((xr1 - xr1mn)/(xr1mx-xr1mn))*1000.) ;
 
                s4 = new Scrollbar(Scrollbar.HORIZONTAL,i4,10,0,1000);
                s5 = new Scrollbar(Scrollbar.HORIZONTAL,i5,10,0,1000);
                s6 = new Scrollbar(Scrollbar.HORIZONTAL,i6,10,0,1000);

                add(s6) ;
                add(s4) ;
                add(s5) ;
                add(new Label(" ", Label.CENTER)) ;
                add(new Label(" ", Label.CENTER)) ;
              }

              public boolean handleEvent(Event evt) {
                if(evt.id == Event.SCROLL_ABSOLUTE) {
                   this.handleBar(evt) ;
                   return true ;
                }
                if(evt.id == Event.SCROLL_LINE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_LINE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                else return false ;
             }

             public void handleBar(Event evt) {
                int i1, i3, i4, i5, i6 ;
                double v1, v2, v3, v4, v5, v6 ;
                float fl1, fl3, fl4, fl5, fl6 ;

                i4 = s4.getValue() ;
                i5 = s5.getValue() ;
                i6 = s6.getValue() ;

                v4 = i4 * (cwymx - cwymn)/ 1000. + cwymn ;
                fl4 = (float) v4 ;
                inleft.f4.setText(String.valueOf(fl4)) ;
                cwly = v4 ;

                v5 = i5 * (cwxmx - cwxmn)/ 1000. + cwxmn ;
                fl5 = (float) v5 ;
                inleft.f5.setText(String.valueOf(fl5)) ;
                cwlx = v5 ;

                v6 = i6 * (xr1mx - xr1mn)/ 1000. + xr1mn ;
                fl6 = (float) v6 ;
                inleft.f6.setText(String.valueOf(fl6)) ;
                xr1 = v6 ;

                comPute() ;    
             }
           }

           class Inleft extends Panel {
             Moc outerparent ;
             TextField f3, f4, f5, f6;
             Label l1,l2,l3,l4 ;
   
             Inleft (Moc target) {
            
               outerparent = target ;
               setLayout(new GridLayout(5,2,5,5)) ;

               f4 = new TextField(String.valueOf(cwly),5) ;
               f5 = new TextField(String.valueOf(cwlx),5) ;
               f6 = new TextField(String.valueOf(xr1),5) ;

               add(new Label(" Throat ", Label.CENTER)) ;
               add(f6) ;

               add(new Label(" Cowl Height", Label.CENTER)) ;
               add(f4) ;

               add(new Label(" Cowl X-Loc", Label.CENTER)) ;
               add(f5) ;

               add(new Label(" ", Label.CENTER)) ;
               add(new Label(" ", Label.CENTER)) ;

               add(new Label(" ", Label.CENTER)) ;
               add(new Label(" ", Label.CENTER)) ;
             }
  
             public boolean handleEvent(Event evt) {
               if(evt.id == Event.ACTION_EVENT) {
                  this.handleText(evt) ;
                  return true ;
               }
               else return false ;
             }

             public void handleText(Event evt) {
                Double V1,V2,V3,V4,V5 ;
                double v1,v2,v3,v4,v5 ;
                float fl1,fl2 ;
                int i1,i3,i4,i5 ;

         // ramp length # 1  for prob = 3
                V1 = Double.valueOf(f6.getText()) ;
                v1 = V1.doubleValue() ;
                if(v1 < xr1mn) {
                  v1 = xr1mn ;
                  fl1 = (float) v1 ;
                  f6.setText(String.valueOf(fl1)) ;
                }
                if(v1 > xr1mx) {
                   v1 =  xr1mx ;
                   fl1 = (float) v1 ;
                   f6.setText(String.valueOf(fl1)) ;
                }
        
                i1 = (int) (((v1 - xr1mn)/(xr1mx-xr1mn))*1000.) ;
                inright.s6.setValue(i1) ;
                xr1 = v1 ;

         // cowl height for prob = 3
                V4 = Double.valueOf(f4.getText()) ;
                v4 = V4.doubleValue() ;
                if(v4 < cwymn) {
                  v4 = cwymn ;
                  fl1 = (float) v4 ;
                  f4.setText(String.valueOf(fl1)) ;
                }
                if(v4 > cwymx) {
                   v4 =  cwymx ;
                   fl1 = (float) v4 ;
                   f4.setText(String.valueOf(fl1)) ;
                }
        
                i4 = (int) (((v4 - cwymn)/(cwymx-cwymn))*1000.) ;
                inright.s4.setValue(i4) ;
                cwly = v4 ;

         // cowl x location for prob = 3
                V5 = Double.valueOf(f5.getText()) ;
                v5 = V5.doubleValue() ;
                if(v5 < cwxmn) {
                  v5 = cwxmn ;
                  fl1 = (float) v5 ;
                  f5.setText(String.valueOf(fl1)) ;
                }
                if(v5 > cwxmx) {
                   v5 =  cwxmx ;
                   fl1 = (float) v5 ;
                   f5.setText(String.valueOf(fl1)) ;
                }
        
                i5 = (int) (((v5 - cwxmn)/(cwxmx-cwxmn))*1000.) ;
                inright.s4.setValue(i5) ;
                cwlx = v5 ;

                comPute() ;
             } // end handle Text
           }// end Inleft
        } // end extwdg

        class Isnwdg extends Panel {
             Moc outerparent ;
             Inright inright ;
             Inleft inleft ;

             Isnwdg (Moc target) {

               outerparent = target ;
               setLayout(new GridLayout(1,2,2,2)) ;

               inleft = new Inleft(outerparent) ; 
               inright = new Inright(outerparent) ;
 
               add(inleft) ;
               add(inright) ;
            }

            class Inright extends Panel {
              Moc outerparent ;
              Scrollbar s3,s2,s1,s6 ;

              Inright (Moc target) {
                int i1, i2, i3, i6 ;
 
                outerparent = target ;
                setLayout(new GridLayout(5,1,5,5)) ;

                i2 = (int) (((cwly - cwymn)/(cwymx-cwymn))*1000.) ;
                i1 = (int) (((macl - mlaslo)/(mlashi-mlaslo))*1000.) ;
                i6 = (int) (((xr1 - xr1mn)/(xr1mx-xr1mn))*1000.) ;
 
                s2 = new Scrollbar(Scrollbar.HORIZONTAL,i2,10,0,1000);
                s1 = new Scrollbar(Scrollbar.HORIZONTAL,i1,10,0,1000);
                s6 = new Scrollbar(Scrollbar.HORIZONTAL,i6,10,0,1000);

                add(new Label(" ", Label.CENTER)) ;
                add(s2) ;
                add(s1) ;
                add(s6) ;
                add(new Label(" ", Label.CENTER)) ;
              }

              public boolean handleEvent(Event evt) {
                if(evt.id == Event.SCROLL_ABSOLUTE) {
                   this.handleBar(evt) ;
                   return true ;
                }
                if(evt.id == Event.SCROLL_LINE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_LINE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                else return false ;
             }

             public void handleBar(Event evt) {
                int i1, i2, i3, i4, i5, i6 ;
                Double V2 ;
                double v1, v2, v3, v4, v5, v6 ;
                float fl1, fl2, fl3, fl4, fl5, fl6 ;

                i2 = s2.getValue() ;
                i1 = s1.getValue() ;
                i6 = s6.getValue() ;

                v2 = i2 * (cwymx - cwymn)/ 1000. + cwymn ;
                v1 = i1 * (mlashi - mlaslo)/ 1000. + mlaslo ;
                v6 = i6 * (xr1mx - xr1mn)/ 1000. + xr1mn ;

                fl2 = (float) v2 ;
                fl1 = (float) v1 ;
                fl6 = (float) v6 ;

                inleft.f2.setText(String.valueOf(fl2)) ;
                inleft.f1.setText(String.valueOf(fl1)) ;
                inleft.f6.setText(String.valueOf(fl6)) ;

                cwly = v2 ;
                macl = v1 ;
                xr1 = v6 ;

                comPute() ;    
             }
           }

           class Inleft extends Panel {
             Moc outerparent ;
             TextField f3,f2,f1,f4,f6;
             Label l1,l2,l3,l4 ;
   
             Inleft (Moc target) {
            
               outerparent = target ;
               setLayout(new GridLayout(5,2,5,5)) ;
  
               f2 = new TextField(String.valueOf(cwly),5) ;
               f1 = new TextField(String.valueOf(macl),5) ;
               f4 = new TextField(String.valueOf(numray),5) ;
               f6 = new TextField(String.valueOf(xr1),5) ;

               add(new Label(" ", Label.CENTER)) ;
               add(new Label(" ", Label.CENTER)) ;

               add(new Label(" Cowl Height", Label.CENTER)) ;
               add(f2) ;

               add(new Label(" Last Mach", Label.CENTER)) ;
               add(f1) ;

               add(new Label(" Throat ", Label.CENTER)) ;
               add(f6) ;

               add(new Label(" # of Rays", Label.CENTER)) ;
               add(f4) ;
             }
  
             public boolean handleEvent(Event evt) {
               if(evt.id == Event.ACTION_EVENT) {
                  this.handleText(evt) ;
                  return true ;
               }
               else return false ;
             }

             public void handleText(Event evt) {
                Double V1,V2,V3,V4,V5,V6 ;
                double v1,v2,v3,v4,v5,v6 ;
                float fl1,fl2 ;
                int i1,i2,i3,i4,i5,i6 ;

         // cowl height  for prob = 5
                V2 = Double.valueOf(f2.getText()) ;
                v2 = V2.doubleValue() ;
                if(v2 < cwymn) {
                  v2 = cwymn ;
                  fl1 = (float) v2 ;
                  f2.setText(String.valueOf(fl1)) ;
                }
                if(v2 > cwymx) {
                   v2 =  cwymx ;
                   fl1 = (float) v2 ;
                   f2.setText(String.valueOf(fl1)) ;
                }
         // last mach  for prob = 5
                V1 = Double.valueOf(f1.getText()) ;
                v1 = V1.doubleValue() ;
                if(v1 < mlaslo) {
                  v1 = mlaslo ;
                  fl1 = (float) v1 ;
                  f1.setText(String.valueOf(fl1)) ;
                }
                if(v1 > mlashi) {
                   v1 = mlashi ;
                   fl1 = (float) v1 ;
                   f1.setText(String.valueOf(fl1)) ;
                }

         // number of rays  for prob = 5
                V4 = Double.valueOf(f4.getText()) ;
                v4 = V4.doubleValue() ;
                if(v4 < 2.0) {
                  v4 = 2.0 ;
                  fl1 = (float) v4 ;
                  f4.setText(String.valueOf(fl1)) ;
                }
                if(v4 > 20.0) {
                   v4 = 20.0 ;
                   fl1 = (float) v4 ;
                   f4.setText(String.valueOf(fl1)) ;
                }
                numray = (int) v4 ;
 
         // location of throat  for prob = 5
                V6 = Double.valueOf(f6.getText()) ;
                v6 = V6.doubleValue() ;
                if(v6 < xr1mn) {
                  v6 = xr1mn ;
                  fl1 = (float) v6 ;
                  f6.setText(String.valueOf(fl1)) ;
                }
                if(v6 > xr1mx) {
                   v6 = xr1mx ;
                   fl1 = (float) v6 ;
                   f6.setText(String.valueOf(fl1)) ;
                }
      
                i2 = (int) (((v2 - cwymn)/(cwymx-cwymn))*1000.) ;
                i1 = (int) (((v1 - mlaslo)/(mlashi-mlaslo))*1000.) ;
                i6 = (int) (((v6 - xr1mn)/(xr1mx-xr1mn))*1000.) ;

                inright.s2.setValue(i2) ;
                inright.s1.setValue(i1) ;
                inright.s6.setValue(i6) ;

                cwly = v2 ;
                macl = v1 ;
                xr1 = v6 ;

                comPute() ;
             } // end handle Text
           }// end Inleft
        } // end Isnwdg

        class Jetinp extends Panel {
             Moc outerparent ;
             Inright inright ;
             Inleft inleft ;

             Jetinp (Moc target) {

               outerparent = target ;
               setLayout(new GridLayout(1,2,2,2)) ;

               inleft = new Inleft(outerparent) ; 
               inright = new Inright(outerparent) ;
 
               add(inleft) ;
               add(inright) ;
            }

            class Inright extends Panel {
              Moc outerparent ;
              Scrollbar s3,s4 ;

              Inright (Moc target) {
                int i1, i3,i4 ;
 
                outerparent = target ;
                setLayout(new GridLayout(5,1,5,5)) ;

                i4 = (int) (((nzlg - nzlglo)/(nzlghi - nzlglo))*1000.) ;
                i3 = (int) (((nprat - pratlo)/(prathi - pratlo))*1000.) ;
 
                s4 = new Scrollbar(Scrollbar.HORIZONTAL,i4,10,0,1000);
                s3 = new Scrollbar(Scrollbar.HORIZONTAL,i3,10,0,1000);

                add(s4) ;
                add(s3) ;
                add(new Label(" ", Label.CENTER)) ;
                add(new Label(" ", Label.CENTER)) ;
                add(new Label(" ", Label.CENTER)) ;
              }

              public boolean handleEvent(Event evt) {
                if(evt.id == Event.SCROLL_ABSOLUTE) {
                   this.handleBar(evt) ;
                   return true ;
                }
                if(evt.id == Event.SCROLL_LINE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_LINE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                else return false ;
             }

             public void handleBar(Event evt) {
                int i1, i3, i4, i5 ;
                Double V2 ;
                double v1, v2, v3, v4, v5 ;
                float fl1, fl3, fl4, fl5 ;

                i4 = s4.getValue() ;
                v4 = i4 * (nzlghi - nzlglo)/ 1000. + nzlglo ;
                fl4 = (float) v4 ;
                inleft.f4.setText(String.valueOf(fl4)) ;
                nzlg = v4 ;

                i3 = s3.getValue() ;
                v3 = i3 * (prathi - pratlo)/ 1000. + pratlo ;
                fl4 = (float) v3 ;
                inleft.f3.setText(String.valueOf(fl4)) ;
                nprat = v3 ;

                comPute() ;    
             }
           }

           class Inleft extends Panel {
             Moc outerparent ;
             TextField f3,f4,f5;
             Label l1,l2,l3,l4 ;
   
             Inleft (Moc target) {
            
               outerparent = target ;
               setLayout(new GridLayout(5,2,5,5)) ;
  
               f4 = new TextField(String.valueOf(nzlg),5) ;
               f3 = new TextField(String.valueOf(nprat),5) ;
               f5 = new TextField(String.valueOf(ncycle),5) ;

               add(new Label(" Noz-Length ", Label.CENTER)) ;
               add(f4) ;

               add(new Label(" NPR (P1/Pexit)", Label.CENTER)) ;
               add(f3) ;

               add(new Label("Number of Cycles ", Label.CENTER)) ;
               add(f5) ;

               add(new Label(" ", Label.CENTER)) ;
               add(new Label(" ", Label.CENTER)) ;

               add(new Label(" ", Label.CENTER)) ;
               add(new Label(" ", Label.CENTER)) ;
             }
  
             public boolean handleEvent(Event evt) {
               if(evt.id == Event.ACTION_EVENT) {
                  this.handleText(evt) ;
                  return true ;
               }
               else return false ;
             }

             public void handleText(Event evt) {
                Double V1,V2,V3,V4,V5 ;
                double v1,v2,v3,v4,v5 ;
                float fl1,fl2 ;
                int i1,i3,i4,i5 ;

 
         // nozzle length  for prob = 6
                V4 = Double.valueOf(f4.getText()) ;
                v4 = V4.doubleValue() ;
                if(v4 < nzlglo) {
                  v4 = nzlglo ;
                  fl1 = (float) v4 ;
                  f4.setText(String.valueOf(fl1)) ;
                }
                if(v4 > nzlghi) {
                   v4 =  nzlghi ;
                   fl1 = (float) v4 ;
                   f4.setText(String.valueOf(fl1)) ;
                } 
       
                i4 = (int) (((v4 - nzlglo)/(nzlghi-nzlglo))*1000.) ;
                inright.s4.setValue(i4) ;
                nzlg = v4 ;

         // static pressure ratio  for prob = 6
                V3 = Double.valueOf(f3.getText()) ;
                v3 = V3.doubleValue() ;
                if(v3 < pratlo) {
                  v3 = pratlo ;
                  fl1 = (float) v3 ;
                  f3.setText(String.valueOf(fl1)) ;
                }
                if(v3 > prathi) {
                   v3 =  prathi ;
                   fl1 = (float) v3 ;
                   f3.setText(String.valueOf(fl1)) ;
                } 
       
                i3 = (int) (((v3 - pratlo)/(prathi-pratlo))*1000.) ;
                inright.s3.setValue(i3) ;
                nprat = v3 ;

                V5 = Double.valueOf(f5.getText()) ;
                v5 = V5.doubleValue() ;
                ncycle = (int) v5 ;

                comPute() ;
             } // end handle Text
           }// end Inleft
        } // end Jetinp

        class Nzexp extends Panel {
             Moc outerparent ;
             Inright inright ;
             Inleft inleft ;

             Nzexp (Moc target) {

               outerparent = target ;
               setLayout(new GridLayout(1,2,2,2)) ;

               inleft = new Inleft(outerparent) ; 
               inright = new Inright(outerparent) ;
 
               add(inleft) ;
               add(inright) ;
            }

            class Inright extends Panel {
              Moc outerparent ;
              Scrollbar s3,s4 ;

              Inright (Moc target) {
                int i1, i3,i4 ;
 
                outerparent = target ;
                setLayout(new GridLayout(5,1,5,5)) ;

                i4 = (int) (((nzlg - nzlglo)/(nzlghi - nzlglo))*1000.) ;
                i3 = (int) (((delx - delxlo)/(delxhi - delxlo))*1000.) ;
 
                s4 = new Scrollbar(Scrollbar.HORIZONTAL,i4,10,0,1000);
                s3 = new Scrollbar(Scrollbar.HORIZONTAL,i3,10,0,1000);

                add(s4) ;
                add(s3) ;
                add(new Label(" ", Label.CENTER)) ;
                add(new Label(" ", Label.CENTER)) ;
                add(new Label(" ", Label.CENTER)) ;
              }

              public boolean handleEvent(Event evt) {
                if(evt.id == Event.SCROLL_ABSOLUTE) {
                   this.handleBar(evt) ;
                   return true ;
                }
                if(evt.id == Event.SCROLL_LINE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_LINE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_DOWN) {
                    this.handleBar(evt) ;
                    return true ;
                }
                if(evt.id == Event.SCROLL_PAGE_UP) {
                    this.handleBar(evt) ;
                    return true ;
                }
                else return false ;
             }

             public void handleBar(Event evt) {
                int i1, i3, i4, i5 ;
                Double V2 ;
                double v1, v2, v3, v4, v5 ;
                float fl1, fl3, fl4, fl5 ;

                i4 = s4.getValue() ;
                v4 = i4 * (nzlghi - nzlglo)/ 1000. + nzlglo ;
                fl4 = (float) v4 ;
//                inleft.f4.setText(String.valueOf(fl4)) ;
                inleft.f4.setText(String.valueOf(filter3(fl4))) ;
                nzlg = v4 ;

                i3 = s3.getValue() ;
                v3 = i3 * (delxhi - delxlo)/ 1000. + delxlo ;
                fl3 = (float) v3 ;
//                inleft.f3.setText(String.valueOf(fl3)) ;
                inleft.f3.setText(String.valueOf(filter3(fl3))) ;
                delx = v3 ;

                comPute() ;    
             }
           }

           class Inleft extends Panel {
             Moc outerparent ;
             TextField f3,f4,f5;
             Label l1,l2,l3,l4 ;
             Button bt1,bt2 ;
   
             Inleft (Moc target) {
            
               outerparent = target ;
               setLayout(new GridLayout(5,2,5,5)) ;
  
               f4 = new TextField(String.valueOf(nzlg),5) ;
               f3 = new TextField(String.valueOf(delx),5) ;
               f5 = new TextField(String.valueOf(numray),5) ;

               bt1 = new Button("Points") ;
               bt1.setBackground(Color.yellow) ;
               bt1.setForeground(Color.black) ;
               bt2 = new Button("Geometry") ;
               bt2.setBackground(Color.white) ;
               bt2.setForeground(Color.blue) ;

               add(new Label(" Noz-UpLngth ", Label.CENTER)) ;
               add(f4) ;

               add(new Label(" Delx ", Label.CENTER)) ;
               add(f3) ;

               add(new Label(" # of Rays ", Label.CENTER)) ;
               add(f5) ;

               add(bt1) ;
               add(bt2) ;

               add(new Label(" ", Label.CENTER)) ;
               add(new Label(" ", Label.CENTER)) ;
             }

              public boolean action(Event evt, Object arg) {
                 if(evt.target instanceof TextField) {
                    this.handleText(evt) ;
                    return true ;
                 }
                 if(evt.target instanceof Button) {
                     this.handleBut(evt,arg)  ;
                     return true ;
                 }
                 else return false ;
              }
  
              public void handleBut(Event evt, Object arg) {
                 String label = (String)arg ;
                 int i, j, k ;

                 comPute() ;

                 if(label.equals("Points")) {
                   if (prob == 7) layout.show(out, "sixth")  ;
                   if (prob == 8) layout.show(out, "seventh")  ;
                   if (prob == 9) layout.show(out, "seventh")  ;
                   bt1.setBackground(Color.yellow) ;
                   bt1.setForeground(Color.black) ;
                   bt2.setBackground(Color.white) ;
                   bt2.setForeground(Color.blue) ;
                 }
                 if(label.equals("Geometry")) {
                   layout.show(out, "eigth")  ;
                   bt1.setBackground(Color.white) ;
                   bt1.setForeground(Color.blue) ;
                   bt2.setBackground(Color.yellow) ;
                   bt2.setForeground(Color.black) ; 

                   num.out.nozg.prnt.appendText("\n ") ;
                   num.out.nozg.prnt.appendText("\n  X  \t Y \t X  \t Y \t X \t Y \t X \t Y") ;
                   j = numray / 4 ;
                   for(k=1; k<=j; ++k) {
                       num.out.nozg.prnt.appendText("\n"  
                          + filter3(wxbgn[(k-1)*4 + 1]) + "\t" + filter3(wybgn[(k-1)*4 + 1]) + "\t"
                          + filter3(wxbgn[(k-1)*4 + 2]) + "\t" + filter3(wybgn[(k-1)*4 + 2]) + "\t"
                          + filter3(wxbgn[(k-1)*4 + 3]) + "\t" + filter3(wybgn[(k-1)*4 + 3]) + "\t"
                          + filter3(wxbgn[(k-1)*4 + 4]) + "\t" + filter3(wybgn[(k-1)*4 + 4]) ) ;
                   }
                   i = (numray + 1) - (4 * j) ;
                   num.out.nozg.prnt.appendText("\n ") ;
                   for(k=1; k<=i; ++k) {
                       num.out.nozg.prnt.appendText( 
                          + filter3(wxbgn[4*j + k]) + "\t" + filter3(wybgn[4*j + k]) + "\t" ) ;
                   }
                 }

                 return ;
             }  
 
             public void handleText(Event evt) {
                Double V1,V2,V3,V4,V5 ;
                double v1,v2,v3,v4,v5 ;
                float fl1,fl2 ;
                int i1,i3,i4,i5 ;

         // nozzle length  for graphics
                V4 = Double.valueOf(f4.getText()) ;
                v4 = V4.doubleValue() ;
                if(v4 < nzlglo) {
                  v4 = nzlglo ;
                  fl1 = (float) v4 ;
                  f4.setText(String.valueOf(fl1)) ;
                }
                if(v4 > nzlghi) {
                   v4 =  nzlghi ;
                   fl1 = (float) v4 ;
                   f4.setText(String.valueOf(fl1)) ;
                } 
       
                i4 = (int) (((v4 - nzlglo)/(nzlghi-nzlglo))*1000.) ;
                inright.s4.setValue(i4) ;
                nzlg = v4 ;

         // delx 
                V3 = Double.valueOf(f3.getText()) ;
                v3 = V3.doubleValue() ;
                if(v3 < delxlo) {
                  v3 = delxlo ;
                  fl1 = (float) v3 ;
                  f3.setText(String.valueOf(fl1)) ;
                }
                if(v3 > delxhi) {
                   v3 =  delxhi ;
                   fl1 = (float) v3 ;
                   f3.setText(String.valueOf(fl1)) ;
                } 
       
                i3 = (int) (((v3 - delxlo)/(delxhi-delxlo))*1000.) ;
                inright.s3.setValue(i3) ;
                delx = v3 ;

         // number of rays  -- must be an even number
                V5 = Double.valueOf(f5.getText()) ;
                v5 = V5.doubleValue() ;
 
                numray = (int) v5 ;
                if(numray / 2 * 2 < numray) {
                    numray  = numray + 1 ;
                    f5.setText(String.valueOf(numray)) ;
                }

                comPute() ;
             } // end handle Text
           }// end Inleft
        } // end Nzexp
      } // end Dwn
     } // end Inp

     class Out extends Panel {
       Moc outerparent ;
       Wdg wdg;
       Con con ;
       Dub dub ;
       Diag diag;
       Isn isn ;
       Jet jet ;
       Noz noz ;
       Nozp nozp ;
       Nozg nozg ;

       Out (Moc target) { 
         outerparent = target ;
         layout = new CardLayout() ;
         setLayout(layout) ;

         wdg = new Wdg(outerparent) ;
         con = new Con(outerparent) ;
         dub = new Dub(outerparent) ;
         diag = new Diag(outerparent) ;
         isn = new Isn(outerparent) ;
         jet = new Jet(outerparent) ;
         noz = new Noz(outerparent) ;
         nozp = new Nozp(outerparent) ;
         nozg = new Nozg(outerparent) ;

         add ("first", wdg) ;
         add ("second", con) ;
         add ("last", diag) ;
         add ("third", dub) ;
         add ("fourth", isn) ;
         add ("fifth", jet) ;
         add ("sixth", noz) ;
         add ("seventh", nozp) ;
         add ("eigth", nozg) ;

       }

       class Wdg extends Panel {
          TextField o1, o2, o3, o4, o6, o7 ;
          Label lo1,lo2,lo3,lo4 ;

          Wdg (Moc target) {
             outerparent = target ;
             setLayout(new GridLayout(7,4,2,2)) ;
 
             o1 = new TextField() ;
             o1.setBackground(Color.black) ;
             o1.setForeground(Color.green) ;
             o2 = new TextField() ;
             o2.setBackground(Color.black) ;
             o2.setForeground(Color.green) ;
             o3 = new TextField() ;
             o3.setBackground(Color.black) ;
             o3.setForeground(Color.green) ;
             o4 = new TextField() ;
             o4.setBackground(Color.black) ;
             o4.setForeground(Color.green) ;
             o6 = new TextField() ;
             o6.setBackground(Color.black) ;
             o6.setForeground(Color.green) ;
             o7 = new TextField() ;
             o7.setBackground(Color.black) ;
             o7.setForeground(Color.green) ;

             lo3 = new Label(" ", Label.CENTER) ;
             lo3.setForeground(Color.blue) ;
             lo4 = new Label("Zone 1", Label.CENTER) ;
             lo4.setForeground(Color.blue) ;

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 

             add(lo3) ;
             add(lo4) ;
             add(new Label("p / p0", Label.CENTER)) ;  
             add(o2) ;

             add(new Label("Mach", Label.CENTER)) ;  
             add(o1) ;  
             add(new Label("T / T0", Label.CENTER)) ; 
             add(o6) ;

             add(new Label("Shock Angle", Label.CENTER)) ; 
             add(o3) ; 
             add(new Label("rho / rho0", Label.CENTER)) ;  
             add(o7) ; 

             add(new Label("pt / pt0", Label.CENTER)) ;  
             add(o4) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
          }
        } // end Wdg

       class Con extends Panel {
          TextField o1, o2, o3, o4, o6, o7, o8, o9, o10 ;
          Label lo1,lo2,lo3,lo4 ;
          Choice raych,draych ;

          Con (Moc target) {
             outerparent = target ;
             setLayout(new GridLayout(7,4,2,2)) ;
 
             o1 = new TextField() ;
             o1.setBackground(Color.black) ;
             o1.setForeground(Color.green) ;
             o2 = new TextField() ;
             o2.setBackground(Color.black) ;
             o2.setForeground(Color.green) ;
             o3 = new TextField() ;
             o3.setBackground(Color.black) ;
             o3.setForeground(Color.green) ;
             o4 = new TextField() ;
             o4.setBackground(Color.black) ;
             o4.setForeground(Color.green) ;
             o6 = new TextField() ;
             o6.setBackground(Color.black) ;
             o6.setForeground(Color.green) ;
             o7 = new TextField() ;
             o7.setBackground(Color.black) ;
             o7.setForeground(Color.green) ;
             o8 = new TextField() ;
             o8.setBackground(Color.black) ;
             o8.setForeground(Color.green) ;
             o9 = new TextField() ;
             o9.setBackground(Color.black) ;
             o9.setForeground(Color.green) ;
             o10 = new TextField() ;
             o10.setBackground(Color.black) ;
             o10.setForeground(Color.green) ;

             raych = new Choice() ;
             raych.addItem("Surface") ;
             raych.addItem("Ray 2");
             raych.addItem("Ray 3");
             raych.addItem("Ray 4");
             raych.addItem("Ray 5");
             raych.addItem("Ray 6");
             raych.addItem("Ray 7");
             raych.addItem("Ray 8");
             raych.addItem("Ray 9");
             raych.addItem("Shock") ;
             raych.setBackground(Color.white) ;
             raych.setForeground(Color.blue) ;
             raych.select(0) ;

             draych = new Choice() ;
             draych.addItem("Hide");
             draych.addItem("Show") ;
             draych.setBackground(Color.white) ;
             draych.setForeground(Color.blue) ;
             draych.select(0) ;

             lo3 = new Label("Ray:", Label.RIGHT) ;
             lo3.setForeground(Color.black) ;

             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;  
             add(lo3) ;
             add(raych) ;

             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label("Mach", Label.CENTER)) ;  
             add(o8) ;

             add(new Label("Surface Mach", Label.CENTER)) ;  
             add(o1) ; 
             add(new Label("p / p0", Label.CENTER)) ;  
             add(o2) ; 

             add(new Label("Shock Angle", Label.CENTER)) ; 
             add(o3) ; 
             add(new Label("T / T0", Label.CENTER)) ; 
             add(o6) ;

             add(new Label("pt / pt0", Label.CENTER)) ; 
             add(o4) ;
             add(new Label("rho / rho0", Label.CENTER)) ;  
             add(o7) ; 

             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label("Turning ", Label.CENTER)) ; 
             add(o9) ; 

             add(new Label("Ray Plot:", Label.CENTER)) ; 
             add(draych) ; 
             add(new Label("Angle", Label.CENTER)) ; 
             add(o10) ; 
          }

          public boolean handleEvent(Event evt) {
            if(evt.id == Event.ACTION_EVENT) {
               outzn  = 1 + raych.getSelectedIndex() ;
               drawray  = draych.getSelectedIndex() ;

               comPute() ;
               return true ;
            }
            else return false ;
          }    
        } // end Con

       class Dub extends Panel {
          TextField o1, o2, o3, o4, o5, o6, o7, o8, o9, o10 ;
          TextField o11, o12, o13, o18, o19 ;
          Label lo1,lo2,lo3,lo4 ;
          Label lab1, lab2, lab3 ;
          Choice zonch ;

          Dub (Moc target) {
             outerparent = target ;
             setLayout(new GridLayout(8,6,2,2)) ;
 
             o1 = new TextField() ;
             o1.setBackground(Color.black) ;
             o1.setForeground(Color.green) ;
             o2 = new TextField() ;
             o2.setBackground(Color.black) ;
             o2.setForeground(Color.green) ;
             o3 = new TextField() ;
             o3.setBackground(Color.black) ;
             o3.setForeground(Color.green) ;
             o4 = new TextField() ;
             o4.setBackground(Color.black) ;
             o4.setForeground(Color.green) ;
             o5 = new TextField() ;
             o5.setBackground(Color.black) ;
             o5.setForeground(Color.green) ;
             o6 = new TextField() ;
             o6.setBackground(Color.black) ;
             o6.setForeground(Color.green) ;
             o7 = new TextField() ;
             o7.setBackground(Color.black) ;
             o7.setForeground(Color.green) ;
             o8 = new TextField() ;
             o8.setBackground(Color.black) ;
             o8.setForeground(Color.green) ;
             o9 = new TextField() ;
             o9.setBackground(Color.black) ;
             o9.setForeground(Color.green) ;
             o10 = new TextField() ;
             o10.setBackground(Color.black) ;
             o10.setForeground(Color.green) ;
             o11 = new TextField() ;
             o11.setBackground(Color.black) ;
             o11.setForeground(Color.green) ;
             o12 = new TextField() ;
             o12.setBackground(Color.black) ;
             o12.setForeground(Color.green) ;
             o13 = new TextField() ;
             o13.setBackground(Color.black) ;
             o13.setForeground(Color.green) ;

             lab1 = new Label("Total", Label.RIGHT) ;
             lab1.setForeground(Color.red) ;
             lab2 = new Label("Turning", Label.RIGHT) ;
             lab2.setForeground(Color.red) ;
             o18 = new TextField() ;
             o18.setBackground(Color.black) ;
             o18.setForeground(Color.cyan) ;
             lab3 = new Label("Recovery", Label.RIGHT) ;
             lab3.setForeground(Color.red) ;
             o19 = new TextField() ;
             o19.setBackground(Color.black) ;
             o19.setForeground(Color.cyan) ;

             zonch = new Choice() ;
             zonch.setBackground(Color.white) ;
             zonch.setForeground(Color.red) ;
             zonch.addItem("Zone 0") ;
             zonch.addItem("Zone 1");
             zonch.addItem("Zone 2");
             zonch.addItem("Zone 3");
             zonch.addItem("Zone 4");
             zonch.addItem("Zone 5");
             zonch.addItem("Zone 6");
             zonch.select(1) ;

             lo3 = new Label("Zone:", Label.RIGHT) ;
             lo3.setForeground(Color.black) ;

             add(new Label(" ", Label.RIGHT)) ;  
             add(new Label(" ", Label.LEFT)) ;  
             add(lo3) ;
             add(zonch) ;
             add(new Label("Free ", Label.RIGHT)) ;  
             add(new Label(" Stream", Label.LEFT)) ;  

             add(new Label("Mach-up", Label.CENTER)) ;  
             add(o13) ;  
             add(new Label("Mach-dwn", Label.CENTER)) ;  
             add(o1) ;  
             add(new Label("Shock", Label.CENTER)) ;  
             add(o3) ;  

             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label("Angle", Label.CENTER)) ; 
             add(o2) ;  
             add(new Label("Turning", Label.CENTER)) ; 
             add(o4) ; 
   
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label("p/p(up)", Label.CENTER)) ;  
             add(o5) ;
             add(new Label("p/p0", Label.CENTER)) ;  
             add(o6) ;

             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label("pt/pt(up)", Label.CENTER)) ;  
             add(o7) ; 
             add(new Label("pt/pt0", Label.CENTER)) ;  
             add(o8) ; 
 
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;    
             add(new Label("T/T(up)", Label.CENTER)) ; 
             add(o9) ;
             add(new Label("T/T0", Label.CENTER)) ; 
             add(o10) ;

             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;     
             add(new Label("r/r(up)", Label.CENTER)) ;  
             add(o11) ; 
             add(new Label("r/r0", Label.CENTER)) ;  
             add(o12) ; 

             add(new Label(" ", Label.CENTER)) ;  
             add(lab1) ;     
             add(lab2) ;  
             add(o18) ; 
             add(lab3) ;  
             add(o19) ; 
          }

          public boolean handleEvent(Event evt) {
            if(evt.id == Event.ACTION_EVENT) {
               outzn  = 1 + zonch.getSelectedIndex() ;

               comPute() ;
               return true ;
            }
            else return false ;
          }    
        } // end Dub

       class Isn extends Panel {
          TextField o1, o2, o3, o4, o5, o6, o7, o8, o9, o10 ;
          TextField o11, o12, o13, o14, o15, o16, o17, o18, o19 ;
          Label lo1,lo2,lo3,lo4 ;
          Label lab1,lab2,lab3;
          Button bt1, bt2 ;
          Choice zonch ;

          Isn (Moc target) {
             outerparent = target ;
             setLayout(new GridLayout(8,6,2,2)) ;
 
             o1 = new TextField() ;
             o1.setBackground(Color.black) ;
             o1.setForeground(Color.green) ;
             o2 = new TextField() ;
             o2.setBackground(Color.black) ;
             o2.setForeground(Color.green) ;
             o3 = new TextField() ;
             o3.setBackground(Color.black) ;
             o3.setForeground(Color.green) ;
             o4 = new TextField() ;
             o4.setBackground(Color.black) ;
             o4.setForeground(Color.green) ;
             o5 = new TextField() ;
             o5.setBackground(Color.black) ;
             o5.setForeground(Color.green) ;
             o6 = new TextField() ;
             o6.setBackground(Color.black) ;
             o6.setForeground(Color.green) ;
             o7 = new TextField() ;
             o7.setBackground(Color.black) ;
             o7.setForeground(Color.green) ;
             o8 = new TextField() ;
             o8.setBackground(Color.black) ;
             o8.setForeground(Color.green) ;
             o9 = new TextField() ;
             o9.setBackground(Color.black) ;
             o9.setForeground(Color.green) ;
             o10 = new TextField() ;
             o10.setBackground(Color.black) ;
             o10.setForeground(Color.green) ;
             o11 = new TextField() ;
             o11.setBackground(Color.black) ;
             o11.setForeground(Color.green) ;
             o12 = new TextField() ;
             o12.setBackground(Color.black) ;
             o12.setForeground(Color.green) ;
             o13 = new TextField() ;
             o13.setBackground(Color.black) ;
             o13.setForeground(Color.green) ;
             o14 = new TextField() ;
             o14.setBackground(Color.black) ;
             o14.setForeground(Color.green) ;
             o15 = new TextField() ;
             o15.setBackground(Color.black) ;
             o15.setForeground(Color.green) ;
             o16 = new TextField() ;
             o16.setBackground(Color.black) ;
             o16.setForeground(Color.green) ;
             o17 = new TextField() ;
             o17.setBackground(Color.black) ;
             o17.setForeground(Color.green) ;

             lab1 = new Label("Total", Label.RIGHT) ;
             lab1.setForeground(Color.red) ;
             lab2 = new Label("Turning", Label.RIGHT) ;
             lab2.setForeground(Color.red) ;
             o18 = new TextField() ;
             o18.setBackground(Color.black) ;
             o18.setForeground(Color.cyan) ;
             lab3 = new Label("Recovery", Label.RIGHT) ;
             lab3.setForeground(Color.red) ;
             o19 = new TextField() ;
             o19.setBackground(Color.black) ;
             o19.setForeground(Color.cyan) ;

             bt1 = new Button("Design") ;
             bt1.setBackground(Color.yellow) ;
             bt1.setForeground(Color.black) ;
             bt2 = new Button("Off-Design") ;
             bt2.setBackground(Color.white) ;
             bt2.setForeground(Color.blue) ;

             zonch = new Choice() ;
             zonch.setBackground(Color.white) ;
             zonch.setForeground(Color.red) ;
             zonch.addItem("Zone 0") ;
             zonch.addItem("Zone 1");
             zonch.addItem("Zone 2");
             zonch.addItem("Zone 3");
             zonch.addItem("Zone 4");
             zonch.addItem("Zone 5");
             zonch.addItem("Zone 6");
             zonch.addItem("Zone 7");
             zonch.addItem("Zone 8");
             zonch.addItem("Zone 9");
             zonch.addItem("Zone 10");
             zonch.addItem("Zone 11");
             zonch.addItem("Zone 12");
             zonch.addItem("Zone 13");
             zonch.addItem("Zone 14");
             zonch.addItem("Zone 15");
             zonch.addItem("Zone 16");
             zonch.addItem("Zone 17");
             zonch.addItem("Zone 18");
             zonch.addItem("Zone 19");

             zonch.select(1) ;

             lo3 = new Label("Zone:", Label.RIGHT) ;
             lo3.setForeground(Color.black) ;

             add(new Label(" ", Label.RIGHT)) ;  
             add(new Label(" ", Label.LEFT)) ;  
             add(lo3) ;
             add(zonch) ;
             add(new Label("Free ", Label.RIGHT)) ;  
             add(new Label(" Stream", Label.LEFT)) ;  

             add(new Label("Mach-up", Label.CENTER)) ;  
             add(o13) ;  
             add(new Label("Mach-dwn", Label.CENTER)) ;  
             add(o1) ;  
             add(new Label("Shock", Label.CENTER)) ;  
             add(o3) ;  

             add(new Label("Mach-Angle ", Label.CENTER)) ;  
             add(o14) ;  
             add(new Label("Deflect", Label.CENTER)) ; 
             add(o2) ;  
             add(new Label("Turning", Label.CENTER)) ; 
             add(o4) ; 
   
             add(new Label("P-M Angle ", Label.CENTER)) ;  
             add(o15) ;  
             add(new Label("p/p(up)", Label.CENTER)) ;  
             add(o5) ;
             add(new Label("p/p0", Label.CENTER)) ;  
             add(o6) ;

             add(new Label(" X-Begin ", Label.CENTER)) ;  
             add(o16) ;  
             add(new Label("pt/pt(up)", Label.CENTER)) ;  
             add(o7) ; 
             add(new Label("pt/pt0", Label.CENTER)) ;  
             add(o8) ; 
 
             add(new Label(" Y-Begin ", Label.CENTER)) ;  
             add(o17) ;    
             add(new Label("T/T(up)", Label.CENTER)) ; 
             add(o9) ;
             add(new Label("T/T0", Label.CENTER)) ; 
             add(o10) ;

             add(bt1) ;  
             add(bt2) ;     
             add(new Label("r/r(up)", Label.CENTER)) ;  
             add(o11) ; 
             add(new Label("r/r0", Label.CENTER)) ;  
             add(o12) ; 

             add(new Label(" ", Label.CENTER)) ;  
             add(lab1) ;  
             add(lab2) ;  
             add(o18) ;  
             add(lab3) ;  
             add(o19) ;  
          }

          public boolean action(Event evt, Object arg) {
             if(evt.target instanceof Choice) {
                outzn  = 1 + zonch.getSelectedIndex() ;
 
                comPute() ;
                return true ;
             }
             if(evt.target instanceof Button) {
                handleBut (evt,arg) ;
                return true ;
             }
             else return false ;
          }

          public void handleBut(Event evt, Object arg) {
              String label = (String)arg ;

              if(label.equals("Design")) {
                 bt1.setBackground(Color.yellow) ;
                 bt1.setForeground(Color.black) ;
                 bt2.setBackground(Color.white) ;
                 bt2.setForeground(Color.blue) ;
                 layin.show(inp.dwn, "fourth")  ;
                 desmod = 1;
              }
              if(label.equals("Off-Design")) {
                 bt1.setBackground(Color.white) ;
                 bt1.setForeground(Color.blue) ;
                 bt2.setBackground(Color.yellow) ;
                 bt2.setForeground(Color.black) ;
                 layin.show(inp.dwn, "blank")  ;
                 desmod = 0 ;
              }

              comPute() ;
              return ;
          }  
        } // end Isn

       class Jet extends Panel {
          TextField o1, o2, o3, o4, o5, o6, o7, o8, o9, o10 ;
          TextField o11, o12, o13, o14, o15, o16, o17, o18, o19 ;
          Label lo1,lo2,lo3,lo4 ;
          Label lab1,lab2,lab3;
          Choice zonch ;

          Jet (Moc target) {
             outerparent = target ;
             setLayout(new GridLayout(8,6,2,2)) ;
 
             o1 = new TextField() ;
             o1.setBackground(Color.black) ;
             o1.setForeground(Color.green) ;
             o2 = new TextField() ;
             o2.setBackground(Color.black) ;
             o2.setForeground(Color.green) ;
             o3 = new TextField() ;
             o3.setBackground(Color.black) ;
             o3.setForeground(Color.green) ;
             o4 = new TextField() ;
             o4.setBackground(Color.black) ;
             o4.setForeground(Color.green) ;
             o5 = new TextField() ;
             o5.setBackground(Color.black) ;
             o5.setForeground(Color.green) ;
             o6 = new TextField() ;
             o6.setBackground(Color.black) ;
             o6.setForeground(Color.green) ;
             o7 = new TextField() ;
             o7.setBackground(Color.black) ;
             o7.setForeground(Color.green) ;
             o8 = new TextField() ;
             o8.setBackground(Color.black) ;
             o8.setForeground(Color.green) ;
             o9 = new TextField() ;
             o9.setBackground(Color.black) ;
             o9.setForeground(Color.green) ;
             o10 = new TextField() ;
             o10.setBackground(Color.black) ;
             o10.setForeground(Color.green) ;
             o11 = new TextField() ;
             o11.setBackground(Color.black) ;
             o11.setForeground(Color.green) ;
             o12 = new TextField() ;
             o12.setBackground(Color.black) ;
             o12.setForeground(Color.green) ;
             o13 = new TextField() ;
             o13.setBackground(Color.black) ;
             o13.setForeground(Color.green) ;
             o14 = new TextField() ;
             o14.setBackground(Color.black) ;
             o14.setForeground(Color.green) ;
             o15 = new TextField() ;
             o15.setBackground(Color.black) ;
             o15.setForeground(Color.green) ;
             o16 = new TextField() ;
             o16.setBackground(Color.black) ;
             o16.setForeground(Color.green) ;
             o17 = new TextField() ;
             o17.setBackground(Color.black) ;
             o17.setForeground(Color.green) ;

             zonch = new Choice() ;
             zonch.setBackground(Color.white) ;
             zonch.setForeground(Color.red) ;
             zonch.addItem("Zone 0") ;
             zonch.addItem("Zone 1");
             zonch.addItem("Zone 2");
             zonch.addItem("Zone 3");
             zonch.addItem("Zone 4");
             zonch.addItem("Zone 5");
             zonch.addItem("Zone 6");
             zonch.addItem("Zone 7");

             zonch.select(1) ;

             lo3 = new Label("Zone:", Label.RIGHT) ;
             lo3.setForeground(Color.black) ;

             add(new Label(" ", Label.RIGHT)) ;  
             add(new Label(" ", Label.LEFT)) ;  
             add(lo3) ;
             add(zonch) ;
             add(new Label(" Ratio to ", Label.RIGHT)) ;  
             add(new Label(" Throat", Label.LEFT)) ;  

             add(new Label("Mach-up", Label.CENTER)) ;  
             add(o13) ;  
             add(new Label("Mach-dwn", Label.CENTER)) ;  
             add(o1) ;  
             add(new Label("Deflect", Label.CENTER)) ; 
             add(o2) ;  

             add(new Label("P-M Angle ", Label.CENTER)) ;  
             add(o15) ;  
             add(new Label("p/p(up)", Label.CENTER)) ;  
             add(o5) ;
             add(new Label("p/p0", Label.CENTER)) ;  
             add(o6) ;

             add(new Label("Mach-Angle ", Label.CENTER)) ;  
             add(o14) ;  
             add(new Label("pt/pt(up)", Label.CENTER)) ;  
             add(o7) ; 
             add(new Label("pt/pt0", Label.CENTER)) ;  
             add(o8) ; 
 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ;  
 //            add(o16) ;    
             add(new Label("T/T(up)", Label.CENTER)) ; 
             add(o9) ;
             add(new Label("T/T0", Label.CENTER)) ; 
             add(o10) ;

             add(new Label(" ", Label.CENTER)) ;
             add(new Label(" ", Label.CENTER)) ;    
  //           add(o17) ;       
             add(new Label("r/r(up)", Label.CENTER)) ;  
             add(o11) ; 
             add(new Label("r/r0", Label.CENTER)) ;  
             add(o12) ; 

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ; 

             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;   
          }

          public boolean action(Event evt, Object arg) {
             if(evt.target instanceof Choice) {
                outzn  = 1 + zonch.getSelectedIndex() ;
 
                comPute() ;
                return true ;
             }
             else return false ;
          }

        } // end Jet

       class Noz extends Panel {
          TextField o1, o2, o3, o4, o5, o6, o7, o8, o9, o10 ;
          TextField o11, o12, o13, o14, o15, o16, o17, o18, o19, o20 ;
          Label lo1,lo2,lo3,lo4 ;
          Label lab1,lab2,lab3;
          TextField rowch,colch ;
          Button rowp,rowm,colp,colm ;

          Noz (Moc target) {
             outerparent = target ;
             setLayout(new GridLayout(10,6,10,2)) ;
 
             o1 = new TextField() ;
             o1.setBackground(Color.black) ;
             o1.setForeground(Color.green) ;
             o2 = new TextField() ;
             o2.setBackground(Color.black) ;
             o2.setForeground(Color.green) ;
             o13 = new TextField() ;
             o13.setBackground(Color.black) ;
             o13.setForeground(Color.green) ;

             o5 = new TextField() ;
             o5.setBackground(Color.black) ;
             o5.setForeground(Color.green) ;
             o6 = new TextField() ;
             o6.setBackground(Color.black) ;
             o6.setForeground(Color.green) ;
             o15 = new TextField() ;
             o15.setBackground(Color.black) ;
             o15.setForeground(Color.green) ;

             o9 = new TextField() ;
             o9.setBackground(Color.black) ;
             o9.setForeground(Color.green) ;
             o10 = new TextField() ;
             o10.setBackground(Color.black) ;
             o10.setForeground(Color.green) ;
             o14 = new TextField() ;
             o14.setBackground(Color.black) ;
             o14.setForeground(Color.green) ;

             o11 = new TextField() ;
             o11.setBackground(Color.black) ;
             o11.setForeground(Color.green) ;
             o12 = new TextField() ;
             o12.setBackground(Color.black) ;
             o12.setForeground(Color.green) ;
             o16 = new TextField() ;
             o16.setBackground(Color.black) ;
             o16.setForeground(Color.green) ;

             o3 = new TextField() ;
             o3.setBackground(Color.black) ;
             o3.setForeground(Color.cyan) ;
             o4 = new TextField() ;
             o4.setBackground(Color.black) ;
             o4.setForeground(Color.cyan) ;

             o7 = new TextField() ;
             o7.setBackground(Color.black) ;
             o7.setForeground(Color.cyan) ;
             o8 = new TextField() ;
             o8.setBackground(Color.black) ;
             o8.setForeground(Color.cyan) ;

             o17 = new TextField() ;
             o17.setBackground(Color.black) ;
             o17.setForeground(Color.cyan) ;
             o18 = new TextField() ;
             o18.setBackground(Color.black) ;
             o18.setForeground(Color.cyan) ;

             o19 = new TextField() ;
             o19.setBackground(Color.black) ;
             o19.setForeground(Color.cyan) ;
             o20 = new TextField() ;
             o20.setBackground(Color.black) ;
             o20.setForeground(Color.cyan) ;

             rowch = new TextField(String.valueOf(row),5) ;
             rowch.setBackground(Color.white) ;
             rowch.setForeground(Color.black) ;
             colch = new TextField(String.valueOf(col),5) ;
             colch.setBackground(Color.white) ;
             colch.setForeground(Color.black) ;

             rowp = new Button("Row UP") ;
             rowp.setBackground(Color.white) ;
             rowp.setForeground(Color.blue) ;
             rowm = new Button("Row DOWN") ;
             rowm.setBackground(Color.white) ;
             rowm.setForeground(Color.blue) ;

             colp = new Button("Col UP") ;
             colp.setBackground(Color.white) ;
             colp.setForeground(Color.red) ;
             colm = new Button("Col DOWN") ;
             colm.setBackground(Color.white) ;
             colm.setForeground(Color.red) ;

             lo3 = new Label("Row:", Label.RIGHT) ;
             lo3.setForeground(Color.black) ;
             lo4 = new Label("Column:", Label.RIGHT) ;
             lo4.setForeground(Color.black) ;

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" Geometry ", Label.CENTER)) ;  
             add(rowp) ; 
             add(rowch) ;
             add(rowm) ;

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ;
             add(colp) ;
             add(colch) ;  
             add(colm) ;  

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label("X low left ", Label.CENTER)) ; 
             add(o3) ;
             add(new Label("Y low left ", Label.CENTER)) ; 
             add(o4) ;

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ;
             add(new Label("X up left ", Label.CENTER)) ;  
             add(o7) ;  
             add(new Label("Y up left ", Label.CENTER)) ;  
             add(o8) ; 

             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label("X up right ", Label.CENTER)) ;  
             add(o17) ; 
             add(new Label("Y up right ", Label.CENTER)) ;  
             add(o18) ;   

             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label("X low right ", Label.CENTER)) ;  
             add(o19) ; 
             add(new Label("Y low right ", Label.CENTER)) ;  
             add(o20) ;   

             add(new Label("Mach", Label.CENTER)) ;  
             add(o1) ;  
             add(new Label("Deflect", Label.CENTER)) ; 
             add(o2) ;  
             add(new Label("Turning", Label.CENTER)) ;  
             add(o13) ;  

             add(new Label("P-M Angle ", Label.CENTER)) ;  
             add(o15) ;   
             add(new Label("p/p(up)", Label.CENTER)) ;  
             add(o5) ;
             add(new Label("p/p0", Label.CENTER)) ;  
             add(o6) ;

             add(new Label("Mach-Angle ", Label.CENTER)) ;  
             add(o14) ;  
             add(new Label("T/T(up)", Label.CENTER)) ; 
             add(o9) ;
             add(new Label("T/T0", Label.CENTER)) ; 
             add(o10) ;

             add(new Label("A/A*", Label.CENTER)) ;
             add(o16) ;  
             add(new Label("r/r(up)", Label.CENTER)) ;  
             add(o11) ; 
             add(new Label("r/r0", Label.CENTER)) ;  
             add(o12) ; 

          }

          public boolean action(Event evt, Object arg) {
              if(evt.target instanceof TextField) {
                    this.handleText(evt) ;
                    return true ;
              }
              if(evt.target instanceof Button) {
                   this.handleBut(evt,arg)  ;
                   return true ;
              }
              else return false ;
           }

          public void handleText(Event evt) {
             Double V1,V2,V3,V4,V5 ;
             double v1,v2,v3,v4,v5 ;
             float fl1,fl2 ;
             int i1,i3,i4,i5 ;

         // row index
             V1 = Double.valueOf(rowch.getText()) ;
             v1 = V1.doubleValue() ;
             if(v1 < rowlo) {
                v1 = rowlo ;
                fl1 = (float) v1 ;
                rowch.setText(String.valueOf(fl1)) ;
             }
             if(v1 > rowhi) {
                v1 = rowhi ;
                fl1 = (float) v1 ;
                rowch.setText(String.valueOf(fl1)) ;
             }
             row = (int) v1 ;

         // column index
             V2 = Double.valueOf(colch.getText()) ;
             v2 = V2.doubleValue() ;
             if(v2 < collo) {
                v2 = collo ;
                fl1 = (float) v2 ;
                colch.setText(String.valueOf(fl1)) ;
             }
             if(v2 > colhi) {
                v2 = colhi ;
                fl1 = (float) v2 ;
                colch.setText(String.valueOf(fl1)) ;
             }
             col = (int) v2 ;

             comPute() ;
          }

          public void handleBut(Event evt, Object arg) {
              String label = (String)arg ;

              if(label.equals("Row UP")) {
                 row = row + 1 ;
                 if(row > rowhi) row = rowhi ;
                 rowch.setText(String.valueOf(row)) ;
              }
              if(label.equals("Row DOWN")) {
                 row = row - 1 ;
                 if(row < rowlo) row = rowlo ;
                 rowch.setText(String.valueOf(row)) ;
              }
              if(label.equals("Col UP")) {
                 col = col + 1 ;
                 if(col > colhi) col = colhi ;
                 colch.setText(String.valueOf(col)) ;
              }
              if(label.equals("Col DOWN")) {
                 col = col - 1 ;
                 if(col < collo) col = collo ;
                 colch.setText(String.valueOf(col)) ;
              }
              comPute() ;
              return ;
          }  

        } // end Noz

       class Nozp extends Panel {
          TextField o1, o2, o3, o4, o5, o6, o7, o8, o9, o10 ;
          TextField o11, o12, o13, o14, o15, o16, o17, o18, o19, o20 ;
          Label lo1,lo2,lo3,lo4 ;
          Label lab1,lab2,lab3;
          TextField rowch,colch ;
          Button rowp,rowm,colp,colm ;

          Nozp (Moc target) {
             outerparent = target ;
             setLayout(new GridLayout(10,6,10,2)) ;
 
             o1 = new TextField() ;
             o1.setBackground(Color.black) ;
             o1.setForeground(Color.green) ;
             o2 = new TextField() ;
             o2.setBackground(Color.black) ;
             o2.setForeground(Color.green) ;
             o13 = new TextField() ;
             o13.setBackground(Color.black) ;
             o13.setForeground(Color.green) ;

             o5 = new TextField() ;
             o5.setBackground(Color.black) ;
             o5.setForeground(Color.green) ;
             o6 = new TextField() ;
             o6.setBackground(Color.black) ;
             o6.setForeground(Color.green) ;
             o15 = new TextField() ;
             o15.setBackground(Color.black) ;
             o15.setForeground(Color.green) ;

             o9 = new TextField() ;
             o9.setBackground(Color.black) ;
             o9.setForeground(Color.green) ;
             o10 = new TextField() ;
             o10.setBackground(Color.black) ;
             o10.setForeground(Color.green) ;
             o14 = new TextField() ;
             o14.setBackground(Color.black) ;
             o14.setForeground(Color.green) ;

             o11 = new TextField() ;
             o11.setBackground(Color.black) ;
             o11.setForeground(Color.green) ;
             o12 = new TextField() ;
             o12.setBackground(Color.black) ;
             o12.setForeground(Color.green) ;
             o16 = new TextField() ;
             o16.setBackground(Color.black) ;
             o16.setForeground(Color.green) ;

             o3 = new TextField() ;
             o3.setBackground(Color.black) ;
             o3.setForeground(Color.cyan) ;
             o4 = new TextField() ;
             o4.setBackground(Color.black) ;
             o4.setForeground(Color.cyan) ;

             o7 = new TextField() ;
             o7.setBackground(Color.black) ;
             o7.setForeground(Color.cyan) ;
             o8 = new TextField() ;
             o8.setBackground(Color.black) ;
             o8.setForeground(Color.cyan) ;

             o17 = new TextField() ;
             o17.setBackground(Color.black) ;
             o17.setForeground(Color.cyan) ;
             o18 = new TextField() ;
             o18.setBackground(Color.black) ;
             o18.setForeground(Color.cyan) ;

             o19 = new TextField() ;
             o19.setBackground(Color.black) ;
             o19.setForeground(Color.cyan) ;
             o20 = new TextField() ;
             o20.setBackground(Color.black) ;
             o20.setForeground(Color.cyan) ;

             rowch = new TextField(String.valueOf(row),5) ;
             rowch.setBackground(Color.white) ;
             rowch.setForeground(Color.black) ;
             colch = new TextField(String.valueOf(col),5) ;
             colch.setBackground(Color.white) ;
             colch.setForeground(Color.black) ;

             rowp = new Button("Row UP") ;
             rowp.setBackground(Color.white) ;
             rowp.setForeground(Color.blue) ;
             rowm = new Button("Row DOWN") ;
             rowm.setBackground(Color.white) ;
             rowm.setForeground(Color.blue) ;

             colp = new Button("Col UP") ;
             colp.setBackground(Color.white) ;
             colp.setForeground(Color.red) ;
             colm = new Button("Col DOWN") ;
             colm.setBackground(Color.white) ;
             colm.setForeground(Color.red) ;

             lo3 = new Label("Vary Q:", Label.RIGHT) ;
             lo3.setForeground(Color.black) ;
             lo4 = new Label("Vary R:", Label.RIGHT) ;
             lo4.setForeground(Color.black) ;

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" Geometry", Label.CENTER)) ;  
             add(rowp) ; 
             add(rowch) ;
             add(rowm) ;

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(colp) ;
             add(colch) ;  
             add(colm) ; 

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ;
             add(new Label(" R ", Label.CENTER)) ;  
             add(o8) ; 
             add(new Label(" Q ", Label.CENTER)) ;  
             add(o7) ;  

             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label("X ", Label.CENTER)) ; 
             add(o3) ;
             add(new Label("Y ", Label.CENTER)) ; 
             add(o4) ;

             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" alpha ", Label.CENTER)) ;  
             add(o17) ; 
             add(new Label(" beta ", Label.CENTER)) ;  
             add(o18) ;   

             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ; 
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label(" ", Label.CENTER)) ;  
             add(new Label("Mach", Label.CENTER)) ;  
             add(o1) ;  
             add(new Label("Deflect", Label.CENTER)) ; 
             add(o2) ;  
             add(new Label("Turning", Label.CENTER)) ;  
             add(o13) ;  

             add(new Label("P-M Angle ", Label.CENTER)) ;  
             add(o15) ;  
             add(new Label("p/p(up)", Label.CENTER)) ;  
             add(o5) ;
             add(new Label("p/p0", Label.CENTER)) ;  
             add(o6) ;

             add(new Label("Mach-Angle ", Label.CENTER)) ;  
             add(o14) ;  
             add(new Label("T/T(up)", Label.CENTER)) ; 
             add(o9) ;
             add(new Label("T/T0", Label.CENTER)) ; 
             add(o10) ;
 
             add(new Label("A/A*", Label.CENTER)) ;
             add(o16) ; 
             add(new Label("r/r(up)", Label.CENTER)) ;  
             add(o11) ; 
             add(new Label("r/r0", Label.CENTER)) ;  
             add(o12) ; 

/*
             add(new Label("X lower right ", Label.CENTER)) ;  
             add(o19) ; 
             add(new Label("Y lower right ", Label.CENTER)) ;  
             add(o20) ; 
*/  
          }

          public boolean action(Event evt, Object arg) {
              if(evt.target instanceof TextField) {
                    this.handleText(evt) ;
                    return true ;
              }
              if(evt.target instanceof Button) {
                   this.handleBut(evt,arg)  ;
                   return true ;
              }
              else return false ;
           }

          public void handleText(Event evt) {
             Double V1,V2,V3,V4,V5 ;
             double v1,v2,v3,v4,v5 ;
             float fl1,fl2 ;
             int i1,i3,i4,i5 ;

         // row index
             V1 = Double.valueOf(rowch.getText()) ;
             v1 = V1.doubleValue() ;
             if(v1 < rowlo) {
                v1 = rowlo ;
                fl1 = (float) v1 ;
                rowch.setText(String.valueOf(fl1)) ;
             }
             if(v1 > rowhi) {
                v1 = rowhi ;
                fl1 = (float) v1 ;
                rowch.setText(String.valueOf(fl1)) ;
             }
             row = (int) v1 ;

         // column index
             V2 = Double.valueOf(colch.getText()) ;
             v2 = V2.doubleValue() ;
             if(v2 < collo) {
                v2 = collo ;
                fl1 = (float) v2 ;
                colch.setText(String.valueOf(fl1)) ;
             }
             if(v2 > colhi) {
                v2 = colhi ;
                fl1 = (float) v2 ;
                colch.setText(String.valueOf(fl1)) ;
             }
             col = (int) v2 ;

             comPute() ;
          }

          public void handleBut(Event evt, Object arg) {
              String label = (String)arg ;

              if(label.equals("Row UP")) {
                 row = row + 1 ;
                 if(row > rowhi) row = rowhi ;
                 rowch.setText(String.valueOf(row)) ;
              }
              if(label.equals("Row DOWN")) {
                 row = row - 1 ;
                 if(row < rowlo) row = rowlo ;
                 rowch.setText(String.valueOf(row)) ;
              }
              if(label.equals("Col UP")) {
                 col = col + 1 ;
                 if(col > colhi) col = colhi ;
                 colch.setText(String.valueOf(col)) ;
              }
              if(label.equals("Col DOWN")) {
                 col = col - 1 ;
                 if(col < collo) col = collo ;
                 colch.setText(String.valueOf(col)) ;
              }

              comPute() ;
              return ;
          }  

        } // end Nozp

        class Nozg extends Panel {
          TextArea prnt ;

          Nozg (Moc target) {
            outerparent = target ;
            setLayout(new GridLayout(1,1,10,10)) ;

            prnt = new TextArea() ;
            prnt.setEditable(false) ;

            prnt.appendText("MOC Version 1.6g ") ;
            add(prnt) ;
          }
        } //Nozg

        class Diag extends Panel {
          TextArea prnt ;

          Diag (Moc target) {
            outerparent = target ;
            setLayout(new GridLayout(1,1,10,10)) ;

            prnt = new TextArea() ;
            prnt.setEditable(false) ;

            prnt.appendText("MOC Version 1.6g ") ;
            add(prnt) ;
          }
        } //Diag
     } // end Out
  } // end Num

  class Viewer extends Canvas 
         implements Runnable{
     Moc outerparent ;
     Thread runner ;
     Point locate,anchor;

     Viewer (Moc target) {
         setBackground(Color.white) ;
         runner = null ;
     }

     public boolean mouseDown(Event evt, int x, int y) {
        anchor = new Point(x,y) ;
        return true;
     }

     public boolean mouseUp(Event evt, int x, int y) {
        handleb(x,y) ;
        return true;
     }
 
     public boolean mouseDrag(Event evt, int x, int y) {
        handle(x,y) ;
        return true;
     }

     public void handle(int x, int y) {
         // determine location
         if (x >= 76 && x <= 700) {
            if (y >= 21 && y <= 300) {
               locate = new Point(x,y) ;
               yt =  yt + (int) (.2*(locate.y - anchor.y)) ;
               xt =  xt + (int) (.4*(locate.x - anchor.x))  ;
            }
         }
         if (x <= 30 ) {  
            if (y >= 20 && y <= 240 ) {   // zoom widget
              sldloc = y ;
              if (sldloc < 20) sldloc = 20;
              if (sldloc > 240) sldloc = 240;
              fact = fac1 + (240 - sldloc)*fac2 ;
            }
         }
/*
 num.out.diag.prnt.appendText("\n xt = " + String.valueOf(xt) +
                              "\n yt = " + String.valueOf(yt) +
                              "\n sldloc = " + String.valueOf(sldloc) +
                              "\n fact = " + String.valueOf(filter3(fact)) +
                              "\n fac1 = " + String.valueOf(filter3(fac1)) +
                              "\n fac2 = " + String.valueOf(filter3(fac2)) ) ;
*/
     }
 
     public void handleb(int x, int y) {
         if (y >= 245 && y <= 260 ) { 

            if (x >= 850 && x <= 900 ) {   // find button
             xt = 80; yt = 20; sldloc = 140;
             fac1 = .625; fac2 = .015;
             fact = 1.9 ;
            }

            if (x >= 40 && x <= 55) {    // decrease factor
              fac1 = fac1 / 2.0 ;
              fac2 = fac2 / 2.0 ;
              fact = fac1 + (240 - sldloc)*fac2 ;
            }
            if (x >= 60 && x <= 75) {    // increase factor
              fac1 = fac1 * 2.0 ;
              fac2 = fac2 * 2.0 ;
              fact = fac1 + (240 - sldloc)*fac2 ;
            }
         }
         view.repaint() ;
     }

     public void start() {
        if (runner == null) {
           runner = new Thread(this) ;
           runner.start() ;
        }
     }

     public void run() {
       int timer ;

       timer = 100 ;
       while (true) {
          try { Thread.sleep(timer); }
          catch (InterruptedException e) {}
          view.repaint() ;
       }
     }
 
     public void update(Graphics g) {
         view.paint(g) ;
     }

     public void paint(Graphics g) {
       int i,j,k ;
       int exes[] = new int[10] ;
       int whys[] = new int[10] ;
       int yorgn = 200 ;
       int xorgn = 50 ;

       offsGg.setColor(Color.white) ;
       offsGg.fillRect(0,0,900,500) ;
       offsGg.setColor(Color.blue) ;
       offsGg.drawString("Upstream", 30, 20) ;
       offsGg.drawString("Flow", xorgn-10, 35) ;
       offsGg.setColor(Color.black) ;
       offsGg.fillRect(xorgn-10,42,30,9) ;
       exes[0] = xorgn + 30 ;
       whys[0] = 46;
       exes[1] = xorgn + 20;
       whys[1] = 56;
       exes[2] = xorgn + 20;
       whys[2] = 36;
       Polygon poly1 = new Polygon(exes,whys,3) ;
       offsGg.fillPolygon(poly1) ;
       offsGg.setColor(Color.red) ;
       offsGg.drawString("0",xorgn + 5,70) ;
          // draw geometry
       offsGg.setColor(Color.red) ;
       offsGg.drawString("1", exes[1]-70, whys[1]-10) ;

//  basic flow problems - cone - ramp - double ramp 
       if (prob <= 3) {
// draw ramps
         for (i = 1; i <= nramps; ++i) {
           exes[0] = xorgn + (int) (.5*fact*wxbgn[i]) + xt ;
           whys[0] = yorgn - (int) (.5*fact*(wxbgn[i]*wslope[i] + winter[i])) + yt ;
           exes[1] = xorgn + (int) (.5*fact*wxnd[i]) + xt ;
           whys[1] = yorgn - (int) (.5*fact*(wxnd[i]*wslope[i] + winter[i])) + yt ;
           exes[2] = exes[1] ;
           whys[2] = whys[0] ;
           offsGg.setColor(Color.red) ;
           Polygon poly = new Polygon(exes,whys,3) ;
           offsGg.fillPolygon(poly) ;
           if (prob == 0 || prob == 3) {
             offsGg.drawString("1", exes[0]+20, whys[0]-10) ;
           }
           if (prob == 2 && i == 1) {
             offsGg.drawString("1", exes[0]+20, whys[0]-10) ;
           }
           if (prob == 2 && i == 2) {
             offsGg.drawString("2", exes[0]+20, whys[0]-10) ;
           }
           if (i > 1) {
               whys[1] = whys[0] ;
               whys[2] = yorgn + yt ;
               whys[3] = whys[2] ;
               exes[3] = exes[0] ;
               offsGg.fillPolygon(exes,whys,4) ;
           }
           if (prob <=2) {
             exes[0] = xorgn + (int) (.5*fact*wxbgn[i]) + xt ;
             whys[0] = yorgn - (int) (.5*fact*(wxbgn[i]*wslope[i] + winter[i])) + yt ;
             exes[1] = xorgn + (int) (.5*fact*xlong) + xt ;
             whys[1] = yorgn - (int) (.5*fact*(xlong*wslope[i] + winter[i])) + yt ;
             offsGg.setColor(Color.black) ;
             offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             offsGg.setColor(Color.white) ;
             offsGg.drawString("Angle", exes[1]-50, whys[0]-5) ;
           }
           if ( prob == 3) {
             exes[0] = xorgn + (int) (.5*fact*wxnd[1]) + xt ;
             whys[0] = yorgn - (int) (.5*fact*(wxnd[1]*wslope[1] + winter[1])) + yt ;
             exes[1] = xorgn + (int) (.5*fact*xlong) + xt ;
             whys[1] = yorgn - 40 + yt;
             exes[2] = exes[1] ;
             whys[2] = yorgn + yt ;
             exes[3] = exes[0] ;
             whys[3] = whys[2] ;
             offsGg.setColor(Color.red) ;
             offsGg.fillPolygon(exes,whys,4) ;
           }
         }
 // draw cowl
         if(prob == 3 || prob == 5) {
           exes[0] = xorgn + (int) (.5*fact*cwlx) + xt ;
           whys[0] = yorgn - (int) (.5*fact*cwly) + yt ;
           exes[1] = xorgn + (int) (.5*fact*(cwlx + xlong)) + xt ;
           whys[1] = yorgn - (int) (.5*fact*(cwly + (cwlx + xlong)*Math.tan(convdr*5.0))) + yt ;
           exes[2] = exes[1] ;
           whys[2] = whys[0] ;
           offsGg.setColor(Color.red) ;
           Polygon poly = new Polygon(exes,whys,3) ;
           offsGg.fillPolygon(poly) ;
         }
 // draw shock waves
         for(i=1; i <= nshocks; ++i) {
            exes[0] = xorgn + (int) (.5*fact*sxbgn[i]) + xt ;
            whys[0] = yorgn - (int) (.5*fact*sybgn[i]) + yt ;
            exes[1] = xorgn + (int) (.5*fact*sxnd[i]) + xt ;
            whys[1] = yorgn - (int) (.5*fact*synd[i]) + yt ;
            if (detach[i]) {
              offsGg.setColor(Color.magenta) ; 
              offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
            }
            else {
              offsGg.setColor(Color.blue) ; 
              offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;  
            }
            if(i==2) {
              offsGg.setColor(Color.red) ; 
              if (prob == 3) {
                offsGg.drawString("2", exes[1]+20, whys[1]-10) ;
              }
            }
            if(i==3) {
              offsGg.setColor(Color.red) ; 
              offsGg.drawString("3", exes[0]+20, whys[0]-10) ;
            }
            if(i==4) {
              offsGg.setColor(Color.red) ; 
              offsGg.drawString("4", exes[1]+20, whys[1]-10) ;
            }
          }
       }

//   centered prandtl-meyer expansion
       if(prob == 4) {
           exes[0] = xorgn + (int) (.5*fact*(wxbgn[1] - 100.)) + xt ;
           whys[0] = yorgn - (int) (.5*fact*(wybgn[1])) + yt ;
           exes[1] = xorgn + (int) (.5*fact*(wxbgn[1])) + xt ;
           whys[1] = whys[0] ;
           exes[2] = exes[1] ;
           whys[2] = yorgn - (int) (.5*fact*(xlong*wslope[1] + winter[1])) + yt  ;
           exes[3] = exes[0] ;
           whys[3] = whys[2] ;
           offsGg.setColor(Color.red) ;
           Polygon poly = new Polygon(exes,whys,4) ;
           offsGg.fillPolygon(poly) ;
           exes[0] = xorgn + (int) (.5*fact*wxbgn[1]) + xt ;
           whys[0] = yorgn - (int) (.5*fact*(wybgn[1])) + yt ;
           exes[1] = xorgn + (int) (.5*fact*xlong) + xt ;
           whys[1] = yorgn - (int) (.5*fact*(xlong*wslope[1] + winter[1])) + yt ;
           exes[2] = exes[0] ;
           whys[2] = whys[1] ;
           offsGg.setColor(Color.red) ;
           poly = new Polygon(exes,whys,3) ;
           offsGg.fillPolygon(poly) ;
           offsGg.drawString("1", (exes[0]+exes[1])/2, (whys[0]+whys[1])/2 - 10) ;

// draw expansion fan
          exes[0] = xorgn + (int) (.5*fact* sxbgn[1]) + xt ;
          whys[0] = yorgn + (int) (.5*fact*-sybgn[1]) + yt ;
          exes[1] = xorgn + (int) (.5*fact* sxnd[1]) + xt ;
          whys[1] = yorgn + (int) (.5*fact*-synd[1]) + yt ;
          exes[2] = xorgn + (int) (.5*fact* exnd[1]) + xt ;
          whys[2] = yorgn + (int) (.5*fact*-eynd[1]) + yt ;
          offsGg.setColor(Color.black) ; 
          offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;  
          offsGg.drawLine(exes[0],whys[0],exes[2],whys[2]) ;  
       }

 // isentropic compression inlet
       if (prob == 5) {
 // draw ramps
         for (i = 1; i <= nramps; ++i) {
           exes[0] = xorgn + (int) (.5*fact*wxbgn[i]) + xt ;
           whys[0] = yorgn - (int) (.5*fact*wybgn[i]) + yt ;
           exes[1] = xorgn + (int) (.5*fact*wxnd[i]) + xt ;
           whys[1] = yorgn - (int) (.5*fact*wynd[i]) + yt ;
           exes[2] = exes[1] ;
           whys[2] = whys[0] ;
           offsGg.setColor(Color.red) ;
           Polygon poly = new Polygon(exes,whys,3) ;
           offsGg.fillPolygon(poly) ;
           if (i > 1) {
               whys[1] = whys[0] ;
               whys[2] = yorgn + yt ;
               whys[3] = whys[2] ;
               exes[3] = exes[0] ;
               offsGg.fillPolygon(exes,whys,4) ;
           }
           if (i == 1) offsGg.drawString("1", exes[0]+20, whys[0]-10) ;
           if (i == nramps) offsGg.drawString(String.valueOf(nramps-1), exes[0]+20, whys[0]-10) ;
           if (i == nramps) {
             exes[0] = xorgn + (int) (.5*fact*wxnd[i]) + xt ;
             whys[0] = yorgn - (int) (.5*fact*wynd[i]) + yt ;
             exes[1] = xorgn + (int) (.5*fact*xlong) + xt ;
             whys[1] = yorgn - 40 + yt;
             exes[2] = exes[1] ;
             whys[2] = yorgn + yt ;
             exes[3] = exes[0] ;
             whys[3] = whys[2] ;
             offsGg.setColor(Color.red) ;
             offsGg.fillPolygon(exes,whys,4) ;
           }
         }
 // draw cowl
         exes[0] = xorgn + (int) (.5*fact*cwlx) + xt ;
         whys[0] = yorgn - (int) (.5*fact*cwly) + yt ;
         exes[1] = xorgn + (int) (.5*fact*(cwlx + xlong)) + xt ;
         whys[1] = yorgn - (int) (.5*fact*(cwly + (cwlx + xlong)*Math.tan(convdr*5.0))) + yt ;
         exes[2] = exes[1] ;
         whys[2] = whys[0] ;
         offsGg.setColor(Color.red) ;
         Polygon poly = new Polygon(exes,whys,3) ;
         offsGg.fillPolygon(poly) ;
 // draw shock waves
         for(i=1; i <= nshocks; ++i) {
            exes[0] = xorgn + (int) (.5*fact*sxbgn[i]) + xt ;
            whys[0] = yorgn - (int) (.5*fact*sybgn[i]) + yt ;
            exes[1] = xorgn + (int) (.5*fact*sxnd[i]) + xt ;
            whys[1] = yorgn - (int) (.5*fact*synd[i]) + yt ;
            offsGg.setColor(Color.red) ;
            if (i == nshocks) offsGg.drawString(String.valueOf(nramps), exes[1]+20, whys[1]-30) ;
            if (detach[i]) {
              offsGg.setColor(Color.magenta) ; 
              offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
            }
            else {
              offsGg.setColor(Color.blue) ; 
              offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;  
            }
          }
 //   draw mach lines for problem 5
         for (i = 2; i <= numray+1; ++i) {
            exes[0] = xorgn + (int) (.5*fact*(exbgn[i])) + xt ;
            whys[0] = yorgn - (int) (.5*fact*(eybgn[i])) + yt ;
            exes[1] = xorgn + (int) (.5*fact*(exnd[i])) + xt ;
            whys[1] = yorgn - (int) (.5*fact*(eynd[i])) + yt ;
            offsGg.setColor(Color.black) ; 
            offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
         }
         if(outzn >= 2 && outzn <= numray +1) { 
           exes[0] = xorgn + (int) (.5*fact*(exbgn[outzn])) + xt ;
           whys[0] = yorgn - (int) (.5*fact*(eybgn[outzn])) + yt ;
           exes[1] = xorgn + (int) (.5*fact*(exnd[outzn])) + xt ;
           whys[1] = yorgn - (int) (.5*fact*(eynd[outzn])) + yt ;
           offsGg.setColor(Color.red) ; 
           offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ; 
         }
       }

       if (prob == 6) {  // over or under expanded jet
 // draw geom
           for(i=1; i<=2; ++i) {
             exes[0] = xorgn + (int) (.5*fact*nzht*wxbgn[i]) + xt ;
             whys[0] = yorgn - (int) (.5*fact*nzht*wybgn[i]) + yt ;
             exes[1] = xorgn + (int) (.5*fact*nzht*wxnd[i]) + xt ;
             whys[1] = yorgn - (int) (.5*fact*nzht*wynd[i]) + yt ;
             exes[2] = exes[1]-5 ;
             whys[2] = whys[1]  ;
             exes[3] = exes[0]-5 ;
             whys[3] = whys[0]  ;
             offsGg.setColor(Color.red) ;
             offsGg.fillPolygon(exes,whys,4) ;
           }
           for(i=3; i<=4; ++i) {
             exes[0] = xorgn + (int) (.5*fact*nzht*wxbgn[i]) + xt ;
             whys[0] = yorgn - (int) (.5*fact*nzht*wybgn[i]) + yt ;
             exes[1] = xorgn + (int) (.5*fact*nzht*wxnd[i]) + xt ;
             whys[1] = yorgn - (int) (.5*fact*nzht*wynd[i]) + yt ;
             exes[2] = exes[1] ;
             whys[2] = whys[1] + 5 ;
             exes[3] = exes[0] ;
             whys[3] = whys[0] + 5 ;
             offsGg.setColor(Color.red) ;
             offsGg.fillPolygon(exes,whys,4) ;
           }
 // draw flow
           if(nprat > 1.0) {  // under expanded
  // waves
             for (j=1 ; j<=ncycle; ++j) { 
               for (i = 1; i <= 8; ++i) {
                 exes[0] = xorgn + (int) (.5*fact*(exbgn[i]
                             +(j-1)*(strx[4][1]-strx[1][1])) ) + xt ;   
                 whys[0] = yorgn - (int) (.5*fact*eybgn[i]) + yt ;
                 exes[1] = xorgn + (int) (.5*fact*(exnd[i]
                            +(j-1)*(strx[4][1]-strx[1][1])) ) + xt ;
                 whys[1] = yorgn - (int) (.5*fact*eynd[i]) + yt ;
                 if(efamily[i] == 1) offsGg.setColor(Color.blue) ; 
                 if(efamily[i] == 2) offsGg.setColor(Color.red) ; 
                 offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
               }
             }
 // streamlines
             offsGg.setColor(Color.black) ; 
             exes[3] = xorgn + (int) (.5*fact*strx[1][3]) + xt ;   
             whys[3] = yorgn - (int) (.5*fact*.5 *(stry[1][3]+stry[1][1])) + yt ;
             offsGg.drawString("1", exes[3]-10, whys[3]) ;

             for (j=1 ; j<=ncycle; ++j) {
               for (i = 1; i <= 3; ++i) {
                 exes[0] = xorgn + (int) (.5*fact*(strx[1][i] 
                             +(j-1)*(strx[4][1]-strx[1][1])) ) + xt ;   
                 whys[0] = yorgn - (int) (.5*fact*stry[1][i]) + yt ;
                 exes[1] = xorgn + (int) (.5*fact*(strx[2][i]
                             +(j-1)*(strx[4][1]-strx[1][1])) ) + xt ;
                 whys[1] = yorgn - (int) (.5*fact*stry[2][i]) + yt ;
                 offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                 exes[4] = (int) (.5*(exes[0] + exes[1])) ;
                 whys[4] = whys[0]+10 ;
                 if(i==3) {
                    offsGg.drawString("2", exes[4], whys[4]) ;
                 }   

                 exes[0] = xorgn + (int) (.5*fact*(strx[2][i]
                            +(j-1)*(strx[4][1]-strx[1][1])) ) + xt ;   
                 whys[0] = yorgn - (int) (.5*fact*stry[2][i]) + yt ;
                 exes[1] = xorgn + (int) (.5*fact*(strx[3][i]
                            +(j-1)*(strx[4][1]-strx[1][1])) ) + xt ;
                 whys[1] = yorgn - (int) (.5*fact*stry[3][i]) + yt ;
                 offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                 exes[3] = (int) (.5*(exes[0] + exes[1])) ;
                 if(i==3) {
                    offsGg.drawString("3", exes[3], whys[3]) ;
                 }   

                 exes[0] = xorgn + (int) (.5*fact*(strx[3][i]
                           +(j-1)*(strx[4][1]-strx[1][1])) ) + xt ;   
                 whys[0] = yorgn - (int) (.5*fact*stry[3][i]) + yt ;
                 exes[1] = xorgn + (int) (.5*fact*(strx[4][i]
                           +(j-1)*(strx[4][1]-strx[1][1])) ) + xt ;
                 whys[1] = yorgn - (int) (.5*fact*stry[4][i]) + yt ;
                 offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                 exes[4] = (int) (.5*(exes[0] + exes[1])) ;
                 exes[3] = exes[1] - 10 ;
                 if(i==3) {
                    offsGg.drawString("4", exes[4], whys[4]) ;
                    offsGg.drawString("5", exes[3], whys[3]) ;
                 }   
               }
             }
           }

           if (nprat <= 1.0) { // over expanded
  // waves
               exes[0] = xorgn + (int) (.5*fact*exbgn[1]) + xt ;   
               whys[0] = yorgn - (int) (.5*fact*eybgn[1]) + yt ;
               exes[1] = xorgn + (int) (.5*fact*exnd[1]) + xt ;
               whys[1] = yorgn - (int) (.5*fact*eynd[1]) + yt ;
               if(efamily[1] == 1) offsGg.setColor(Color.blue) ; 
               if(efamily[1] == 2) offsGg.setColor(Color.red) ; 
               offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
               exes[0] = xorgn + (int) (.5*fact*exbgn[2]) + xt ;   
               whys[0] = yorgn - (int) (.5*fact*eybgn[2]) + yt ;
               exes[1] = xorgn + (int) (.5*fact*exnd[2]) + xt ;
               whys[1] = yorgn - (int) (.5*fact*eynd[2]) + yt ;
               if(efamily[2] == 1) offsGg.setColor(Color.blue) ; 
               if(efamily[2] == 2) offsGg.setColor(Color.red) ; 
               offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
            for (j=1 ; j<=ncycle; ++j) {
             for (i = 3; i<= 10; ++i) {
               exes[0] = xorgn + (int) (.5*fact*(exbgn[i]+(j-1)*(strx[5][1]-strx[2][1]))) + xt ;   
               whys[0] = yorgn - (int) (.5*fact*eybgn[i]) + yt ;
               exes[1] = xorgn + (int) (.5*fact*(exnd[i]+(j-1)*(strx[5][1]-strx[2][1]))) + xt ;
               whys[1] = yorgn - (int) (.5*fact*eynd[i]) + yt ;
               if(efamily[i] == 1) offsGg.setColor(Color.blue) ; 
               if(efamily[i] == 2) offsGg.setColor(Color.red) ; 
               offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             }
            }
 // streamlines
             offsGg.setColor(Color.black) ; 
             exes[0] = xorgn + (int) (.5*fact*strx[1][1]) + xt ;   
             whys[0] = yorgn - (int) (.5*fact*stry[1][1]) + yt ;
             exes[1] = xorgn + (int) (.5*fact*strx[2][1]) + xt ;
             whys[1] = yorgn - (int) (.5*fact*stry[2][1]) + yt ;
             offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[3] = xorgn + (int) (.5*fact*strx[1][3]) + xt ;   
             whys[3] = yorgn - (int) (.5*fact*.5 *(stry[1][3]+stry[1][1])) + yt ;
             offsGg.drawString("1", exes[3]-10, whys[3]) ;

             exes[0] = xorgn + (int) (.5*fact*strx[1][2]) + xt ;   
             whys[0] = yorgn - (int) (.5*fact*stry[1][2]) + yt ;
             exes[1] = xorgn + (int) (.5*fact*strx[2][2]) + xt ;
             whys[1] = yorgn - (int) (.5*fact*stry[2][2]) + yt ;
             offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = xorgn + (int) (.5*fact*strx[1][3]) + xt ;   
             whys[0] = yorgn - (int) (.5*fact*stry[1][3]) + yt ;
             exes[1] = xorgn + (int) (.5*fact*strx[2][3]) + xt ;
             whys[1] = yorgn - (int) (.5*fact*stry[2][3]) + yt ;
             offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[4] = (int) (.5*(exes[0] + exes[1])) ;
             whys[4] = whys[1]+10 ;
             offsGg.drawString("2", exes[4], whys[4]) ;
             exes[3] = exes[1] - 10 ;
             offsGg.drawString("3", exes[3], whys[3]) ;

            for (j=1 ; j<=ncycle; ++j) {
             for (i = 1; i <= 3; ++i) {
               exes[0] = xorgn + (int) (.5*fact*(strx[2][i] +(j-1)*(strx[5][1]-strx[2][1]))) + xt ;
               whys[0] = yorgn - (int) (.5*fact*stry[2][i]) + yt ;    
               exes[1] = xorgn + (int) (.5*fact*(strx[3][i]+(j-1)*(strx[5][1]-strx[2][1]))) + xt ;
               whys[1] = yorgn - (int) (.5*fact*stry[3][i]) + yt ;
               offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
               exes[4] = (int) (.5*(exes[0] + exes[1])) ;
               if(i==3) {
                  offsGg.drawString("4", exes[4], whys[4]) ;
               }   

               exes[0] = xorgn + (int) (.5*fact*(strx[3][i]+(j-1)*(strx[5][1]-strx[2][1]))) + xt ;   
               whys[0] = yorgn - (int) (.5*fact*stry[3][i]) + yt ;
               exes[1] = xorgn + (int) (.5*fact*(strx[4][i]+(j-1)*(strx[5][1]-strx[2][1]))) + xt ;
               whys[1] = yorgn - (int) (.5*fact*stry[4][i]) + yt ;
               offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
               exes[3] = (int) (.5*(exes[0] + exes[1])) ;
               if(i==3) {
                  offsGg.drawString("5", exes[3], whys[3]) ;
               }   

               exes[0] = xorgn + (int) (.5*fact*(strx[4][i]+(j-1)*(strx[5][1]-strx[2][1]))) + xt ;   
               whys[0] = yorgn - (int) (.5*fact*stry[4][i]) + yt ;
               exes[1] = xorgn + (int) (.5*fact*(strx[5][i]+(j-1)*(strx[5][1]-strx[2][1]))) + xt ;
               whys[1] = yorgn - (int) (.5*fact*stry[5][i]) + yt ;
               offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
               exes[4] = (int) (.5*(exes[0] + exes[1])) ;
               exes[3] = exes[1] - 10 ;
               if(i==3) {
                  offsGg.drawString("6", exes[4], whys[4]) ;
                  offsGg.drawString("7", exes[3], whys[3]) ;
               }   
             }
            }
          }
       }

       if (prob == 7) {  // moc nozzle by planes
 // draw geom
           for(i=1; i<=numray; ++i) {
             exes[0] = xorgn + (int) (.5*fact*nzht*wxbgn[i]) + xt ;
             whys[0] = yorgn - (int) (.5*fact*nzht*wybgn[i]) + yt ;
             exes[1] = xorgn + (int) (.5*fact*nzht*wxnd[i]) + xt ;
             whys[1] = yorgn - (int) (.5*fact*nzht*wynd[i]) + yt ;
             exes[2] = exes[1] ;
             whys[2] = whys[1] - 5 ;
             exes[3] = exes[0] ;
             whys[3] = whys[0] - 5 ;
             offsGg.setColor(Color.red) ;
             offsGg.fillPolygon(exes,whys,4) ;
           }
           exes[0] = xorgn + (int) (.5*fact*0.0) + xt ;   
           whys[0] = yorgn - (int) (.5*fact*0.0) + yt ;
           exes[1] = xorgn + (int) (.5*fact*500.0) + xt ;
           whys[1] = yorgn - (int) (.5*fact*0.0) + yt ;
           offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

           exes[0] = xorgn + (int) (.5*fact*nzht*mcxll[1][1]) + xt ;
           whys[0] = yorgn - (int) (.5*fact*nzht*mcyll[1][1]) + yt ;
           exes[1] = xorgn + (int) (.5*fact*nzht*mcxul[1][1]) + xt ;
           whys[1] = yorgn - (int) (.5*fact*nzht*mcyul[1][1]) + yt ;
           exes[2] = xorgn + (int) (.5*fact*nzht*mcxur[1][1]) + xt ;
           whys[2] = yorgn - (int) (.5*fact*nzht*mcyur[1][1]) + yt ;
           exes[3] = xorgn + (int) (.5*fact*nzht*mcxlr[1][1]) + xt ;
           whys[3] = yorgn - (int) (.5*fact*nzht*mcylr[1][1]) + yt ;
           offsGg.setColor(Color.blue) ;
           offsGg.drawPolygon(exes,whys,4) ;

           for(i=2; i<=numray/2+1 ; ++i) {
             exes[0] = xorgn + (int) (.5*fact*nzht*mcxll[1][i]) + xt ;
             whys[0] = yorgn - (int) (.5*fact*nzht*mcyll[1][i]) + yt ;
             exes[1] = xorgn + (int) (.5*fact*nzht*mcxul[1][i]) + xt ;
             whys[1] = yorgn - (int) (.5*fact*nzht*mcyul[1][i]) + yt ;
             exes[2] = xorgn + (int) (.5*fact*nzht*mcxur[1][i]) + xt ;
             whys[2] = yorgn - (int) (.5*fact*nzht*mcyur[1][i]) + yt ;
             exes[3] = xorgn + (int) (.5*fact*nzht*mcxlr[1][i]) + xt ;
             whys[3] = yorgn - (int) (.5*fact*nzht*mcylr[1][i]) + yt ;
             offsGg.setColor(Color.blue) ;
             offsGg.drawPolygon(exes,whys,4) ;
           }

           exes[0] = xorgn + (int) (.5*fact*nzht*mcxll[2][2]) + xt ;
           whys[0] = yorgn - (int) (.5*fact*nzht*mcyll[2][2]) + yt ;
           exes[1] = xorgn + (int) (.5*fact*nzht*mcxul[2][2]) + xt ;
           whys[1] = yorgn - (int) (.5*fact*nzht*mcyul[2][2]) + yt ;
           exes[2] = xorgn + (int) (.5*fact*nzht*mcxur[2][2]) + xt ;
           whys[2] = yorgn - (int) (.5*fact*nzht*mcyur[2][2]) + yt ;
           exes[3] = xorgn + (int) (.5*fact*nzht*mcxlr[2][2]) + xt ;
           whys[3] = yorgn - (int) (.5*fact*nzht*mcylr[2][2]) + yt ;
           offsGg.setColor(Color.blue) ;
           offsGg.drawPolygon(exes,whys,4) ;

           for (k=3; k<= numray/2+1; ++k) {
              for(i=k; i<=numray/2+1 ; ++i) {
                exes[0] = xorgn + (int) (.5*fact*nzht*mcxll[k-1][i]) + xt ;
                whys[0] = yorgn - (int) (.5*fact*nzht*mcyll[k-1][i]) + yt ;
                exes[1] = xorgn + (int) (.5*fact*nzht*mcxul[k-1][i]) + xt ;
                whys[1] = yorgn - (int) (.5*fact*nzht*mcyul[k-1][i]) + yt ;
                exes[2] = xorgn + (int) (.5*fact*nzht*mcxur[k-1][i]) + xt ;
                whys[2] = yorgn - (int) (.5*fact*nzht*mcyur[k-1][i]) + yt ;
                exes[3] = xorgn + (int) (.5*fact*nzht*mcxlr[k-1][i]) + xt ;
                whys[3] = yorgn - (int) (.5*fact*nzht*mcylr[k-1][i]) + yt ;
                offsGg.setColor(Color.blue) ;
                offsGg.drawPolygon(exes,whys,4) ;
              }
              exes[0] = xorgn + (int) (.5*fact*nzht*mcxll[k][k]) + xt ;
              whys[0] = yorgn - (int) (.5*fact*nzht*mcyll[k][k]) + yt ;
              exes[1] = xorgn + (int) (.5*fact*nzht*mcxul[k][k]) + xt ;
              whys[1] = yorgn - (int) (.5*fact*nzht*mcyul[k][k]) + yt ;
              exes[2] = xorgn + (int) (.5*fact*nzht*mcxur[k][k]) + xt ;
              whys[2] = yorgn - (int) (.5*fact*nzht*mcyur[k][k]) + yt ;
              exes[3] = xorgn + (int) (.5*fact*nzht*mcxlr[k][k]) + xt ;
              whys[3] = yorgn - (int) (.5*fact*nzht*mcylr[k][k]) + yt ;
              offsGg.setColor(Color.blue) ;
              offsGg.drawPolygon(exes,whys,4) ;
          }
// chosen zone
          exes[0] = xorgn + (int) (.5*fact*nzht*mcxll[row][col]) + xt ;
          whys[0] = yorgn - (int) (.5*fact*nzht*mcyll[row][col]) + yt ;
          exes[1] = xorgn + (int) (.5*fact*nzht*mcxul[row][col]) + xt ;
          whys[1] = yorgn - (int) (.5*fact*nzht*mcyul[row][col]) + yt ;
          exes[2] = xorgn + (int) (.5*fact*nzht*mcxur[row][col]) + xt ;
          whys[2] = yorgn - (int) (.5*fact*nzht*mcyur[row][col]) + yt ;
          exes[3] = xorgn + (int) (.5*fact*nzht*mcxlr[row][col]) + xt ;
          whys[3] = yorgn - (int) (.5*fact*nzht*mcylr[row][col]) + yt ;
          offsGg.setColor(Color.red) ;
          offsGg.drawPolygon(exes,whys,4) ;
       }

       if (prob == 8) {  // moc nozzle by points
 // draw geom
           for(i=1; i<=numray+1; ++i) {
             exes[0] = xorgn + (int) (.5*fact*nzht*wxbgn[i]) + xt ;
             whys[0] = yorgn - (int) (.5*fact*nzht*wybgn[i]) + yt ;
             exes[1] = xorgn + (int) (.5*fact*nzht*wxnd[i]) + xt ;
             whys[1] = yorgn - (int) (.5*fact*nzht*wynd[i]) + yt ;
             exes[2] = exes[1] ;
             whys[2] = whys[1] - 5 ;
             exes[3] = exes[0] ;
             whys[3] = whys[0] - 5 ;
             offsGg.setColor(Color.red) ;
             offsGg.fillPolygon(exes,whys,4) ;
           }
           exes[0] = xorgn + (int) (.5*fact*0.0) + xt ;   
           whys[0] = yorgn - (int) (.5*fact*0.0) + yt ;
           exes[1] = xorgn + (int) (.5*fact*500.0) + xt ;
           whys[1] = yorgn - (int) (.5*fact*0.0) + yt ;
           offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

           offsGg.setColor(Color.blue) ;
           for(i=0; i<=numray/2; ++i) {
              exes[0] = xorgn + (int) (.5*fact*nzht*mcx[i][1]) + xt ;
              whys[0] = yorgn - (int) (.5*fact*nzht*mcy[i][1]) + yt ;
              exes[1] = xorgn + (int) (.5*fact*nzht*mcx[i+1][1]) + xt ;
              whys[1] = yorgn - (int) (.5*fact*nzht*mcy[i+1][1]) + yt ;
              offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
           }
           k = numray/2 + 1 ;
           for(i=1; i<=numray/2; ++i) {
              exes[0] = xorgn + (int) (.5*fact*nzht*mcx[k][i]) + xt ;
              whys[0] = yorgn - (int) (.5*fact*nzht*mcy[k][i]) + yt ;
              exes[1] = xorgn + (int) (.5*fact*nzht*mcx[k][i+1]) + xt ;
              whys[1] = yorgn - (int) (.5*fact*nzht*mcy[k][i+1]) + yt ;
              offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
           }
           for(i=1; i<=numray/2; ++i) {
              for (k=1; k<=i; ++k) {
                 exes[0] = xorgn + (int) (.5*fact*nzht*mcx[i][k]) + xt ;
                 whys[0] = yorgn - (int) (.5*fact*nzht*mcy[i][k]) + yt ;
                 exes[1] = xorgn + (int) (.5*fact*nzht*mcx[i][k+1]) + xt ;
                 whys[1] = yorgn - (int) (.5*fact*nzht*mcy[i][k+1]) + yt ;
                 offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
              }
           }
           for(i=2; i<=numray/2+1; ++i) {
              for (k=i-1; k<=numray/2; ++k) {
                 exes[0] = xorgn + (int) (.5*fact*nzht*mcx[k][i]) + xt ;
                 whys[0] = yorgn - (int) (.5*fact*nzht*mcy[k][i]) + yt ;
                 exes[1] = xorgn + (int) (.5*fact*nzht*mcx[k+1][i]) + xt ;
                 whys[1] = yorgn - (int) (.5*fact*nzht*mcy[k+1][i]) + yt ;
                 offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
              }
           }

// chosen point
          exes[0] = xorgn + (int) (.5*fact*nzht*mcx[row][col]) + xt ;
          whys[0] = yorgn - (int) (.5*fact*nzht*mcy[row][col]) + yt ;
          offsGg.setColor(Color.red) ;
          offsGg.drawString("o",exes[0]-3,whys[0]+5) ;
       }

       if (prob == 9) {  // axi moc nozzle by points
 // draw geom
           for(i=1; i<=numray+1; ++i) {
             exes[0] = xorgn + (int) (.5*fact*nzht*wxbgn[i]) + xt ;
             whys[0] = yorgn - (int) (.5*fact*nzht*wybgn[i]) + yt ;
             exes[1] = xorgn + (int) (.5*fact*nzht*wxnd[i]) + xt ;
             whys[1] = yorgn - (int) (.5*fact*nzht*wynd[i]) + yt ;
             exes[2] = exes[1] ;
             whys[2] = whys[1] - 5 ;
             exes[3] = exes[0] ;
             whys[3] = whys[0] - 5 ;
             offsGg.setColor(Color.red) ;
             offsGg.fillPolygon(exes,whys,4) ;
           }

           exes[0] = xorgn + (int) (.5*fact*0.0) + xt ;   
           whys[0] = yorgn - (int) (.5*fact*0.0) + yt ;
           exes[1] = xorgn + (int) (.5*fact*500.0) + xt ;
           whys[1] = yorgn - (int) (.5*fact*0.0) + yt ;
           offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

           offsGg.setColor(Color.blue) ;
           for(i=0; i<=numray/2; ++i) {
              exes[0] = xorgn + (int) (.5*fact*nzht*mcx[i][1]) + xt ;
              whys[0] = yorgn - (int) (.5*fact*nzht*mcy[i][1]) + yt ;
              exes[1] = xorgn + (int) (.5*fact*nzht*mcx[i+1][1]) + xt ;
              whys[1] = yorgn - (int) (.5*fact*nzht*mcy[i+1][1]) + yt ;
              offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
           }
           k = numray/2 + 1 ;
           for(i=1; i<=numray/2; ++i) {
              exes[0] = xorgn + (int) (.5*fact*nzht*mcx[k][i]) + xt ;
              whys[0] = yorgn - (int) (.5*fact*nzht*mcy[k][i]) + yt ;
              exes[1] = xorgn + (int) (.5*fact*nzht*mcx[k][i+1]) + xt ;
              whys[1] = yorgn - (int) (.5*fact*nzht*mcy[k][i+1]) + yt ;
              offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
           }
           for(i=1; i<=numray/2; ++i) {
              for (k=1; k<=i; ++k) {
                 exes[0] = xorgn + (int) (.5*fact*nzht*mcx[i][k]) + xt ;
                 whys[0] = yorgn - (int) (.5*fact*nzht*mcy[i][k]) + yt ;
                 exes[1] = xorgn + (int) (.5*fact*nzht*mcx[i][k+1]) + xt ;
                 whys[1] = yorgn - (int) (.5*fact*nzht*mcy[i][k+1]) + yt ;
                 offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
              }
           }
           for(i=2; i<=numray/2+1; ++i) {
              for (k=i-1; k<=numray/2; ++k) {
                 exes[0] = xorgn + (int) (.5*fact*nzht*mcx[k][i]) + xt ;
                 whys[0] = yorgn - (int) (.5*fact*nzht*mcy[k][i]) + yt ;
                 exes[1] = xorgn + (int) (.5*fact*nzht*mcx[k+1][i]) + xt ;
                 whys[1] = yorgn - (int) (.5*fact*nzht*mcy[k+1][i]) + yt ;
                 offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
              }
           }

// chosen point
          exes[0] = xorgn + (int) (.5*fact*nzht*mcx[row][col]) + xt ;
          whys[0] = yorgn - (int) (.5*fact*nzht*mcy[row][col]) + yt ;
          offsGg.setColor(Color.red) ;
          offsGg.drawString("o",exes[0]-3,whys[0]+5) ;
       }

 // draw rays for prob = 1
       if(drawray == 1) {
          if(prob == 1 && detach[1] != true) {
             offsGg.setColor(Color.blue) ; 
             exes[0] = xorgn + (int) (.5*fact*sxbgn[1]) + xt ;
             whys[0] = yorgn - (int) (.5*fact*sybgn[1]) + yt ;
             for(i=1; i<=numray; ++ i) {
                exes[1] = xorgn + (int) (.5*fact*(sxbgn[1] + xlong*Math.cos(rthet[i]))) + xt ;
                whys[1] = yorgn - (int) (.5*fact*(sybgn[1] + xlong*Math.sin(rthet[i]))) + yt ;
                offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             }
             exes[1] = xorgn + (int) (.5*fact*(sxbgn[1] + xlong*Math.cos(rthet[outzn]))) + xt ;
             whys[1] = yorgn - (int) (.5*fact*(sybgn[1] + xlong*Math.sin(rthet[outzn]))) + yt ;
             offsGg.setColor(Color.red) ; 
             offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
          }
       }
          // draw streamlines
    if(prob == 0 || detach[1] == true) {
       offsGg.setColor(Color.black) ; 
        for (i = 1; i <= 4; ++i) {
          exes[0] = (int) (.5*fact*strx[1][i]) + xt + xorgn ;
          whys[0] = (int) (.5*fact*stry[1][i]) + yt + yorgn ;
          for(k=2; k<=4; ++k) {
            exes[1] = (int) (.5*fact*strx[k][i]) + xt + xorgn ;
            whys[1] = (int) (.5*fact*stry[k][i]) + yt + yorgn ;
            offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
            exes[0] = exes[1] ;
            whys[0] = whys[1] ;
          }
        }
     }

     if(prob == 1 && detach[1] != true) {
        offsGg.setColor(Color.black) ; 
        for (i = 1; i <= 4; ++i) {
          exes[0] = (int) (.5*fact*strx[1][i]) + xt + xorgn ;
          whys[0] = (int) (.5*fact*stry[1][i]) + yt + yorgn ;
          for(k=2; k<=12; ++k) {
            exes[1] = (int) (.5*fact*strx[k][i]) + xt + xorgn ;
            whys[1] = (int) (.5*fact*stry[k][i]) + yt + yorgn ;
            offsGg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
            exes[0] = exes[1] ;
            whys[0] = whys[1] ;
          }
        }
      }
 // zoom widget
        offsGg.setColor(Color.white) ;
        offsGg.fillRect(0,0,30,275) ;
        offsGg.setColor(Color.black) ;
        offsGg.drawLine(15,20,15,240) ;
        offsGg.fillRect(5,sldloc,20,5) ;
 // bottom labels
        offsGg.setColor(Color.white) ;
        offsGg.fillRect(0,245,900,30) ;

        offsGg.setColor(Color.black) ;
        offsGg.drawString("Zoom",2,260) ;
        offsGg.drawRect(40,245,15,15) ;
        offsGg.drawString("-",45,257) ;
        offsGg.drawRect(60,245,15,15) ;
        offsGg.drawString("+",65,257) ;

        offsGg.setColor(Color.blue) ;
        offsGg.drawString("Input-Upstream",200,260) ;
        offsGg.setColor(Color.red) ;
        offsGg.drawString("Supersonic Flows - Version 1.6g",400,260) ;
        offsGg.setColor(Color.blue) ;
        offsGg.drawString("Output-Downstream",700,260) ;

 // buttons
        offsGg.setColor(Color.blue) ;
        offsGg.fillRect(845,245,50,15) ;
        offsGg.setColor(Color.white) ;
        offsGg.drawString("Find",850,256) ;

       g.drawImage(offscreenImg,0,0,this) ;
    }
  }
}
