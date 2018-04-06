/*
                   Wright 1901 Wind Tunnel Simulation
                        angle of attack input
                  special graphics for tunnel lift balance
                    made from FoilSim II  - Airfoil  mode
   
                           A Java Applet
               to perform Kutta-Joukowski Airfoil analysis

                     Version 1.4a   - 5 Sep 02

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
 
  New test -
            * add all of the Wright airfoil models
            * add the drag balance
            * calibrate against Wright data
            * smooth transitions on graphics
 
                                           TJB  9 Sep 02

*/
package Tunl;

import java.awt.*;
import java.lang.Math ;

public class Foil extends java.applet.Applet {
 
   static double convdr = 3.1415926/180. ;
   static double pid2 = 3.1415926/2.0 ;
   static double rval,ycval,xcval,gamval,alfval,thkval,camval,chrd,clift ;
   static double thkinpt,caminpt ;                 /* MODS 10 Sep 99 */
   static double leg,teg,lem,tem;
   static double usq,vsq,alt,altmax,area,armax,armin ;
   static double chord,span,aspr,arold,chrdold,spnold ; /* Mod 13 Jan 00 */
   static double q0,ps0,pt0,ts0,rho,rlhum,temf,presm ;
   static double lyg,lrg,lthg,lxgt,lygt,lrgt,lthgt;/* MOD 20 Jul */
   static double lxm,lym,lxmt,lymt,vxdir;/* MOD 20 Jul */
   static double deltb,xflow ;             /* MODS  20 Jul 99 */
   static double delx,delt,vfsd,spin,spindr,yoff ;
   static double vel,pres,lift,side,omega,radcrv,relsy,angr ;
   static double weight,aread,sped ;
   static double angout,glidang,angdraw,angdraw2 ;

   static double rg[][]  = new double[25][50] ; 
   static double thg[][] = new double[25][50] ; 
   static double xg[][]  = new double[25][50] ; 
   static double yg[][]  = new double[25][50] ; 
   static double xm[][]  = new double[25][50] ; 
   static double ym[][]  = new double[25][50] ; 
   static double plp[]   = new double[50] ;
   static double plv[]   = new double[50] ;

   int nptc,npt2,nlnc,nln2,rdflag,browflag,probflag,anflag;
   int foil,flflag,lunits,lftout,planet,plane ;
   int conflag,displ,dispp,antim,ancol;   /* MODS  2 APR 99 - 22 APR -27 JUL */
   int varflg ;
   int stepper,balance,group,model,ysl,pegged,trans1,trans2 ;
       /* units data */
   static double vmn,almn,angmn,vmx,almx,angmx ;
   static double camn,thkmn,camx,thkmx ;
   static double chrdmn,spnmn,armn,chrdmx,spnmx,armx ;
   static double vconv,vmaxa,vmaxb ;
   static double pconv,pmax,pmin,lconv,fconv,fmax,fmaxb;
   int lflag,gflag,plscale,nond;
       /*  plot & probe data */
   static double fact,xpval,ypval,pbval,factp;
   static double prg,pthg,pxg,pyg,pxm,pym ;
   int pboflag,xt,yt,ntikx,ntiky,npt,xtp,ytp ;
   int lines,nord,nabs,ntr ;
   static double begx,endx,begy,endy ;
   static String labx,labxu,laby,labyu ;
   static double pltx[][]  = new double[3][50] ;
   static double plty[][]  = new double[3][50] ;

   Solver solve ;
   Sq sq ;
   Lq lq ;
   CardLayout layin ;
   Image offImg1 ;
   Graphics off1Gg ;
   Image offImg2 ;
   Graphics off2Gg ;
   Image offImg3 ;
   Graphics off3Gg ;
   Image offImg4 ;
   Graphics off4Gg ;
   Image offImg5 ;
   Graphics off5Gg ;
   Image partimg ;

   public void init() {
     int i;
     Foil a = new Foil() ;
     solve = new Solver() ;

     offImg1 = createImage(this.size().width,
                      this.size().height) ;
     off1Gg = offImg1.getGraphics() ;
     offImg2 = createImage(this.size().width,
                      this.size().height) ;
     off2Gg = offImg2.getGraphics() ;
     offImg3 = createImage(this.size().width,
                      this.size().height) ;
     off3Gg = offImg3.getGraphics() ;
     offImg4 = createImage(this.size().width,
                      this.size().height) ;
     off4Gg = offImg4.getGraphics() ;
     offImg5 = createImage(this.size().width,
                      this.size().height) ;
     off5Gg = offImg5.getGraphics() ;

     setLayout(new GridLayout(2,1,5,5)) ;

     solve.setDefaults () ;
 
     sq = new Sq(this) ;
     lq = new Lq(this) ;

     add(lq) ;
     add(sq) ;

     solve.getFreeStream ();
     computeFlow () ;
     partimg = getImage(getCodeBase(),"tunnel.gif");
     lq.view.start() ;
     sq.inl.start() ;
     sq.inc.start() ;
     sq.inr.repaint() ;
  }
 
  public Insets insets() {
     return new Insets(5,5,5,5) ;
  }

  public void computeFlow() { 

     if (flflag == 1) {
         solve.getCirc ();                   /* get circulation */
         solve.genFlow () ;
         solve.getFreeStream () ;
     }

     solve.getProbe() ;
 
     loadOut() ;
     lq.ob.inright.l.repaint() ;
     sq.inl.repaint() ;
     sq.inc.repaint() ;
  }

  public int filter0(double inumbr) {
        //  output only to .
       int number ;
       int intermed ;
 
       number = (int) (inumbr);
       return number ;
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
 
  public void setUnits() {   // Switching Units
       double ovs,chords,spans,aros,chos,spos ;
       double alts,ares,ards,spds ;

       alts = alt / lconv ;
       chords = chord / lconv ;
       spans = span / lconv ;
       ares = area /lconv/lconv ;
       aros = arold /lconv/lconv ;
       chos = chrdold / lconv ;
       spos = spnold / lconv ;
       ovs = vfsd / vconv ;
       ards = aread /lconv/lconv ;
       spds = sped / vconv ;

       switch (lunits) {
          case 0: {                             /* English */
            lconv = 1.;                      /*  feet    */
            vconv = .6818; vmaxa = 50.; vmaxb = 10. ;  /*  mph  */
            fconv = 1.0; fmax = 100000.; fmaxb = .5;  /* pounds   */
            pconv = 14.7  ;                   /* lb/sq in */
            break;
          }
          case 1: {                             /* Metric */
            lconv = .3048;                    /* meters */
            vconv = 1.097;  vmaxa = 80.; vmaxb = 17.;   /* km/hr  */
            fconv = 4.448 ; fmax = 500000.; fmaxb = 2.5; /* newtons */
            pconv = 101.3 ;               /* kilo-pascals */
            break ;
          }
       }
 
       alt = alts * lconv ;
       chord = chords * lconv ;
       span = spans * lconv ;
       area = ares * lconv * lconv ;
       arold = aros * lconv * lconv ;
       chrdold = chos * lconv ;
       spnold = spos * lconv ;
       vfsd  = ovs * vconv;
       aread = ards * lconv * lconv ;
       sped = spds * vconv ;

       return ;
  }

  public void loadOut() {   // output routine
     double stfact,cdrag,corfact ;
     double af,bf,cf ;
     double y0,y1,y2 ;
     Double V1 ;
     double v1 ;
     float fl1 ;
//  determine output angle for dial
//     lift balance
     if (balance == 1) {
                          // stall model
        stfact = 1.0 ;
        pegged = 0 ;
        if (anflag == 1) {
            if (alfval > 10.0 ) {
               stfact = .5 + .1 * alfval - .005 * alfval * alfval ;
               if (model == 1) stfact = 1.5 ;
               if (model == 2) stfact = 1.01 - .01 * alfval  ;
               if (model == 3) {
                    if (alfval <= 16.0) stfact = 1.1 - .015 * alfval  ;
                    if (alfval > 16.0) stfact = .86 - .05 * (alfval-16.)  ;
               }
               if (model == 4) {
                    if (alfval < 15.0) stfact = 1. + .2 * alfval ;
                    if (alfval >= 15.0) stfact = 1. + .35 * alfval ;
               }
               if (model == 5) {
                    if (alfval <= 16.0) stfact = 1. + .1 * alfval ;
                    if (alfval > 16.0) stfact = 1. + .12 * alfval ;
               }
               if (model == 6) {
                    if (alfval <= 16.0) stfact = 1. + .1 * alfval ;
                    if (alfval > 16.0) stfact = 1. + .12 * alfval ;
               }
               if (model == 7) {
                    if (alfval <= 14.0) stfact = .8 ;
                    if ((alfval > 14.0) && (alfval <=17) ) stfact = .8 - (alfval - 14.0)*.07 ;
                    if (alfval > 17.0) stfact = .6 - (alfval - 17.0)*.04 ;
               }
               if (model == 8) {
                    if (alfval <= 14.0) stfact = .8 ;
                    if ((alfval > 14.0) && (alfval <=17) ) stfact = .8 - (alfval - 14.0)*.07 ;
                    if (alfval > 17.0) stfact = .6 - (alfval - 17.0)*.04 ;
               }
               if (model == 9) {
                    if (alfval < 14.0) stfact = .85 ;
                    if ((alfval >= 14.0) && (alfval <=17) ) stfact = .75 - (alfval - 14.0)*.07 ;
                    if (alfval > 17.0) stfact = .55 - (alfval - 17.0)*.04 ;
               }
               if (model == 10) stfact = .95 - (alfval-10.0) * .02;
               if (model == 11) stfact = .95 - (alfval-10.0) * .02;
               if (model == 12) stfact = .95 - (alfval-10.0) * .02;
               if (model == 15) {
                    if (alfval < 15.0) stfact = 1. + .02 * alfval ;
                    if (alfval >= 15.0) stfact = 1.5 + .2 * (alfval - 15.0) ;
               }
               if (model == 16) {
                    if (alfval < 15.0) stfact = 1. + .02 * alfval ;
                    if (alfval >= 15.0) stfact = 1.5 + .2 * (alfval - 15.0) ;
               }
               if (model == 17) {
                    if (alfval < 15.0) stfact = 1. + .02 * alfval ;
                    if (alfval >= 15.0) stfact = 1.5 + .18 * (alfval - 15.0) ;
               }
               if (model == 18) {
                    if (alfval <= 15.0) stfact = 1.0 ;
                    if (alfval > 15.0) stfact = 1.0 - (alfval - 15.0)*.04 ;
               }
               if (model == 19) {
                    if (alfval < 15.0) stfact = 1.0 ;
                    if (alfval >= 15.0) stfact = .95 - (alfval - 15.0)*.045 ;
               }
               if (model == 20) {
                    if (alfval < 14.0) stfact = 1.0 ;
                    if (alfval >= 14.0) stfact = 1.0 + (alfval - 14.0)*.1 ;
               }
               if (model == 21) {
                    if (alfval < 15.0) stfact = .95 - (alfval - 10.0)*.02  ;
                    if (alfval >= 15.0) stfact = .8 - (alfval - 15.0)*.05 ;
               }
               if (model == 23) {
                    if (alfval < 15.0) stfact = .95 - (alfval-10.0) * .02;
                    if (alfval >= 15.0) stfact = .8 - (alfval-15.0) * .05;
               }
               if (model == 24) {
                    if (alfval < 15.0) stfact = .90 - (alfval-10.0) * .03;
                    if (alfval >= 15.0) stfact = .77 - (alfval-15.0) * .05;
               }
               if (model == 25) stfact = .95 - (alfval - 10.0)*.02  ;
               if (model == 27) {
                    if (alfval < 15.0) stfact = .95 - (alfval-10.0) * .03;
                    if (alfval >= 15.0) stfact = .77 - (alfval-15.0) * .04;
               }
               if (model == 30) stfact = 1.0 ;
               if (model == 31) {
                     if (alfval < 15.0) stfact =  1.0 - (alfval - 10.0)*.04 ;
                     if (alfval >= 15.0) stfact =  .8 - (alfval - 15.0)*.03 ;
               }
               if (model == 33) stfact = .95 - (alfval - 10.0)*.025  ;
               if (model == 35) stfact = .9 - (alfval - 10.0)*.04;
               if (model == 40) {
                    if (alfval <= 13.0) stfact = 1.0 - (alfval - 10.0)*.02 ;
               }
               if (model == 41) {
                    if (alfval <= 13.0) stfact = 1.0 - (alfval - 10.0)*.02 ;
               }
               if (model == 43) {
                    if (alfval <= 16.0) stfact = 1.0 - (alfval -10.0)*.02  ;
                    if (alfval > 16.0) stfact = .88 - (alfval - 16.0)*.03 ;
               }
               if (model == 51) stfact = .9 - (alfval - 10.0)*.04;
            }
            clift = clift * stfact ;
        }
   //     correction for aspect ratio
   
        clift = clift /(1.0 + clift/(3.14159*aspr)) ;
        cdrag = 1.28 ;
 
        corfact = .5 ;
        if(model == 1) corfact = .3 ;
        if(model == 2) corfact = .5 ;
        if(model == 3) corfact = .53 ;
        if(model == 4) { 
              if (alfval <= 10.0) corfact = .3 + alfval * .006 ;
              if (alfval  > 10.0) corfact = .25 + (alfval - 10.0) * .004 ;
        }
        if(model == 5){
              if (alfval <= 10.0) corfact = .3 + alfval * .006 ;
              if (alfval  > 10.0) corfact = .25 + (alfval - 10.0) * .008 ;
        } 
        if(model == 6) {
              if (alfval <= 10.0) corfact = .3 + alfval * .006 ;
              if (alfval  > 10.0) corfact = .25 + (alfval - 10.0) * .008 ;
        }
        if(model == 7) corfact = .20 + .03*alfval ;
        if(model == 8) corfact = .22 + .03*alfval ;
        if(model == 9) corfact = .25 + .03*alfval ;
        if(model == 10) {
               if (alfval <= 6.0)  corfact = .22 + .025*alfval ;
               if (alfval > 6.0)  corfact = .37  ;
        }
        if(model == 11) {
               if (alfval <= 6.0)  corfact = .26 + .025*alfval ;
               if (alfval > 6.0)  corfact = .41  ;
        }
        if(model == 12) {
               if (alfval <= 6.0)  corfact = .28 + .025*alfval ;
               if (alfval > 6.0 && alfval < 15.0)  corfact = .43 + .005*alfval  ;
               if (alfval > 15.0)  corfact = .5 - (alfval - 15.)*.01  ;
        }
        if(model == 15) {
               if (alfval <= 10.0)  corfact = .15 + .02*alfval ;
               if (alfval > 10.0)  corfact = .35  ;
        }
        if(model == 16) {
               if (alfval <= 10.0)  corfact = .15 + .02*alfval ;
               if (alfval > 10.0)  corfact = .35  ;
        }
        if(model == 17) {
               if (alfval <= 10.0)  corfact = .15 + .02*alfval ;
               if (alfval > 10.0)  corfact = .35  ;
        }
        if(model == 18) { 
              if (alfval <= 10.0) corfact = .3 + alfval * .015 ;
              if (alfval  > 10.0) corfact = .45 + (alfval - 10.0) * .005 ;
        }
        if(model == 19) { 
              if (alfval <= 10.0) corfact = .3 + alfval * .015 ;
              if (alfval  > 10.0) corfact = .45 + (alfval - 10.0) * .005 ;
        }
        if(model == 20) { 
              if (alfval < 5.0) corfact = .2 + alfval * .015 ;
              if (alfval >= 5.0) corfact = .25 - (alfval - 5.0) * .007 ;
        }
        if(model == 21) { 
              if (alfval <= 8.0) corfact = .2 + alfval * .04 ;
              if (alfval  > 8.0) corfact = .5 + (alfval - 8.0) * .005 ;
        }
        if(model == 23) {
               if (alfval <= 10.0)  corfact = .10 + .008*alfval ;
               if (alfval > 10.0)   corfact = .2 + .005*(alfval-10.0)  ;
        }
        if(model == 24) {
               if (alfval <= 10.0)  corfact = .15 + .02*alfval ;
               if (alfval > 10.0)   corfact = .4 + .01*(alfval-10.0)  ;
        }
        if(model == 25) { 
              if (alfval < 5.0) corfact = .2 + alfval * .04 ;
              if (alfval >= 5.0) corfact = .35 - (alfval - 5.0) * .01 ;
        }
        if(model == 27) {
               if (alfval < 10.0)  corfact = .3 + .025*alfval ;
               if (alfval >= 10.0)   corfact = .52 + .02*(alfval-10.0)  ;
        }
        if(model == 30) corfact = .07 + alfval * .004 ; 
        if(model == 31) {
              if (alfval <= 6.0) corfact = .35 + alfval * .04 ; 
              if (alfval > 6.0) corfact = .60 + (alfval - 6.0) * .003 ; 
        }
        if(model == 33) { 
              if (alfval < 5.0) corfact = .4 + alfval * .04 ;
              if (alfval >= 5.0) corfact = .6 - (alfval - 5.0) * .01 ;
        }
        if(model == 35) {
               if (alfval <= 5.0)  corfact = .2 + .04*alfval ;
               if (alfval > 5.0)  corfact = .4 + .01*(alfval - 5.0)  ;
        }
        if(model == 40) { 
              if (alfval <= 10.0) corfact = .55 + alfval * .015 ;
              if (alfval  > 10.0) corfact = .7 + (alfval - 10.0) * .005 ;
              if (alfval > 13.0) pegged = 1 ;
        }
        if(model == 41) { 
              if (alfval <= 5.0) corfact = .55 + alfval * .04 ;
              if (alfval > 5.0) corfact = .55 + alfval * .025 ;
              if (alfval > 10.0) pegged = 1 ;
        }
        if(model == 43) { 
              if (alfval <= 10.0) corfact = .48 + alfval * .01 ;
              if (alfval  > 10.0) corfact = .6 + (alfval - 10.0) * .005 ;
        }
        if(model == 51) {
               if (alfval <= 6.0)  corfact = .2 + .02*alfval ;
               if (alfval > 6.0)  corfact = .3 + .01*(alfval - 6.0)  ;
        }

        clift = corfact * clift ;
 
        if (stepper == 3) {
           angout = .9 * Math.asin(clift/cdrag)/convdr ;
           if (pegged == 1) angout = 90.0 ;
        }
        if (stepper == 4) {
           angout =  Math.asin(clift/cdrag)/convdr ;
           if (pegged == 1) angout = 90.0 ;
        }
     }
//    drag balance -  curve fits
     if (balance == 2) {
        y0 = 19.5 ;
        y1 = 2.75 ;
        y2 = -2.0 ;
        if(model == 4) { 
           y0 = 19.5 ;  y1 = 2.75;  y2 = -2.0 ;
        }
        if(model == 5) { 
           y0 = 16.0 ;  y1 = 1.0;  y2 = -2.5 ;
        }
        if(model == 6) { 
           y0 = 13.5 ;  y1 = 0.75;  y2 = -2.0 ;
        }
        if(model == 7) { 
           y0 = 18.0 ;  y1 = -1.6;  y2 = -5.0 ;
        }
        if(model == 8) { 
           y0 = 12.75 ;  y1 = -3.0;  y2 = -3.75 ;
        }
        if(model == 9) { 
           y0 = 11.0 ;  y1 = -3.75;  y2 = -3.0 ;
        }
        if(model == 10) { 
           y0 = 17.5 ;  y1 = 0.0;  y2 = -2.25 ;
        }
        if(model == 11) { 
           y0 = 15.5 ;  y1 = -1.825;  y2 = -2.75 ;
        }
        if(model == 12) { 
           y0 = 14.75 ;  y1 = -3.25;  y2 = -2.25 ;
        }
        if(model == 15) { 
           y0 = 19.75 ;  y1 = 3.0;  y2 = -3.125 ;
        }
        if(model == 16) { 
           y0 = 19.75 ;  y1 = 2.0;  y2 = -4.0 ;
        }
        if(model == 17) { 
           y0 = 20.0 ;  y1 = 1.125;  y2 = -4.5 ;
        }
        if(model == 18) { 
           y0 = 12.0 ;  y1 = -2.75;  y2 = -2.125 ;
        }
        if(model == 19) { 
           y0 = 15.0 ;  y1 = -2.825;  y2 = -3.25 ;
        }
        if(model == 20) { 
           y0 = 13.0 ;  y1 = 9.0;  y2 = 6.78 ;
        }
        if(model == 21) { 
           y0 = 20.5 ;  y1 = -3.25;  y2 = -1.25 ;
        }
        if(model == 23) { 
           y0 = 20. ;  y1 = -0.125;  y2 = -1.0 ;
        }
        if(model == 24) { 
           y0 = 20. ;  y1 = 0.5;  y2 = -1.75 ;
        }
        if(model == 25) { 
           y0 = 15.5 ;  y1 = 2.75;  y2 = 4.25 ;
        }
        if(model == 27) { 
           y0 = 23.0 ;  y1 = 1.25;  y2 = -1.5 ;
        }
        if(model == 30) { 
           y0 = 38.0 ;  y1 = 18.0;  y2 = 9.5 ;
        }
        if(model == 31) { 
           y0 = 11.0 ;  y1 = -2.3;  y2 = -4.5 ;
        }
        if(model == 33) { 
           y0 = 21.0 ;  y1 = 2.75;  y2 = 3.0 ;
        }
        if(model == 35) { 
           y0 = 20.5 ;  y1 = -1.0;  y2 = 1.125 ;
        }
        if(model == 40) { 
           y0 = 15.0 ;  y1 = -1.0;  y2 = -4.25 ;
        }
        if(model == 41) { 
           y0 = 15.0 ;  y1 = -1.125;  y2 = -3.5 ;
        }
        if(model == 43) { 
           y0 = 16.75 ;  y1 = -0.5;  y2 = -5.0 ;
        }
        if(model == 51) { 
           y0 = 14.5 ;  y1 = 3.0;  y2 = -0.125 ;
        }

        cf = y0;
        bf = (4.0*y1 - y2 - 3.0*y0)/20.0 ;
        af = (y2 - 2.0*y1 + y0)/200.0 ;
        
//        v1 = bf ;
//        fl1 = (float) v1 ;
//        lq.ob.inleft.o1.setText(String.valueOf(fl1)) ;
//        v1 = af ;
//        fl1 = (float) v1 ;
//        lq.ob.inleft.o2.setText(String.valueOf(fl1)) ;

        if (stepper >= 3) {
           angout = af * alfval * alfval + bf * alfval + cf ;
              //  special for flat plate
           stfact = 1.0 ;
           if (model < 4) {
              if (anflag == 1) {
                  if (alfval > 10.0 ) {
                     if (model == 1) stfact = 1.5 ;
                     if (model == 2) stfact = 1.01 - .01 * alfval  ;
                     if (model == 3) {
                          if (alfval <= 16.0) stfact = 1.1 - .015 * alfval  ;
                          if (alfval > 16.0) stfact = .86 - .05 * (alfval-16.)  ;
                     }
                  }
                  clift = clift * stfact ;
              }
         //     correction for aspect ratio
   
              clift = clift /(1.0 + clift/(3.14159*aspr)) ;
 
              corfact = .5 ;
              if(model == 1) corfact = .3 ;
              if(model == 2) corfact = .5 ;
              if(model == 3) corfact = .53 ;
              clift = corfact * clift ;
//              if (clift < .001) clift = .001;
              angout = 90.0 - Math.atan(1.0 / clift)/convdr - alfval ;
           }
        }
     }

     return ;
  }
 
  public void loadProbe() {   // probe output routine

     pbval = 0.0 ;
     if (pboflag == 1) pbval = vel * vfsd ;           // velocity
     if (pboflag == 2) pbval = ((ps0 + pres * q0)/2116.) * pconv ; // pressure

     return ;
  }

  class Solver {
 
     Solver () {
     }

     public void setDefaults() {

        trans2 = 0 ;
        trans1 = 0 ;
        pegged = 0 ;
        balance = 1;
        group = 1 ;
        model = 1 ;
        stepper = 1 ;
        ysl = 80 ;
        varflg = 0 ;
        planet = 0 ;
        lunits = 0 ;
        lftout = 0 ;
        nlnc = 21 ;
        nln2 = nlnc/2 + 1 ;
        nptc = 45 ;
        npt2 = nptc/2 + 1 ;
        deltb = .5 ;
        foil = 3 ;
        flflag = 1;
        thkval = .040 ;
        thkinpt = 1.0 ;                   /* MODS 10 SEP 99 */
        camval = 0.0 ;
        caminpt = 0.0 ;
        alfval = 0.0 ;
        gamval = 0.0 ;
        spin = 0.0 ;
        spindr = 1.0 ;
        rval = 1.0 ;
        ycval = 0.0 ;
        xcval = 0.0 ;
        conflag = 1 ;                             /* MODS  2 Apr 99 */
        displ = 0 ;                              /* MODS  22 Apr 99 */
        dispp = 2 ;
        lift = 53. ;
        weight = 50. ;
        aread = 6. ;
        sped  = 15. ;
        plane = 0 ;
 
        xpval = 2.1;
        ypval = -.5 ;
        pboflag = 0 ;
        xflow = -12.0;                             /* MODS  20 Jul 99 */

        pconv = 14.7;
        pmin = .5 ;
        pmax = 1.0 ;
        fconv = 1.0 ;
        fmax = 100000. ;
        fmaxb = .50 ;
        vconv = .6818 ;
        vfsd = 0.0 ;
        vmaxa = 50. ;
        vmaxb = 10. ;
        lconv = 1.0 ;

        alt = 0.0 ;
        altmax = 50000. ;
        chrdold = chord = 5.0 ;
        spnold = span = 17.0 ;
        aspr = 6.0 ;
        arold = area = 6.0 ;
        armax = 510.01 ;
        armin = .01 ;                 /* MODS 9 SEP 99 */
 
        xt = 130;  yt = 80; fact = 15.0 ;
        xtp = 95; ytp = 130; factp = 25.0 ;
 
        probflag = 2 ;
        anflag = 1 ;
        vmn = 0.0;     vmx = 50.0 ;
        almn = 0.0;    almx = 50000.0 ;
        angmn = -5.0; angmx = 20.0 ;
        camn = -25.0;  camx = 25.0 ;
        thkmn = 1.0; thkmx = 26.0 ;
        chrdmn = .1 ;  chrdmx = 10.1 ;
        spnmn = .1 ;  spnmx = 100.1 ;
        armn = .01 ;  armx = 510.01 ;

        angout = 0.0 ;
        return ;
     }

     public void getFreeStream() {    //  free stream conditions
       double hite,pvap,rgas,gama ;       /* MODS  19 Jan 00  whole routine*/

       rgas = 1718. ;                /* ft2/sec2 R */
       gama = 1.4 ;
       hite = alt/lconv ;
       if (planet == 0) {    // Earth  standard day
         if (conflag == 1) {
           if (hite <= 36152.) {           // Troposphere
              ts0 = 518.6 - 3.56 * hite/1000. ;
              ps0 = 2116. * Math.pow(ts0/518.6,5.256) ;
           }
           if (hite >= 36152. && hite <= 82345.) {   // Stratosphere
              ts0 = 389.98 ;
              ps0 = 2116. * .2236 *
                 Math.exp((36000.-hite)/(53.35*389.98)) ;
           }
           if (hite >= 82345.) {
              ts0 = 389.98 + 1.645 * (hite-82345)/1000. ;
              ps0 = 2116. *.02456 * Math.pow(ts0/389.98,-11.388) ;
           }
           rlhum = 0.0 ;
           temf = ts0 - 459.6 ;
           if (temf <= 0.0) temf = 0.0 ;                    
           presm = ps0 * 29.92 / 2116. ;
         }
         if (conflag == 2) {
            ts0 = temf + 459.6 ;
            if (temf < 0.0) {
                  temf = 0.0 ;
                  rlhum = 0.0 ;
            }
             ps0 = presm * 2116. / 29.92 ;
         }
         pvap = rlhum*(2.685+.00353*Math.pow(temf,2.245));/* Eq 1:6A  Domasch */
         rho = (ps0 - .379*pvap)/(rgas * ts0) ;  /* effect of humidty */
         rho = ps0/(53.3 * 32.17 * ts0) ;
       }

       if (planet == 1) {   // Mars - curve fit of orbiter data
         rgas = 1149. ;                /* ft2/sec2 R */
         gama = 1.29 ;

         if (hite <= 22960.) {
            ts0 = 434.02 - .548 * hite/1000. ;
            ps0 = 14.62 * Math.pow(2.71828,-.00003 * hite) ;
         }
         if (hite > 22960.) {
            ts0 = 449.36 - 1.217 * hite/1000. ;
            ps0 = 14.62 * Math.pow(2.71828,-.00003 * hite) ;
         }
         rho = ps0/(rgas*ts0) ;
       }

       q0  = .5 * rho * vfsd * vfsd / (vconv * vconv) ;
       pt0 = ps0 + q0 ;

       return ;
     }

     public void getCirc() {   // circulation from Kutta condition
       double thet,rdm,thtm ;
       double beta,rball;
       int index;

       xcval = 0.0 ;
       switch (foil)  {
          case 0: {         /* get circulation from spin for baseball */
              rball = .1 ;         /* baseball radius = .1 ft = 1.2 in */
              gamval = 4.0 * 3.1415926 * 3.1415926 *spin * rball * rball
                                 / (vfsd/vconv) ;
              gamval = gamval * spindr ;
              ycval = .0001 ;
              break ;
          }
          case 1:  {                  /* Juokowski geometry*/
              ycval = camval / 2.0 ;
              rval = thkval/4.0 +Math.sqrt(thkval*thkval/16.0+ycval*ycval +1.0);
              xcval = 1.0 - Math.sqrt(rval*rval - ycval*ycval) ;
              beta = Math.asin(ycval/rval)/convdr ;     /* Kutta condition */
              gamval = 2.0*rval*Math.sin((alfval+beta)*convdr) ;
              break ;
          }
          case 2:  {                  /* Elliptical geometry*/
              ycval = camval / 2.0 ;
              rval = thkval/4.0 + Math.sqrt(thkval*thkval/16.0+ycval*ycval+1.0);
              beta = Math.asin(ycval/rval)/convdr ;    /* Kutta condition */
              gamval = 2.0*rval*Math.sin((alfval+beta)*convdr) ;
              break ;
          }
          case 3:  {                  /* Plate geometry*/
              ycval = camval / 2.0 ;
              rval = Math.sqrt(ycval*ycval+1.0);
              beta = Math.asin(ycval/rval)/convdr ;    /* Kutta condition */
              gamval = 2.0*rval*Math.sin((alfval+beta)*convdr) ;
              break ;
          }
       }
                                                   /* geometry */
       for (index =1; index <= nptc; ++index) {
           thet = (index -1)*360./(nptc-1) ;
           xg[0][index] = rval * Math.cos(convdr * thet) + xcval ;
           yg[0][index] = rval * Math.sin(convdr * thet) + ycval ;
           rg[0][index] = Math.sqrt(xg[0][index]*xg[0][index] +
                                    yg[0][index]*yg[0][index])  ;
           thg[0][index] = Math.atan2(yg[0][index],xg[0][index])/convdr;
           xm[0][index] = (rg[0][index] + 1.0/rg[0][index])*
                        Math.cos(convdr*thg[0][index]) ;
           ym[0][index] = (rg[0][index] - 1.0/rg[0][index])*
                        Math.sin(convdr*thg[0][index]) ;
           rdm = Math.sqrt(xm[0][index]*xm[0][index] +
                           ym[0][index]*ym[0][index])  ;
           thtm = Math.atan2(ym[0][index],xm[0][index])/convdr;
           xm[0][index] = rdm * Math.cos((thtm - alfval)*convdr);
           ym[0][index] = rdm * Math.sin((thtm - alfval)*convdr);
           getVel(rval,thet) ;
           plp[index] = ((ps0 + pres * q0)/2116.) * pconv ;
           plv[index] = vel * vfsd ;
       }

       return ;
     }

     public void genFlow() {   // generate flowfield
       double rnew,thet,psv,fxg;
       int k,index;
                              /* all lines of flow  except stagnation line*/
       for (k=1; k<=nlnc; ++k) {
         psv = -.5*(nln2-1) + .5*(k-1) ;
         fxg = xflow ;
         for (index =1; index <=nptc; ++ index) {
           solve.getPoints (fxg,psv) ;
           xg[k][index]  = lxgt ;
           yg[k][index]  = lygt ;
           rg[k][index]  = lrgt ;
           thg[k][index] = lthgt ;
           xm[k][index]  = lxmt ;
           ym[k][index]  = lymt ;
           if (anflag == 1) {           // stall model
              if (alfval > 10.0 && psv > 0.0) {
                   if (xm[k][index] > 0.0) {
                      ym[k][index] = ym[k][index -1] ;
                   }
              }
              if (alfval < -10.0 && psv < 0.0) {
                   if (xm[k][index] > 0.0) {
                      ym[k][index] = ym[k][index -1] ;
                   }
              }
           }
           solve.getVel(lrg,lthg) ;
           fxg = fxg + vxdir*deltb ;
         }
       }
                                              /*  stagnation line */
       k = nln2 ;
       psv = 0.0 ;
                                              /*  incoming flow */
       for (index =1; index <= npt2; ++ index) {
           rnew = 10.0 - (10.0 - rval)*Math.sin(pid2*(index-1)/(npt2-1)) ;
           thet = Math.asin(.999*(psv - gamval*Math.log(rnew/rval))/
                                   (rnew - rval*rval/rnew)) ;
           fxg =  - rnew * Math.cos(thet) ;
           solve.getPoints (fxg,psv) ;
           xg[k][index]  = lxgt ;
           yg[k][index]  = lygt ;
           rg[k][index]  = lrgt ;
           thg[k][index] = lthgt ;
           xm[k][index]  = lxmt ;
           ym[k][index]  = lymt ;
       }
                                              /*  downstream flow */
       for (index = 1; index <= npt2; ++ index) {
           rnew = 10.0 + .01 - (10.0 - rval)*Math.cos(pid2*(index-1)/(npt2-1)) ;
           thet = Math.asin(.999*(psv - gamval*Math.log(rnew/rval))/
                                      (rnew - rval*rval/rnew)) ;
           fxg =   rnew * Math.cos(thet) ;
           solve.getPoints (fxg,psv) ;
           xg[k][npt2+index]  = lxgt ;
           yg[k][npt2+index]  = lygt ;
           rg[k][npt2+index]  = lrgt ;
           thg[k][npt2+index] = lthgt ;
           xm[k][npt2+index]  = lxmt ;
           ym[k][npt2+index]  = lymt ;
       }
                                              /*  stagnation point */
       xg[k][npt2]  = xcval ;
       yg[k][npt2]  = ycval ;
       rg[k][npt2]  = Math.sqrt(xcval*xcval+ycval*ycval) ;
       thg[k][npt2] = Math.atan2(ycval,xcval)/convdr ;
       xm[k][npt2]  = (xm[k][npt2+1] + xm[k][npt2-1])/2.0 ;
       ym[k][npt2]  = (ym[0][nptc/4+1] + ym[0][nptc/4*3+1])/2.0 ;
                                /*  compute lift coefficient */
       leg = xcval - Math.sqrt(rval*rval - ycval*ycval) ;
       teg = xcval + Math.sqrt(rval*rval - ycval*ycval) ;
       lem = leg + 1.0/leg ;
       tem = teg + 1.0/teg ;
       chrd = tem - lem ;
       clift = gamval*4.0*3.1415926/chrd ;

       return ;
     }

     public void getPoints(double fxg, double psv) {   // flow in x-psi
       double radm,thetm ;                /* MODS  20 Jul 99  whole routine*/
       double fnew,ynew,yold,rfac,deriv ;
       double xold,xnew,thet ;
       double rmin,rmax ;
       int iter,isign;
                       /* get variables in the generating plane */
                           /* iterate to find value of yg */
       ynew = 10.0 ;
       yold = 10.0 ;
       if (psv < 0.0) ynew = -10.0 ;
       if (Math.abs(psv) < .001 && alfval < 0.0) ynew = rval ;
       if (Math.abs(psv) < .001 && alfval >= 0.0) ynew = -rval ;
       fnew = 0.1 ;
       iter = 1 ;
       while (Math.abs(fnew) >= .00001 && iter < 25) {
           ++iter ;
           rfac = fxg*fxg + ynew*ynew ;
           if (rfac < rval*rval) rfac = rval*rval + .01 ;
           fnew = psv - ynew*(1.0 - rval*rval/rfac)
                  - gamval*Math.log(Math.sqrt(rfac)/rval) ;
           deriv = - (1.0 - rval*rval/rfac)
               - 2.0 * ynew*ynew*rval*rval/(rfac*rfac)
               - gamval * ynew / rfac ;
           yold = ynew ;
           ynew = yold  - .5*fnew/deriv ;
       }
       lyg = yold ;
                                     /* rotate for angle of attack */
       lrg = Math.sqrt(fxg*fxg + lyg*lyg) ;
       lthg = Math.atan2(lyg,fxg)/convdr ;
       lxgt = lrg * Math.cos(convdr*(lthg + alfval)) ;
       lygt = lrg * Math.sin(convdr*(lthg + alfval)) ;
                              /* translate cylinder to generate airfoil */
       lxgt = lxgt + xcval ;
       lygt = lygt + ycval ;
       lrgt = Math.sqrt(lxgt*lxgt + lygt*lygt) ;
       lthgt = Math.atan2(lygt,lxgt)/convdr ;
                               /*  Kutta-Joukowski mapping */
       lxm = (lrgt + 1.0/lrgt)*Math.cos(convdr*lthgt) ;
       lym = (lrgt - 1.0/lrgt)*Math.sin(convdr*lthgt) ;
                              /* tranforms for view fixed with free stream */
                /* take out rotation for angle of attack mapped and cylinder */
       radm = Math.sqrt(lxm*lxm+lym*lym) ;
       thetm = Math.atan2(lym,lxm)/convdr ;
       lxmt = radm*Math.cos(convdr*(thetm-alfval)) ;
       lymt = radm*Math.sin(convdr*(thetm-alfval)) ;

       lxgt = lxgt - xcval ;
       lygt = lygt - ycval ;
       lrgt = Math.sqrt(lxgt*lxgt + lygt*lygt)  ;
       lthgt = Math.atan2(lygt,lxgt)/convdr;
       lxgt = lrgt * Math.cos((lthgt - alfval)*convdr);
       lygt = lrgt * Math.sin((lthgt - alfval)*convdr);

       return ;
     }
 
     public void getVel(double radius, double theta) {  //velocity and pressure 
      double ur,uth,jake1,jake2,jakesq ;
      double xloc,yloc,thrad,alfrad ;

      thrad = convdr * theta ;
      alfrad = convdr * alfval ;
                                /* get x, y location in cylinder plane */
      xloc = radius * Math.cos(thrad) ;
      yloc = radius * Math.sin(thrad) ;
                                /* velocity in cylinder plane */
      ur  = Math.cos(thrad-alfrad)*(1.0-(rval*rval)/(radius*radius)) ;
      uth = -Math.sin(thrad-alfrad)*(1.0+(rval*rval)/(radius*radius))
                            - gamval/radius;
      usq = ur*ur + uth*uth ;
      vxdir = ur * Math.cos(thrad) - uth * Math.sin(thrad) ; // MODS  20 Jul 99 
                                /* translate to generate airfoil  */
      xloc = xloc + xcval ;
      yloc = yloc + ycval ;
                                   /* compute new radius-theta  */
      radius = Math.sqrt(xloc*xloc + yloc*yloc) ;
      thrad  = Math.atan2(yloc,xloc) ;
                                   /* compute Joukowski Jacobian  */
      jake1 = 1.0 - Math.cos(2.0*thrad)/(radius*radius) ;
      jake2 = Math.sin(2.0*thrad)/(radius*radius) ;
      jakesq = jake1*jake1 + jake2*jake2 ;
      if (Math.abs(jakesq) <= .01) jakesq = .01 ;  /* protection */
      vsq = usq / jakesq ;
          /* vel is velocity ratio - pres is coefficient  (p-p0)/q0   */
      if (foil > 0) {
           vel = Math.sqrt(vsq) ;
           pres = 1.0 - vsq ;
      }
      if (foil == 0) {
           vel = Math.sqrt(usq) ;
           pres = 1.0 - usq ;
      }
      return ;
    }

    public void getProbe () { /* all of the information needed for the probe */
      double prxg;
      int index;
                       /* get variables in the generating plane */
      if (Math.abs(ypval) < .01) ypval = .05 ;
      solve.getPoints (xpval,ypval) ;

      solve.getVel(lrg,lthg) ;
      loadProbe() ;

      pxg = lxgt ;
      pyg = lygt ;
      prg = lrgt ;
      pthg = lthgt ;
      pxm = lxmt ;
      pym = lymt ;
                                    /* smoke */
      if (pboflag == 3 ) {
        prxg = xpval ;
        for (index =1; index <=nptc; ++ index) {
          solve.getPoints (prxg,ypval) ;
          xg[19][index] = lxgt ;
          yg[19][index] = lygt ;
          rg[19][index] = lrgt ;
          thg[19][index] = lthgt ;
          xm[19][index] = lxmt ;
          ym[19][index] = lymt ;
          if (anflag == 1) {           // stall model
             if (xpval > 0.0) {
                if (alfval > 10.0 && ypval > 0.0) {
                   ym[19][index] = ym[19][1] ;
                }
                if (alfval < -10.0 && ypval < 0.0) {
                     ym[19][index] = ym[19][1] ;
                }
             }
             if (xpval < 0.0) {
                if (alfval > 10.0 && ypval > 0.0) {
                   if (xm[19][index] > 0.0) {
                       ym[19][index] = ym[19][index-1] ;
                   }
                }
                if (alfval < -10.0 && ypval < 0.0) {
                   if (xm[19][index] > 0.0) {
                     ym[19][index] = ym[19][index-1] ;
                   }
                }
             }
          }
          solve.getVel(lrg,lthg) ;
          prxg = prxg + vxdir*deltb ;
        }
      }
      return ;
    }
  }
 
  class Lq extends Panel {
     Foil outerparent ;
     Viewer view ;
     Ob ob ;

     Lq (Foil target) { 
         outerparent = target ;
         setLayout(new GridLayout(1,2,5,5)) ;

         view  = new Viewer(outerparent) ;
         ob = new Ob(outerparent) ;

         add(view) ;
         add(ob) ;
     }

     class Ob extends Panel {
        Foil outerparent ;
        Inleft inleft ;
        Inright inright ;

        Ob (Foil target) {

           outerparent = target ;
           setLayout(new GridLayout(1,2,5,5)) ;

           inleft = new Inleft(outerparent) ;
           inright = new Inright(outerparent) ;

           add(inright) ;
           add(inleft) ;
        }

        class Inleft extends Panel {
           Foil outerparent ;
           Ipt ipt ;
           Label l1,l2,l3 ;
           TextField o1,o2 ;
           Button bt1,bt2,bt3,bt4 ;

           Inleft (Foil target) {
            outerparent = target ;
            setLayout(new GridLayout(8,1,2,2)) ;

            bt1 = new Button("Step 1: Select Model") ;
            bt1.setBackground(Color.white) ;
            bt1.setForeground(Color.red) ;
            bt2 = new Button("Step 2:Set Angle of Attack") ;
            bt2.setBackground(Color.red) ;
            bt2.setForeground(Color.white) ;
            bt3 = new Button("Step 3: Start Tunnel") ;
            bt3.setBackground(Color.red) ;
            bt3.setForeground(Color.white) ;
            bt4 = new Button("Step 4: Adjust for Drag") ;
            bt4.setBackground(Color.red) ;
            bt4.setForeground(Color.white) ;
 
            o1 = new TextField("0.0",5) ;
            l1 = new Label("Step 5: Record Data", Label.CENTER) ;
 
            o2 = new TextField("0.0",5) ;
            l3 = new Label("Step 6: Reduce Data", Label.CENTER) ;
 
            l2 = new Label("Procedure", Label.CENTER) ;
            l2.setBackground(Color.white) ;
            l2.setForeground(Color.red) ;
 
            ipt = new Ipt(outerparent) ;

            add(l2) ;
            add(bt1) ;
            add(bt2) ;
            add(ipt);
            add(bt3) ;
            add(bt4) ;
            add(l1) ;
            add(l3) ;
//            add(o1) ;
//            add(o2) ;
          }

          public boolean action(Event evt, Object arg) {
             Double V3 ;
             double v3 ;
             float fl1 ;

            if(evt.target instanceof Button) {
               String label = (String)arg ;
               if(label.equals("Step 1: Select Model")) {
                   stepper = 1 ;
                   angout = 0.0 ;
                   vfsd = 0.0 ;
                   alfval = v3 = 0.0  ;
                     fl1 = (float) v3 ;
                     ipt.f3.setText(String.valueOf(fl1)) ;
                   ysl = (int) (75 + ((alfval - angmn)/(angmx-angmn))*25.) ;
                   bt1.setBackground(Color.white) ;
                   bt1.setForeground(Color.red) ;
                   bt2.setBackground(Color.red) ;
                   bt2.setForeground(Color.white) ;
                   bt3.setBackground(Color.red) ;
                   bt3.setForeground(Color.white) ;
                   bt4.setBackground(Color.red) ;
                   bt4.setForeground(Color.white) ;
                   ipt.f3.setBackground(Color.red) ;
                   ipt.f3.setForeground(Color.white) ;
                   if (balance == 1) {
                      partimg = getImage(getCodeBase(),"wrightbal.gif");
                   }
                   if (balance == 2) {
                      partimg = getImage(getCodeBase(),"dbalance.gif");
                   }
               }
               if(label.equals("Step 2:Set Angle of Attack")) {
                   stepper = 2 ;
                   angout = 0.0 ;
                   vfsd = 0.0 ;
                   bt1.setBackground(Color.red) ;
                   bt1.setForeground(Color.white) ;
                   bt2.setBackground(Color.white) ;
                   bt2.setForeground(Color.red) ;
                   bt3.setBackground(Color.red) ;
                   bt3.setForeground(Color.white) ;
                   bt4.setBackground(Color.red) ;
                   bt4.setForeground(Color.white) ;
                   ipt.f3.setBackground(Color.white) ;
                   ipt.f3.setForeground(Color.red) ;
                   if (balance == 1) {
                      partimg = getImage(getCodeBase(),"ang.gif");
                   }
                   if (balance == 2) {
                      partimg = getImage(getCodeBase(),"dangle.gif");
                   }
               }
               if(label.equals("Step 3: Start Tunnel")) {
                   stepper = 3 ;
                   trans1 = 0 ;
                   V3 = Double.valueOf(ipt.f3.getText()) ;
                   v3 = V3.doubleValue() ;

                   alfval = v3 ;
                   if(v3 < angmn) {
                     alfval = v3 = angmn  ;
                     fl1 = (float) v3 ;
                     ipt.f3.setText(String.valueOf(fl1)) ;
                   }
                   if(v3 > angmx) {
                      alfval = v3 = angmx ;
                      fl1 = (float) v3 ;
                      ipt.f3.setText(String.valueOf(fl1)) ;
                   }
 
                   ysl = (int) (75 + ((alfval - angmn)/(angmx-angmn))*25.) ;
      
                   vfsd = 20.0 ;
                   bt1.setBackground(Color.red) ;
                   bt1.setForeground(Color.white) ;
                   bt2.setBackground(Color.red) ;
                   bt2.setForeground(Color.white) ;
                   bt3.setBackground(Color.white) ;
                   bt3.setForeground(Color.red) ;
                   bt4.setBackground(Color.red) ;
                   bt4.setForeground(Color.white) ;
                   ipt.f3.setBackground(Color.red) ;
                   ipt.f3.setForeground(Color.white) ;
                   if (balance == 1) {
                      partimg = getImage(getCodeBase(),"bal1.gif");
                   }
                   if (balance == 2) {
                      partimg = getImage(getCodeBase(),"bald1.gif");
                   }
               }
               if(label.equals("Step 4: Adjust for Drag")) {
                   stepper = 4 ;
                   trans2 = 0 ;
                   vfsd = 20.0 ;
                   bt1.setBackground(Color.red) ;
                   bt1.setForeground(Color.white) ;
                   bt2.setBackground(Color.red) ;
                   bt2.setForeground(Color.white) ;
                   bt3.setBackground(Color.red) ;
                   bt3.setForeground(Color.white) ;
                   bt4.setBackground(Color.white) ;
                   bt4.setForeground(Color.red) ;
                   ipt.f3.setBackground(Color.red) ;
                   ipt.f3.setForeground(Color.white) ;
                   if (balance == 1) {
                      partimg = getImage(getCodeBase(),"bal2.gif");
                   }
                   if (balance == 2) {
                      partimg = getImage(getCodeBase(),"bald2.gif");
                   }
               }

               sq.inr.repaint() ;
               computeFlow() ;
               return true ;
            }
            else return false ;
          } // Handler

          class Ipt extends Panel {
             Foil outerparent ;
             TextField f3,o1,o2 ;
             Label l3,l2 ;
      
             Ipt (Foil target) {
      
               outerparent = target ;
               setLayout(new GridLayout(1,3,2,2)) ;

               l3 = new Label("Angle of", Label.RIGHT) ;
               l2 = new Label(" Attack", Label.LEFT) ;
               f3 = new TextField("0.0",5) ;
               f3.setBackground(Color.red) ;
               f3.setForeground(Color.white) ;

               add(l3) ;
               add(l2) ;
               add(f3) ;
             }

             public boolean action(Event evt, Object arg) {
                 if(evt.id == Event.ACTION_EVENT) {
                    this.handleText(evt) ;
                    return true ;
                 }
                 else return false ;
             }

             public void  handleText(Event evt) {
               Double V3 ;
               double v3 ;
               float fl1,fl3 ;

               if (stepper == 2) {
                  V3 = Double.valueOf(f3.getText()) ;
                  v3 = V3.doubleValue() ;

                  alfval = v3 ;
                  if(v3 < angmn) {
                    alfval = v3 = angmn  ;
                    fl1 = (float) v3 ;
                    f3.setText(String.valueOf(fl1)) ;
                  }
                  if(v3 > angmx) {
                     alfval = v3 = angmx ;
                     fl1 = (float) v3 ;
                     f3.setText(String.valueOf(fl1)) ;
                   }

                   ysl = (int) (75 + ((alfval - angmn)/(angmx-angmn))*25.) ;
      
                   computeFlow() ;
                }
                else {
                   fl3 = (float) alfval ;
                   f3.setText(String.valueOf(fl3)) ;
                }
             } // Text Handler
           }  // Ipt
        }  // Inleft

        class Inright extends Panel {
            Foil outerparent ;
            L l;
            Label lab1 ;
            Btp btp ;

            Inright (Foil target) {

             outerparent = target ;
             setLayout(new BorderLayout(5,5)) ;

             l = new L(outerparent) ;
             btp = new Btp(outerparent) ;
             lab1 = new Label("Detailed Dial", Label.CENTER) ;

             add("South",btp) ;
             add("Center",l) ;
             add("North",lab1) ;
           }

           class Btp extends Panel {
              Foil outerparent ;
              Bta bta ;
              Label l1 ;

              Btp (Foil target) {
               outerparent = target ;
               setLayout(new GridLayout(1,2,2,2)) ;

               bta = new Bta(outerparent) ;

               l1 = new Label("Model Group", Label.CENTER) ;
               l1.setBackground(Color.blue) ;
               l1.setForeground(Color.white) ;

               add(l1) ;
               add(bta) ;
             }

             class Bta extends Panel {
                Foil outerparent ;
                Button bt1,bt2,bt3 ;

                Bta (Foil target) {
                 outerparent = target ;
                 setLayout(new GridLayout(1,3,2,2)) ;
  
                 bt1 = new Button("1") ;
                 bt1.setBackground(Color.white) ;
                 bt1.setForeground(Color.blue) ;
                 bt2 = new Button("2") ;
                 bt2.setBackground(Color.blue) ;
                 bt2.setForeground(Color.white) ;
                 bt3 = new Button("3") ;
                 bt3.setBackground(Color.blue) ;
                 bt3.setForeground(Color.white) ;
  
                 add(bt1) ;
                 add(bt2) ;
                 add(bt3) ;
               }
  
               public boolean action(Event evt, Object arg) {
  
                 if(evt.target instanceof Button) {
                    String label = (String)arg ;
                    if(label.equals("1")) {
                        group = 1;
                        bt1.setBackground(Color.white) ;
                        bt1.setForeground(Color.blue) ;
                        bt2.setBackground(Color.blue) ;
                        bt2.setForeground(Color.white) ;
                        bt3.setBackground(Color.blue) ;
                        bt3.setForeground(Color.white) ;
                    }
                    if(label.equals("2")) {
                        group = 2;
                        bt1.setBackground(Color.blue) ;
                        bt1.setForeground(Color.white) ;
                        bt2.setBackground(Color.white) ;
                        bt2.setForeground(Color.blue) ;
                        bt3.setBackground(Color.blue) ;
                        bt3.setForeground(Color.white) ;
                    }
                    if(label.equals("3")) {
                        group = 3;
                        bt1.setBackground(Color.blue) ;
                        bt1.setForeground(Color.white) ;
                        bt2.setBackground(Color.blue) ;
                        bt2.setForeground(Color.white) ;
                        bt3.setBackground(Color.white) ;
                        bt3.setForeground(Color.blue) ;
                    }
                    sq.inl.repaint() ;
                    sq.inr.repaint() ;
                    computeFlow() ;
                    return true ;
                 }
                 else return false ;
               } // Handler
             } // end Bta
           } // end Btp

           class L extends Canvas  {
              Foil outerparent ;

              L (Foil target) {
                setBackground(Color.black) ;
              }

              public void update(Graphics g) {
                ob.inright.l.paint(g) ;
              }

              public void paint(Graphics g) {
                int ex,ey,index ;
                double tick ;
                Color brown ;
   
                brown = new Color((float) (.81),(float) (.59),(float) (.40)) ;

                off3Gg.setColor(brown) ;
                off3Gg.fillRect(0,0,150,150) ;

                off3Gg.setColor(Color.gray) ;
                off3Gg.fillArc(-75,20,200,200,0,90) ;

         // tick marks
                off3Gg.setColor(Color.yellow) ;
                for (index = 0; index <= 8; ++index) {
                   tick = 5. + index*10. ;
                   ex = 25 + (int) (100.0 * Math.sin(convdr * tick)) ;
                   ey = 120 - (int) (100.0 * Math.cos(convdr* tick)) ;
                   off3Gg.drawLine(25,120,ex,ey) ;
                }

                off3Gg.setColor(Color.gray) ;
                off3Gg.fillArc(-65,30,180,180,0,90) ;

                off3Gg.setColor(Color.yellow) ;
                for (index = 0; index <= 9; ++index) {
                   tick = index*10. ;
                   ex = 25 + (int) (100.0 * Math.sin(convdr * tick)) ;
                   ey = 120 - (int) (100.0 * Math.cos(convdr* tick)) ;
                   off3Gg.drawLine(25,120,ex,ey) ;
                }

                off3Gg.setColor(Color.gray) ;
                off3Gg.fillArc(-60,35,170,170,0,90) ;

                off3Gg.setColor(Color.yellow) ;
                if (balance == 1) off3Gg.drawString("0",25,20) ;
                if (balance == 2) off3Gg.drawString("25",25,20) ;
                off3Gg.setColor(Color.yellow) ;
                ex = 25 + (int) (100.0 * Math.sin(convdr *10.0)) ;
                ey = 120 - (int) (100.0 * Math.cos(convdr *10.0)) ;
                if(balance == 1) off3Gg.drawString("10",ex,ey) ;
                if(balance == 2) off3Gg.drawString("20",ex,ey) ;
                ex = 25 + (int) (100.0 * Math.sin(convdr *20.0)) ;
                ey = 120 - (int) (100.0 * Math.cos(convdr *20.0)) ;
                if(balance == 1) off3Gg.drawString("20",ex,ey) ;
                if(balance == 2) off3Gg.drawString("15",ex,ey) ;
                ex = 25 + (int) (100.0 * Math.sin(convdr *30.0)) ;
                ey = 120 - (int) (100.0 * Math.cos(convdr *30.0)) ;
                if(balance == 1) off3Gg.drawString("30",ex,ey) ;
                if(balance == 2) off3Gg.drawString("10",ex,ey) ;
                ex = 25 + (int) (100.0 * Math.sin(convdr *40.0)) ;
                ey = 120 - (int) (100.0 * Math.cos(convdr *40.0)) ;
                if(balance == 1) off3Gg.drawString("40",ex,ey) ;
                if(balance == 2) off3Gg.drawString("5",ex,ey) ;
                ex = 25 + (int) (100.0 * Math.sin(convdr *50.0)) ;
                ey = 120 - (int) (100.0 * Math.cos(convdr *50.0)) ;
                if(balance == 1) off3Gg.drawString("50",ex,ey) ;
                if(balance == 2) off3Gg.drawString("0",ex,ey) ;
                ex = 25 + (int) (100.0 * Math.sin(convdr *60.0)) ;
                ey = 120 - (int) (100.0 * Math.cos(convdr *60.0)) ;
                if(balance == 1) off3Gg.drawString("60",ex,ey) ;
                if(balance == 2) off3Gg.drawString("-5",ex,ey) ;
                ex = 25 + (int) (100.0 * Math.sin(convdr *70.0)) ;
                ey = 120 - (int) (100.0 * Math.cos(convdr *70.0)) ;
                if(balance == 1) off3Gg.drawString("70",ex,ey) ;
                if(balance == 2) off3Gg.drawString("-10",ex,ey) ;
                ex = 25 + (int) (100.0 * Math.sin(convdr *80.0)) ;
                ey = 120 - (int) (100.0 * Math.cos(convdr *80.0)) ;
                if(balance == 1) off3Gg.drawString("80",ex,ey) ;
                if(balance == 2) off3Gg.drawString("-15",ex,ey) ;
                if(balance == 1) off3Gg.drawString("90",125,120) ;
                if(balance == 2) off3Gg.drawString("-20",125,120) ;

                off3Gg.setColor(Color.black) ;
                if (balance == 1) {
                   ex = 25 + (int) (110.0 * Math.sin(convdr *
                           angdraw)) ;
                   ey = 120 - (int) (110.0 * Math.cos(convdr *
                           angdraw)) ;
                }
                if (balance == 2) {
                   ex = 25 + (int) (110.0 * Math.sin(convdr *
                           (50.0 - 2.0*angdraw))) ;
                   ey = 120 - (int) (110.0 * Math.cos(convdr *
                           (50.0 - 2.0*angdraw))) ;
                }
                off3Gg.drawLine(25,120,ex,ey) ;
 
                off3Gg.fillOval(22,117,6,6);

                g.drawImage(offImg3,0,0,this) ;
             }
           } //Lor2a
        }  // Inright
     }  // Outprb

     class Viewer extends Canvas  
         implements Runnable{
        Foil outerparent ;
        Thread runner ;
        Point locate,anchor;
   
        Viewer (Foil target) {
            setBackground(Color.black) ;
            runner = null ;
        }    

        public Insets insets() {
           return new Insets(0,10,0,10) ;
        }
 

        public boolean mouseDrag(Event evt, int x, int y) {
           handle(x,y) ;
           return true;
        }
 
        public void handle(int x, int y) {
            float fl3 ;
         // determine location
            if (stepper == 2) {
               if (x < 105 ) {   
                  if (y >= 75  && y <= 100) {
                     locate = new Point(x,y) ;
                     ysl = locate.y ;
                     alfval = angmn + (double) ((ysl - 75)*(angmx - angmn)/25.) ; 
                     fl3 = (float) alfval;
                     lq.ob.inleft.ipt.f3.setText(String.valueOf(fl3)) ;
                     computeFlow();
                  }
               }
            }
        }

        public boolean mouseUp(Event evt, int x, int y) {
           handleb(x,y) ;
           return true;
        }

        public void handleb(int x, int y) {
          if (stepper == 1) {
             if (y >= 10 && y < 20) {  // tunnel image
                if (x >= 20 && x < 100) {  
                   partimg = getImage(getCodeBase(),"tunnel.gif");
                }
                if (x >= 190 && x < 300) {  // lift balance
                   balance = 1;
                   partimg = getImage(getCodeBase(),"wrightbal.gif");
                }
             }
             if (y >= 21 && y < 35) {  
                if (x >= 190 && x < 300) { //drag balance
                   balance = 2;
                   partimg = getImage(getCodeBase(),"dbalance.gif");
                }
             }
             sq.inr.repaint() ;
             computeFlow() ;
          }
        }

        public void start() {
           if (runner == null) {
              runner = new Thread(this) ;
              runner.start() ;
           }
           antim = 0 ;                              /* MODS  21 JUL 99 */
           ancol = 1 ;                              /* MODS  27 JUL 99 */
        }
   
        public void run() {
          int timer ;
    
          timer = 100 ;
          while (true) {
             ++ trans2 ;
             ++ trans1 ;
             ++ antim ;
             try { Thread.sleep(timer); }
             catch (InterruptedException e) {}
             lq.view.repaint() ;
             if (antim == 3) {
                antim = 0;
                ancol = - ancol ;               /* MODS 27 JUL 99 */
             }
             timer = 135 - (int) (.227 *vfsd/vconv) ;
             angdraw = trans1 * angout / 30.0 ;
             if (trans1 > 30) angdraw = angout ;
             lq.ob.inright.l.repaint() ;
             angdraw2 = angout - trans2 * angout / 30.0 ;
             if (trans2 > 30) angdraw2 = 0.0 ;
          }
        }
   
        public void update(Graphics g) {
           lq.view.paint(g) ;
        }
    
        public void paint(Graphics g) {
           int i,j,k,n ;
           int xlabel,ylabel,ind,inmax,inmin ;
           int xcen,ycen,lstrut,wide,tall ;
           int xtm,ytm ;
           int xh1,yh1,xh2,yh2 ;
           int xh3,yh3,xh4,yh4 ;
           int xh5,yh5,xh6,yh6 ;
           int deg1,deg2 ;
           int exes[] = new int[8] ;
           int whys[] = new int[8] ;
           double offx,scalex,offy,scaley,waste,incy,incx;
           double xtrans,ytrans,xl,yl;
           int camx[] = new int[19] ;
           int camy[] = new int[19] ;
           Color col,brown ;
   
           col = new Color(0,0,0) ;
           if(planet == 0) col = Color.cyan ;
           brown = new Color((float) (.81),(float) (.59),(float) (.40)) ;
           off1Gg.setColor(brown) ;
           off1Gg.fillRect(0,0,300,300) ;

           xtm = xt ;
           ytm = yt ;
// lift balance
           if (balance == 1) {
// mounting bracket
              xcen = 170 ;
              ycen = 165 ;
              lstrut = 40 ;
   
              wide = 16 ;
              tall = 160 ;
              exes[0] = xcen - wide/2 ;
              whys[0] = ycen + wide/2 ;
              exes[1] = exes[0];
              whys[1] = whys[0] - tall - wide/2 ;
              exes[2] = exes[1] + wide ;
              whys[2] = whys[1] ;
              exes[3] = exes[2];
              whys[3] = whys[0] ;
              off1Gg.setColor(Color.darkGray) ;
              off1Gg.fillPolygon(exes,whys,4) ;
// dial
              off1Gg.setColor(Color.gray) ;
              off1Gg.fillArc(120,120,100,100,0,90) ;
              off1Gg.setColor(Color.black) ;
              off1Gg.drawLine(xcen,ycen,
                   (int) (xcen + 40 * Math.sin(convdr*angdraw)),
                   (int) (ycen - 40 * Math.cos(convdr*angdraw))) ;
// lower arms 
              wide = 12 ;
              tall = 4 ;
              exes[0] = (int) (xcen + tall*Math.sin(convdr*(45.0 - angdraw))) ;
              whys[0] = (int) (ycen  -116 + tall*Math.cos(convdr*(45.0 - angdraw))) ;
              exes[1] = (int) (xcen - lstrut*Math.cos(convdr*angdraw) 
                          - tall*Math.sin(convdr*(45.0 + angdraw))) ;
              whys[1] = (int) (ycen - lstrut*Math.sin(convdr*angdraw) -116
                          + tall*Math.sin(convdr*(45.0 + angdraw))) ;
              exes[2] = (int) (xcen - lstrut*Math.cos(convdr*angdraw) 
                          - tall*Math.sin(convdr*(45.0 - angdraw))) ;
              whys[2] = (int) (ycen - lstrut*Math.sin(convdr*angdraw) -116
                          - tall*Math.sin(convdr*(45.0 - angdraw))) ;
              exes[3] = (int) (xcen + tall*Math.sin(convdr*(45.0 + angdraw))) ;
              whys[3] = (int) (ycen -116 - tall*Math.cos(convdr*(45.0 + angdraw))) ; ;
              off1Gg.setColor(Color.gray) ;
              off1Gg.fillPolygon(exes,whys,4) ;
   
              wide = 12 ;
              tall = 4 ;
              exes[0] = (int) (xcen + tall*Math.sin(convdr*(45.0 - angdraw))) ;
              whys[0] = (int) (ycen + tall*Math.cos(convdr*(45.0 - angdraw))) ;
              exes[1] = (int) (xcen - lstrut*Math.cos(convdr*angdraw) 
                             - tall*Math.sin(convdr*(45.0 + angdraw))) ;
              whys[1] = (int) (ycen - lstrut*Math.sin(convdr*angdraw)
                          + tall*Math.sin(convdr*(45.0 + angdraw))) ;
              exes[2] = (int) (xcen - lstrut*Math.cos(convdr*angdraw) 
                          - tall*Math.sin(convdr*(45.0 - angdraw))) ;
              whys[2] = (int) (ycen - lstrut*Math.sin(convdr*angdraw)
                          - tall*Math.sin(convdr*(45.0 - angdraw))) ;
              exes[3] = (int) (xcen + tall*Math.sin(convdr*(45.0 + angdraw))) ;
              whys[3] = (int) (ycen - tall*Math.cos(convdr*(45.0 + angdraw))) ; ;
              off1Gg.setColor(Color.gray) ;
              off1Gg.fillPolygon(exes,whys,4) ;
// lower cross beam 
              wide = 8 ;
              tall = 120 ;
              exes[0] = (int) (xcen - lstrut*Math.cos(convdr*angdraw) -wide/2) ;
              whys[0] = (int) (ycen - lstrut*Math.sin(convdr*angdraw) +wide/2) ;
              exes[1] = exes[0];
              whys[1] = whys[0] - tall -wide/2;
              exes[2] = exes[1] + wide ;
              whys[2] = whys[1] ;
              exes[3] = exes[2];
              whys[3] = whys[0] ;
              off1Gg.setColor(Color.gray) ;
              off1Gg.fillPolygon(exes,whys,4) ;
              off1Gg.setColor(Color.black) ;
              off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
              off1Gg.drawLine(exes[1],whys[1],exes[2],whys[2]) ;
              off1Gg.drawLine(exes[2],whys[2],exes[3],whys[3]) ;
              off1Gg.drawLine(exes[3],whys[3],exes[0],whys[0]) ;
  // plates 
              exes[0] = exes[3];
              whys[0] = whys[2] + 5;
              exes[1] = exes[0] + 5;
              whys[1] = whys[0] + 10;
              exes[2] = exes[0];
              whys[2] = whys[1] + 10;
              off1Gg.fillPolygon(exes,whys,3) ;
   
              exes[0] = exes[3];
              whys[0] = whys[2] + 10;
              exes[1] = exes[0] + 3;
              whys[1] = whys[0] + 10;
              exes[2] = exes[0];
              whys[2] = whys[1] + 10;
              off1Gg.fillPolygon(exes,whys,3) ;
   
              exes[0] = exes[3];
              whys[0] = whys[2] + 10;
              exes[1] = exes[0] + 3;
              whys[1] = whys[0] + 10;
              exes[2] = exes[0];
              whys[2] = whys[1] + 10;
              off1Gg.fillPolygon(exes,whys,3) ;
   
              exes[0] = exes[3];
              whys[0] = whys[2] + 10;
              exes[1] = exes[0] + 5;
              whys[1] = whys[0] + 10;
              exes[2] = exes[0];
              whys[2] = whys[1] + 10;
              off1Gg.fillPolygon(exes,whys,3) ;
   
//   upper cross beam
              if (stepper == 4 && pegged == 0) {
                 xtm = (int) (xt + lstrut*(1.0 - Math.cos(convdr*angdraw2))) ;
                 ytm = (int) (yt - lstrut*Math.sin(convdr*angdraw2)) ;
                 wide = 12 ;
                 tall = 6 ;
                 exes[0] = (int) (xcen + tall*Math.sin(convdr*(45.0 - angdraw2))) ;
                 whys[0] = (int) (ycen  -116 + tall*Math.cos(convdr*(45.0 - angdraw2))) ;
                 exes[1] = (int) (xcen - lstrut*Math.cos(convdr*angdraw2) 
                                - tall*Math.sin(convdr*(45.0 + angdraw2))) ;
                 whys[1] = (int) (ycen - lstrut*Math.sin(convdr*angdraw2) -116
                                + tall*Math.sin(convdr*(45.0 + angdraw2))) ;
                 exes[2] = (int) (xcen - lstrut*Math.cos(convdr*angdraw2) 
                          - tall*Math.sin(convdr*(45.0 - angdraw2))) ;
                 whys[2] = (int) (ycen - lstrut*Math.sin(convdr*angdraw2) -116
                                - tall*Math.sin(convdr*(45.0 - angdraw2))) ;
                 exes[3] = (int) (xcen + tall*Math.sin(convdr*(45.0 + angdraw2))) ;
                 whys[3] = (int) (ycen -116 - tall*Math.cos(convdr*(45.0 + angdraw2))) ; ;
                 off1Gg.setColor(Color.lightGray) ;
                 off1Gg.fillPolygon(exes,whys,4) ;
         
                 wide = 12 ;
                 tall = 6 ;
                 exes[0] = (int) (xcen + tall*Math.sin(convdr*(45.0 - angdraw2))) ;
                 whys[0] = (int) (ycen + tall*Math.cos(convdr*(45.0 - angdraw2))) ;
                 exes[1] = (int) (xcen - lstrut*Math.cos(convdr*angdraw2) 
                                - tall*Math.sin(convdr*(45.0 + angdraw2))) ;
                 whys[1] = (int) (ycen - lstrut*Math.sin(convdr*angdraw2)
                                + tall*Math.sin(convdr*(45.0 + angdraw2))) ;
                 exes[2] = (int) (xcen - lstrut*Math.cos(convdr*angdraw2) 
                                - tall*Math.sin(convdr*(45.0 - angdraw2))) ;
                 whys[2] = (int) (ycen - lstrut*Math.sin(convdr*angdraw2)
                                - tall*Math.sin(convdr*(45.0 - angdraw2))) ;
                 exes[3] = (int) (xcen + tall*Math.sin(convdr*(45.0 + angdraw2))) ;
                 whys[3] = (int) (ycen - tall*Math.cos(convdr*(45.0 + angdraw2))) ; ;
                 off1Gg.setColor(Color.lightGray) ;
                 off1Gg.fillPolygon(exes,whys,4) ;
      
                 wide = 10 ;
                 tall = 125 ;
                 exes[0] = (int) (xcen - lstrut*Math.cos(convdr*angdraw2) -wide/2) ;
                 whys[0] = (int) (ycen - lstrut*Math.sin(convdr*angdraw2) +wide/2) ;
                 exes[1] = exes[0];
                 whys[1] = whys[0] - tall -wide/2;
                 exes[2] = exes[1] + wide ;
                 whys[2] = whys[1] ;
                 exes[3] = exes[2];
                 whys[3] = whys[0] ;
                 off1Gg.setColor(Color.lightGray) ;
                 off1Gg.fillPolygon(exes,whys,4) ;
                 off1Gg.setColor(Color.black) ;
                 off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                 off1Gg.drawLine(exes[1],whys[1],exes[2],whys[2]) ;
                 off1Gg.drawLine(exes[2],whys[2],exes[3],whys[3]) ;
                 off1Gg.drawLine(exes[3],whys[3],exes[0],whys[0]) ;
   
                 off1Gg.setColor(Color.black) ;
                 off1Gg.fillOval(xcen-2,ycen-2,4,4);
                 off1Gg.fillOval(xcen-2,ycen-118,4,4);
                 off1Gg.fillOval((int) (xcen-lstrut*Math.cos(convdr*angdraw2)-2),
                          (int) (ycen-lstrut*Math.sin(convdr*angdraw2)-2),4,4);
                 off1Gg.fillOval((int) (xcen-lstrut*Math.cos(convdr*angdraw2)-2),
                             (int) (ycen-lstrut*Math.sin(convdr*angdraw2)-118),4,4);
              }
              else {
                 xtm = (int) (xt + lstrut*(1.0 - Math.cos(convdr*angdraw))) ;
                 ytm = (int) (yt - lstrut*Math.sin(convdr*angdraw)) ;
   
                 wide = 12 ;
                 tall = 6 ;
                 exes[0] = (int) (xcen + tall*Math.sin(convdr*(45.0 - angdraw))) ;
                 whys[0] = (int) (ycen  -116 + tall*Math.cos(convdr*(45.0 - angdraw))) ;
                 exes[1] = (int) (xcen - lstrut*Math.cos(convdr*angdraw) 
                                - tall*Math.sin(convdr*(45.0 + angdraw))) ;
                 whys[1] = (int) (ycen - lstrut*Math.sin(convdr*angdraw) -116
                                + tall*Math.sin(convdr*(45.0 + angdraw))) ;
                 exes[2] = (int) (xcen - lstrut*Math.cos(convdr*angdraw) 
                          - tall*Math.sin(convdr*(45.0 - angdraw))) ;
                 whys[2] = (int) (ycen - lstrut*Math.sin(convdr*angdraw) -116
                                - tall*Math.sin(convdr*(45.0 - angdraw))) ;
                 exes[3] = (int) (xcen + tall*Math.sin(convdr*(45.0 + angdraw))) ;
                 whys[3] = (int) (ycen -116 - tall*Math.cos(convdr*(45.0 + angdraw))) ; ;
                 off1Gg.setColor(Color.lightGray) ;
                 off1Gg.fillPolygon(exes,whys,4) ;
         
                 wide = 12 ;
                 tall = 6 ;
                 exes[0] = (int) (xcen + tall*Math.sin(convdr*(45.0 - angdraw))) ;
                 whys[0] = (int) (ycen + tall*Math.cos(convdr*(45.0 - angdraw))) ;
                 exes[1] = (int) (xcen - lstrut*Math.cos(convdr*angdraw) 
                                - tall*Math.sin(convdr*(45.0 + angdraw))) ;
                 whys[1] = (int) (ycen - lstrut*Math.sin(convdr*angdraw)
                                + tall*Math.sin(convdr*(45.0 + angdraw))) ;
                 exes[2] = (int) (xcen - lstrut*Math.cos(convdr*angdraw) 
                                - tall*Math.sin(convdr*(45.0 - angdraw))) ;
                 whys[2] = (int) (ycen - lstrut*Math.sin(convdr*angdraw)
                                - tall*Math.sin(convdr*(45.0 - angdraw))) ;
                 exes[3] = (int) (xcen + tall*Math.sin(convdr*(45.0 + angdraw))) ;
                 whys[3] = (int) (ycen - tall*Math.cos(convdr*(45.0 + angdraw))) ; ;
                 off1Gg.setColor(Color.lightGray) ;
                 off1Gg.fillPolygon(exes,whys,4) ;
      
                 wide = 10 ;
                 tall = 125 ;
                 exes[0] = (int) (xcen - lstrut*Math.cos(convdr*angdraw) -wide/2) ;
                 whys[0] = (int) (ycen - lstrut*Math.sin(convdr*angdraw) +wide/2) ;
                 exes[1] = exes[0];
                 whys[1] = whys[0] - tall -wide/2;
                 exes[2] = exes[1] + wide ;
                 whys[2] = whys[1] ;
                 exes[3] = exes[2];
                 whys[3] = whys[0] ;
                 off1Gg.setColor(Color.lightGray) ;
                 off1Gg.fillPolygon(exes,whys,4) ;
                 off1Gg.setColor(Color.black) ;
                 off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                 off1Gg.drawLine(exes[1],whys[1],exes[2],whys[2]) ;
                 off1Gg.drawLine(exes[2],whys[2],exes[3],whys[3]) ;
                 off1Gg.drawLine(exes[3],whys[3],exes[0],whys[0]) ;
   
                 off1Gg.setColor(Color.black) ;
                 off1Gg.fillOval(xcen-2,ycen-2,4,4);
                 off1Gg.fillOval(xcen-2,ycen-118,4,4);
                 off1Gg.fillOval((int) (xcen-lstrut*Math.cos(convdr*angdraw)-2),
                          (int) (ycen-lstrut*Math.sin(convdr*angdraw)-2),4,4);
                 off1Gg.fillOval((int) (xcen-lstrut*Math.cos(convdr*angdraw)-2),
                             (int) (ycen-lstrut*Math.sin(convdr*angdraw)-118),4,4);
              }
           }
// drag balance
           if (balance == 2) {
              xcen = 130 ;
              ycen = 120 ;
              
              glidang = alfval + angdraw ;
              xtm = (int) (xcen - 25 * Math.sin(convdr*glidang)) ;
              ytm = (int) (ycen - 25 * Math.cos(convdr*glidang)) ;

              lstrut = 80 ;
              wide = 16 ;

              xh1 = (int) (xcen + lstrut * Math.cos(convdr*alfval)) ;
              yh1 = (int) (ycen - lstrut * Math.sin(convdr*alfval)) ;
              xh2 = (int) (xcen - lstrut * Math.cos(convdr*alfval)) ;
              yh2 = (int) (ycen + lstrut * Math.sin(convdr*alfval)) ;
// mounting bracket
              exes[0] = (int) (xh1 + wide/2 * Math.cos(convdr*(45.0 + alfval)));
              whys[0] = (int) (yh1 - wide/2 * Math.sin(convdr*(45.0 + alfval)));
              exes[1] = (int) (xh1 + wide/2 * Math.cos(convdr*(45.0 - alfval)));
              whys[1] = (int) (yh1 + wide/2 * Math.sin(convdr*(45.0 - alfval)));
              exes[2] = (int) (xh2 - wide/2 * Math.cos(convdr*(45.0 + alfval)));
              whys[2] = (int) (yh2 + wide/2 * Math.sin(convdr*(45.0 - alfval)));
              exes[3] = (int) (xh2 - wide/2 * Math.cos(convdr*(45.0 - alfval)));
              whys[3] = (int) (yh2 - wide/2 * Math.sin(convdr*(45.0 + alfval)));
              off1Gg.setColor(Color.darkGray) ;
              off1Gg.fillPolygon(exes,whys,4) ;
// dial
              deg1 = (int) (135 + alfval) ;
              off1Gg.setColor(Color.gray) ;
              off1Gg.fillArc(xh1-40,yh1-40,80,80,deg1,90) ;
              off1Gg.setColor(Color.black) ;
              off1Gg.drawLine(xh1,yh1,
                 (int) (xh1 - 40 * Math.cos(convdr*glidang)),
                 (int) (yh1 + 40 * Math.sin(convdr*glidang))) ;
              tall = 45 ;
              xh3 = (int) (xh1 - tall * Math.sin(convdr*glidang)) ;
              yh3 = (int) (yh1 - tall * Math.cos(convdr*glidang)) ;
              xh4 = (int) (xh1 + tall * Math.sin(convdr*glidang)) ;
              yh4 = (int) (yh1 + tall * Math.cos(convdr*glidang)) ;
              xh5 = (int) (xh2 - tall * Math.sin(convdr*glidang)) ;
              yh5 = (int) (yh2 - tall * Math.cos(convdr*glidang)) ;
              xh6 = (int) (xh2 + tall * Math.sin(convdr*glidang)) ;
              yh6 = (int) (yh2 + tall * Math.cos(convdr*glidang)) ;
// arms 
              wide = 12 ;

              exes[0] = (int) (xh3 + wide/2 * Math.cos(convdr*(45.0 +glidang)));
              whys[0] = (int) (yh3 - wide/2 * Math.sin(convdr*(45.0 +glidang)));
              exes[1] = (int) (xh4 + wide/2 * Math.cos(convdr*(45.0 -glidang)));
              whys[1] = (int) (yh4 + wide/2 * Math.sin(convdr*(45.0 -glidang)));
              exes[2] = (int) (xh4 - wide/2 * Math.cos(convdr*(45.0 +glidang)));
              whys[2] = (int) (yh4 + wide/2 * Math.sin(convdr*(45.0 +glidang)));
              exes[3] = (int) (xh3 - wide/2 * Math.cos(convdr*(45.0 -glidang)));
              whys[3] = (int) (yh3 - wide/2 * Math.sin(convdr*(45.0 -glidang)));
              off1Gg.setColor(Color.lightGray) ;
              off1Gg.fillPolygon(exes,whys,4) ;

              exes[0] = (int) (xh5 + wide/2 * Math.cos(convdr*(45.0 +glidang)));
              whys[0] = (int) (yh5 - wide/2 * Math.sin(convdr*(45.0 +glidang)));
              exes[1] = (int) (xh6 + wide/2 * Math.cos(convdr*(45.0 -glidang)));
              whys[1] = (int) (yh6 + wide/2 * Math.sin(convdr*(45.0 -glidang)));
              exes[2] = (int) (xh6 - wide/2 * Math.cos(convdr*(45.0 +glidang)));
              whys[2] = (int) (yh6 + wide/2 * Math.sin(convdr*(45.0 +glidang)));
              exes[3] = (int) (xh5 - wide/2 * Math.cos(convdr*(45.0 -glidang)));
              whys[3] = (int) (yh5 - wide/2 * Math.sin(convdr*(45.0 -glidang)));
              off1Gg.setColor(Color.lightGray) ;
              off1Gg.fillPolygon(exes,whys,4) ;
// beams 
              wide = 8 ;

              exes[0] = (int) (xh3 + wide/2 * Math.cos(convdr*(45.0 + alfval)));
              whys[0] = (int) (yh3 - wide/2 * Math.sin(convdr*(45.0 + alfval)));
              exes[1] = (int) (xh3 + wide/2 * Math.cos(convdr*(45.0 - alfval)));
              whys[1] = (int) (yh3 + wide/2 * Math.sin(convdr*(45.0 - alfval)));
              exes[2] = (int) (xh5 - wide/2 * Math.cos(convdr*(45.0 + alfval)));
              whys[2] = (int) (yh5 + wide/2 * Math.sin(convdr*(45.0 - alfval)));
              exes[3] = (int) (xh5 - wide/2 * Math.cos(convdr*(45.0 - alfval)));
              whys[3] = (int) (yh5 - wide/2 * Math.sin(convdr*(45.0 + alfval)));
              off1Gg.setColor(Color.lightGray) ;
              off1Gg.fillPolygon(exes,whys,4) ;

              exes[0] = (int) (xh4 + wide/2 * Math.cos(convdr*(45.0 + alfval)));
              whys[0] = (int) (yh4 - wide/2 * Math.sin(convdr*(45.0 + alfval)));
              exes[1] = (int) (xh4 + wide/2 * Math.cos(convdr*(45.0 - alfval)));
              whys[1] = (int) (yh4 + wide/2 * Math.sin(convdr*(45.0 - alfval)));
              exes[2] = (int) (xh6 - wide/2 * Math.cos(convdr*(45.0 + alfval)));
              whys[2] = (int) (yh6 + wide/2 * Math.sin(convdr*(45.0 - alfval)));
              exes[3] = (int) (xh6 - wide/2 * Math.cos(convdr*(45.0 - alfval)));
              whys[3] = (int) (yh6 - wide/2 * Math.sin(convdr*(45.0 + alfval)));
              off1Gg.setColor(Color.lightGray) ;
              off1Gg.fillPolygon(exes,whys,4) ;

              off1Gg.setColor(Color.black) ;
              off1Gg.fillOval(xh1-2,yh1-2,4,4);
              off1Gg.fillOval(xh2-2,yh2-2,4,4);
              off1Gg.fillOval(xh3-2,yh3-2,4,4);
              off1Gg.fillOval(xh4-2,yh4-2,4,4);
              off1Gg.fillOval(xh5-2,yh5-2,4,4);
              off1Gg.fillOval(xh6-2,yh6-2,4,4);
           }

// flowfield 
           if (displ <= 3) {  // Side View
            if (vfsd > .01) {
                                               /* plot airfoil flowfield */
             for (j=1; j<=nln2-1; ++j) {           /* lower half */
                exes[1] = (int) (fact*(-xm[j][1])) + xtm ;
                whys[1] = (int) (fact*(-ym[j][1])) + ytm ;
                for (i=2 ; i<= nptc; ++i) {
                   exes[0] = exes[1] ;
                   whys[0] = whys[1] ;
                   exes[1] = (int) (fact*(-xm[j][i])) + xtm ;
                   whys[1] = (int) (fact*(-ym[j][i])) + ytm ;
                   if (displ == 2) {                   /* MODS  21 JUL 99 */
                     off1Gg.setColor(Color.yellow) ;
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                   if (displ == 1 && (i/3*3 == i) ) {
                     off1Gg.setColor(col) ;
                     for (n=1 ; n <= 4 ; ++n) {
                        if(i == 6 + (n-1)*9) off1Gg.setColor(Color.red);
                     }
                     if(i/9*9 == i) off1Gg.setColor(Color.white);
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                   if (displ == 0 && ((i-antim)/3*3 == (i-antim)) ) {
                     if (ancol == -1) {          /* MODS  27 JUL 99 */
                       if((i-antim)/6*6 == (i-antim))off1Gg.setColor(col);
                       if((i-antim)/6*6 != (i-antim))off1Gg.setColor(Color.white);
                     }
                     if (ancol == 1) {          /* MODS  27 JUL 99 */
                       if((i-antim)/6*6 == (i-antim))off1Gg.setColor(Color.white);
                       if((i-antim)/6*6 != (i-antim))off1Gg.setColor(col);
                     }
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                }
             }
             for (j=nln2+1; j<=nlnc; ++j) {          /* upper half */
                exes[1] = (int) (fact*(-xm[j][1])) + xtm ;
                whys[1] = (int) (fact*(-ym[j][1])) + ytm ;
                for (i=2 ; i<= nptc; ++i) {
                   exes[0] = exes[1] ;
                   whys[0] = whys[1] ;
                   exes[1] = (int) (fact*(-xm[j][i])) + xtm ;
                   whys[1] = (int) (fact*(-ym[j][i])) + ytm ;
                   if (displ == 2) {                     /* MODS  21 JUL 99 */
                     off1Gg.setColor(col) ;
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                   if (displ == 1 && (i/3*3 == i) ) {
                     off1Gg.setColor(col);   /* MODS  27 JUL 99 */
                     for (n=1 ; n <= 4 ; ++n) {
                        if(i == 6 + (n-1)*9) off1Gg.setColor(Color.red);
                     }
                     if(i/9*9 == i) off1Gg.setColor(Color.white);
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                   if (displ == 0 && ((i-antim)/3*3 == (i-antim)) ) {
                     if (ancol == -1) {          /* MODS  27 JUL 99 */
                       if((i-antim)/6*6 == (i-antim))off1Gg.setColor(col);
                       if((i-antim)/6*6 != (i-antim))off1Gg.setColor(Color.white);
                     }
                     if (ancol == 1) {          /* MODS  27 JUL 99 */
                       if((i-antim)/6*6 == (i-antim))off1Gg.setColor(Color.white);
                       if((i-antim)/6*6 != (i-antim))off1Gg.setColor(col);
                     }
                     off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                   }
                }
             }
            }
     // draw the airfoil geometry
             exes[1] = (int) (fact*(-xm[0][npt2])) + xtm ;
             whys[1] = (int) (fact*(-ym[0][npt2])) + ytm ;
             for (i=1 ; i<= npt2-1; ++i) {
                exes[0] = exes[1] ;
                whys[0] = whys[1] ;
                exes[1] = (int) (fact*(-xm[0][npt2-i])) + xtm ;
                whys[1] = (int) (fact*(-ym[0][npt2-i])) + ytm ;
                off1Gg.setColor(Color.black) ;
                off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             }
             if (model == 24 || model == 27 || model == 33 || model == 40){
                exes[1] = (int) (fact*(-xm[0][npt2])) + xtm ;
                whys[1] = (int) (fact*(-ym[0][npt2])) + ytm + 30 ;
                for (i=1 ; i<= npt2-1; ++i) {
                   exes[0] = exes[1] ;
                   whys[0] = whys[1] ;
                   exes[1] = (int) (fact*(-xm[0][npt2-i])) + xtm ;
                   whys[1] = (int) (fact*(-ym[0][npt2-i])) + ytm + 30 ;
                   off1Gg.setColor(Color.black) ;
                   off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                }
             }
             if (model == 27 || model == 43) {
                exes[1] = (int) (fact*(-xm[0][npt2])) + xtm ;
                whys[1] = (int) (fact*(-ym[0][npt2])) + ytm + 15 ;
                for (i=1 ; i<= npt2-1; ++i) {
                   exes[0] = exes[1] ;
                   whys[0] = whys[1] ;
                   exes[1] = (int) (fact*(-xm[0][npt2-i])) + xtm ;
                   whys[1] = (int) (fact*(-ym[0][npt2-i])) + ytm + 15 ;
                   off1Gg.setColor(Color.black) ;
                   off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                }
             }
             if (model == 41) {
                exes[1] = (int) (fact*(-xm[0][npt2])) + xtm ;
                whys[1] = (int) (fact*(-ym[0][npt2])) + ytm + 45 ;
                for (i=1 ; i<= npt2-1; ++i) {
                   exes[0] = exes[1] ;
                   whys[0] = whys[1] ;
                   exes[1] = (int) (fact*(-xm[0][npt2-i])) + xtm ;
                   whys[1] = (int) (fact*(-ym[0][npt2-i])) + ytm + 45 ;
                   off1Gg.setColor(Color.black) ;
                   off1Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
                }
             }
          }

   //   some notes
          if (stepper == 1) {
               off1Gg.setColor(Color.white) ;
               off1Gg.drawString("1901 Wind Tunnel",20,20) ;
               off1Gg.setColor(Color.black) ;
               off1Gg.drawString("Test Section",30,35) ;
               off1Gg.drawString("(Overhead View)",22,50) ;

               off1Gg.setColor(Color.black) ;
               off1Gg.drawString("Direction of",220,65) ;
               off1Gg.drawString("Air Flow",240,80) ;
                off1Gg.drawLine(255,90,285,90) ;
                exes[0] = 255 ;  exes[1] = 245; exes[2] = 255 ;
                whys[0] = 85 ;   whys[1] = 90;  whys[2] = 95 ;
                off1Gg.fillPolygon(exes,whys,3) ;

               off1Gg.setColor(Color.yellow) ;
               off1Gg.drawString("Select Step 2.",210,190) ;

               off1Gg.setColor(Color.yellow) ;
                off1Gg.drawLine(200,140,300,140) ;
                exes[0] = 280 ;  exes[1] = 293; exes[2] = 280 ;
                whys[0] = 135 ;  whys[1] = 140; whys[2] = 145 ;
                off1Gg.fillPolygon(exes,whys,3) ;
               off1Gg.setColor(brown) ;
               off1Gg.fillRect(230,130,30,30) ;
               off1Gg.setColor(Color.yellow) ;
               off1Gg.drawString("Dial",235,142) ;

               off1Gg.setColor(Color.white) ;
               if (balance == 1)  off1Gg.setColor(Color.yellow) ;
               off1Gg.drawString("Lift Balance",190,20) ;
               off1Gg.setColor(Color.white) ;
               if (balance == 2)  off1Gg.setColor(Color.yellow) ;
               off1Gg.drawString("Drag Balance",190,35) ;

               off1Gg.setColor(Color.yellow) ;
               off1Gg.drawString("Select Model",50,180) ;
                off1Gg.drawLine(40,165,40,200) ;
                exes[0] = 35 ;  exes[1] = 40; exes[2] = 45 ;
                whys[0] = 182 ;  whys[1] = 192; whys[2] = 182 ;
                off1Gg.fillPolygon(exes,whys,3) ;
                off1Gg.setColor(Color.yellow) ;
                if (model == 1)  off1Gg.drawString("Model #1",25,65) ;
                if (model == 2)  off1Gg.drawString("Model #2",25,65) ;
                if (model == 3)  off1Gg.drawString("Model #3",25,65) ;
                if (model == 4)  off1Gg.drawString("Model #4",25,65) ;
                if (model == 5)  off1Gg.drawString("Model #5",25,65) ;
                if (model == 6)  off1Gg.drawString("Model #6",25,65) ;
                if (model == 7)  off1Gg.drawString("Model #7",25,65) ;
                if (model == 8)  off1Gg.drawString("Model #8",25,65) ;
                if (model == 9)  off1Gg.drawString("Model #9",25,65) ;
                if (model == 10)  off1Gg.drawString("Model #10",25,65) ;
                if (model == 11)  off1Gg.drawString("Model #11",25,65) ;
                if (model == 12)  off1Gg.drawString("Model #12",25,65) ;
                if (model == 15)  off1Gg.drawString("Model #15",25,65) ;
                if (model == 16)  off1Gg.drawString("Model #16",25,65) ;
                if (model == 17)  off1Gg.drawString("Model #17",25,65) ;
                if (model == 18)  off1Gg.drawString("Model #18",25,65) ;
                if (model == 19)  off1Gg.drawString("Model #19",25,65) ;
                if (model == 20)  off1Gg.drawString("Model #20",25,65) ;
                if (model == 21)  off1Gg.drawString("Model #21",25,65) ;
                if (model == 23)  off1Gg.drawString("Model #23",25,65) ;
                if (model == 24)  off1Gg.drawString("Model #24",25,65) ;
                if (model == 25)  off1Gg.drawString("Model #25",25,65) ;
                if (model == 27)  off1Gg.drawString("Model #27",25,65) ;
                if (model == 30)  off1Gg.drawString("Model #30",25,65) ;
                if (model == 31)  off1Gg.drawString("Model #31",25,65) ;
                if (model == 33)  off1Gg.drawString("Model #33",25,65) ;
                if (model == 35)  off1Gg.drawString("Model #35",25,65) ;
                if (model == 40)  off1Gg.drawString("Model #40",25,65) ;
                if (model == 41)  off1Gg.drawString("Model #41",25,65) ;
                if (model == 43)  off1Gg.drawString("Model #43",25,65) ;
                if (model == 51)  off1Gg.drawString("Model #51",25,65) ;
 //               off1Gg.drawString("Mounted Here",10,125) ;
                off1Gg.drawLine(85,60,100,80) ;
                exes[0] = 90 ;  exes[1] = 100; exes[2] = 102 ;
                whys[0] = 75 ;   whys[1] = 80;  whys[2] = 70 ;
                off1Gg.fillPolygon(exes,whys,3) ;
          }
          if (stepper == 2) {
               off1Gg.setColor(Color.yellow) ;
               off1Gg.drawString("Set Angle of Attack",30,20) ;
               off1Gg.drawString("Use Mouse to Move Arm",20,35) ;
               off1Gg.drawString("Or Type into",190,35) ;
               off1Gg.drawString(" Input Box",200,50) ;
               off1Gg.drawString(" then Hit Return",190,65) ;
               off1Gg.drawString("Select Step 3.",210,190) ;
 //  Wilbur's arm
               off1Gg.setColor(Color.black) ;
               off1Gg.fillRect(0,ysl-15,80,30) ;
               off1Gg.setColor(Color.white) ;
               off1Gg.fillRect(80,ysl-8,25,12) ;
               off1Gg.setColor(Color.black) ;
               off1Gg.drawLine(95,ysl-4,105,ysl-4) ;
               off1Gg.drawLine(97,ysl,105,ysl) ;
               off1Gg.setColor(Color.white) ;
               off1Gg.fillRect(80,ysl-12,15,4) ;
               off1Gg.fillRect(80,ysl+4,22,4) ;
               off1Gg.setColor(Color.black) ;
               off1Gg.drawLine(85,ysl-9,95,ysl-9) ;
               off1Gg.drawLine(97,ysl+4,102,ysl+4) ;
          }
          if (stepper == 3) {
               off1Gg.setColor(Color.yellow) ;
               if (balance == 1) off1Gg.drawString("Select Step 4.",210,190) ;
               if (balance == 2) off1Gg.drawString("Record Data.",210,190) ;
          }
          if (stepper == 4) {
               if (balance ==1 && trans2 <= 30) {
 //  Wilbur's arm
                 off1Gg.setColor(Color.black) ;
                 off1Gg.fillRect(0,ytm-15,xtm-10,30) ;
                 off1Gg.setColor(Color.white) ;
                 off1Gg.fillRect(xtm-10,ytm-8,25,12) ;
                 off1Gg.setColor(Color.black) ;
                 off1Gg.drawLine(xtm+5,ytm-4,xtm+15,ytm-4) ;
                 off1Gg.drawLine(xtm+7,ytm,xtm+15,ytm) ;
                 off1Gg.setColor(Color.white) ;
                 off1Gg.fillRect(xtm-10,ytm-12,15,4) ;
                 off1Gg.fillRect(xtm-10,ytm+4,22,4) ;
                 off1Gg.setColor(Color.black) ;
                 off1Gg.drawLine(xtm-5,ytm-9,xtm+5,ytm-9) ;
                 off1Gg.drawLine(xtm+7,ytm+4,xtm+12,ytm+4) ;
               }
               if (trans2 > 30 && pegged == 1) {
                 off1Gg.setColor(Color.yellow) ;
                 off1Gg.drawString("No Adjustment Possible",10,160) ;
                 off1Gg.drawString("at this angle",20,175) ;
               }
               off1Gg.setColor(Color.yellow) ;
               off1Gg.drawString("Record Data.",210,190) ;
          }
          g.drawImage(offImg1,0,0,this) ;   
       }
    } // end Viewer
  }  // end lokpnl
 
  class Sq extends Panel {
     Foil outerparent ;
     Inl inl ;
     Inc inc ;
     Inr inr ;

     Sq (Foil target) { 
         outerparent = target ;
         setLayout(new GridLayout(1,3,5,5)) ;

         inl = new Inl(outerparent) ;
         inc = new Inc(outerparent) ;
         inr = new Inr(outerparent) ;

         add(inl) ;
         add(inr) ;
         add(inc) ;
     }

     class Inl extends Canvas  
         implements Runnable{
        Foil outerparent ;
        Thread runner ;

        Inl (Foil target) {
          setBackground(Color.black) ;
          runner = null ;
        }

        public boolean mouseUp(Event evt, int x, int y) {
           handleb(x,y) ;
           return true;
        }

        public void handleb(int x, int y) {
          if (stepper == 1) {
             if (group == 1) {
                if (y >= 15 && y < 45) {  // model 1
                   model = 1 ;
                   foil = 3 ;
                   camval = 0.0 ;
                   aspr = 1.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod1.gif");
                }
                if (y >= 45 && y < 75) {  // model 2
                   model = 2 ;
                   foil = 3 ;
                   camval = 0.0 ;
                   aspr = 4.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod2.gif");
                }
                if (y >= 75 && y < 105) {  // model 3
                   model = 3 ;
                   foil = 3 ;
                   camval = 0.0 ;
                   aspr = 6.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod3.gif");
                }
                if (y >= 105 && y < 135) {  // model 4
                   model = 4 ;
                   foil = 3 ;
                   camval = 0.33333 ;
                   aspr = 1.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod4.gif");
                }
                if (y >= 135 && y < 165) {  // model 5
                   model = 5 ;
                   foil = 3 ;
                   camval = 0.25 ;
                   aspr = 1.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod5.gif");
                }
                if (y >= 165 && y < 195) {  // model 6
                   model = 6 ;
                   foil = 3 ;
                   camval = 0.2 ;
                   aspr = 1.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod6.gif");
                }
             }
             if (group == 2) {
                if (y >= 15 && y < 45) {  // model 15
                   model = 15 ;
                   foil = 1 ;
                   camval = 0.333 ;
                   thkval = .5 ;
                   aspr = 1.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod15.gif");
                }
                if (y >= 45 && y < 75) {  // model 16
                   model = 16 ;
                   foil = 1 ;
                   camval = 0.25 ;
                   thkval = .5 ;
                   aspr = 1.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod16.gif");
                }
                if (y >= 75 && y < 105) {  // model 17
                   model = 17 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 1.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"nonav.gif");
                }
                if (y >= 105 && y < 135) {  // model 18
                   model = 18 ;
                   foil = 3 ;
                   camval = 0.2 ;
                   aspr = 4.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod18.gif");
                }
                if (y >= 135 && y < 165) {  // model 19
                   model = 19 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 4.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod19.gif");
                }
                if (y >= 165 && y < 195) {  // model 20
                   model = 20 ;
                   foil = 1 ;
                   camval = 0.7 ;
                   thkval = .5 ;
                   aspr = 4.7 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod20.gif");
                }
             }
             if (group == 3) {
                if (y >= 15 && y < 75) {  // model 24
                   model = 24 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 6.7 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod23.gif");
                }
                if (y >= 75 && y < 135) {  // model 27
                   model = 27 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 6.7 ;
                   area = 9.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod23.gif");
                }
                if (y >= 135 && y < 195) {  // model 33
                   model = 33 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 8.0 ;
                   area = 8.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod25.gif");
                }
             }
             inl.repaint() ;
             inr.repaint() ;
             computeFlow() ;
          }
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
             inl.repaint() ;
          }
        }

        public void update(Graphics g) {
          inl.paint(g) ;
        }

        public void paint(Graphics g) {
          int ex,ey,index ;
          int exes[] = new int[8] ;
          int whys[] = new int[8] ;
   
          off2Gg.setColor(Color.blue) ;
          off2Gg.fillRect(0,0,300,200) ;

          off2Gg.setColor(Color.white) ;
          off2Gg.drawString("#",4,12) ;
          off2Gg.drawString("Shape",20,12) ;
          off2Gg.drawString("Area-sq in",70,12) ;
          off2Gg.drawString("Cam",140,12) ;
          off2Gg.drawString("AR",175,12) ;
 // Group 1
          if (group == 1) {
             off2Gg.setColor(Color.white) ;
             if (model == 1) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,15,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,15,198,30) ;
             off2Gg.drawString("1",2,35) ;
             exes[0] = 20 ;
             whys[0] = 30 ;
             exes[1] = 60 ;
             whys[1] = 30 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,20,20,20) ;
             off2Gg.drawString("6",80,35) ;
             off2Gg.drawString("1",180,35) ;
    
             off2Gg.setColor(Color.white) ;
             if (model == 2) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,45,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,45,198,30) ;
             off2Gg.drawString("2",2,65) ;
             exes[0] = 20 ;
             whys[0] = 60 ;
             exes[1] = 60 ;
             whys[1] = 60 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,55,35,15) ;
             off2Gg.drawString("6",85,68) ;
             off2Gg.drawString("4",180,65) ;
   
             off2Gg.setColor(Color.white) ;
             if (model == 3) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,75,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,75,198,30) ;
             off2Gg.drawString("3",2,95) ;
             exes[0] = 20 ;
             whys[0] = 90 ;
             exes[1] = 60 ;
             whys[1] = 90 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,84,60,12) ;
             off2Gg.drawString("6",100,95) ;
             off2Gg.drawString("6",180,95) ;
   
             off2Gg.setColor(Color.white) ;
             if (model == 4) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,105,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,105,200,30) ;
             off2Gg.drawString("4",2,125) ;
             exes[0] = 20 ;
             whys[0] = 121 ;
             exes[1] = 30 ;
             whys[1] = 119 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 119 ;
             exes[1] = 40 ;
             whys[1] = 118 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 118 ;
             exes[1] = 50 ;
             whys[1] = 119 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 119 ;
             exes[1] = 60 ;
             whys[1] = 121 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,110,20,20) ;
             off2Gg.drawString("6",80,125) ;
             off2Gg.drawString("1/12",140,125) ;
             off2Gg.drawString("1",180,125) ;
   
             off2Gg.setColor(Color.white) ;
             if (model == 5) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,135,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,135,200,30) ;
             off2Gg.drawString("5",2,155) ;
             exes[0] = 20 ;
             whys[0] = 150 ;
             exes[1] = 30 ;
             whys[1] = 149 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 149 ;
             exes[1] = 40 ;
             whys[1] = 148 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 148 ;
             exes[1] = 50 ;
             whys[1] = 149 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 149 ;
             exes[1] = 60 ;
             whys[1] = 150 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,140,20,20) ;
             off2Gg.drawString("6",80,155) ;
             off2Gg.drawString("1/16",140,155) ;
             off2Gg.drawString("1",180,155) ;
   
             off2Gg.setColor(Color.white) ;
             if (model == 6) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,165,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,165,200,30) ;
             off2Gg.drawString("6",2,185) ;
             exes[0] = 20 ;
             whys[0] = 180 ;
             exes[1] = 30 ;
             whys[1] = 179 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 179 ;
             exes[1] = 40 ;
             whys[1] = 178 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 178 ;
             exes[1] = 50 ;
             whys[1] = 179 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 179 ;
             exes[1] = 60 ;
             whys[1] = 180 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,170,20,20) ;
             off2Gg.drawString("6",80,185) ;
             off2Gg.drawString("1/20",140,185) ;
             off2Gg.drawString("1",180,185) ;
          }
 // Group 2
          if (group == 2) {
             off2Gg.setColor(Color.white) ;
             if (model == 15) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,15,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,15,198,30) ;
             off2Gg.drawString("15",2,35) ;
             exes[0] = 20 ;
             whys[0] = 31 ;
             exes[1] = 30 ;
             whys[1] = 30 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 30 ;
             exes[1] = 40 ;
             whys[1] = 29 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 29 ;
             exes[1] = 50 ;
             whys[1] = 28 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 28 ;
             exes[1] = 60 ;
             whys[1] = 31 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,20,20,20) ;
             off2Gg.drawString("6",80,35) ;
             off2Gg.drawString("1/12",140,35) ;
             off2Gg.drawString("1",180,35) ;
    
             off2Gg.setColor(Color.white) ;
             if (model == 16) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,45,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,45,198,30) ;
             off2Gg.drawString("16",2,65) ;
             exes[0] = 20 ;
             whys[0] = 60 ;
             exes[1] = 30 ;
             whys[1] = 59 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 59 ;
             exes[1] = 40 ;
             whys[1] = 59 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 59 ;
             exes[1] = 50 ;
             whys[1] = 58 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 58 ;
             exes[1] = 60 ;
             whys[1] = 60 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,50,20,20) ;
             off2Gg.drawString("6",80,65) ;
             off2Gg.drawString("1/16",140,65) ;
             off2Gg.drawString("1",180,65) ;
   
             off2Gg.setColor(Color.white) ;
             if (model == 17) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,75,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,75,198,30) ;
             off2Gg.drawString("17",2,95) ;
             exes[0] = 20 ;
             whys[0] = 90 ;
             exes[1] = 30 ;
             whys[1] = 89 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 89 ;
             exes[1] = 40 ;
             whys[1] = 89 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 89 ;
             exes[1] = 50 ;
             whys[1] = 88 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 88 ;
             exes[1] = 60 ;
             whys[1] = 90 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,80,20,20) ;
             off2Gg.drawString("6",80,95) ;
             off2Gg.drawString("1/20",140,95) ;
             off2Gg.drawString("1",180,95) ;
   
             off2Gg.setColor(Color.white) ;
             if (model == 18) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,105,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,105,200,30) ;
             off2Gg.drawString("18",2,125) ;
             exes[0] = 20 ;
             whys[0] = 120 ;
             exes[1] = 30 ;
             whys[1] = 119 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 119 ;
             exes[1] = 40 ;
             whys[1] = 118 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 118 ;
             exes[1] = 50 ;
             whys[1] = 119 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 119 ;
             exes[1] = 60 ;
             whys[1] = 120 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,115,35,15) ;
             off2Gg.drawString("6",85,128) ;
             off2Gg.drawString("1/20",140,125) ;
             off2Gg.drawString("4",180,125) ;
   
             off2Gg.setColor(Color.white) ;
             if (model == 19) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,135,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,135,200,30) ;
             off2Gg.drawString("19",2,155) ;
             exes[0] = 20 ;
             whys[0] = 150 ;
             exes[1] = 30 ;
             whys[1] = 149 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 149 ;
             exes[1] = 40 ;
             whys[1] = 149 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 149 ;
             exes[1] = 50 ;
             whys[1] = 148 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 148 ;
             exes[1] = 60 ;
             whys[1] = 150 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,145,35,15) ;
             off2Gg.drawString("6",85,158) ;
             off2Gg.drawString("1/20",140,155) ;
             off2Gg.drawString("4",180,155) ;
   
             off2Gg.setColor(Color.white) ;
             if (model == 20) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,165,200,30) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,165,200,30) ;
             off2Gg.drawString("20",2,185) ;
             exes[0] = 20 ;
             whys[0] = 180 ;
             exes[1] = 30 ;
             whys[1] = 179 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 179 ;
             exes[1] = 40 ;
             whys[1] = 178 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 178 ;
             exes[1] = 50 ;
             whys[1] = 179 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 179 ;
             exes[1] = 60 ;
             whys[1] = 180 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
          //  off2Gg.drawRect(70,170,20,20) ;
             exes[0] = 70 ;
             whys[0] = 175 ;
             exes[1] = 130 ;
             whys[1] = 175 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 126 ;
             whys[0] = 179 ;
             off2Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 122 ;
             whys[1] = 182 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 114 ;
             whys[0] = 185 ;
             off2Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 100 ;
             whys[1] = 187 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 86 ;
             whys[0] = 185 ;
             off2Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 78 ;
             whys[1] = 182 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 74 ;
             whys[0] = 179 ;
             off2Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 70 ;
             whys[1] = 175 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawString("6",100,185) ;
             off2Gg.drawString("4.7",175,185) ;
          }
 // Group 3
          if (group == 3) {
             off2Gg.setColor(Color.white) ;
             if (model == 24) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,15,200,60) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,15,198,60) ;
             off2Gg.drawString("24",2,50) ;
             exes[0] = 20 ;
             whys[0] = 30 ;
             exes[1] = 30 ;
             whys[1] = 29 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 29 ;
             exes[1] = 40 ;
             whys[1] = 29 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 29 ;
             exes[1] = 50 ;
             whys[1] = 28 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 28 ;
             exes[1] = 60 ;
             whys[1] = 30 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

             exes[0] = 20 ;
             whys[0] = 60 ;
             exes[1] = 30 ;
             whys[1] = 59 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 59 ;
             exes[1] = 40 ;
             whys[1] = 59 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 59 ;
             exes[1] = 50 ;
             whys[1] = 58 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 58 ;
             exes[1] = 60 ;
             whys[1] = 60 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,40,60,12) ;
             off2Gg.drawString("6",100,51) ;

             off2Gg.drawString("1/20",140,35) ;
             off2Gg.drawString("6.7",175,35) ;
             off2Gg.drawString("Gap",150,50) ;
             off2Gg.drawString("5/8",150,65) ;
             off2Gg.drawString("Foils 22 & 23",65,70) ;
     // 27
             off2Gg.setColor(Color.white) ;
             if (model == 27) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,75,200,60) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,75,198,60) ;
             off2Gg.drawString("27",2,110) ;
             exes[0] = 20 ;
             whys[0] = 90 ;
             exes[1] = 30 ;
             whys[1] = 89 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 89 ;
             exes[1] = 40 ;
             whys[1] = 89 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 89 ;
             exes[1] = 50 ;
             whys[1] = 88 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 88 ;
             exes[1] = 60 ;
             whys[1] = 90 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

             exes[0] = 20 ;
             whys[0] = 105 ;
             exes[1] = 30 ;
             whys[1] = 104 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 104 ;
             exes[1] = 40 ;
             whys[1] = 104 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 104 ;
             exes[1] = 50 ;
             whys[1] = 103 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 103 ;
             exes[1] = 60 ;
             whys[1] = 105 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

             exes[0] = 20 ;
             whys[0] = 120 ;
             exes[1] = 30 ;
             whys[1] = 119 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 119 ;
             exes[1] = 40 ;
             whys[1] = 119 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 119 ;
             exes[1] = 50 ;
             whys[1] = 118 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 118 ;
             exes[1] = 60 ;
             whys[1] = 120 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off2Gg.drawRect(70,100,60,12) ;
             off2Gg.drawString("9",100,111) ;

             off2Gg.drawString("1/20",140,95) ;
             off2Gg.drawString("6.7",175,95) ;
             off2Gg.drawString("Gap",150,110) ;
             off2Gg.drawString("5/8",150,125) ;
             off2Gg.drawString("Foils 22, 23, 26",60,130) ;
   // 33
             off2Gg.setColor(Color.white) ;
             if (model == 33) {
                off2Gg.setColor(Color.white) ;
                off2Gg.fillRect(0,135,200,60) ;
                off2Gg.setColor(Color.blue) ;
             }
             off2Gg.drawRect(0,135,200,60) ;
             off2Gg.drawString("33",2,170) ;
             exes[0] = 20 ;
             whys[0] = 150 ;
             exes[1] = 30 ;
             whys[1] = 149 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 149 ;
             exes[1] = 40 ;
             whys[1] = 149 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 149 ;
             exes[1] = 50 ;
             whys[1] = 148 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 148 ;
             exes[1] = 60 ;
             whys[1] = 150 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

             exes[0] = 20 ;
             whys[0] = 180 ;
             exes[1] = 30 ;
             whys[1] = 179 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 179 ;
             exes[1] = 40 ;
             whys[1] = 179 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 179 ;
             exes[1] = 50 ;
             whys[1] = 178 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 178 ;
             exes[1] = 60 ;
             whys[1] = 180 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

             exes[0] = 70 ;
             whys[0] = 160 ;
             exes[1] = 130 ;
             whys[1] = 160 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 128 ;
             whys[0] = 163 ;
             off2Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 123 ;
             whys[1] = 167 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 120 ;
             whys[0] = 170 ;
             off2Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 115 ;
             whys[1] = 172 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 85 ;
             whys[0] = 172 ;
             off2Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 80 ;
             whys[1] = 170 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 77 ;
             whys[0] = 167 ;
             off2Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 72 ;
             whys[1] = 163 ;
             off2Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 70 ;
             whys[0] = 160 ;
             off2Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             off2Gg.drawString("8",100,170) ;

             off2Gg.drawString("1/20",140,155) ;
             off2Gg.drawString("8",178,155) ;
             off2Gg.drawString("Gap",150,170) ;
             off2Gg.drawString("11/16",150,185) ;
             off2Gg.drawString("Foils 25 & 32",65,190) ;
          }

          g.drawImage(offImg2,0,0,this) ;
       }
     } //Inleft

     class Inc extends Canvas
         implements Runnable{
        Foil outerparent ;
        Thread runner ;

        Inc (Foil target) {
          setBackground(Color.black) ;
          runner = null ;
        }

        public boolean mouseUp(Event evt, int x, int y) {
           handleb(x,y) ;
           return true;
        }

        public void handleb(int x, int y) {
          if (stepper == 1) {
             if (group == 1) {
                if (y >= 15 && y < 45) {  // model 7
                   model = 7 ;
                   foil = 3 ;
                   camval = 0.333 ;
                   aspr = 6.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod7.gif");
                }
                if (y >= 45 && y < 75) {  // model 8
                   model = 8 ;
                   foil = 3 ;
                   camval = 0.25 ;
                   aspr = 6.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod8.gif");
                }
                if (y >= 75 && y < 105) {  // model 9
                   model = 9 ;
                   foil = 3 ;
                   camval = 0.2 ;
                   aspr = 6.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod9.gif");
                }
                if (y >= 105 && y < 135) {  // model 10
                   model = 10 ;
                   foil = 1 ;
                   camval = 0.33333 ;
                   thkval = .5 ;
                   aspr = 6.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod10.gif");
                }
                if (y >= 135 && y < 165) {  // model 11
                   model = 11 ;
                   foil = 1 ;
                   camval = 0.25 ;
                   thkval = .5 ;
                   aspr = 6.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod11.gif");
                }
                if (y >= 165 && y < 195) {  // model 12
                   model = 12 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 6.0 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod12.gif");
                }
             }
             if (group == 2) {
                if (y >= 15 && y < 45) {  // model 21
                   model = 21 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 4.7 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod21.gif");
                }
                if (y >= 45 && y < 75) {  // model 23
                   model = 23 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 6.75 ;
                   area = 6.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod23.gif");
                }
                if (y >= 75 && y < 105) {  // model 25
                   model = 25 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 8.0 ;
                   angout = 0.0 ;
                   area = 4.0 ;
                   partimg = getImage(getCodeBase(),"mod25.gif");
                }
                if (y >= 105 && y < 135) {  // model 30
                   model = 30 ;
                   foil = 1 ;
                   camval = 0.7 ;
                   thkval = 1.0 ;
                   aspr = 4.1 ;
                   area = 5.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod30.gif");
                }
                if (y >= 135 && y < 165) {  // model 31
                   model = 31 ;
                   foil = 1 ;
                   camval = 0.33 ;
                   thkval = .5 ;
                   aspr = 4.64 ;
                   area = 8.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod31.gif");
                }
                if (y >= 165 && y < 195) {  // model 35
                   model = 35 ;
                   foil = 1 ;
                   camval = 0.37 ;
                   thkval = .5 ;
                   aspr = 6.0 ;
                   area = 6.00 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod35.gif");
                }
             }
             if (group == 3) {
                if (y >= 15 && y < 65) {  // model 40
                   model = 40 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 4. ;
                   area = 12.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod18.gif");
                }
                if (y >= 65 && y < 115) {  // model 41 
                   model = 41 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 4.0 ;
                   area = 12.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod18.gif");
                }
                if (y >= 115 && y < 165) {  // model 43
                   model = 43 ;
                   foil = 1 ;
                   camval = 0.2 ;
                   thkval = .5 ;
                   aspr = 4.0 ;
                   area = 12.0 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod18.gif");
                }
                if (y >= 165 && y < 195) {  // model 51
                   model = 51 ;
                   foil = 1 ;
                   camval = 0.37 ;
                   thkval = .5 ;
                   aspr = 4.0 ;
                   area = 6.25 ;
                   angout = 0.0 ;
                   partimg = getImage(getCodeBase(),"mod51.gif");
                }
             }
             inc.repaint() ;
             inr.repaint() ;
             computeFlow() ;
          }
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
             inc.repaint() ;
          }
        }

        public void update(Graphics g) {
          inc.paint(g) ;
        }

        public void paint(Graphics g) {
          int ex,ey,index ;
          int exes[] = new int[8] ;
          int whys[] = new int[8] ;
   
          off4Gg.setColor(Color.blue) ;
          off4Gg.fillRect(0,0,300,200) ;

          off4Gg.setColor(Color.white) ;
          off4Gg.drawString("#",4,12) ;
          off4Gg.drawString("Shape",20,12) ;
          off4Gg.drawString("Area-sq in",70,12) ;
          off4Gg.drawString("Cam",140,12) ;
          off4Gg.drawString("AR",175,12) ;
  // Group 1
          if (group == 1) {
             off4Gg.setColor(Color.white) ;
             if (model == 7) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,15,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,15,198,30) ;
             off4Gg.drawString("7",2,35) ;
             exes[0] = 20 ;
             whys[0] = 31 ;
             exes[1] = 30 ;
             whys[1] = 29 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 29 ;
             exes[1] = 40 ;
             whys[1] = 28 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 28 ;
             exes[1] = 50 ;
             whys[1] = 29 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 29 ;
             exes[1] = 60 ;
             whys[1] = 31 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,24,60,12) ;
             off4Gg.drawString("6",100,35) ;
             off4Gg.drawString("1/12",140,35) ;
             off4Gg.drawString("6",180,35) ;
 
             off4Gg.setColor(Color.white) ;
             if (model == 8) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,45,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,45,198,30) ;
             off4Gg.drawString("8",2,65) ;
             exes[0] = 20 ;
             whys[0] = 61 ;
             exes[1] = 30 ;
             whys[1] = 60 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 60 ;
             exes[1] = 40 ;
             whys[1] = 59 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 59 ;
             exes[1] = 50 ;
             whys[1] = 60 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 60 ;
             exes[1] = 60 ;
             whys[1] = 61 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,54,60,12) ;
             off4Gg.drawString("6",100,65) ;
             off4Gg.drawString("1/16",140,65) ;
             off4Gg.drawString("6",180,65) ;
   
             off4Gg.setColor(Color.white) ;
             if (model == 9) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,75,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,75,198,30) ;
             off4Gg.drawString("9",2,95) ;
             exes[0] = 20 ;
             whys[0] = 91 ;
             exes[1] = 30 ;
             whys[1] = 90 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 90 ;
             exes[1] = 40 ;
             whys[1] = 89 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 89 ;
             exes[1] = 50 ;
             whys[1] = 90 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 90 ;
             exes[1] = 60 ;
             whys[1] = 91 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,84,60,12) ;
             off4Gg.drawString("6",100,95) ;
             off4Gg.drawString("1/20",140,95) ;
             off4Gg.drawString("6",180,95) ;
   
             off4Gg.setColor(Color.white) ;
             if (model == 10) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,105,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,105,200,30) ;
             off4Gg.drawString("10",2,125) ;
             exes[0] = 20 ;
             whys[0] = 121 ;
             exes[1] = 30 ;
             whys[1] = 120 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 120 ;
             exes[1] = 40 ;
             whys[1] = 119 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 119 ;
             exes[1] = 50 ;
             whys[1] = 118 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 118 ;
             exes[1] = 60 ;
             whys[1] = 121 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,114,60,12) ;
             off4Gg.drawString("6",100,125) ;
             off4Gg.drawString("1/12",140,125) ;
             off4Gg.drawString("6",180,125) ;

             off4Gg.setColor(Color.white) ;
             if (model == 11) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,135,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,135,200,30) ;
             off4Gg.drawString("11",2,155) ;
             exes[0] = 20 ;
             whys[0] = 150 ;
             exes[1] = 30 ;
             whys[1] = 149 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 149 ;
             exes[1] = 40 ;
             whys[1] = 149 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 149 ;
             exes[1] = 50 ;
             whys[1] = 148 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 148 ;
             exes[1] = 60 ;
             whys[1] = 150 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,144,60,12) ;
             off4Gg.drawString("6",100,155) ;
             off4Gg.drawString("1/16",140,155) ;
             off4Gg.drawString("6",180,155) ;
   
             off4Gg.setColor(Color.white) ;
             if (model == 12) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,165,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,165,200,30) ;
             off4Gg.drawString("12",2,185) ;
             exes[0] = 20 ;
             whys[0] = 180 ;
             exes[1] = 30 ;
             whys[1] = 179 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 179 ;
             exes[1] = 40 ;
             whys[1] = 179 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 179 ;
             exes[1] = 50 ;
             whys[1] = 178 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 178 ;
             exes[1] = 60 ;
             whys[1] = 180 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,174,60,12) ;
             off4Gg.drawString("6",100,185) ;
             off4Gg.drawString("1/20",140,185) ;
             off4Gg.drawString("6",180,185) ;
          }
  // Group 2
          if (group == 2) {
             off4Gg.setColor(Color.white) ;
             if (model == 21) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,15,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,15,198,30) ;
             off4Gg.drawString("21",2,35) ;
             exes[0] = 20 ;
             whys[0] = 30 ;
             exes[1] = 30 ;
             whys[1] = 29 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 29 ;
             exes[1] = 40 ;
             whys[1] = 29 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 29 ;
             exes[1] = 50 ;
             whys[1] = 28 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 28 ;
             exes[1] = 60 ;
             whys[1] = 30 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
          //   off4Gg.drawRect(70,24,60,12) ;
             exes[0] = 70 ;
             whys[0] = 24 ;
             exes[1] = 130 ;
             whys[1] = 24 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 126 ;
             whys[0] = 28 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 122 ;
             whys[1] = 31 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 114 ;
             whys[0] = 34 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 100 ;
             whys[1] = 36 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 86 ;
             whys[0] = 34 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 78 ;
             whys[1] = 31 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 74 ;
             whys[0] = 28 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 70 ;
             whys[1] = 24 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

             off4Gg.drawString("6",100,35) ;
             off4Gg.drawString("1/20",140,35) ;
             off4Gg.drawString("4.7",175,35) ;
 
             off4Gg.setColor(Color.white) ;
             if (model == 23) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,45,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,45,198,30) ;
             off4Gg.drawString("23",2,65) ;
             exes[0] = 20 ;
             whys[0] = 60 ;
             exes[1] = 30 ;
             whys[1] = 59 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 59 ;
             exes[1] = 40 ;
             whys[1] = 59 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 59 ;
             exes[1] = 50 ;
             whys[1] = 58 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 58 ;
             exes[1] = 60 ;
             whys[1] = 60 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,54,60,12) ;
             off4Gg.drawString("3",100,65) ;
             off4Gg.drawString("1/20",140,65) ;
             off4Gg.drawString("6.7",175,65) ;
   
             off4Gg.setColor(Color.white) ;
             if (model == 25) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,75,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,75,198,30) ;
             off4Gg.drawString("25",2,95) ;
             exes[0] = 20 ;
             whys[0] = 90 ;
             exes[1] = 30 ;
             whys[1] = 89 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 89 ;
             exes[1] = 40 ;
             whys[1] = 89 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 89 ;
             exes[1] = 50 ;
             whys[1] = 88 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 88 ;
             exes[1] = 60 ;
             whys[1] = 90 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
        //    off4Gg.drawRect(70,84,60,12) ;
             exes[0] = 70 ;
             whys[0] = 84 ;
             exes[1] = 130 ;
             whys[1] = 84 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 128 ;
             whys[0] = 87 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 123 ;
             whys[1] = 91 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 120 ;
             whys[0] = 94 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 115 ;
             whys[1] = 96 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 85 ;
             whys[0] = 96 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 80 ;
             whys[1] = 94 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 77 ;
             whys[0] = 91 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 72 ;
             whys[1] = 87 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 70 ;
             whys[0] = 84 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             off4Gg.drawString("4",100,95) ;
             off4Gg.drawString("1/20",140,95) ;
             off4Gg.drawString("8",180,95) ;
   
             off4Gg.setColor(Color.white) ;
             if (model == 30) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,105,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,105,200,30) ;
             off4Gg.drawString("30",2,125) ;
             exes[0] = 20 ;
             whys[0] = 120;
             exes[1] = 30 ;
             whys[1] = 119 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 119 ;
             exes[1] = 40 ;
             whys[1] = 118 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 118 ;
             exes[1] = 50 ;
             whys[1] = 117 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 117 ;
             exes[1] = 60 ;
             whys[1] = 122 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,114,60,12) ;
             off4Gg.drawString("5",100,125) ;
             off4Gg.drawString("1/6",140,125) ;
             off4Gg.drawString("4.1",175,125) ;

             off4Gg.setColor(Color.white) ;
             if (model == 31) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,135,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,135,200,30) ;
             off4Gg.drawString("31",2,155) ;
             exes[0] = 20 ;
             whys[0] = 151 ;
             exes[1] = 30 ;
             whys[1] = 150 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 150 ;
             exes[1] = 40 ;
             whys[1] = 149 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 149 ;
             exes[1] = 50 ;
             whys[1] = 148 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 148 ;
             exes[1] = 60 ;
             whys[1] = 151 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
        //    off4Gg.drawRect(70,144,60,12) ;
             exes[1] = 130 ;
             whys[1] = 150 ;
             exes[0] = 126 ;
             whys[0] = 148 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 122 ;
             whys[1] = 147 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 114 ;
             whys[0] = 145 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 100 ;
             whys[1] = 144 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 86 ;
             whys[0] = 145 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 78 ;
             whys[1] = 147 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 74 ;
             whys[0] = 148 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 70 ;
             whys[1] = 150 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 74 ;
             whys[0] = 152 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 78 ;
             whys[1] = 153 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 86 ;
             whys[0] = 155 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 100 ;
             whys[1] = 156 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 114 ;
             whys[0] = 155 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 122 ;
             whys[1] = 153 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 126 ;
             whys[0] = 152 ;
             off4Gg.drawLine(exes[1],whys[1],exes[0],whys[0]) ;
             exes[1] = 130 ;
             whys[1] = 150 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawString("8",100,155) ;
             off4Gg.drawString("1/12",140,155) ;
             off4Gg.drawString("4.6",175,155) ;
   
             off4Gg.setColor(Color.white) ;
             if (model == 35) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,165,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,165,200,30) ;
             off4Gg.drawString("35",2,185) ;
             exes[0] = 20 ;
             whys[0] = 180 ;
             exes[1] = 30 ;
             whys[1] = 179 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 179 ;
             exes[1] = 40 ;
             whys[1] = 179 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 179 ;
             exes[1] = 50 ;
             whys[1] = 178 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 178 ;
             exes[1] = 60 ;
             whys[1] = 180 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,174,60,12) ;
             off4Gg.drawString("6",100,185) ;
             off4Gg.drawString("1/10",140,185) ;
             off4Gg.drawString("6",180,185) ;
          }
 // Group 3
          if (group == 3) {
             off4Gg.setColor(Color.white) ;
             if (model == 40) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,15,200,50) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,15,198,50) ;
             off4Gg.drawString("40",2,45) ;
             exes[0] = 20 ;
             whys[0] = 25 ;
             exes[1] = 30 ;
             whys[1] = 24 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 24 ;
             exes[1] = 40 ;
             whys[1] = 24 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 24 ;
             exes[1] = 50 ;
             whys[1] = 23 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 23 ;
             exes[1] = 60 ;
             whys[1] = 25 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

             exes[0] = 20 ;
             whys[0] = 50 ;
             exes[1] = 30 ;
             whys[1] = 49 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 49 ;
             exes[1] = 40 ;
             whys[1] = 49 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 49 ;
             exes[1] = 50 ;
             whys[1] = 48 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 48 ;
             exes[1] = 60 ;
             whys[1] = 50 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,30,60,12) ;
             off4Gg.drawString("12",100,40) ;

             off4Gg.drawString("1/20",140,28) ;
             off4Gg.drawString("4",178,28) ;
             off4Gg.drawString("Gap",150,43) ;
             off4Gg.drawString("11/16",150,58) ;
             off4Gg.drawString("Foils 18 & 19",65,63) ;
     //41 
             off4Gg.setColor(Color.white) ;
             if (model == 41) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,65,200,50) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,65,198,50) ;
             off4Gg.drawString("41",2,95) ;
             exes[0] = 20 ;
             whys[0] = 75 ;
             exes[1] = 30 ;
             whys[1] = 74 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 74 ;
             exes[1] = 40 ;
             whys[1] = 74 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 74 ;
             exes[1] = 50 ;
             whys[1] = 73 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 73 ;
             exes[1] = 60 ;
             whys[1] = 75 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

             exes[0] = 20 ;
             whys[0] = 100 ;
             exes[1] = 30 ;
             whys[1] = 99 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 99 ;
             exes[1] = 40 ;
             whys[1] = 99 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 99 ;
             exes[1] = 50 ;
             whys[1] = 98 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 98 ;
             exes[1] = 60 ;
             whys[1] = 100 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,80,60,12) ;
             off4Gg.drawString("12",100,90) ;

             off4Gg.drawString("1/20",140,78) ;
             off4Gg.drawString("4",178,78) ;
             off4Gg.drawString("Gap",150,93) ;
             off4Gg.drawString("18/16",150,108) ;
             off4Gg.drawString("Foils 18 & 19",65,113) ;
   // 43
             off4Gg.setColor(Color.white) ;
             if (model == 43) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,115,200,50) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,115,200,50) ;
             off4Gg.drawString("43",2,145) ;
             exes[0] = 20 ;
             whys[0] = 125 ;
             exes[1] = 30 ;
             whys[1] = 124 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 124 ;
             exes[1] = 40 ;
             whys[1] = 124 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 124 ;
             exes[1] = 50 ;
             whys[1] = 123 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 123 ;
             exes[1] = 60 ;
             whys[1] = 125 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;

             exes[0] = 20 ;
             whys[0] = 150 ;
             exes[1] = 30 ;
             whys[1] = 149 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 149 ;
             exes[1] = 40 ;
             whys[1] = 149 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 149 ;
             exes[1] = 50 ;
             whys[1] = 148 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 148 ;
             exes[1] = 60 ;
             whys[1] = 150 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,130,60,12) ;
             off4Gg.drawString("12",100,140) ;

             off4Gg.drawString("1/20",140,128) ;
             off4Gg.drawString("4",178,128) ;
             off4Gg.drawString("Gap",150,143) ;
             off4Gg.drawString("1/4",150,158) ;
             off4Gg.drawString("Foils 18 & 19",65,163) ;
    // 51
             off4Gg.setColor(Color.white) ;
             if (model == 51) {
                off4Gg.setColor(Color.white) ;
                off4Gg.fillRect(0,165,200,30) ;
                off4Gg.setColor(Color.blue) ;
             }
             off4Gg.drawRect(0,165,200,30) ;
             off4Gg.drawString("51",2,185) ;
             exes[0] = 20 ;
             whys[0] = 180 ;
             exes[1] = 30 ;
             whys[1] = 179 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 30 ;
             whys[0] = 179 ;
             exes[1] = 40 ;
             whys[1] = 179 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 40 ;
             whys[0] = 179 ;
             exes[1] = 50 ;
             whys[1] = 178 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             exes[0] = 50 ;
             whys[0] = 178 ;
             exes[1] = 60 ;
             whys[1] = 180 ;
             off4Gg.drawLine(exes[0],whys[0],exes[1],whys[1]) ;
             off4Gg.drawRect(70,174,60,12) ;
             off4Gg.drawString("6.25",95,185) ;
             off4Gg.drawString("1/11",140,185) ;
             off4Gg.drawString("4",180,185) ;
          }

          g.drawImage(offImg4,0,0,this) ;
       }
      }  // Incenter

     class Inr extends Canvas { 
         Foil outerparent ;

         Inr (Foil target) {
          setBackground(Color.black) ;
        }

        public void update(Graphics g) {
           inr.paint(g) ;
        }
 
        public void paint(Graphics g) {

          off5Gg.setColor(Color.white) ;
          off5Gg.fillRect(0,0,300,200) ;

           off5Gg.drawImage(partimg,0,0,this) ;
           g.drawImage(offImg5,0,0,this) ;
        }
      }  // Inright 
  } //Slvpnl 
}
