/*
                Interactive Isentropic Flow  Program

     Program to perform one dimensional isentropic flow analysis 
                        from NACA 1135

                     Version 1.2a   - 30 May 06

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


      New Test :  * Derive from Nozzle program and vucalc
                  * Two panels -- input (Up) - output (Dn)
                    add compute button
                    move input variables


                                                     TJB 30 May 06
*/
package Isentrop;

import java.awt.*;
import java.lang.Math ;

public class Isentrop extends java.applet.Applet {

   final double convdr = 3.14515926/180.;
   
   static double gama, gm1, gp1 ;
   static double mach1, mmin, mmax ;
   static double arat, amin, amax ;
   static double prat, pmin, pmax ;
   static double trat, tmin, tmax ;
   static double rhorat, rhomin, rhomax ;
   static double qrat, qmin, qmax ;
   static double wcor, wmin, wmax ;
   static double nu, numin, numax ;
   static double mu, mumin, mumax ;
   static double inptval ;
   static int inopt, isub ;

   In in ;

   public void init() {

     setLayout(new GridLayout(1,2,0,0)) ;

     setDefaults () ;

     in = new In(this) ;

     add(in) ;

     comPute() ;
  }

  public Insets insets() {
     return new Insets(10,10,10,10) ;
  }

  public void setDefaults() {

       inopt = 0 ;
       isub = 0 ;
       inptval = 2.0 ;
       gama = 1.4 ;
 
       mmin = .01;  mmax = 3.0 ;
       pmin = .000001;  pmax = 1.0 ;
       tmin = .005;  tmax = 1.0 ;
       rhomin = .0005;  rhomax = 1.0 ;
       qmin = .0005;  qmax =  400. ;
       amin = 1.0 ;  amax = 60.0  ;
       wmin = .00002  ;  wmax = .343 ;
       mumin = 2.0  ;  mumax = 90.0 ;
       numin = 2.0  ;  numax = 120. ;
  }

  public void comPute() {

    gm1 = gama -1.0 ;
    gp1 = gama +1.0 ;

    switch (inopt) {
        case 0: {                               /* mach number given */
             mach1 = inptval ;
             break ;
        }
        case 1: {                               /* pressure ratio given */
             prat = inptval ;
             mach1 = Math.sqrt(2.0*(Math.pow(prat,-gm1/gama) - 1.0)/gm1) ;
             break ;
        }
        case 2: {                               /* temperature ratio given */
             trat = inptval ;
             mach1 = Math.sqrt(2.0*(Math.pow(trat,-1.0) - 1.0)/gm1) ;
             break ;
        }
        case 3: {                               /* density ratio given */
             rhorat = inptval ;
             mach1 = Math.sqrt(2.0*(Math.pow(rhorat,-gm1) - 1.0)/gm1) ;
             break ;
        }
        case 4: {                               /* area ratio given */
             arat = inptval ;
             getMachArat() ;
             break ;
        }
        case 5: {                               /* q ratio given */
             qrat = inptval ;
             mach1 = Math.sqrt(2.0*qrat/gama) ;
             break ;
        }
        case 6: {                               /* Corrected Airflow given */
             wcor = inptval ;
             mach1 = getMach(isub, wcor) ;
             break ;
        }
        case 7: {                               /* Mach angle given */
             mu = inptval ;
             mach1 = 1.0 / Math.sin(mu * convdr) ;
             break ;
        }
        case 8: {                               /* Prandtl-Meyer angle given */
             nu = inptval ;
             getMachpm() ;
             break ;
        }
    }
    getIsen () ;
    return ; 
  }   
 
  public void getIsen () {                /* isentropic relations */
     double mach1s,msm1,fac,fac1 ;            /* mach is given        */
                                              /* prat and trat are ratiod */
     mach1s = mach1*mach1 ;                   /* to total conditions - */
     msm1 = mach1s - 1.0;
     fac = 1.0 + .5*gm1*mach1s;

     prat = Math.pow(1.0/fac,gama/gm1) ;                  /* EQ 44 */
     trat = 1.0 / fac ;                                   /* EQ 43 */
     rhorat = Math.pow(1.0/fac,1.0/gm1) ;                 /* EQ 45 */
     fac1 = gp1/(2.0*gm1) ;
     arat = mach1 * Math.pow(fac,-fac1) * Math.pow(gp1/2.0,fac1) ; /* EQ 80 */
     arat = 1.0/arat ;
     qrat = gama*mach1s*.5 ;                              /* EQ 47 */
     wcor = getAir(mach1) ;
     mu   = (Math.asin(1.0/mach1))/convdr ;
     nu   = Math.sqrt(gp1/gm1)*Math.atan(Math.sqrt(gm1*msm1/gp1)) 
             - Math.atan(Math.sqrt(msm1)) ;
     nu   = nu / convdr;

     loadOut() ;
     return;
  }

  public double getAir(double mach) {
/* Utility to get the corrected airflow per area given the Mach number */
       double number,fac1,fac2;
       fac2 = (gama+1.0)/(2.0*(gama-1.0)) ;
       fac1 = Math.pow((1.0+.5*(gama-1.0)*mach*mach),fac2);
       number =  .50161*Math.sqrt(gama) * mach/ fac1 ;

       return(number) ;
  }

  public double getMach (int sub, double corair) {
/* Utility to get the Mach number given the corrected airflow per area */
         double number,chokair;              /* iterate for mach number */
         double deriv,machn,macho,airo,airn;
         int iter ;

         chokair = getAir(1.0) ;
         if (corair > chokair) {
           number = 1.0 ;
           return (number) ;
         }
         else {
           airo = .25618 ;                 /* initial guess */
           if (sub == 2) macho = 1.0 ;   /* sonic */
           else {
              if (sub == 0) macho = 1.703 ; /* supersonic */
              else macho = .5;                /* subsonic */
              iter = 1 ;
              machn = macho - .2  ;
              while (Math.abs(corair - airo) > .0001 && iter < 20) {
                 airn =  getAir(machn) ;
                 deriv = (airn-airo)/(machn-macho) ;
                 airo = airn ;
                 macho = machn ;
                 machn = macho + (corair - airo)/deriv ;
                 ++ iter ;
              }
           }
           number = macho ;
         }
         return(number) ;
   }

   public void getMachArat()  {                 /*  get the Mach number */
                                            /* given the area ratio A/Astar */
      double deriv,machn,macho,aro,arn,fac,fac1;

      fac1 = gp1/(2.0*gm1) ;

      aro = 2.0 ;                 /* initial guess */
      macho = 2.2 ;                      /* supersonic */
      if (isub == 1) macho = .30;        /* subsonic */
      machn = macho + .05  ;
      while (Math.abs(arat - aro) > .0001) {
         fac = 1.0 + .5*gm1*machn*machn;
         arn = 1.0/(machn * Math.pow(fac,-fac1) * Math.pow(gp1/2.0,fac1))  ;
         deriv = (arn-aro)/(machn-macho) ;
         aro = arn ;
         macho = machn ;
         machn = macho + (arat - aro)/deriv ;
      }
      mach1 = macho ;
      return ;
   }

   public void getMachpm ()  {                 /* get the Mach number */
                                             /* given the Prandtl-meyer angle */
      double msm1o,msm1n;
      double nur,nuo,nun,deriv ;

      nur = nu*convdr ;

      msm1o = 2.0;                                  /* iterate for mach */
      nuo = Math.sqrt(gp1/gm1)*Math.atan(Math.sqrt(gm1*msm1o/gp1)) 
             - Math.atan(Math.sqrt(msm1o));
      msm1n = msm1o+.01 ;
      while (Math.abs(nur - nuo) > .0001) {
         nun = Math.sqrt(gp1/gm1)*Math.atan(Math.sqrt(gm1*msm1n/gp1)) 
                - Math.atan(Math.sqrt(msm1n));
         deriv = (nun-nuo)/(msm1n-msm1o) ;
         nuo = nun ;
         msm1o = msm1n ;
         msm1n = msm1o + (nur-nuo)/deriv ;
       }
      mach1 = Math.sqrt(msm1o + 1.0);
      return ;
   }

   public void loadOut() {

      in.dn.o1.setText(String.valueOf(filter3(mach1))) ;
      in.dn.o2.setText(String.valueOf(filter5(prat)));
      in.dn.o3.setText(String.valueOf(filter3(trat)));
      in.dn.o4.setText(String.valueOf(filter5(rhorat)));
      in.dn.o5.setText(String.valueOf(filter3(arat)));
      in.dn.o6.setText(String.valueOf(filter3(qrat)));
      in.dn.o7.setText(String.valueOf(filter5(wcor)));
      in.dn.o8.setText(String.valueOf(filter3(mu)));
      in.dn.o9.setText(String.valueOf(filter3(nu)));
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

   class In extends Panel {
        Isentrop outerparent ;
        Up up ;
        Dn dn ;

        In (Isentrop target) {

           outerparent = target ;
           setLayout(new GridLayout(2,1,5,5)) ;

           up = new Up(outerparent) ;
           dn = new Dn(outerparent) ;

           add(up) ;
           add(dn) ;
        }

        class Up extends Panel {
             Isentrop outerparent ;
             Label l1,l2,l3 ;
             TextField f1,f2;
             Choice inch,sub ;
             Button bt1 ;

             Up (Isentrop target) {

                outerparent = target ;
                setLayout(new GridLayout(4,3,5,5)) ;

                l1 = new Label("Isentropic Flow", Label.RIGHT) ;
                l1.setForeground(Color.red) ;
                l2 = new Label(" Calculator", Label.LEFT) ;
                l2.setForeground(Color.red) ;

                l3 = new Label("Input:", Label.LEFT) ;
                l3.setForeground(Color.red) ;

                inch = new Choice() ;
                inch.addItem("Mach Number - M =") ;
                inch.addItem("Pressure Ratio  p/pt =");
                inch.addItem("Temperature Ratio  T/Tt =");
                inch.addItem("Density Ratio  rho/rhot =");
                inch.addItem("Area Ratio  A/A* =");
                inch.addItem("Dynamic Press Ratio  q/p =");
                inch.addItem("Flow per Area Wcor/A =");
                inch.addItem("Mach Angle mu =");
                inch.addItem("Prandtl-Meyer Angle nu =");
                inch.select(0) ;
 
                sub = new Choice() ;
                sub.addItem("Supersonic");
                sub.addItem("Subsonic") ;
                sub.select(0) ;

                f1 = new TextField(String.valueOf(inptval),5) ;
                f2 = new TextField(String.valueOf(gama),5) ;

                bt1 = new Button("COMPUTE") ;
                bt1.setBackground(Color.red) ;
                bt1.setForeground(Color.white) ;

                add(l1) ;
                add(l2) ;
                add(sub) ;

                add(l3) ;
                add(new Label("Gamma", Label.RIGHT)) ;
                add(f2) ;

                add(new Label("Input Variable", Label.RIGHT)) ;
                add(inch) ;
                add(f1) ;

                add(new Label(" ", Label.RIGHT)) ;
                add(bt1) ;
                add(new Label(" ", Label.RIGHT)) ;
             }
 
             public boolean handleEvent(Event evt) {
                Double V1, V2 ;
                double v1, v2, varmin, varmax ;
                float fl1 ;

                if(evt.id == Event.ACTION_EVENT) {
                  inopt = inch.getSelectedIndex() ;
                  isub = sub.getSelectedIndex() ;

                  V2 = Double.valueOf(f2.getText()) ;
                  v2 = V2.doubleValue() ;
                  if(v2 < 1.0) {
                       v2 = 1.0 ;
                       fl1 = (float) v2 ;
                       f2.setText(String.valueOf(fl1)) ;
                  }
                  if(v2 > 1.6) {
                       v2 = 1.6 ;
                       fl1 = (float) v2 ;
                       f2.setText(String.valueOf(fl1)) ;
                  }
                  gama = v2 ;

                  varmin = 0.0; varmax = 1.0 ;
                  switch (inopt) {
                     case 0: {
                         varmin = mmin ; varmax = mmax;
                         break ;
                     }
                     case 1: {
                         varmin = pmin ; varmax = pmax;
                         break ;
                     }
                     case 2: {
                         varmin = tmin ; varmax = tmax;
                         break ;
                     }
                     case 3: {
                         varmin = rhomin ; varmax = rhomax;
                         break ;
                     }
                     case 4: {
                         varmin = amin ; varmax = amax;
                         break ;
                     }
                     case 5: {
                         varmin = qmin ; varmax = qmax;
                         break ;
                     }
                     case 6: {
                         varmin = wmin ; varmax = wmax;
                         break ;
                     }
                     case 7: {
                         varmin = mumin ; varmax = mumax;
                         break ;
                     }
                     case 8: {
                         varmin = numin ; varmax = numax;
                         break ;
                     }
                  }
                  V1 = Double.valueOf(f1.getText()) ;
                  v1 = V1.doubleValue() ;
                  if(v1 < varmin) {
                      v1 = varmin ;
                      fl1 = (float) v1 ;
                      f1.setText(String.valueOf(fl1)) ;
                  }
                  if(v1 > varmax) {
                      v1 = varmax ;
                      fl1 = (float) v1 ;
                      f1.setText(String.valueOf(fl1)) ;
                  }
                  inptval = v1 ;
           // set sub flag  -- if applicable
                  if (inopt == 0 && inptval <= 1.0) {
                      isub = 1 ;
                      sub.select(1) ;
                  }
                  if (inopt == 0 && inptval >= 1.0) {
                      isub = 0 ;
                      sub.select(0) ;
                  }
                  if (inopt == 1 && inptval >= .5283) {
                      isub = 1 ;
                      sub.select(1) ;
                  }
                  if (inopt == 1 && inptval <= .5283) {
                      isub = 0 ;
                      sub.select(0) ;
                  }
                  if (inopt == 2 && inptval >= .8333) {
                      isub = 1 ;
                      sub.select(1) ;
                  }
                  if (inopt == 2 && inptval <= .8333) {
                      isub = 0 ;
                      sub.select(0) ;
                  }
                  if (inopt == 3 && inptval >= .6339) {
                      isub = 1 ;
                      sub.select(1) ;
                  }
                  if (inopt == 3 && inptval <= .6339) {
                      isub = 0 ;
                      sub.select(0) ;
                  }
 
                  comPute () ;

                  return true ;
                }
                else return false ;
             }
        } // end Up

        class Dn extends Panel {
             Isentrop outerparent ;
             Label l1 ;
             TextField o1,o2,o3,o4,o5,o6,o7,o8,o9 ;
             Label lo1,lo2,lo3,lo4,lo5,lo6,lo7,lo8,lo9 ;

             Dn (Isentrop target) {

                outerparent = target ;
                setLayout(new GridLayout(4,6,5,5)) ;

                l1 = new Label("Output", Label.LEFT) ;
                l1.setForeground(Color.blue) ;
 
                lo1 = new Label("Mach", Label.RIGHT) ;
                o1 = new TextField() ;
                o1.setBackground(Color.black) ;
                o1.setForeground(Color.yellow) ;

                lo2 = new Label("p/pt", Label.RIGHT) ;
                o2 = new TextField() ;
                o2.setBackground(Color.black) ;
                o2.setForeground(Color.yellow) ;

                lo3 = new Label("T/Tt", Label.RIGHT) ;
                o3 = new TextField() ;
                o3.setBackground(Color.black) ;
                o3.setForeground(Color.yellow) ;

                lo4 = new Label("rho/rhot", Label.RIGHT) ;
                o4 = new TextField() ;
                o4.setBackground(Color.black) ;
                o4.setForeground(Color.yellow) ;

                lo5 = new Label("A/A*", Label.RIGHT) ;
                o5 = new TextField() ;
                o5.setBackground(Color.black) ;
                o5.setForeground(Color.yellow) ;

                lo6 = new Label("q/p", Label.RIGHT) ;
                o6 = new TextField() ;
                o6.setBackground(Color.black) ;
                o6.setForeground(Color.yellow) ;

                lo7 = new Label("Wcor/A", Label.RIGHT) ;
                o7 = new TextField() ;
                o7.setBackground(Color.black) ;
                o7.setForeground(Color.yellow) ;

                lo8 = new Label("Mach Ang", Label.RIGHT) ;
                o8 = new TextField() ;
                o8.setBackground(Color.black) ;
                o8.setForeground(Color.yellow) ;

                lo9 = new Label("P-M Ang", Label.RIGHT) ;
                o9 = new TextField() ;
                o9.setBackground(Color.black) ;
                o9.setForeground(Color.yellow) ;

                add(l1) ;
                add(new Label(" ", Label.LEFT)) ;
                add(new Label(" ", Label.LEFT)) ;
                add(new Label(" ", Label.LEFT)) ;
                add(new Label(" ", Label.LEFT)) ;
                add(new Label(" ", Label.LEFT)) ;

                add(lo1) ;
                add(o1) ;
                add(lo2) ;
                add(o2) ;
                add(lo6) ;
                add(o6) ;

                add(lo8) ;
                add(o8) ;
                add(lo3) ;
                add(o3) ;
                add(lo5) ;
                add(o5) ;

                add(lo9) ;
                add(o9) ;
                add(lo4) ;
                add(o4) ;
                add(lo7) ;
                add(o7) ;
             }
        } // end Dn
   } // end In
}

