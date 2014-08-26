/*
 *  =================================================================
 *  Complex.java: A library of methods for complex number arithmetic
 *  The algoritms for arithmetic are taken from the text "Numerical
 *  Recipies in C."
 *
 *  Written By : Mark Austin                                 May 1998
 *  =================================================================
 *  source = http://www.isr.umd.edu/~austin/ence688r.d/java-basics/appl-complex1.d/Complex.java
 */
package fftpolymultrec;
import java.lang.Math;
import java.text.DecimalFormat;

public class Complex {
   protected double dReal, dImaginary;

   // Constructor methods ....

   public Complex() {}

   public Complex( double dReal, double dImaginary ) {
      this.dReal      = dReal;
      this.dImaginary = dImaginary;
   }

   // Convert complex number to a string ...

   public String toString() {
      if (dImaginary >= 0)
          return dReal + "+" +  dImaginary + "i";
      else
          return dReal + "-" + -dImaginary + "i";
   }

   // ==============================
   // Complex number arithmetic ...
   // ==============================

   // Complex negative of complex number ...

   public Complex Negate() {
   Complex negate = new Complex();

       negate.dReal      = - dReal;
       negate.dImaginary = - dImaginary;

       return (negate);
   }

   // Compute sum of two complex numbers cA + cB.....

   public Complex Add( Complex cB ) {
   Complex sum = new Complex();

       sum.dReal      = dReal      + cB.dReal;
       sum.dImaginary = dImaginary + cB.dImaginary;

       return (sum);
   }

   // Compute difference of two complex numbers cA - cB.....

   public Complex Sub( Complex cB ) {
   Complex diff = new Complex();

       diff.dReal      = dReal      - cB.dReal;
       diff.dImaginary = dImaginary - cB.dImaginary;

       return (diff);
   }

   // Compute product of two complex numbers cA * cB.....

   public Complex Mult( Complex cB ) {
   Complex prod = new Complex();
   DecimalFormat oneDigit = new DecimalFormat("#.######");
       prod.dReal      = dReal*cB.dReal - dImaginary*cB.dImaginary;
       prod.dImaginary = dImaginary*cB.dReal + dReal*cB.dImaginary;
       //System.out.println("inside mult");
       //System.out.println(this.toString()+" * "+ cB.toString()+" = "+prod.toString());
       return (prod);
   }
   // Author of this function is Rohan D. Shah and not Mark Austin
   public Complex threeMult(Complex cB)
   {
       Complex prod = new Complex();
       prod.dReal = dReal*cB.dReal - dImaginary*cB.dImaginary;
       double aPlusB = dReal+dImaginary;
       double cPlusD = cB.dImaginary+cB.dReal;
       
       prod.dImaginary = (aPlusB*cPlusD) - (dReal*cB.dReal)-(dImaginary*cB.dImaginary);
       return prod;
       
   }
   // Compute divisor of two complex numbers cA / cB.....

   public Complex Div( Complex cB ) {
   Complex div = new Complex();
   double dR, dDen;
 
       if(Math.abs( cB.dReal ) >= Math.abs( cB.dImaginary )) {
          dR   = cB.dImaginary/cB.dReal;
          dDen = cB.dReal + dR*cB.dImaginary;
          div.dReal      = (dReal      + dR*dImaginary)/dDen;
          div.dImaginary = (dImaginary - dR*dReal)/dDen;
       } else {
          dR   = cB.dReal/cB.dImaginary;
          dDen = cB.dImaginary + dR*cB.dReal;
          div.dReal      = (dR*dReal      + dImaginary)/dDen;
          div.dImaginary = (dR*dImaginary - dReal)/dDen;
       }

       return (div);
   }

   // Scale complex number by double precision number.....

   public Complex Scale( double dFactor ) {
   Complex scale = new Complex();

      scale.dReal      = dFactor*dReal;
      scale.dImaginary = dFactor*dImaginary;

      return (scale);
   }

   // Compute complex number conjugate....

   public Complex Conjugate() {
   Complex conj = new Complex();

      conj.dReal      =  dReal;
      conj.dImaginary = -dImaginary;

      return (conj);
   }

   // Compute absolute value of complex number ....

   public double Abs() {
   double dX, dY, dTemp, dAnswer;

      dX = Math.abs( dReal );
      dY = Math.abs( dImaginary );
      if(dX == 0)
         dAnswer = dY;
      else if(dY == 0)
         dAnswer = dX;
      else if(dX > dY) {
         dTemp = dY/dX;
         dAnswer = dX*Math.sqrt((double) (1.0 + dTemp*dTemp));
      } else {
         dTemp = dX/dY;
         dAnswer = dY*Math.sqrt((double) (1.0 + dTemp*dTemp));
      } 

      return ( dAnswer );
   }

   // Compute square root of complex number ....

   public Complex Sqrt() {
   Complex csqrt = new Complex();
   double dX, dY, dW, dR;

      if((dReal == 0) && (dImaginary == 0.0)) {
          csqrt.dReal      = 0.0;
          csqrt.dImaginary = 0.0;
          return (csqrt);
      } 

      dX = Math.abs(dReal);
      dY = Math.abs(dImaginary);

      if( dX >= dY ) {
          dR = dY/dX;
          dW = Math.sqrt(dX)*Math.sqrt(0.5*(1.0 + Math.sqrt(1+dR*dR)));
      } else {
          dR = dX/dY;
          dW = Math.sqrt(dY)*Math.sqrt(0.5*(dR + Math.sqrt(1+dR*dR)));
      }

      if(dReal >= 0.0) {
         csqrt.dReal = dW;
         csqrt.dImaginary = dImaginary/(2.0*dW);
      } else {
         csqrt.dImaginary = (dImaginary > 0.0) ? dW : -dW;
         csqrt.dReal      =  dImaginary/(2.0*csqrt.dImaginary);
      }

      return (csqrt);
   }
}