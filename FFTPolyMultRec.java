/**
 * Name: Rohan D. Shah A01943549
 * Language: JAVA
 * IDE: NetBeans IDE 7.4
 * Program: Assignment 4 CS 5050 Polynomial Multiplication using Fast Fourier transform.
 * 
 * Description:
 * This program takes two polynomials, computes the fourier transforms of the two polynomials, multiplies point
 * to point and then takes the inverse transform of the multiplied array, to get the actual multiplication answer.
 * This is faster than other approaches and uses complex numbers to compute the fourier transforms.
 * The fourier transform in this file is a recursive function.
 * 
 * Algorithm:
 * 1) Take input two polynomials
 * 2) pad the polynomials with zero upto places double of the number of coefficients
 * 3) take fourier transforms of the polynomials
 *   3.1) split the array into two arrays with even indices in one array and odd indices in the other
 *   3.2) Multiply with the 8th root of unity by various powers of omega depending on the level
 *   3.3) Combine the answer by the formula, A = Aeven + X.Aodd
 * 4) multiply the The two transforms point by point
 * 5) take inverse transforms of the multiplication by multiplying with omega inverse
 * 6) reshuffle the values by dividing the values by 2*number of coefficients
 * 
 * Input: 
 * P = 0,1,2,3
 * Q = 10,11,12,13
 * 
 * Output:
 * P*Q = 0,10,30.9999999999,63.99999999999,70,62,39,0
 * 
 */

/**
 *
 * @author Rohan D. Shah
 */

package fftpolymultrec;

import java.text.DecimalFormat;
import java.util.Random;
import java.lang.*;

public class FFTPolyMultRec
{

    /**
     * @param args the command line arguments
     */
    // order of the polynomials
    static int noOfCoefficients = 65536;
    
    // storing first polynomial, no of coefficients = order+1 hence size = order + 1
    static Complex[] p = new Complex[2*noOfCoefficients];
    
    // storing second polynomial, no of coefficients = order+1 hence size = order + 1
    static Complex[] q = new Complex[2*noOfCoefficients];
    
    static Complex [] omega = new Complex[2*noOfCoefficients];
            //new Complex(Math.cos(Math.PI  / (noOfCoefficients)), Math.sin(Math.PI / (noOfCoefficients)));
    static Complex [] omegaInverse = new Complex[2*noOfCoefficients];
            //new Complex(Math.cos(Math.PI  / (noOfCoefficients)), -(Math.sin(Math.PI  / (noOfCoefficients))));
    
    // splits the polynomial into two arrays based on odd and even indices
    public static Complex[] split(Complex[] a, int start)
    {    
        Complex[] splitted = new Complex[a.length/2];
        int i = 0;
        while(start < a.length)
        {
            splitted[i] = a[start];
            start+=2;
            i++;
        }
      
        return splitted;
    }
    // The function that computes the fourier transform
    public static Complex[] FFTmult(Complex[] a, int length, int pow)
    {
        if(length == 1)
        {
            return a;
        }
        
        Complex [] aEven = split(a, 0);
        Complex [] aOdd = split(a, 1);
        
        // recursion
        Complex [] fEven = FFTmult(aEven, length/2, 2*pow);
        Complex [] fOdd = FFTmult(aOdd, length/2, 2*pow);
        Complex [] f = new Complex[length];
        // combining the solution       
        for(int i=0 ; i<length/2 ; i++)
        {
            Complex temp = omega[i*pow].Mult(fOdd[i]);
            f[i] = fEven[i].Add(temp);
            f[i+(length/2)] = fEven[i].Sub(temp);
           
        }
        return f;
    }
    public static void main(String[] args) 
    {
        double forDistribution; 
        double temp;
        //assign sizes and values randomly
        Random randomGenerator = new Random();
        
         //format to 5 decimal place
        DecimalFormat oneDigit = new DecimalFormat("#.#####");
        // Fill in size and values of items with random numbers between -1 and 1 inclusive
        for(int i=0 ; i<noOfCoefficients ; i++)
        {   
            // the nextdouble method gives values from 0 to 1 but we want to have values from -1 to 1
            // hence we multiply by 2 and then subtract 1 from it to get values between that range
            forDistribution = randomGenerator.nextDouble() * 2 - 1;
            temp = Double.valueOf(oneDigit.format(forDistribution));
            p[i] = new Complex(temp,0);

            forDistribution = randomGenerator.nextDouble() * 2 - 1;            
            temp = Double.valueOf(oneDigit.format(forDistribution));
            q[i] = new Complex(temp,0);
        }
        // padding with zeros
        for(int t = noOfCoefficients ; t<2*noOfCoefficients ; t++)
        {
            p[t] = new Complex(0, 0);
            q[t] = new Complex(0, 0);
        }
        /*System.out.println("------------------------------------");
        System.out.println("The first polynomial=");
        System.out.println("------------------------------------");
        for(int a=0 ; a<2*noOfCoefficients ; a++)
            System.out.println(p[a].toString());
        
        
        System.out.println("------------------------------------");
         System.out.println("The second polynomial=");
         System.out.println("------------------------------------");
         for(int a=0 ; a<2*noOfCoefficients ; a++)
            System.out.println(q[a].toString());        
        System.out.println("------------------------------------");
        
        // stores the complex 8th root of unity
        */for(int a=0 ; a<omega.length ; a++)
        {
            omega[a] = new Complex(Math.cos(Math.PI*a/noOfCoefficients), Math.sin(Math.PI*a/noOfCoefficients));
        }
        
        // stores the inverse of the complex 8th root of unity
        for(int a=0 ; a<omega.length ; a++)
        {
            omegaInverse[a] = new Complex(omega[a].dReal, -omega[a].dImaginary);
        }
        
        //Computing the transforms
        Complex[] pFFT = FFTmult(p, p.length, 1);
        Complex[] qFFT = FFTmult(q, q.length, 1);
        /*System.out.println("===========PFFT==========");
        for(int t=0;t<pFFT.length;t++)
        {
            System.out.println((pFFT[t].toString()));
        }
        System.out.println("=========QFFT============");
        
        for(int t=0;t<qFFT.length;t++)
        {
            System.out.println((qFFT[t].toString()));
        }*/
        // stores the point by point multiplication values of the transforms 
        Complex[] afterMult = new Complex[2*noOfCoefficients];
        
        // point by point multiplication
        for(int t=0 ; t<2*noOfCoefficients ; t++)
            afterMult[t] = pFFT[t].Mult(qFFT[t]);        
        
       // using the same array for the inverse values
        for(int a=0 ; a<omega.length ; a++)
        {
            omega[a] = new Complex(Math.cos(Math.PI*a/noOfCoefficients), -Math.sin(Math.PI*a/noOfCoefficients));
        }
        Complex[] inverseFFT = FFTmult(afterMult, afterMult.length, 1);
        Complex[] answer = new Complex[inverseFFT.length];
        System.out.println("actual output");
        for(int t=0;t<answer.length;t++)
        {
            System.out.println((inverseFFT[t].dReal/p.length)+" X^"+t);
        }
        
    }    
}
