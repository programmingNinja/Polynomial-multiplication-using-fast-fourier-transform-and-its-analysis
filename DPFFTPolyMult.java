
package fftpolymultrec;
import java.text.DecimalFormat;
import java.util.Random;

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
 * The fourier transform in this file is not a recursive function but dynamic programming and reverse bit shuffling
 * is used to save space.
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
 * 6) reshuffle the values by dividing the values by 2*number of coeeficients
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
 * @author Rohan
 */
public class DPFFTPolyMult 
{
    /**
     * @param args the command line arguments
     */
    // order of the polynomials
    static int noOfCoefficients = 4;
    
    // storing first polynomial, no of coefficients = order+1 hence size = order + 1
    static Complex[] p = new Complex[2*noOfCoefficients];
    
    // storing second polynomial, no of coefficients = order+1 hence size = order + 1
    static Complex[] q = new Complex[2*noOfCoefficients];
    
    static Complex [] omega = new Complex[2*noOfCoefficients];
                
    static int[] rbsArray = new int[2*noOfCoefficients];
            //new Complex(Math.cos(Math.PI  / (noOfCoefficients)), -(Math.sin(Math.PI  / (noOfCoefficients))));
    static int logN =  (int)( Math.log(2*noOfCoefficients)/Math.log(2));
    // initializing the array which stores the reverse bit shuffled values
    public static void RBSinit()
    {
        for(int a=0;a<rbsArray.length;a++)
        {
            rbsArray[a] = rbs(a, logN);
        }
    }
    // function that shuffles the bits
    public static int rbs(int i, int k)
    {
        if(k==0)
            return i;
        if(i%2 == 1)
            return (int) (Math.pow(2,k-1)+rbs(i/2,k-1));
        else
            return rbs(i/2,k-1);
    }
    // FFT by dynamic programming
    public static Complex[] DPFFTmult(Complex[] a, int length)
    {          
        Complex[][] sol = new Complex[logN+1][length];
        
        for(int q=0;q<length;q++)
        {
            sol[0][rbsArray[q]] = new Complex(a[q].dReal, a[q].dImaginary);
        }
        
        // going up each level
        int pow = length/2;
        int size = 2;
        Complex odd = new Complex();
        for(int q=1 ; q<=logN ; q++)
        {
            for(int w=0 ; w<length ; w+=size)
            {
                for(int e=0 ; e<size/2 ; e++)
                {
                    odd = omega[e*pow].threeMult(sol[q-1][w+e+size/2]);
                    sol[q][e+w]=sol[q-1][e+w].Add(odd);
                    sol[q][w+e+size/2]=sol[q-1][e+w].Sub(odd);
                }
            }
            pow=pow/2;
            size=size*2;
        }
        // sol[logn] contains the answer
        return sol[logN];        
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
       /* System.out.println("------------------------------------");
        System.out.println("The first polynomial=");
        System.out.println("------------------------------------");
        for(int a=0 ; a<2*noOfCoefficients ; a++)
            System.out.println(p[a].toString());
        
        
        System.out.println("------------------------------------");
         System.out.println("The second polynomial=");
         System.out.println("------------------------------------");
         for(int a=0 ; a<2*noOfCoefficients ; a++)
            System.out.println(q[a].toString());        
        System.out.println("------------------------------------");*/
        // initializing the omega array with complex 8th root of unity
        for(int a=0 ; a<omega.length ; a++)
        {
            omega[a] = new Complex(Math.cos(Math.PI*a/noOfCoefficients), Math.sin(Math.PI*a/noOfCoefficients));
        }
        
        // intializing the array which stores the reverse bits
        RBSinit();
        
        // computing the fourier transforms
        Complex[] pFFT = DPFFTmult(p, p.length);
        Complex[] qFFT = DPFFTmult(q, q.length);
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
        
        // using threeMult function. you can change to Mult(Complex cB) for four multiplication function
        for(int t=0 ; t<2*noOfCoefficients ; t++)
            afterMult[t] = pFFT[t].threeMult(qFFT[t]);        
        
        // using the same array for the inverse values
        for(int a=0 ; a<omega.length ; a++)
        {
            omega[a] = new Complex(Math.cos(Math.PI*a/noOfCoefficients), -Math.sin(Math.PI*a/noOfCoefficients));
        }
        
        Complex[] inverseFFT = DPFFTmult(afterMult, afterMult.length);
        
        System.out.println("actual output");
        for(int t=0;t<inverseFFT.length;t++)
        {
            System.out.println((inverseFFT[t].dReal/p.length)+" X^"+t);
        }
        
    }    
}
