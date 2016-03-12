package teaspoon.adaptation;

import java.util.Random;
import java.util.Calendar;

public class Samplers {
  private static Random rng = new Random(
      Calendar.getInstance().getTimeInMillis() +
      Thread.currentThread().getId());
  double[] p;
  
  public Samplers(double[] observations){
	  this.p=observations;
  }
  
  
  public static double sampleGamma(double k, double theta) {
    boolean accept = false;
    if (k < 1) {
 // Weibull algorithm
 double c = (1 / k);
 double d = ((1 - k) * Math.pow(k, (k / (1 - k))));
 double u, v, z, e, x;
 do {
  u = rng.nextDouble();
  v = rng.nextDouble();
  z = -Math.log(u);
  e = -Math.log(v);
  x = Math.pow(z, c);
  if ((z + e) >= (d + x)) {
   accept = true;
  }
 } while (!accept);
 return (x * theta);
    } else {
 // Cheng's algorithm
 double b = (k - Math.log(4));
 double c = (k + Math.sqrt(2 * k - 1));
 double lam = Math.sqrt(2 * k - 1);
 double cheng = (1 + Math.log(4.5));
 double u, v, x, y, z, r;
 do {
  u = rng.nextDouble();
  v = rng.nextDouble();
  y = ((1 / lam) * Math.log(v / (1 - v)));
  x = (k * Math.exp(y));
  z = (u * v * v);
  r = (b + (c * y) - x);
  if ((r >= ((4.5 * z) - cheng)) ||
                    (r >= Math.log(z))) {
   accept = true;
  }
 } while (!accept);
 return (x * theta);
    }
  }

  public double[] Dirichlet(){
	  double total = 0;
	  double[] D = new double[4];
	   D[0] = sampleGamma(p[0],1);total+=D[0];
	   D[1] = sampleGamma(p[1],1);total+=D[1];
	   D[2] = sampleGamma(p[2],1);total+=D[2];
	   D[3] = sampleGamma(p[3],1);total+=D[3];
	   
	   D[0]=D[0]/total;
	   D[1]=D[1]/total;
	   D[2]=D[2]/total;
	   D[3]=D[3]/total;
	   return D;  
	  
  }

 public double Beta() {

     double total = 0;
     double[] D = new double[2];

     D[0] = sampleGamma(p[0],1);total+=D[0];
     D[1] = sampleGamma(p[1],1);total+=D[1];

     D[0] = D[0]/total;
     D[1] = D[1]/total;

     return D[0];

 }
}