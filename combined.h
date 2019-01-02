#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <conio.h>
#include <time.h>

#define	D_max 200  // Max number of dimensions of the search space
#define	S_max 200 // Max swarm size
#define R_max 200 // Max number of runs

int D;
int LocalLimit; // Spider Monkey
int GlobalLimit; // Spider Monkey
int limit;// ABC
double Population[S_max][D_max]; 
double fun_val[S_max ];  
double fitness[S_max ]; 
double prob[S_max ]; 
double new_position[D_max]; 
double ObjValSol; 
double FitnessSol; 
int neighbour, param2change; 
double GlobalMin; 
double GlobalLeaderPosition[D_max]; 
double LocalMin[S_max /2]; 
double LocalLeaderPosition[S_max /2][D_max]; 
int LocalLimitCount[S_max /2];
double GlobalMins[R_max]; 
int GlobalLimitCount;
int gpoint[S_max ][2];
double r,r1,r2; 
FILE * f_run,*f_run1;
int Pr;
double part;
double acc_err;
double lb[D_max],ub[D_max];
double feval;
double mean_feval,total_feval,mean_error;
int lo,group_size,g,iter;
double obj_val;
double cr;
double Foods[S_max ][D_max]; /*Foods is the population of food sources. Each row of Foods matrix is a vector holding D parameters to be optimized. The number of rows of Foods matrix equals to the FoodNumber*/
double f[S_max ];  /*f is a vector holding objective function values associated with food sources */
double fitness_abc[S_max ]; /*fitness is a vector holding fitness (quality) values associated with food sources*/
double trial[S_max ]; /*trial is a vector holding trial numbers through which solutions can not be improved*/
double prob_abc[S_max ]; /*prob is a vector holding probabilities of food sources (solutions) to be chosen*/
double solution [D_max]; /*New solution (neighbour) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomlu chosen solution different from i*/
double GlobalParams[D_max]; /*Parameters of the optimum solution*/
double GlobalMins_abc[R_max]; /*GlobalMins holds the GlobalMin of each run in multiple runs*/
//double pi;
//double E;
long double E;			// exp(1). Useful for some test functions
struct landscape funct; // When (x,f(x)) is read from a file
//int iter;
struct archive memPos;
double nEval;			// Number of fitness evaluations
long double pi;			// Useful for some test functions
int run;

double fun(double sol[]);
double distance (double x1[], double x2[],double L) 
{  // Distance between two positions
	// L = 2 => Euclidean	
	int d;     
	double n;

	n = 0;

	for (d = 0; d < 3; d++)
		n = n + pow (fabs(x1[d] - x2[d]), L);

	n = pow (n, 1/L);
	return n;    
}
int sign (double x) 
{    // Be careful: not the classical one 
	if (x <= 0)	return -1;    
	return 1;   
}
//************* optimization fun *************//(-32,32)
double fun(double x[])
   {
     feval++;
      
       /////////////////////////////////
        double c;
    int d, d1;
    int i, j, k;
    double f=0, f1, p, xd, x1, x2, x3, x4;
    double min;
    double sum1, sum2;
    double t0, tt, t1;
     // For polynomial fitting problem
    int const M = 60;
    double py, y = -1, dx = ( double )M;

    switch (Pr)
    {
      case 0: // Parabola (Sphere)
		  {
        f = 0;
        p = 0; // Shift
        for ( d = 0; d < D; d++ )
        {
          xd = x[d] - p;
          f = f + xd * xd;
        }
      break;
		  }
       /* case 1: // De Jong's f4
		  { f = 0;
        p = 0; // Shift
        for ( d = 0; d < D; d++ )
        {
          xd = x[d] - p;
          f = f + (d+1)*xd*xd*xd*xd;
        }
		break;}


    case 2: // Griewank
		  { f = 0;
        p = 1;
        for ( d = 0; d < D; d++ )
        {
          xd = x[d];
          f = f + xd * xd;
          p = p * cos( xd / sqrt( d + 1 ) );
        }
        f = f / 4000 - p + 1;
      break;
		  }
      case 3: // Rosenbrock
		  { f = 0;
        t0 = x[0];
        for ( d = 1; d < D; d++ )
        {
          t1 = x[d];
          tt = 1 - t0;
          f += tt * tt;
          tt = t1 - t0 * t0;
          f += 100 * tt * tt;
          t0 = t1;
        }
      break;
		  }*/
      case 1: // Rastrigin. Minimum value 0. Solution (0,0 ...0) (in list 4)
		  {k = 10;
        f = 0;
        for ( d = 0; d < D; d++ )
        {
          xd = x[d];
          f += xd * xd - k * cos( 2 * pi * xd );
        }
        f += D * k;
      break;
		  }
      /*case 5: // Ackley
		  {sum1 = 0;
        sum2 = 0;
        for ( d = 0; d < D; d++ )
        {
          xd = x[d];
          sum1 += xd * xd;
          sum2 += cos( 2 * pi * xd );
        }
        y = D;
        f = ( -20 * exp( -0.2 * sqrt( sum1 / y ) ) - exp( sum2 / y ) + 20 + E );
      break;
		  }
   //===========================================================
  case 6: //Alpine
 	 {


	f=0;
	for(d=0;d<D;d++)
 f+=fabs(x[d]*sin(x[d])+0.1*x[d]);
  break;
}
case 7: //Michalewicz
{
	/* f=0;
for(d=0;d<n;d++)
f+=Math.pow((x[d]-d-1),2);
f=0;
for(d=0;d<D;d++)
	f+=-sin(x[d])*pow(sin((d+1)*pow(x[d],2)/pi),20);


 break;
 }
 case 8:  //Cosine Mixture [-1,1] f(0,0,...0)=-D*0.1

		 {
			 double s1=0;                                  
        for(d=0;d<D;d++)
		{ s1+=x[d]*x[d];}
        double p1=0;
        for(d=0;d < D;d++)
		{ p1+=cos(5*pi*x[d]);}
    f=s1-(0.1)*p1;
	break;
		 
		 }
 case 9: //Exponential [-1,1] f(0,0,...0)= -1
		 { double s2=0;                  
        for(d=0;d < D;d++)
		{ s2+=x[d]*x[d];}
         f = -(exp(-0.5*s2));
		 break; }



 case 10:     //Zakharov's [-5.12,5.12] f(0,0,0,...,0)=0

		 {
double s1=0;                             
double s2=0;
for(d=0;d < D;d++)
s1+=x[d]*x[d];
for(d=0;d < D;d++)
s2+=((d+1)*x[d])/2;
f = s1+s2*s2+s2*s2*s2*s2;		 
		 break;
		 }*/

 case 2: //Cigar [-10,10]  f(0,0,0,...,0)=0

		 {
double s=0;          //Cigar (madam paper)
for(d=1;d < D;d++)
s+=x[d]*x[d];
f=x[0]*x[0]+100000*s;

break;
}

/* case 12:  //brown3 [-1,4]  f(0,0,0.....,0)=0
		 {
		 
		 f=0;                  
       for(d=0;d < D-1;d++)
   f+=pow(x[d]*x[d],x[d+1]*x[d+1]+1)+pow(x[d+1]*x[d+1],x[d]*x[d]+1);		 
		 	   
    break;
		 }*/



case 3:  //Schewel prob 2.22
{
double s1=0;                                    
double s2=1;
for(d=0;d < D;d++)
s1+=fabs(x[d]);
for(d=0;d < D;d++)
s2=s2*fabs(x[d]);
f=s1+s2;
break;
}
/*case 14: //Salomon Problem (SAL)
{
f=0;  
double s1 =0,s2;   
for(d=0;d<D;d++)
s1+=pow(x[d],2);
s2=sqrt(s1);
f=1-cos(2*3.14*s2)+0.1*s2;
		 
   break;
}*/
case 4: // Axis parallel hyperellipsoid
{
        f = 0;
        for ( d = 0; d < D; d++ )
        {
           f = f + (d+1)*x[d] * x[d];
        }
      break;
 }
 /*case 16: // Pathological
{
        f = 0;
        for ( d = 0; d < D-1; d++ )
        {
           f = f + (0.5+(pow(sin(sqrt(100*x[d]*x[d]+x[d+1]*x[d+1])),2)-0.5)/(1+0.001*pow((x[d]*x[d]-2*x[d]*x[d+1]+x[d+1]*x[d+1]),2)));
        }
      break;
 }
 case 17: //Sum of different powers
 {

f=0;
for(d=0;d < D;d++){
f+=pow(fabs(x[d]),d+1);

}
break;
} */
case 5: //step function [-100, 100] f(-0.5<=x<=0.5)=0
 {

f=0;
for(d=0;d < D;d++){
f+=pow(floor(x[d]+0.5),2);

}
break;
}   
/*case 19: //Quartic function, i.e., noise [-1.28, 1.28] f(0000..00)=0
 {

f=0;
for(d=0;d < D;d++){
f+=(d+1)*pow(x[d],4)+ ((double)rand()/((double)(RAND_MAX)+(double)(1)) );

}
break;
}  
case 20: //Inverted cosine wave function (Masters) //Inverted cosine wave function (Masters) [-5, 5] f(000..0)=-D+1
{
double s1=0,s2=0;      
f=0;            
for(d=0;d < (D-1);d++)
{ 
          f+=-(exp(-(x[d]*x[d]+x[d+1]*x[d+1]+0.5*x[d]*x[d+1])/8.0)*cos(4*sqrt(x[d]*x[d]+x[d+1]*x[d+1]+0.5*x[d]*x[d+1])));
}

        break; 
}
case 21: //Neumaier 3 Problem (NF3) (Neumaier, 2003b)
 {

f=0;
double s1=0,s2=0;
for(d=0;d<D;d++)
   s1+=pow((x[d]-1),2);
for(d=1;d<D;d++)
    s2+=x[d]*x[d-1];
f=s1-s2;

break;
}  
case 22: //Rotated hyper-ellipsoid function
 {

f=0;
int j=0;
for(d=0;d<D;d++)
{
for(j=0;j<=d;j++)
{
    f+=x[j]*x[j];             
}                
}

break;
}  */
case 6: //Levy montalvo 1 
 {	                                  
   
double s =0;
for(d=0;d<D-1;d++)
s+=((1+0.25*(x[d]+1))-1)*((1+0.25*(x[d]+1))-1)*(1+10*sin(3.14*(1+0.25*(x[d+1]+1)))*sin(3.14*(1+0.25*(x[d+1]+1))));
f = (3.14/D)*(10*sin(3.14*(1+0.25*(x[0]+1)))*sin(3.14*(1+0.25*(x[0]+1)))+s+((1+0.25*(x[D-1]+1))-1)*((1+0.25*(x[D-1]+1))-1));
		 break;
}
  //===========================================================
 case 7: //Levy montalvo 2 
{
		                                                        
double s =0;
for(d=0;d<D-1;d++)
          s+=(x[d]-1)*(x[d]-1)*(1+sin(3*3.14*x[d+1])*sin(3*3.14*x[d+1]));
         f = 0.1*(sin(3*3.14*x[0])*sin(3*3.14*x[0])
            	+s+(x[D-1]-1)*(x[D-1]-1)*(1+sin(2*3.14*x[D-1])*sin(2*3.14*x[D-1])));
 		 break;
		 }

  /* case 25: //Ellipsoidal Ellipsoidal [-D,D] f(1,2,3,...,D)=0
	   {
		 f=0;     
for(d=0;d<D;d++)
f+=pow((x[d]-d-1),2);
		 
        break;
	   }*/
	   case 8:		//Beale function [-4.5,4.5] f(3, 0.5)=0
			
       {
           	x1=x[0];
            x2=x[1];	
	
        
            f = pow((1.5-x1*(1-x2)),2)+pow((2.25-x1*(1-x2*x2)),2)+pow((2.625-x1*(1-x2*x2*x2)),2);
           
        
          
		break;
    }
     /*case 27:		//Colville function [-10,10] f(1111)=0
			
       {
           	x1=x[0];
            x2=x[1];	
	        x3=x[2];
            x4=x[3];
        
            f = 100*pow((x2-x1*x1),2)+pow((1-x1),2)+90*pow((x4-x3*x3),2)+pow((1-x3),2)+10.1*(pow((x2-1),2)+pow((x4-1),2))+19.8*(x2-1)*(x4-1);
           
        
          
		break;
    }*/
    case 9: //Branins\92s function [-5,10][0,15] f(-pi, 12.275)=0.3979		
			
       {
           	double a,b,c,d,e,fav;
           	a=1;
           	b=5.1/(4*pi*pi);
           	c=5/pi;
           	d=6;
           	e=10;
           	fav=1/(8*pi);
            x1=x[0];
            x2=x[1];
            
            f = a*pow((x2-b*x1*x1+c*x1-d),2)+e*(1-fav)*cos(x1)+e;
           
        
          
		break;
    }
     case 26: //Kowalik function [-5,5] f(0.192833, 0.190836, 0.123117, 0.135766)=0.000307486
			
       {
            x1=x[0];
            x2=x[1];	
	        x3=x[2];
            x4=x[3];
           	double a[]={0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627, 0.0456, 0.0342, 0.0323, 0.0235,0.0246};
           	double b[]={4.0, 2.0,1.0, 0.5, 0.25, 0.1667, 0.125, 0.1, 0.0833,  0.07143, 0.0625};
            for(i=0;i<11;i++)
            {
                f = f+pow((a[i]-(x1*(b[i]*b[i]+b[i]*x2))/(b[i]*b[i]+b[i]*x3+x4)),2);
            }
           
        
          
		break;
    }
    // Minibench Mark Problems: 4-Problem Taken from clerc website
   /* case 30:		//// 2D Tripod function [-100,100] f(0, 50)=-50
       {
           	double s11, s12, s21, s22;		
			x1 = x[0] ; 
		x2 = x[1];  
		s11 = (1.0 - sign (x1)) / 2;
		s12 = (1.0 + sign (x1)) / 2; 
		s21 = (1.0 - sign (x2)) / 2;
		s22 = (1.0 + sign (x2)) / 2;

		//f = s21 * (fabs (x1) - x2); // Solution on (0,0)
		f = s21 * (fabs (x1) +fabs(x2+50)); // Solution on (0,-50)  
		f = f + s22 * (s11 * (1 + fabs (x1 + 50) + fabs (x2 - 50)) 
		               + s12 * (2 +fabs (x1 - 50) + fabs (x2 - 50)));	  
		break;
    }
     case 31:	//Shifted CEC 2005 Rosenbrock F6  [-100, 100], solution point is O + (1; : : : ; 1) where f = 390.
	 {
        double o[]={81.0232,-48.395,19.2316,-2.5231,70.4338,47.1774,-7.8358,-86.6693,57.8532,-9.9533};
        
        for (i=1;i<D;i++)
        {
            double zd1=x[i-1] - o[i-1];
            double zd2=x[i] - o[i];
            f = f + ( 100 *(zd1*zd1 - zd2)*(zd1*zd1 - zd2) + ( zd1 - 1)*( zd1 - 1) );
        }
            f=f+390;
      	break;
    }*/
  	case 10: // /// Shifted Parabola/Sphere (CEC 2005 benchmark)	x\81[-100,100] , Global optimum: x* = offset  f(x) = f_bias = - 450
   	{
          	static double offset_0[30] =
            { 
                   -3.9311900e+001, 5.8899900e+001, -4.6322400e+001, -7.4651500e+001, -1.6799700e+001,
                   -8.0544100e+001, -1.0593500e+001, 2.4969400e+001, 8.9838400e+001, 9.1119000e+000, 
                   -1.0744300e+001, -2.7855800e+001, -1.2580600e+001, 7.5930000e+000, 7.4812700e+001,
                   6.8495900e+001, -5.3429300e+001, 7.8854400e+001, -6.8595700e+001, 6.3743200e+001, 
                   3.1347000e+001, -3.7501600e+001, 3.3892900e+001, -8.8804500e+001, -7.8771900e+001, 
                   -6.6494400e+001, 4.4197200e+001, 1.8383600e+001, 2.6521200e+001, 8.4472300e+001
            };		
		f=-450;
	    for (i=0;i<D;i++)
        {
			xd = x[i]-offset_0[i];
			f = f + xd * xd;  
		}
		break;
    }
			
   /* case 33: //Shifted CEC 2005  Rastrigin  x\81[-5,5] , Global optimum x* = offset , f(x*) = f_bias = - 330
   {
        // Shifted Rastrigin (CEC 2005)
        static double offset_3[30] =
        { 
               1.9005000e+000, -1.5644000e+000, -9.7880000e-001, -2.2536000e+000,  2.4990000e+000,
               -3.2853000e+000,  9.7590000e-001, -3.6661000e+000,  9.8500000e-002, -3.2465000e+000,
               3.8060000e+000, -2.6834000e+000, -1.3701000e+000,  4.1821000e+000,  2.4856000e+000, 
               -4.2237000e+000,  3.3653000e+000,  2.1532000e+000, -3.0929000e+000,  4.3105000e+000, 
               -2.9861000e+000,  3.4936000e+000, -2.7289000e+000, -4.1266000e+000, -2.5900000e+000, 
                1.3124000e+000, -1.7990000e+000, -1.1890000e+000, -1.0530000e-001, -3.1074000e+000
        };

	for (i=0;i<D;i++)
		{
			x[i]=x[i]-offset_3[i];
		}
		f=-330;
	    k = 10;  
	
	    for (i=0;i<D;i++) 
	    {     
			xd = x[i];	      
			f =f+ xd * xd - k * cos (2 * pi * xd);	    
		}	  
		f =f+ D * k;  
		break;	
  }*/
case 11: //Shifted CEC 2005 Schwefel [-100,100], Global optimum x* = offset , f(x*) = f_bias = - 450
  {
       // Shifted Schwefel (F2 CEC 2005)
       static double offset_4[30] =
       { 
              3.5626700e+001, -8.2912300e+001, -1.0642300e+001, -8.3581500e+001,  8.3155200e+001,
              4.7048000e+001, -8.9435900e+001, -2.7421900e+001,  7.6144800e+001, -3.9059500e+001,
              4.8885700e+001, -3.9828000e+000, -7.1924300e+001,  6.4194700e+001, -4.7733800e+001,
              -5.9896000e+000 ,-2.6282800e+001, -5.9181100e+001,  1.4602800e+001, -8.5478000e+001,
              -5.0490100e+001,  9.2400000e-001,  3.2397800e+001,  3.0238800e+001, -8.5094900e+001,
              6.0119700e+001, -3.6218300e+001, -8.5883000e+000, -5.1971000e+000,  8.1553100e+001 
       };

		for (i=0;i<D;i++)
		{
			x[i]=x[i]-offset_4[i];
		}

    f = -450;
    for (i=0;i<D;i++)
    {
        sum2 = 0.0;
        for (j=0; j<=i; j++)
        {
            sum2 += x[j];
        }
        f += sum2*sum2;
    }
		break;
  }
 case 12: //Shifted CEC 2005 Griewank. WARNING: in the CEC 2005 benchmark it is rotated
  {
       // Shifted Griewank (CEC 2005)
       static double offset_5[30] =
       { 
         -2.7626840e+002, -1.1911000e+001, -5.7878840e+002, -2.8764860e+002, -8.4385800e+001,
         -2.2867530e+002, -4.5815160e+002, -2.0221450e+002, -1.0586420e+002, -9.6489800e+001,
         -3.9574680e+002, -5.7294980e+002, -2.7036410e+002, -5.6685430e+002, -1.5242040e+002,
         -5.8838190e+002, -2.8288920e+002, -4.8888650e+002, -3.4698170e+002, -4.5304470e+002,
         -5.0658570e+002, -4.7599870e+002, -3.6204920e+002, -2.3323670e+002, -4.9198640e+002,
         -5.4408980e+002, -7.3445600e+001, -5.2690110e+002, -5.0225610e+002, -5.3723530e+002 
       };

	   sum1 = 0.0;
	   sum2 = 1.0;
	   f=-180;
	   for (i=0;i<D;i++)
		 {
		    xd=x[i]-offset_5[i];
				sum1 += xd*xd;
        sum2 *= cos(xd/sqrt(1.0+i));
 		 }
    f =f+ 1.0 + sum1/4000.0 - sum2;
		break;
  }
 case 13:  // Shifted Ackley (CEC 2005) [-32,32], Global optimum x* = offset , f(x*) = f_bias = - 140
  {
  // Shifted Ackley (CEC 2005)
  static double offset_6[30] =
  { 
    -1.6823000e+001,  1.4976900e+001,  6.1690000e+000,  9.5566000e+000,  1.9541700e+001,
    -1.7190000e+001, -1.8824800e+001,  8.5110000e-001, -1.5116200e+001,  1.0793400e+001,
    7.4091000e+000,  8.6171000e+000, -1.6564100e+001, -6.6800000e+000,  1.4543300e+001,
    7.0454000e+000, -1.8621500e+001,  1.4556100e+001, -1.1594200e+001, -1.9153100e+001,
    -4.7372000e+000,  9.2590000e-001,  1.3241200e+001, -5.2947000e+000,  1.8416000e+000,
    4.5618000e+000, -1.8890500e+001,  9.8008000e+000, -1.5426500e+001,  1.2722000e+000
  };
    f=-140;
    sum1 = 0.0;
    sum2 = 0.0;
    for (i=0;i<D;i++)
    {
        xd = x[i]-offset_6[i];
				sum1 += xd*xd;
        sum2 += cos(2.0*pi*xd);
    }
    sum1 = -0.2*sqrt(sum1/D);
    sum2 /= D;
    f = f+ 20.0 + E - 20.0*exp(sum1) - exp(sum2);
		break;
}
 
   case 14:		// Goldstein-Price function  [-2,2] f(0, -1) = 3.
			
       {
           f=0;
           x1=x[0];
           x2=x[1];
           f=(1+pow((x1+x2+1),2)*(19-14*x1+3*x1*x1-14*x2+6*x1*x2+3*x2*x2))* (30+pow((2*x1-3*x2),2)*(18-32*x1+12*x1*x1+48*x2-36*x1*x2+27*x2*x2));
           	break;
       }
       	
 	 case 15:	//Six-hump camel back function  [-5 < x < 5] f(-0.0898,0.7126) = -1.0316.
			
       {
           f=0;
           float x1=x[0];
           float x2=x[1];
           f=4*x1*x1-2.1*pow(x1,4)+(pow(x1,6)/3.0)+x1*x2-4*x2*x2+4*pow(x2,4);
              	break;
       }
        case 16:	//Easom's function  -10<=x(i)<=10, i=1:2. f(x1,x2)=-1; (x1,x2)=(pi,pi).
			
       {
               f=0;
               x1=x[0];
               x2=x[1];
               f=-cos(x1)*cos(x2)*exp((-pow((x1-pi),2)-pow((x2-pi),2)));
                break;
       }

 case 17:	////Dekkers and Aarts Problem (DA)  -20<=x(i)<=20, i=1:2. f(0,15),f(0,-15)=-24777; 
			
       {
               f=0;
               x1=x[0];
              x2=x[1];
               f=1.0e5*x1*x1+x2*x2-pow((x1*x1+x2*x2),2)+1.0e-5*pow((x1*x1+x2*x2),4);
                break;
       }
       /*case 43:	////Hosaki Problem (HSK)   0<=x1<=5, 0<=x2<=6,f(4,2)=-2.3458;
			
       {
               f=0;
               x1=x[0];
               x2=x[1];
               f=(1 - 8*x1 + 7*x1*x1  - (x1*x1*x1*7)/3.0 + (x1*x1*x1*x1)/4)*x2*x2*exp(-x2);
                break;
       }
        case 44:	 // 	McCormick Problem (MC)   -1.5<=x1<=4, -3<=x2<=3,f(4,2)=-2.3458; 
			
       {
               f=0;
               x1=x[0];
               x2=x[1];
               f=sin(x1 + x2) + pow((x1 - x2),2) - 1.5*x1 + 2.5*x2+1;
                break;
       }
      case 45:	 //Meyer and Roth Problem (MR) (Wolfe, 1978)  -10<=x<=10, f(\F03.13; 15.16; 0.78)=0.4*10^(-4); 
			
       {
               f=0;
               x1=x[0];
               x2=x[1];
               x3=x[2];
               double t_arr[]={1.0,2.0,1.0,2.0,0.1};
               double v_arr[]={1.0,1.0,2.0,2.0,0.0};
               double y_arr[]={0.126,0.219,0.076,0.126,0.186};
               for(int i=0;i<5;i++)
               {
                       f=f+pow((((x1*x3*t_arr[i])/(1+x1*t_arr[i]+x2*v_arr[i]))-y_arr[i]),2);
               }
                break;
       }*/
       case 18:	 //Shubert Problem (SBT)  -10<=x<=10, f((7.0835, 4.8580)=-186.7309;
			
       {
               f=1;
              
               
               for(i=0;i<2;i++)
               {
                double temp=0;
                for(j=1;j<=5;j++)
                {
                       temp=temp+j*cos((j+1)*x[i]+j);
                }
                f=f*temp;
               }
                break;
       }
      /* case 47:	 // Sinusoidal Problem (SIN) 0<=x<=180, A= 2.5; B =5; z = 30. f(90 + z; 90 + z; . . . ; 90 + z)= -(A+1)
			
       {
               f=0;
                double temp1=1;  
                 double temp2=1;   
                double aa=2.5,bb=5,zz=30;        
              for(i=0;i<D;i++)
                {
                       temp1=temp1*sin(x[i]-zz);
                       temp2=temp2*sin(bb*(x[i]-zz));
                       
                }
                temp1=temp1*aa;
                f=-(temp1+temp2);
                 break;
               
       }
       case 48: // Moved axis parallel hyper-ellipsoid function [-5.12, 5.12] f(x)=0; x(i)= 5*i, i=1:D.
       {
            f = 0;
            for ( d = 0; d < D; d++ )
            {
                 f = f + 5*(d+1)*x[d] * x[d];
             }
             break;
      }
       case 49: // Pressure vessel (confinement method) 	// (1.125, 0.625, 55.8592, 57.7315) => 7197.729
			// If no granularity => min = 6059.7143
			// We are using function without granularity 
       {
            f = 0;
           	x1=x[0]; // [1.1,12.5] granularity 0.0625 // not using granularity
		    x2=x[1];// [0.6,12.5] granularity 0.0625 // not using granularity
		    x3=x[2]; // [0,240]
		    x4=x[3];// [0,240]

		f=0.6224*x1*x3*x4 + 1.7781*x2*x3*x3 + 3.1611*x1*x1*x4 + 19.84*x1*x1*x3;
		
		double ff[3],c;
		
		ff[0]=0.0193*x[2]-x[0];
		ff[1]=0.00954*x[2]-x[1];
		ff[2]=750*1728-pi*x[2]*x[2]*(x[3]+(4.0/3)*x[2]); 
		
		if (ff[0]>0) {c= 1+pow(10,10)*ff[0];  f=f*c*c; }
		if (ff[1]>0) {c=1+ff[1]; f=f*c*c;  }
		if (ff[2]>0) {c=1+ff[2]; f=f*c*c;  }
	
             break;
      }
        case 50: // lennard_jones (-2,2) for 5 atoms => -9.103852
       {
           	int d;
	        int dim=3;
	        double dist;
	       int i,j;
	        int nPoints=15/dim;
	         double zz;
            double x1[3];
            double x2[3];
	        f=0;
	        for(i=0;i<nPoints-1;i++)
	        {
		                            for(d=0;d<dim;d++)  x1[d]=x[3*i+d];
		                            for(j=i+1;j<nPoints;j++)
		                            {
			                                                for(d=0;d<dim;d++)  x2[d]=x[3*j+d];

			                                                dist=distance(x1,x2,2);

			                                                zz=pow(dist,-6); 
			                                                f=f+zz*(zz-1);
                                     }
            }
	        f=4*f;
             break;
      }// end of case 50
      
       case 51: // Parameter Estimation for Frequency-Modulated (FM) Sound Waves f(x)=0;
       {
            f = 0;
            double q=2*pi/100;
            double a1=x[0];
            double w1=x[1];
            double a2=x[2];
            double w2=x[3];
            double a3=x[4];
            double w3=x[5];
            
            for ( i = 1; i <= 100; i++ )
            {
                 f = f +pow((a1*sin(w1*i*q+a2*sin(w2*i*q+a3*sin(w3*i*q)))- 1.0*sin(5.0*i*q-1.5*sin(4.8*i*q+2.0*sin(4.9*i*q)))),2);
             }
             break;
      }*/
      case 19:  //Schwefel prob 1.2 (in list 1)
      {
       double s1=0;                                    
       double s2=0;
       for(i=0;i < D;i++)
       {
          for(j=0;j < i;j++)
          {
                 s1+=x[j];
          }
          s2=s2+s1*s1;
       }
       f=s2;
       break;
       }
       /*case 53:  //Schwefel prob 2.21
      {
       double s1=0;                                    
       double s2=0;
       double maxdim=fabs(x[0]);
       for(i=1;i < D;i++)
       {
          if(maxdim<fabs(x[i]))
          {
                 maxdim=fabs(x[i]);
          }
          
       }
       f=maxdim;
       break;
       }*/
         case 20:  //Schwefel prob (in list 3)
     {
       double s1=0;                                    
         
       for(i=0;i < D;i++)
       {
          s1=s1+x[i]*sin(sqrt(fabs(x[i])));
          
       }
       f=-s1;
       break;
     }
     
	case 21:          //hartmann function 3
	{ 
		static double p[4][3] =
            { 		0.36890,0.11700,0.26730,
			0.46990,0.43870,0.74700,0.10910,0.87320,0.55470,
			0.03815,0.57430,0.88280
                   
            };	
		static double c[4]=
		{
			1.0,1.2,3.0,3.2

		};
		static double a[4][3] = {
				3.0,0.1,30.0,35.0
			
			};
		//static double x[] = {1.0,2.0,3.0};
			float x1 = x[0];
			float x2 = x[1];
			float x3 = x[2];
		int s= 0, i,j;
		
		for(i=0;i<4;i++){
			int sm = 0;
		for(j=0;j<3;j++){
			sm = sm + a[i][j]*(((x1*x2*x3)-p[i][j])*((x1*x2*x3)-p[i][j]));
		}
		s = s + c[i] * exp(-sm);
		}
			f = -s;	
		//printf("%g \n",f);	
	}

	case 22:          //hartmann function 6
	{ 
		static double p[4][7] =
            { 		0.36890,0.11700,0.26730,
			0.46990,0.43870,0.74700,0.10910,0.87320,0.55470,
			0.03815,0.57430,0.88280
                   
            };	
		static double c[4]=
		{
			1.0,1.2,3.0,3.2

		};
		static double b[4][3] = {
				3.0,0.1,30.0,35.0,25,33
			
			};
		//static double x[] = {1.0,2.0,3.0};
			float x1 = x[0];
			float x2 = x[1];
			float x3 = x[2];
			float x4 = x[3];
			float x5 = x[4];
			float x6 = x[5];
		int o= 0, i,j;
		
		for(i=0;i<4;i++){
			int sm1 = 0;
		for(j=0;j<3;j++){
			sm1 = sm1 + b[i][j]*(((x1*x2*x3*x4*x5*x6)-p[i][j])*((x1*x2*x3*x4*x5*x6)-p[i][j]));
		}
		o = o + c[i] * exp(-sm1);
		}
			f = -o;	
		//printf("%g \n",f);	
	}
	case 23:          //shekel function 5
	{ 
		static double c[5] =
            { 		
		0.1,0.2,0.2,0.4
//,0.3,0.7,0.5,0.5,0.4
                   
            };	
		
		static double a[5][4] = {
				3.0,7.0,2.0,9.0,5.0
//,3.0,8.0,1.0,6.0,2.0
			
			};
		//static double x[] = {1.0,2.0,3.0};
			float x1 = x[0];
			float x2 = x[1];
			float x3 = x[2];
			float x4 = x[3];
			//float x5 = x[4];
			//float x6 = x[5];
		int s= 0, i,j;
		
		for(j=0;j<5;j++){
			int p = 0;
		for(i=0;i<4;i++){
			p = p + ((x1*x2*x3*x4)-a[j][i])*((x1*x2*x3*x4)-a[j][i]);
		}
		s = s + 1/(p+c[j]);
		}
			f = -s;	
		//printf("%g \n",f);	
	}
	case 24:          //Shekel function 7
	{ 
		static double c[7] =
            { 		
		0.1,0.2,0.2,0.4,0.3,0.7,0.5
//,0.5,0.4
                
            };	
		
		static double a[7][4] = {
				3.0,7.0,2.0,9.0,5.0,3.0,8.0
//,1.0,6.0,2.0
			
			};
		//static double x[] = {1.0,2.0,3.0};
			float x1 = x[0];
			float x2 = x[1];
			float x3 = x[2];
			float x4 = x[3];
			//float x5 = x[4];
			//float x6 = x[5];
		int s= 0, i,j;
		
		for(j=0;j<7;j++){
			int p = 0;
		for(i=0;i<4;i++){
			p = p + ((x1*x2*x3*x4)-a[j][i])*((x1*x2*x3*x4)-a[j][i]);
		}
		s = s + 1/(p+c[j]);
		}
			f = -s;	
		//printf("%g \n",f);	
	}
	case 25:          //Shekel function 10
	{ 
		static double c[10] =
            {		
		0.1,0.2,0.2,0.4,0.3,0.7,0.5,0.5,0.4,0.8
                
            };	
		
		static double a[10][4] = {
				3.0,7.0,2.0,9.0,5.0,3.0,8.0,1.0,6.0,2.0
			
			};
		//static double x[] = {1.0,2.0,3.0};
			float x1 = x[0];
			float x2 = x[1];
			float x3 = x[2];
			float x4 = x[3];
			//float x5 = x[4];
			//float x6 = x[5];
		int s= 0, i,j;
		
		for(j=0;j<7;j++){
			int p = 0;
		for(i=0;i<4;i++){
			p = p + ((x1*x2*x3*x4)-a[j][i])*((x1*x2*x3*x4)-a[j][i]);
		}
		s = s + 1/(p+c[j]);
		}
			f = -s;	
		//printf("%g \n",f);	
	}



}
    return f;
       
}// end of the optimization fun



void initilize_params(int Pr)
{
 int d;
 if (Pr==0)// Parabola (Sphere)
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -5.12; ub[d] =5.12;
           }

        }
        /*if (Pr==1)// De Jong's f4
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
            lb[d] = -5.12; ub[d] =5.12;
  
            }
        }
        if (Pr==2)//// Griewank
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -600; ub[d] =600;
           }
        }
        if (Pr==3)//// Rosenbrock
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-2;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -100; ub[d] =100;
           }
        }*/
        if (Pr==1)//// Rastrigin. Minimum value 0. Solution (0,0 ...0)
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
                lb[d] = -5.12; ub[d] = 5.12;
           }
        }
       /* if (Pr==5)//// Ackley
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -30; ub[d] =30;
           }  
        }        
        if (Pr==6)////Alpine
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
                 lb[d] = -10; ub[d] =10;
           }
        }   
        if (Pr==7)//Michalewicz function
        {
           D=10;
           obj_val=-9.66015;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = 0; ub[d] =pi;
           } 
        }
        if (Pr==8)//Cosine Mixture [-1,1] f(0,0,...0)=-D*0.1
        {
           D=30;
           obj_val=-D*0.1;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -1; ub[d] =1;
           }
        } 
        if (Pr==9) //Exponential [-1,1] f(0,0,...0)= -1
        {
           D=30;
           obj_val=-1;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -1; ub[d] =1;
           }
        } 
        if (Pr==10)//Zakharov's [-5.12,5.12] f(0,0,0,...,0)=0
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-2;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -5.12; ub[d] =5.12;
           }
        } */  

        if (Pr==2)//Cigar [-10,10]  f(0,0,0,...,0)=0
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -10; ub[d] =10;
           }
        }
       /* if (Pr==12) //brown3 [-1,4]  f(0,0,0.....,0)=0
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -1; ub[d] =4;
           }
        }*/
        if (Pr==3)//Schewel prob 2.22
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -10; ub[d] =10;
           }
        } 
        /*if (Pr==14)//Salomon Problem (SAL) f(0,0,0,0,0,0)=0
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-1;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -100; ub[d]=100;
           }
        } */
        if (Pr==4)// Axis parallel hyperellipsoid
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -5.12; ub[d] =5.12;
           }
        } 
       /* if (Pr==16)// Pathological
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -100; ub[d] =100;
           }
        } 
        if (Pr==17)//Sum of different powers
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -1; ub[d] =1;
           }
        } */
        if (Pr==5)//step function [-100, 100] f(-0.5<=x<=0.5)=0
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -100; ub[d] =100;
           }
        } 
        /*if (Pr==19) //Quartic function, i.e., noise [-1.28, 1.28] f(0000..00)=0
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -1.28; ub[d] =1.28;
           }
        }
        if (Pr==20)//Inverted cosine wave function (Masters) [-5, 5] f(000..0)=-D+1
        {
           D=10;
           obj_val=-D+1;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -5; ub[d] =5;
           }
        }
        if (Pr==21)//Neumaier 3 Problem (NF3) (Neumaier, 2003b)
        {
           D=10;
           obj_val=-(D*(D+4)*(D-1))/6.0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -900; ub[d] =900;
           }
        }
        if (Pr==22)//Rotated hyper-ellipsoid function
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -65.536; ub[d] =65.536;
           }
        }*/
        if (Pr==6) //Levy montalvo 1 
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -10; ub[d] =10;
           }
        }
        if (Pr==7)//Levy montalvo 2 
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -5; ub[d] =5;
           }
        }
       /* if (Pr==25)//Ellipsoidal Ellipsoidal [-D,D] f(1,2,3,...,D)=0
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -D; ub[d] =D;
           }
        }*/
        if (Pr==8)//Beale function [-4.5,4.5] f(3, 0.5)=0
        {
           D=2;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -4.5; ub[d] =4.5;
           }
        }
       /* if (Pr==27)//Colville function [-10,10] f(1111)=0
        {
           D=4;
           obj_val=0;
           acc_err=1.0e-3;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = -10; ub[d] =10;
           }
        }*/
        if (Pr==9)//Branins\92s function [-5,10][0,15] f(-pi, 12.275)=0.3979
        {
               D=2;
               acc_err=1.0e-5;
               obj_val=0.3979;
               lb[0] = -5; ub[0] =10;
               lb[1] = 0; ub[1] =15;
           
        }
         /*if (Pr==29)//Kowalik function [-5,5] f(0.192833, 0.190836, 0.123117, 0.135766)=0.000307486
        {
               D=4;
               acc_err=1.0e-4;
               obj_val=0.000307486;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -5; ub[d] =5;
                }
           
        }
        // Minibench Mark Problems: 4-Problem Taken from clerc website
         if (Pr==30)//// 2D Tripod function [-100,100] f(0, -50)=0
        {
               D=2;
               acc_err=1.0e-4;
               obj_val=0;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -100; ub[d] =100;
                }
           
        }
        if (Pr==31)//Shifted Rosenbrock F6 CEC 2005 [-100, 100], solution point is O + (1; : : : ; 1) where f = 390.
        {
               D=10;
               acc_err=1.0e-1;
               obj_val=390;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -100; ub[d] =100;
                }
           
        }*/
           if (Pr==10)//// Shifted Parabola/Sphere (CEC 2005 benchmark)	x\81\B8[.100,100] , Global optimum: x* = o , 1( *) 1 F x = f_bias = - 450
        {
                D=10;
               acc_err=1.0e-5;
               obj_val=-450;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -100; ub[d] =100;
                }
           
        }
     /* if (Pr==33)//Shifted CEC 2005  Rastrigin  x\81[-5,5] , Global optimum x* = offset , f(x*) = f_bias = - 330
        {
                D=10;
               acc_err=1.0e-2;
               obj_val=-330;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -5; ub[d] =5;
                }
           
        }*/
         if (Pr==11)//Shifted CEC 2005 Schwefel [-100,100], Global optimum x* = offset , f(x*) = f_bias = - 450
        {
               D=10;
               acc_err=1.0e-5;
               obj_val=-450;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -100; ub[d] =100;
                }
           
        }
         if (Pr==12)//Shifted CEC 2005 Griewank. WARNING: in the CEC 2005 benchmark it is rotated
        {
                D=10;
               acc_err=1.0e-5;
               obj_val=-180;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -600; ub[d] =600;
                }
           
        }
         if (Pr==13)// Shifted Ackley (CEC 2005) [-32,32], Global optimum x* = offset , f(x*) = f_bias = - 140
        {
                D=10;
               acc_err=1.0e-5;
               obj_val=-140;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -32; ub[d] =32;
                }
           
        }
        
        /* if (Pr==37)// // Compression spring  [1,...,70],[0.6,3],[0.207,0.5] f(7; 1:386599591; 0:292) = 2.6254214578.
        {
                D=3;
               acc_err=1.0e-10;
               obj_val=2.6254214578;
               lb[0] = 1; ub[0] =70;
               lb[1] = 0.6; ub[1] =3;
               lb[2] = 0.207; ub[2] =0.5;
        }*/
        /// if (Pr=38) Gear Train Problem
         if (Pr==14)// 	// Goldstein-Price function  [-2,2] f(0, -1) = 3.
        {
                D=2;
               acc_err=1.0e-14;
               obj_val=3;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -2; ub[d] =2;
                }
           
        }
         if (Pr==15)// 	//Six-hump camel back function  [-5 < x < 5] f(-0.0898,0.7126) = -1.0316.
        {
                D=2;
               acc_err=1.0e-5;
               obj_val=-1.0316;
                   lb[0] = -5; ub[0] =5;
                    lb[1] = -5; ub[1] =5;
                
           
        }
          if (Pr==16)// 	//Easom's function  -10<=x(i)<=10, i=1:2. f(x1,x2)=-1; (x1,x2)=(pi,pi).
        {
               D=2;
               acc_err=1.0e-13;
               obj_val=-1;
               for ( d = 0; d < D; d++ )
                {
                    lb[d] = -100; ub[d] =100;
                }
                
           
        }
        if (Pr==17)// 	//Dekkers and Aarts Problem (DA)  -20<=x(i)<=20, i=1:2. f(0,15),f(0,-15)=-24777; 
        {
               D=2;
               acc_err=0.5;
               obj_val=-24777;
               for ( d = 0; d < D; d++ )
                {
                    lb[d] = -20; ub[d] =20;
                }
                
           
        }
         /* if (Pr==43)// 	//Hosaki Problem (HSK)   0<=x1<=5, 0<=x2<=6,f(4,2)=-2.3458; 
        {
               D=2;
               acc_err=1.0e-6;
               obj_val=-2.3458;
                lb[0] = 0; ub[0] =5;
                lb[1] = 0; ub[1] =6;
                
           
        }
         if (Pr==44)// 	McCormick Problem (MC)   -1.5<=x1<=4, -3<=x2<=3,f(-0.547, -1.547)==-1.9133; 
        {
               D=2;
               acc_err=1.0e-4;
               obj_val=-1.9133;
                lb[0] = -1.5; ub[0] =4;
                lb[1] = -3; ub[1] =3;
                
           
        }
        if (Pr==45)//Meyer and Roth Problem (MR) (Wolfe, 1978)  -10<=x<=10, f(\F03.13; 15.16; 0.78)=0.4*10^(-4); 
        {
               D=3;
               acc_err=1.95e-3;
               obj_val=-0.4e-4;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -10; ub[d] =10;
                }
                
           
        }*/
        if (Pr==18)//Shubert Problem (SBT)  -10<=x<=10, f((7.0835, 4.8580)=-186.7309; 
        {
               D=2;
               acc_err=1.0e-5;
               obj_val=-186.7309;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -10; ub[d] =10;
                }
                
           
        }
       /* if (Pr==47)// Sinusoidal Problem (SIN) 0<=x<=180, A= 2.5; B =5; z = 30. f(90 + z; 90 + z; . . . ; 90 + z)= -(A+1)
        {
               D=10;
               acc_err=1.0e-2;
               obj_val=-3.5;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = 0; ub[d] =180;
                }
        }
         if (Pr==48)// Moved axis parallel hyper-ellipsoid function [-5.12, 5.12] f(x)=0; x(i)= 5*i, i=1:D.
        {
               D=30;
               acc_err=1.0e-15;
               obj_val=0;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -5.12; ub[d] =5.12;
                }
        }
        if (Pr==49)// Pressure vessel (confinement method) 	// (1.125, 0.625, 55.8592, 57.7315) => 7197.729
			
        {
               D=4;
               acc_err=1.0e-5;
               obj_val=7197.729;
               lb[0] = 1.125; ub[0] =12.5;
               lb[1] = 0.625; ub[1] =12.5;
               lb[2] = 1.0e-8; ub[2] =240;
               lb[3] = 1.0e-8; ub[3] =240;                
        }
        if (Pr==50)// lennard_jones (-2,2) for 5 atoms => -9.103852
			
        {
               int nAtoms=5; // in {2, ..., 15}
               D=3*nAtoms; 
               acc_err=1.0e-3;
               obj_val=-9.103852;
  	            for ( d = 0; d < D; d++ )
                {
                    lb[d] = -2; ub[d] =2;
                }               
        }
         if (Pr==51)// Parameter Estimation for Frequency-Modulated (FM) Sound Waves f(x)=0;
			
        {
               D=6;
               acc_err=1.0e-2;
               obj_val=0;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -6.4; ub[d] =6.35;
                }
        }*/
          if (Pr==19)  //Schwefel prob 1.2
      {
               D=30;
               acc_err=1.0e-3;
               obj_val=0;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -100; ub[d] =100;
                }
      
       }
        /*if (Pr==20)  //Schwefel prob 2.21
      {
               D=30;
               acc_err=1.0e-3;
               obj_val=0;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -100; ub[d] =100;
                }
      
       }*/
          if (Pr==20)  //Schwefel 
      {
               D=30;
              acc_err=1.0e-3;
               obj_val=-418.9829*D;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -500; ub[d] =500;
                }
      
       }
       /*if (Pr==55)	 // Shekel Foxholes Functions; 
			
       {
               D=2;
               acc_err=1.0e-3;
               obj_val=0.998;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = -65.536; ub[d] =65.536;
                }
              
       }*/
        if (Pr==23)	 // Shekel Functions 5; 
			
       {
               D=4;
               acc_err=1.0e-3;
               obj_val=-10.1532;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = 0; ub[d] =10;
                }
                
       }
     if (Pr==24)	 // Shekel Functions 7; 
			
       {
               D=4;
               acc_err=1.0e-3;
               obj_val=-10.4029;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = 0; ub[d] =10;
                }
               
       }
       if (Pr==25)	 // Shekel Functions 10; 
			
       {
               D=4;
               acc_err=1.0e-3;
               obj_val=-10.5364;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = 0; ub[d] =10;
                }
               
       }
        if (Pr==21)	 // Hartmann Functions 3 dim; 
			
       {
                D=3;
              acc_err=1.0e-3;
               obj_val=-3.86278;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = 0; ub[d] =1;
                }
                
       }
       if (Pr==22)	 // Hartmann Functions 6 dim; 
			
       {
               D=6;
              acc_err=1.0e-3;
               obj_val=-3.32237;
                for ( d = 0; d < D; d++ )
                {
                    lb[d] = 0; ub[d] =1;
                }
                
       }
        
    
}
