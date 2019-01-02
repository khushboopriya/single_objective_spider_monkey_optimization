
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "combined.h"
#include "main.h"
#define Pop_size 50
#define Max_iterations 1000 
#define Total_Run 100
#define max_part 5  /*this represents the maximum group size*/
int try;


double CalculateFitness(double fun)
 {
	 double result=0;
	 if(fun>=0)
	 {result=1/(fun+1);
	 }
	 else
	 {
		 result=1+fabs(fun);
	 }
	 return result;
 }

void create_g()
{
   
int g=0;
  for(lo=0;lo<Pop_size;)
    {

          group_size=(lo+(int)(Pop_size/part)-1);
           gpoint[g][0]=lo;
           gpoint[g][1]=group_size;
           if((Pop_size-group_size)<(Pop_size/part))
                {
   	         gpoint[g][1]=(Pop_size-1);
		break;}
           
    // printf("\n [%d] [%d] \n",lo,group_size);
                g = g+1;
		lo = group_size+1;
    }  
   //printf("\n");
	try=g;
     // printf("this is trial%d \n",g);
	//try = g;
	//printf("thisoooooo is trial%d \n",try);
	//printf("\n %s ","okslll"); 
	//printf("\n %s "," okspp"); 
 }
void initialize()
{
   int i,j,k;
   	for(i=0;i<Pop_size;i++)
	{
        // printf("\nindividual is %d \n",i);
       for (j=0;j<D;j++)
		{
        Population[i][j]=(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*(ub[j]-lb[j])+lb[j];
		new_position[j]=Population[i][j];
                //printf("poplulation is %lf",Population[i][j]);
		}// for j
        
	fun_val[i]=fun(new_position);
	fitness[i]=CalculateFitness(fun_val[i]);
    }// for i
    // Initialize Global Leader Learning
    GlobalMin=fun_val[0];
    for(i=0;i<D;i++)
        GlobalLeaderPosition[i]=Population[0][i];
    GlobalLimitCount=0;
    // Initialize Local Leader Learning
    for(k=0;k<try;k++)
	{
        LocalMin[k]=fun_val[gpoint[k][0]]; 
        LocalLimitCount[k]=0;   
         for(i=0;i<D;i++)
           LocalLeaderPosition[k][i]=Population[gpoint[k][0]][i];            
    }
   //printf("%s","ok"); 
}

void GlobalLearning()
{
   int i,j;
    double G_trial=GlobalMin;
    for(i=0;i<Pop_size;i++)
	{
	if (fun_val[i]<GlobalMin)
		{
        GlobalMin=fun_val[i];
        for(j=0;j<D;j++)
           GlobalLeaderPosition[j]=Population[i][j];
        }//if fun_val[i]
        //printf("\n %d ",i); 
	}// for i
	if(fabs(G_trial-GlobalMin)<acc_err)
	       GlobalLimitCount=GlobalLimitCount+1;
    else
          GlobalLimitCount=0;
	
 }

void LocalLearning()
{
   int i,j,k;
   
    double OldMin[Pop_size/2];
    for(k=0;k<try;k++)
	{
        OldMin[k]=LocalMin[k];
	   //printf("\n %d", k);          
    }// end for k
    
	for(k=0;k<try;k++)
	{
	for(i=gpoint[k][0];i<=gpoint[k][1];i++)
	{
	if (fun_val[i]<LocalMin[k])
		{
        LocalMin[k]=fun_val[i];
        for(j=0;j<D;j++)
           LocalLeaderPosition[k][j]=Population[i][j];
        }//if
   
	}// for i
   }// for k
   
    for(k=0;k<try;k++)
	{
        if (fabs(OldMin[k]-LocalMin[k])<acc_err)
             LocalLimitCount[k]=LocalLimitCount[k]+1;    
        else
              LocalLimitCount[k]=0;        
    }// for k
	//printf("\n %s "," oksss"); 	
 }// end of fun


void LocalLeaderPhase(int k)
{
  int i,j;
  lo=gpoint[k][0];
  group_size=gpoint[k][1];
  for (i=lo;i<=group_size;i++)
  {
        int PopRand;
        do{
           PopRand=(int)((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*(group_size-lo)+lo);
        }while(PopRand==i);
        
        for(j=0;j<D;j++){
        if((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )>=cr)
		{
            new_position[j]=Population[i][j]+(LocalLeaderPosition[k][j]-Population[i][j])*((double)rand() / ((double)(RAND_MAX)+(double)(1)))+
                            (Population[PopRand][j]-Population[i][j])*((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )-0.5)*2;
		     
	}
	  
        else
        {
            new_position[j]=Population[i][j];    
        }
        if (new_position[j]<lb[j])
           new_position[j]=lb[j];
        if (new_position[j]>ub[j])
           new_position[j]=ub[j];
	//printf("\n %d ",i); 
        }// end of for loop
        ObjValSol=fun(new_position);
        FitnessSol=CalculateFitness(ObjValSol);
        if (FitnessSol>fitness[i])
        {
         for(j=0;j<D;j++)
               Population[i][j]=new_position[j];
         fun_val[i]=ObjValSol;
         fitness[i]=FitnessSol;
        }// if fitness
  }// for i
	 
}// end of fun
void GlobalLeaderPhase(int k)
{

   int i,j,l;
  
  lo=gpoint[k][0];
  group_size=gpoint[k][1];
  i=lo;
  l=lo; /* l is the count*/
 while(l<group_size)
  {
         
     
       if((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )<prob[i]) 
        {  
              l++;      
         int PopRand;
        do{
       
        PopRand=(int)((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*(group_size-lo)+lo);
        }while(PopRand==i);
        param2change=(int)((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*D);
        for(j=0;j<D;j++)
            new_position[j]=Population[i][j];
        new_position[param2change]=Population[i][param2change]+(GlobalLeaderPosition[param2change]-Population[i][param2change])*
(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )+(Population[PopRand][param2change]-Population[i][param2change])*((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )-0.5)*2;
        if (new_position[param2change]<lb[param2change])
           new_position[param2change]=lb[param2change];
        if (new_position[param2change]>ub[param2change])
           new_position[param2change]=ub[param2change];
        ObjValSol=fun(new_position);
        FitnessSol=CalculateFitness(ObjValSol);
        if (FitnessSol>fitness[i])
        {
              for(j=0;j<D;j++)
                   Population[i][j]=new_position[j];
              fun_val[i]=ObjValSol;
              fitness[i]=FitnessSol;
        }
    }// if prob[i]
     i++;
     if (i==(group_size))
        i=lo;
 }//while l
}

void GlobalLeaderDecision()
{
if(GlobalLimitCount>GlobalLimit)
{
        
             
             GlobalLimitCount=0;
            
         if(part<max_part) 
        {     
              //printf("Global ... Partition\n");
             part=part+1;
              create_g(); 
              LocalLearning();
        }
        else
           {     //printf("Global ... Partition\n");
             part=1;
              create_g(); 
              LocalLearning();
        }
                        
                              
}
}
void LocalLeaderDecision()
{
     int i,j,k;
for(k=0;k<try;k++)
   {
          if(LocalLimitCount[k]>LocalLimit)
          {
            // printf("Local Limit Cross\n");
             for(i=gpoint[k][0];i<=gpoint[k][1];i++)
              {
                   for(j=0;j<D;j++){
                        if((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )>=cr){
                        Population[i][j]=Population[i][j]+(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*(ub[j]-lb[j]);
                    }// end of if
                    else
                    {
                        Population[i][j]=Population[i][j]+(GlobalLeaderPosition[j]-Population[i][j])*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )+(Population[i][j]-LocalLeaderPosition[k][j])*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );    
                    }// end of else
                     if (Population[i][j]<lb[j])
                       Population[i][j]=lb[j];
                    if (Population[i][j]>ub[j])
                       Population[i][j]=ub[j];
                    }// end of for j
                    
	                fun_val[i]=fun( Population[i]);
	                fitness[i]=CalculateFitness(fun_val[i]);
	           
	           
             }//For i
             
              LocalLimitCount[k]=0;
              //printf("\n %s", "hello");
         }//if
         
          
   }//k
  

}// End of Function


void CalculateProbabilities()
{
     int i;
     double maxfit;
     maxfit=fitness[0];
  for (i=1;i<Pop_size;i++)
        {
           if (fitness[i]>maxfit)
           maxfit=fitness[i];
        }
    for (i=0;i<Pop_size;i++)
        {
         prob[i]=(0.9*(fitness[i]/maxfit))+0.1;
        }

}

int main()
{
int run,j;
double mean;
mean=0;
srand(time(NULL));
 E = exp( 1 );
        pi = acos( -1 );
f_run = fopen( "SpiderMonkey_General1_SBCRor.ods", "w" );
f_run1 = fopen( "SpiderMonkey_General1_SBCR_100or.ods", "w" );    
fprintf(f_run," Prob \t Dim \t   Mean Fun Val   \t SD \t   Error \t Succ Mean Fun Eval \t  RUN \t Succ Rate \t Acceptable Error\n  ");
fprintf(f_run1," Prob  \t No. of RUN  \t Fun Val \t   Error \n  ");
int d;


                  // for combine.h
for ( Pr =0;Pr<27;Pr++)
{
    //if(Pr==37)
      /// Pr=39;

       initilize_params(Pr);
double mean_error=0,error=0,total_feval=0;
double mean=0.0,var=0.0,sd=0.0;
mean_feval=0;
int succ_rate=0;
LocalLimit=D*Pop_size;
GlobalLimit=Pop_size;
for(run=0;run<Total_Run;run++)
{
initialize();
GlobalLearning();
LocalLearning();
feval=0;
part=1;
create_g();
//printf("%s","hello");
//printf("\n%d main loop try \n",try);
iter=0;
error=0;
cr=0.1;
//double cr_gmin[3][2];
//cr=alea(0.1,0.9);
for (iter=0;iter<Max_iterations;iter++)
 {                
    for(int k=0;k<try;k++)
    {
          LocalLeaderPhase(k);
		//printf("%d k value \n ",k);
    }// for k
    CalculateProbabilities();
     for(int k=0;k<try;k++)
    {
          GlobalLeaderPhase(k);
    }// for k
            
    GlobalLearning();
    LocalLearning();
    LocalLeaderDecision();
    GlobalLeaderDecision ();     
                 
   
   
    if(fabs(GlobalMin-obj_val)<=acc_err)
    {
        succ_rate+=1;
        mean_feval=mean_feval+feval;
	//return succ_rate;  
        break;
        
    }//if GlobalMin*/
    //if(iter%((int)(Max_iterations/5.0))==0)
    //{
                  // cr=cr+0.1;
          
    //}
    cr=cr+(0.4/Max_iterations);
   
    //cr=cr+(0.4/1000);
    //printf("cr=%f",cr);
   // cr=alea(0.1,0.9);
}// iter
error=fabs(GlobalMin-obj_val);
mean_error=mean_error+error;
fprintf(f_run1,"%d \t %d  \t %e \t %e \n",Pr,run+1,GlobalMin,error);
printf("Pr %d Dim %d  Run %d  F_val %e Error %e \n",Pr, D, run+1,GlobalMin,error);
GlobalMins[run]=GlobalMin;
mean=mean+GlobalMin;
total_feval=total_feval+feval;
}// end of run
mean=mean/Total_Run;
mean_error=mean_error/Total_Run;
if(succ_rate>0)
 mean_feval=mean_feval/(double)(succ_rate);
total_feval=total_feval/Total_Run;
 for(int k=0;k<Total_Run;k++)
    var=var+pow(GlobalMins[k]-mean,2);
    var=var/Total_Run;
    sd=sqrt(var);
//for(int k=0;k<D;k++)
//fprintf(f_run," %f \t",GlobalLeaderPosition[k]);
printf("Means of %d runs: %e \n",Total_Run,mean);
fprintf(f_run," %d \t %d \t  %e \t %e \t %e \t  %f \t %d \t %d  \t %e \n",Pr,D, mean,sd, mean_error,total_feval,Total_Run, succ_rate,acc_err );
//printf("Means of %d runs: %e \n",Total_Run,mean);
}// end of problems
fclose(f_run);
fclose(f_run1);
//getch();
}

FILE * f_synth;
FILE * f_expl;
FILE * fLandscape;
FILE * f_init;
FILE * f_init_save;

FILE * f_coeff;

FILE *f_swarm; // For information about the variable swarm size



