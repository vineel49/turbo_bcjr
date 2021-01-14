/* Iterative decoding of turbo codes using the logBCJR algorithm */
/* G(D)=[1 (1+D^2)/(1+D+D^2)], QPSK modulation, overall rate is 1/2 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#define columns 512 // num_bit
/* global declaration */
void swap(int *p,int *q);
void logBCJR(double log_Gamma[][columns],int num_bits,double LLR[]);
void logBCJR_END(double log_Gamma[][columns],int num_bits,double LLR[]);
double maxstar (double a, double b); 
/*--------------------------*/
void main()
{
clock_t start, end;
double cpu_time_used;
start = clock();
system("cls"); // clears the screen    
/*--------------------------------------------------------------------------------------------------------------------*/
/*    SIMULATION PARAMETERS */    
int frame_size = 1024; /* frame size */
float SNR_dB = 3; /* SNR per bit (dB) */
long int sim_runs=pow(10,2); /* simulation runs */
int num_bit = 0.5*frame_size; /*number of data bits */
double noise_var_1D; /* 1D AWGN variance */ 
noise_var_1D = (double)2*2/(2*pow(10,0.1*SNR_dB)); /* overall rate is 1/2, Average power of QPSK symbols is 2 */
/*--------------------------------------------------------------------------------------------------------------------*/
/* Interleaving & Deinterleaving mapping (Fisher-Yates shuffle algorithm) */
int Iv_Map[num_bit],Div_Map[num_bit];
for(int i=0;i<num_bit;i++)
Iv_Map[i]=i;

for (int i = num_bit-1; i > 0; i--) 
{ 
int j = rand() % (i+1); 
swap(&Iv_Map[i],&Iv_Map[j]); 
}

for (int i=0;i<num_bit;i++)
Div_Map[Iv_Map[i]]=i;
/*------------------------------------------------------*/
/*       QPSK Constellation                             */
int constell_size = 4; // QPSK constellation size
double re_qpsk_sym[4]={1,1,-1,-1},im_qpsk_sym[4]={1,-1,1,-1};
/*-----------------------------------------------------------*/
long int fr_cnt;
long long int C_BER=0;
for(fr_cnt=0;fr_cnt<sim_runs;fr_cnt++)
{
srand(time(0));
/*----------------------------------------------------------*/
/* source  (0s and 1s generation according to uniform distribution)*/
int a[num_bit];
for(int i=0;i<num_bit;i++)
a[i] = rand()%2;
/*---------------------------------------------------------*/
/* Turbo encoder (outputs QPSK symbols) */
int re_sym1[num_bit],im_sym1[num_bit],re_sym2[num_bit],im_sym2[num_bit];
int D1_1=0,D2_1=0,D1_2=0,D2_2=0,temp;
for(int i=0;i<num_bit;i++)
{
// component encoder 1    
re_sym1[i]=1-2*a[i]; // Re{QPSK symbol} from component encoder 1
temp = a[i]^D1_1; // parity bit from component encoder 1
im_sym1[i]=1-2*temp; // Im{QPSK symbol} from component encoder 1
temp = temp^D2_1;
D2_1=D1_1;
D1_1=temp; // right shift operation

// component encoder 2
re_sym2[i]=1-2*a[Iv_Map[i]]; // Re{QPSK symbol} from component encoder 2
temp = a[Iv_Map[i]]^D1_2; // parity bit from component encoder 2
im_sym2[i]=1-2*temp; // Im{QPSK symbol} from component encoder 2
temp = temp^D2_2;
D2_2=D1_2;
D1_2=temp; // right shift operation

} // for(int i=0;i<num_bit;i++)
/*------------------------------------------------------------------------------*/
/* AWGN (Marsaglia algorithm)*/
/* https://rosettacode.org/wiki/Statistics/Normal_distribution */
double re_chan_op1[num_bit],im_chan_op1[num_bit];
double x1,y1,rsq1,f1;
double re_chan_op2[num_bit],im_chan_op2[num_bit];
double x2,y2,rsq2,f2;
for(int i=0;i<num_bit;i++)
{
do {
x1 = 2.0 * rand() / (double)RAND_MAX - 1.0;
y1 = 2.0 * rand() / (double)RAND_MAX - 1.0;
rsq1 = x1 * x1 + y1 * y1;
}while( rsq1 >= 1. || rsq1 == 0. );
f1 = sqrt( -2.0 * log(rsq1) / rsq1 );
re_chan_op1[i] = re_sym1[i]+sqrt(noise_var_1D)*x1*f1;
im_chan_op1[i] = im_sym1[i]+sqrt(noise_var_1D)*y1*f1;

do {
x2 = 2.0 * rand() / (double)RAND_MAX - 1.0;
y2 = 2.0 * rand() / (double)RAND_MAX - 1.0;
rsq2 = x2 * x2 + y2 * y2;
}while( rsq2 >= 1. || rsq2 == 0. );
f2 = sqrt( -2.0 * log(rsq2) / rsq2 );
re_chan_op2[i] = re_sym2[i]+sqrt(noise_var_1D)*x2*f2;
im_chan_op2[i] = im_sym2[i]+sqrt(noise_var_1D)*y2*f2;
}
/*-------------------------------------------------------------------------*/
/*                       RECEIVER                                          */
/* Gammas for BCJR algorithm */
double log_Gamma1[4][num_bit],log_Gamma2[4][num_bit],temp1,temp2;

for(int u=0;u<constell_size;u++)
{
for(int v=0;v<num_bit;v++)
{
temp1=(re_chan_op1[v]-re_qpsk_sym[u])*(re_chan_op1[v]-re_qpsk_sym[u]);
temp2=(im_chan_op1[v]-im_qpsk_sym[u])*(im_chan_op1[v]-im_qpsk_sym[u]);
log_Gamma1[u][v]=-1*(temp1+temp2)/(2*noise_var_1D);

temp1=(re_chan_op2[v]-re_qpsk_sym[u])*(re_chan_op2[v]-re_qpsk_sym[u]);
temp2=(im_chan_op2[v]-im_qpsk_sym[u])*(im_chan_op2[v]-im_qpsk_sym[u]);
log_Gamma2[u][v]=-1*(temp1+temp2)/(2*noise_var_1D);
}    
}
/*---------------------------------------------------------------*/
// apriori LLR initialization
double LLR[num_bit];
for(int i=0;i<num_bit;i++)
LLR[i]=0;
/*--------------------------------------------------------------*/
void interleave(double LLR[],int num_bit,int Iv_Map[]);/*local declaration of the interleave function */
void deinterleave(double LLR[],int num_bit,int Div_Map[]);/*local declaration of the interleave function */
/*  Iterative MAP decoding using the logBCJR algorithm          */ 


logBCJR(log_Gamma1,num_bit,LLR);
interleave(LLR,num_bit,Iv_Map);
logBCJR(log_Gamma2,num_bit,LLR); // iteration #1

deinterleave(LLR,num_bit,Div_Map);
logBCJR(log_Gamma1,num_bit,LLR);
interleave(LLR,num_bit,Iv_Map);
logBCJR(log_Gamma2,num_bit,LLR); // iteration #2

deinterleave(LLR,num_bit,Div_Map);
logBCJR(log_Gamma1,num_bit,LLR);
interleave(LLR,num_bit,Iv_Map);
logBCJR(log_Gamma2,num_bit,LLR); // iteration #3

deinterleave(LLR,num_bit,Div_Map);
logBCJR(log_Gamma1,num_bit,LLR);
interleave(LLR,num_bit,Iv_Map);
logBCJR(log_Gamma2,num_bit,LLR); // iteration #4

deinterleave(LLR,num_bit,Div_Map);
logBCJR(log_Gamma1,num_bit,LLR);
interleave(LLR,num_bit,Iv_Map);
logBCJR(log_Gamma2,num_bit,LLR); // iteration #5

deinterleave(LLR,num_bit,Div_Map);
logBCJR(log_Gamma1,num_bit,LLR);
interleave(LLR,num_bit,Iv_Map);
logBCJR(log_Gamma2,num_bit,LLR); // iteration #6

deinterleave(LLR,num_bit,Div_Map);
logBCJR(log_Gamma1,num_bit,LLR);
interleave(LLR,num_bit,Iv_Map);
logBCJR(log_Gamma2,num_bit,LLR); // iteration #7

deinterleave(LLR,num_bit,Div_Map);
logBCJR(log_Gamma1,num_bit,LLR);
interleave(LLR,num_bit,Iv_Map);
logBCJR_END(log_Gamma2,num_bit,LLR); // iteration #8 (OUTPUTS APOSTERIORI INFORMATION)
/*-------------------------------------------------------------------------------------*/
/*            HARD DECISION            */
deinterleave(LLR,num_bit,Div_Map);
int a_hat[num_bit];
int errors=0,err;
// calculating bit error rate
for (int u=0;u<num_bit;u++)
{
    a_hat[u]=LLR[u]>0?0:1;
    err=a[u]==a_hat[u]?0:1;
    errors = errors+err;
}

C_BER = C_BER+errors;
} // for fr_cnt

end = clock();
cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
double BER;
BER = (double)C_BER/(num_bit*sim_runs);

printf("\n Simulation runs: %ld",sim_runs);
printf("\n Bit error rate: %lf",BER);
printf("\n Time taken: %lf seconds",cpu_time_used);
} // main()

//        FUNCTIONS
void swap (int *a, int *b) 
{ 
    int temp = *a; 
    *a = *b; 
    *b = temp; 
}
/*------------------------------------------------------------------------------------------*/
void interleave(double LLR[],int num_bit,int Iv_Map[])
{
double LLR_d[num_bit];

for(int i=0;i<num_bit;i++)
LLR_d[i]=LLR[Iv_Map[i]];

for(int i=0;i<num_bit;i++)
LLR[i]=LLR_d[i];
}
/*------------------------------------------------------------------------------*/
void deinterleave(double LLR[],int num_bit,int Div_Map[])
{
double LLR_d[num_bit];

for(int i=0;i<num_bit;i++)
LLR_d[i]=LLR[Div_Map[i]];

for(int i=0;i<num_bit;i++)
LLR[i]=LLR_d[i];
}
/*------------------------------------------------------------------------------*/
double maxstar (double x, double y) 
{ 
    double res;
    if (x>=y)
    {
    res = x + log(1+exp(-x+y)); 
    }
    else
    {
    res = y + log(1+exp(-y+x));
    }
    return res;

}
/*-----------------------------------------------------------------------------------------*/
void logBCJR(double log_Gamma[][columns],int num_bits,double LLR[]) // outputs extrinsic information
{
/* Trellis for G(D)=[1 (1+D^2)/(1+D+D^2)]        */
int Ga_Inx[4][2] = {{0,3},{1,2},{0,3},{1,2}};  // Gamma indices for alpha recursion
int P_State[4][2] = {{0,1},{3,2},{1,0},{2,3}}; // previous state
int P_Ip[4][2] =  {{0,1},{0,1},{0,1},{0,1}} ; // Previous input.
int N_State[4][2] = {{0,2},{2,0},{3,1},{1,3}}; // Next state
int Gb_Inx[4][2] = {{0,3},{0,3},{1,2},{1,2}}; // Gamma indices for beta recursion
int N_Ip[4][2] = {{0,1},{0,1},{0,1},{0,1}}; // Next input
int num_states=4; // number of states in the trellis

double C[num_bits];
for(int i=0;i<num_bits;i++)
C[i] =  exp(-1*LLR[i]/2)/(1+exp(-1*LLR[i])); // Ck

double log_alpha[num_states][num_bits+1],log_beta[num_states][num_bits+1],temp3,temp4;
for(int i=0;i<num_states;i++)
{
log_alpha[i][0]=0;
log_beta[i][num_bits]=0; // assuming that the receiver does not know the starting and ending states    
}
for(int u=0;u<num_bits;u++)
{
for(int v=0;v<num_states;v++)
{
temp3=log_alpha[P_State[v][0]][u] + C[u] +(1-2*P_Ip[v][0])*LLR[u]/2 + log_Gamma[Ga_Inx[v][0]][u];
temp4=log_alpha[P_State[v][1]][u] + C[u] +(1-2*P_Ip[v][1])*LLR[u]/2 + log_Gamma[Ga_Inx[v][1]][u];
log_alpha[v][u+1] = maxstar(temp3,temp4); // Jacobian logarithm

temp3=log_beta[N_State[v][0]][num_bits-u] + C[num_bits-1-u] +(1-2*N_Ip[v][0])*LLR[num_bits-1-u]/2 + log_Gamma[Gb_Inx[v][0]][num_bits-1-u];
temp4=log_beta[N_State[v][1]][num_bits-u] + C[num_bits-1-u] +(1-2*N_Ip[v][1])*LLR[num_bits-1-u]/2 + log_Gamma[Gb_Inx[v][1]][num_bits-1-u];
log_beta[v][num_bits-1-u] = maxstar(temp3,temp4); // Jacobian logarithm
}  //for(int v=0;v<num_states;v++)  
} //for(int u=0;u<num_bit;u++)
/* extrinsic information calculation */
double LLR1_temp[num_states][num_bits],LLR2_temp[num_states][num_bits];
for (int u=0;u<num_bits;u++)
{
for (int v=0;v<num_states;v++)
{
LLR1_temp[v][u]=log_alpha[v][u]+log_Gamma[Gb_Inx[v][0]][u]+log_beta[N_State[v][0]][u+1];
LLR2_temp[v][u]=log_alpha[v][u]+log_Gamma[Gb_Inx[v][1]][u]+log_beta[N_State[v][1]][u+1];
}
temp3 = maxstar(maxstar(LLR1_temp[0][u],LLR1_temp[1][u]),maxstar(LLR1_temp[2][u],LLR1_temp[3][u]));
temp4 = maxstar(maxstar(LLR2_temp[0][u],LLR2_temp[1][u]),maxstar(LLR2_temp[2][u],LLR2_temp[3][u]));
LLR[u]=temp3-temp4;    
}    
}
/*----------------------------------------------------------------------------------*/
void logBCJR_END(double log_Gamma[][columns],int num_bits,double LLR[]) // outputs aposteriori information
{
/* Trellis for G(D)=[1 (1+D^2)/(1+D+D^2)]        */
int Ga_Inx[4][2] = {{0,3},{1,2},{0,3},{1,2}};  // Gamma indices for alpha recursion
int P_State[4][2] = {{0,1},{3,2},{1,0},{2,3}}; // previous state
int P_Ip[4][2] =  {{0,1},{0,1},{0,1},{0,1}} ; // Previous input.
int N_State[4][2] = {{0,2},{2,0},{3,1},{1,3}}; // Next state
int Gb_Inx[4][2] = {{0,3},{0,3},{1,2},{1,2}}; // Gamma indices for beta recursion
int N_Ip[4][2] = {{0,1},{0,1},{0,1},{0,1}}; // Next input
int num_states=4; // number of states in the trellis

double C[num_bits];
for(int i=0;i<num_bits;i++)
C[i] =  exp(-1*LLR[i]/2)/(1+exp(-1*LLR[i])); // Ck

double log_alpha[num_states][num_bits+1],log_beta[num_states][num_bits+1],temp3,temp4;
for(int i=0;i<num_states;i++)
{
log_alpha[i][0]=0;
log_beta[i][num_bits]=0; // assuming that the receiver does not know the starting and ending states    
}
for(int u=0;u<num_bits;u++)
{
for(int v=0;v<num_states;v++)
{
temp3=log_alpha[P_State[v][0]][u] + C[u] +(1-2*P_Ip[v][0])*LLR[u]/2 + log_Gamma[Ga_Inx[v][0]][u];
temp4=log_alpha[P_State[v][1]][u] + C[u] +(1-2*P_Ip[v][1])*LLR[u]/2 + log_Gamma[Ga_Inx[v][1]][u];
log_alpha[v][u+1] = maxstar(temp3,temp4); // Jacobian logarithm

temp3=log_beta[N_State[v][0]][num_bits-u] + C[num_bits-1-u] +(1-2*N_Ip[v][0])*LLR[num_bits-1-u]/2 + log_Gamma[Gb_Inx[v][0]][num_bits-1-u];
temp4=log_beta[N_State[v][1]][num_bits-u] + C[num_bits-1-u] +(1-2*N_Ip[v][1])*LLR[num_bits-1-u]/2 + log_Gamma[Gb_Inx[v][1]][num_bits-1-u];
log_beta[v][num_bits-1-u] = maxstar(temp3,temp4); // Jacobian logarithm
}  //for(int v=0;v<num_states;v++)  
} //for(int u=0;u<num_bit;u++)
/* extrinsic information calculation */
double LLR1_temp[num_states][num_bits],LLR2_temp[num_states][num_bits];
for (int u=0;u<num_bits;u++)
{
for (int v=0;v<num_states;v++)
{
LLR1_temp[v][u]=log_alpha[v][u]+log_Gamma[Gb_Inx[v][0]][u]+LLR[u]/2+log_beta[N_State[v][0]][u+1];
LLR2_temp[v][u]=log_alpha[v][u]+log_Gamma[Gb_Inx[v][1]][u]-LLR[u]/2+log_beta[N_State[v][1]][u+1];
}
temp3 = maxstar(maxstar(LLR1_temp[0][u],LLR1_temp[1][u]),maxstar(LLR1_temp[2][u],LLR1_temp[3][u]));
temp4 = maxstar(maxstar(LLR2_temp[0][u],LLR2_temp[1][u]),maxstar(LLR2_temp[2][u],LLR2_temp[3][u]));
LLR[u]=temp3-temp4;    
}    
}