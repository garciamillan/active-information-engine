//
//  rnt_control.c
//  
//
//  Created by Rosalba Garcia Millan on 02/06/2023.
//  Optimal protocol of run-and-tumble particle in harmonic potential
//  COMPILE: cc -O3 -Wall -o eng rnt_engine.c -lgsl
//  EXECUTE: GSL_RNG_SEED=$RANDOM ./eng
//  code rnt_control.c extended on 2 May 2024 to run the cyclic information engine.

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <string.h>

#define RTP //AOUP //RTP //PASSIVE

//particle
#define U_PART (0.2) //self-propulsion
#define D_PART (0.) //(0.0001) //diffusion constant (thermal diffusion)
//#define ALPHA_PART (1.) //tumbling rate
#define ALPHA_PART (0.382416450134328) //tumbling rate
#define X0_PART (0.) //initial position
#define SIGN0_PART (1.) //initial direction: +-1 for RTP, +1 for AOUP, 0 for passive

//potential
#define K_POTENTIAL (2.)

//protocol
#define LAMBDA0 (0.)
#define LAMBDAF (2.)
#define LAMBDA0 (0.)
#define FINAL_T (10.)

#define DT (0.003)
#define TIME_BEFORE_PROTOCOL (-10.)
#define TOL_ZERO (0.1)//(1.e-5)
#define TOL_STD (0.05)//(4.e-3)
#define RELAX_MAX_STEPS (100)//(4.e-3)
#define CHUNK (1000)
#define NUM_REALISATIONS (1)//(100*CHUNK)
#define NUM_CYCLES (15)

#define SIG (0.001*DT)

#ifndef PARTICLE
#define PARTICLE
struct particle {
    //random variables (updated in time)
    double x;   //position
    double v;   //self-propulsion (either RTP or AOUP)
    double w;   //velocity (self-propulsion + external force)

    //parameters
    double x0;  //initial position
    double v0;  //initial self-propulsion
    double D;   //diffusion constant (thermal diffusion)
    double u;   //RTP self-propulsion
    double alpha;   //tumbling rate
    double tau; //persistence time
    double sigma;   //standard deviation thermal diffusion
    double Dy;  //AOUP diffusion constant in self-propulsion
    double sigma_y; //AOUP self-propulsion standard deviation
};
#endif


#define MALLOC(a,n) if ((a=malloc(sizeof(*a)*(n)))==NULL) { fprintf(stderr, "Not enough memory for %s, requested %i bytes, %i items of size %i. %i::%s\n", #a, (int)(sizeof(*a)*n), n, (int)sizeof(*a), errno, strerror(errno)); exit(EXIT_FAILURE); } else { printf("# Info malloc(3)ed %i bytes (%i items of %i bytes) for %s.\n", (int)(sizeof(*a)*(n)), n, (int)sizeof(*a), #a); }

#define DD fprintf(stderr, "%s::%i\n", __FILE__, __LINE__)

void initialise_particle(struct particle *part, gsl_rng *r, double initial_x0, double initial_v0, double u_part, double D_part, double a_part);
void update_velocity(struct particle *part, gsl_rng *r, double t, double time_f, double lambda_f);
void update_position(struct particle *part, gsl_rng *r);
void print_position(int step, struct particle *part, double t, double time_f, double lambda_f, double shift_lambda, double shift_time);


void print_header(FILE *file, char*const*argv, int argc, gsl_rng *r);
void print_parameters(FILE *file, double u_part, double alpha_part, double D_part, double k, double x0_part, double v0_part, double time_f, double dt, double lambda_f);

int flag_print_POS=1;
int flag_print_WORK=1;
//int flag_distr_WORK=1;
int flag_stationarity=1;
int flag_measure_velocity=1;
int flag_print_ave_TRAJ=0;
//int flag_observables=1;

double k_potential;
double final_t;
double Dt;
//double lambdaf;
double u_part;
double D_part;
double a_part;

double *m0_pos;
double *m1_pos;
double *m0_dwork;
double *m1_dwork;


double protocol(double t, double x0, double v0, double time_f, double lambdaf){ //protocol_jumps
    //protocol with initial and final jumps
    double lambda;
    double tau = 1./a_part; //persistence time is the inverse of the tumbling rate
    double d = (lambdaf - x0 + tau*v0*(exp(-time_f/tau) - 1.)/2.);
    //fprintf(stderr,"  d = %.16G    %.16G  %.16G  %.16G  %.16G  %.16G  %.16G\n",d,lambdaf,x0,tau,v0,time_f,tau);
    if(t<=0) lambda = 0.;
    else if(t<time_f) lambda = x0
                                + (t + 1./k_potential)/(time_f + 2./k_potential)*d
                                + tau*v0*(1. - exp(-t/tau)*(1. + 1./(tau * k_potential)))/2.;
    else lambda = lambdaf;
    return lambda;
}

double dprotocol(double t, double x0, double v0, double time_f, double lambdaf){
    //discretised symmetric time derivative
    return (protocol(t+Dt/2.,x0,v0,time_f,lambdaf)-protocol(t-Dt/2.,x0,v0,time_f,lambdaf));
}

double force_potential(double x, double t, double x0, double v0, double time_f, double lambdaf){
    return k_potential * (x - protocol(t,x0,v0,time_f,lambdaf));
}

double dWork(double x, double t, double x0, double v0, double time_f, double lambdaf){
    //printf("  %.16G  %.16G\n",dprotocol(t),force_potential(x,t));
    return -dprotocol(t,x0,v0,time_f,lambdaf)*force_potential(x,t,x0,v0,time_f,lambdaf);
}

double protocol_smooth(double t, double x0, double v0, double time_f, double lambdaf){
    //smooth protocol (no jumps)
    double lambda, Dlambda0, Dlambdaf, m;
    double tau = 1./a_part; //persistence time is the inverse of the tumbling rate
    double d = (lambdaf - x0 + tau*v0*(exp(-time_f/tau) - 1.)/2.);

    Dlambda0 = d/(k_potential*time_f+2.) + x0 - v0/(2.*k_potential);
    Dlambdaf = d/(k_potential*time_f+2.) + v0/(2.*k_potential)*exp(-time_f/tau);
    m = d/(time_f+2./k_potential);
    
    //steps
    double H;
    if(t<=0) H = 0.;
    else if(t<time_f) H = m*t + v0*(tau + 1./k_potential)*(1. - exp(-t/tau))/2.;
    else H = m*time_f + v0*(tau + 1./k_potential)*(1. - exp(-time_f/tau))/2.;
    
    //smoothening
    lambda = H + Dlambda0*(tanh(t/SIG) + 1.)/2. + Dlambdaf*(tanh((t-time_f)/SIG) + 1.)/2.;

    return lambda;
}

double dWork_smooth(double x, double t, double x0, double v0, double time_f, double lambdaf){
    return -(protocol_smooth(t+Dt/2.,x0,v0,time_f,lambdaf)-protocol_smooth(t-Dt/2.,x0,v0,time_f,lambdaf))*k_potential * (x - protocol_smooth(t,x0,v0,time_f,lambdaf));
}



double theory_work_m1(double x0, double v0, double time_f, double lambdaf){
    double tau = 1./a_part; //persistence time is the inverse of the tumbling rate
    double d = (lambdaf - x0 + tau*v0*(exp(-time_f/tau) - 1.)/2.);
    return d*d/(time_f + 2./k_potential) - k_potential*x0*x0/2. + tau*v0*v0*(exp(-2.*time_f/tau) - 1.)/8.;
}

double theory_work_m2(double x0, double time_f, double lambdaf){
    //double tau = 1./a_part; //persistence time is the inverse of the tumbling rate
    //double d = (lambdaf - initial_x0 + tau*v0*(exp(-time_f/tau) - 1.)/2.);
    double tmp0 = (lambdaf - x0);
    double tmp1 = k_potential * time_f;
    double tmp2 = 1. - exp(-tmp1);
    double tmp3 = (time_f + 2./k_potential);
    
    return (k_potential*lambdaf*lambdaf/2.)*(k_potential*lambdaf*lambdaf/2.) - k_potential*k_potential*lambdaf*lambdaf * tmp0 / tmp3 * (x0*tmp2/k_potential + time_f*time_f*(lambdaf *(tmp1-3.)+x0*(2.*tmp1+9.) )/(6.*tmp3)) + tmp0*tmp0/(tmp3*tmp3) * ( x0*tmp2*x0*tmp2 + D_part * tmp2 * tmp2 * tmp2/k_potential + ((2.*tmp0 - x0*(2.+tmp1))*tmp2 + tmp1 *(lambdaf*(-4.+tmp1) + x0*(8.+tmp1))/2.)*((2.*tmp0 - x0*(2.+tmp1))*tmp2 + tmp1 *(lambdaf*(-4.+tmp1) + x0*(8.+tmp1))/2.)/(tmp3*tmp3));
}


double theory_m1_x(double t, double x0, double v0, double time_f, double lambdaf){
    double tau = 1./a_part; //persistence time is the inverse of the tumbling rate
    double d = (lambdaf - x0 + tau*v0*(exp(-time_f/tau) - 1.)/2.);
    return x0 + d/(time_f + 2./k_potential)*t - tau*v0*(exp(-t/tau) - 1.)/2.;
}

int main(int argc, char *argv[], char *envp[]){
    setlinebuf(stdout);
    
    u_part = U_PART;
    D_part = D_PART;
    a_part = ALPHA_PART;
    k_potential = K_POTENTIAL;
    final_t = FINAL_T;
    double lambdaf = LAMBDAF;

    //double shift_position_particle = 0.;
    double shift_position_trap = 0.;
    double shift_time = 0.;

    int step, n_steps;
    double time;
    
    double initial_x0 = X0_PART;
    double initial_v0;
#ifdef RTP
    initial_v0 = SIGN0_PART*U_PART;
#elifdef AOUP
    initial_v0 = U_PART;
#else
    initial_v0 = 0.;
#endif

    struct particle part;
    

    double obs_work,obs_work_CHUNK;
    double obs_work_smooth,obs_work_CHUNK_smooth,obs_work2_smooth,obs_work_cycle,theor_work_engine;
    double obs_work_engine;

    int ch;
    while ((ch = getopt(argc, argv, "u:d:a:k:l:T:h")) != -1)
            switch (ch) {
                case 'u':
                    u_part=atof(optarg);
                    break;
                case 'd':
                    if ((D_part=atof(optarg))<0) {
                        fprintf(stderr, "diffusion constant D_part must be non-negative but D_part = %lf.\n",D_part);
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'a':
                    if ((a_part=atof(optarg))<0) {
                        fprintf(stderr, "tumbling rate a_part must be non-negative but a_part = %lf.\n",a_part);
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'k':
                    if ((k_potential=atof(optarg))<0) {
                        fprintf(stderr, "spring constant k must be non-negative but k = %lf.\n",k_potential);
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'l':
                    lambdaf=atof(optarg);
                    break;
                case 'T':
                    if ((final_t=atof(optarg))<0) {
                        fprintf(stderr, "final time must be non-negative but final_t = %lf.\n",final_t);
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'h':
                    fprintf(stderr, "Usage: time GSL_RNG_SEED=$RANDOM ./protocol > test.dat &\n");
                    fprintf(stderr, "Input arguments: -u u_part -d D_part -h help\n");
                    return(0);
                    break;
                default:
                    fprintf(stderr, "Flag %c unknown.\n", ch);
                    break;
            }
        
    /* create a generator chosen by the environment variable GSL_RNG_TYPE */
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    n_steps=(int) (final_t/DT);
    Dt = final_t/((double) n_steps);

    print_header(stdout, argv, argc, r);
    print_parameters(stdout,u_part,a_part,D_part,k_potential,initial_x0,initial_v0,final_t,Dt,lambdaf);
    
    int i,N,cycle,Ncycles;
    N = NUM_REALISATIONS;
    Ncycles = NUM_CYCLES;
    
    double m0_work,m1_work,m2_work,std_work;
    double m0_work_CHUNK,m1_work_CHUNK;
    double m0_work_smooth,m1_work_smooth,m2_work_smooth,std_work_smooth;
    double m0_work2_smooth,m1_work2_smooth,m2_work2_smooth,std_work2_smooth;
    double m0_work_CHUNK_smooth,m1_work_CHUNK_smooth,work2_CHUNK_smooth;
    m0_work = m1_work = m2_work = 0.;
    //m0_work_CHUNK = m1_work_CHUNK = 0.;
    m0_work_smooth = m1_work_smooth = m2_work_smooth = 0.;
    m0_work2_smooth = m1_work2_smooth = m2_work2_smooth = 0.;
    //m0_work_CHUNK_smooth = m1_work_CHUNK_smooth = work2_CHUNK_smooth = 0.;
    MALLOC(m0_pos,n_steps);
    MALLOC(m1_pos,n_steps);
    MALLOC(m0_dwork,n_steps);
    MALLOC(m1_dwork,n_steps);
    for(step=0;step<=n_steps;step++) m0_pos[step] = m1_pos[step] = 0.;
    for(step=0;step<=n_steps;step++) m0_dwork[step] = m1_dwork[step] = 0.;

#ifdef RTP
    printf("#Info: RTP (Run-and-Tumble Particle) \n");
    printf("#Info: self-propulsion v = %.16G\n", u_part);
    printf("#Info: tumbling rate alpha = %.16G\n",a_part);
#elifdef AOUP
    printf("#Info: AOUP (Active Ornstein-Uhlenbeck Particle)\n");
    printf("#Info: diffusivity self-propulsion Dy = %.16G\n", u_part * u_part / a_part);
    printf("#Info: persistence time tau = %.16G\n",1./a_part);
#else
    printf("#Info: Passive Particle\n");
#endif
    
    if(flag_stationarity) printf("#Info: initial x0 and v0 measured from stationary state\n");
    else if(flag_measure_velocity) printf("#Info: initial x0 = %.16G, v0 measured from stationary distribution\n",initial_x0);
    else printf("#Info: initial x0 = %.16G , v0 =  %.16G\n",initial_x0,initial_v0);

    {

        for(i=0;i<N;i++){
            
            //default initial state: x_0 and u_0 are given by the parameters X0_PART and SIGN0_PART
            //flag_stationarity: particle is initially at the stationary state, when both position and self-propulsion are measured
            //flag_measure_velocity: initial position x_0 given by the parameters X0_PART, and initial u_0 from stationary velocity distribution
            
            if(flag_stationarity){
                //initialisise particle
                initialise_particle(&part,r,initial_x0,initial_v0,u_part,D_part,a_part);
                
                //to start the protocol at stationarity I run the process until the mean position is zero within error, and its std is small enough
                double mom0_x,mom1_x,mom2_x,std_x;
                mom0_x = mom1_x = mom2_x = 0.;
                double obs, mom0_x_CHUNK, mom1_x_CHUNK;
                int timestep = 0;
                
                double tol_std = TOL_STD;
                double tol_zero = TOL_ZERO;
                time = TIME_BEFORE_PROTOCOL;
                                
                do {
                    mom0_x_CHUNK = mom1_x_CHUNK = 0.;
                    for(step=0;step<CHUNK;step++){
                        update_velocity(&part,r,time,final_t,0.);
                        update_position(&part,r);
                        mom0_x_CHUNK ++;
                        mom1_x_CHUNK += part.x;
                    }
                    obs = mom1_x_CHUNK/mom0_x_CHUNK;
                    mom0_x ++;
                    mom1_x += obs;
                    mom2_x += obs*obs;
                    std_x = sqrt((mom2_x/mom0_x-(mom1_x/mom0_x)*(mom1_x/mom0_x))/(mom0_x-1.));
                    
                    timestep++;
                    
                    //printf("RELAX  %.16G  %.16G  %.16G  %.16G\n",mom0_x,mom0_x/mom0_x,mom1_x/mom0_x,std_x);
                    //printf("RELAX  %.16G  %.16G  %.16G  %.16G\n", fabs(mom1_x/mom0_x),tol_zero,std_x,tol_std);
                } while((fabs(mom1_x/mom0_x) > tol_zero || std_x > tol_std) && timestep < RELAX_MAX_STEPS);
                
                //printf("#Info: system in stationary state.\n");
                printf("#Info: Relaxation timesteps  %d\n",timestep);
                initial_x0 = part.x;
                initial_v0 = part.v;
            } else if (flag_measure_velocity){ //self-propulsion is measured at stationarity
#ifdef RTP
                initial_v0 = (gsl_ran_bernoulli(r,0.5)*2. - 1.)*u_part;
#elifdef AOUP
                initial_v0 = gsl_ran_gaussian(r,sqrt((part.Dy)/(part.tau)));
#else
                initial_v0 = 0.;
#endif
            }
            
            //initialisise process with chosen initial conditions before first protocol
            initialise_particle(&part,r,initial_x0,initial_v0,u_part,D_part,a_part);
            //printf("#Info: initial measurement: x0 =  %.5G, v0 =  %.5G, u_part =  %.5G\n",initial_x0,initial_v0,u_part);
            //printf("#Info: initial measurement: x =  %.5G, v =  %.5G\n",part.x,part.v);
            
            //initially, the engine is centered at the origin, 0
            lambdaf = 0.;
            theor_work_engine = 0.;
            obs_work_engine = 0.;

            //engine starts here
            printf("#Info: switch engine on.\n");
            for(cycle=0;cycle<Ncycles;cycle++){
                obs_work_cycle = 0.;

                //in each cycle, the self-propulsion is measured and the target position of the trap is decided
                //the initial position is taken from the *conditional* average: initial_x0=v0/(k+alpha)
                part.v0 = part.v;
                part.x0 = part.v/(k_potential+a_part);
                part.x -= lambdaf; //position origin is centeret at trap position
                lambdaf = part.x0 + part.v*part.tau*(1.-exp(-final_t*a_part))/2.;
                printf("#Info: cycle %d\n"
                       "#Info: self-propulsion measurement %.16G\n"
                       "#Info: target trap position %.16G\n",cycle,part.v,lambdaf);

                m0_work_CHUNK = m1_work_CHUNK = 0.;
                m0_work_CHUNK_smooth = m1_work_CHUNK_smooth = work2_CHUNK_smooth = 0.;
                
                obs_work_CHUNK = 0.;
                obs_work_CHUNK_smooth = 0.;
                /*part.x = initial_x0;
                 part.sign = initial_sign;
                 part.u = initial_sign * u_part;*/
                
                //simulate trajectories
                for(step=0;step<=n_steps;step++){
                    time = ((double) step)*Dt;
                    
                    
                    obs_work_cycle += dWork_smooth(part.x,time,part.x0,part.v0,final_t,lambdaf);
                    theor_work_engine += dWork_smooth(theory_m1_x(time,part.x0,part.v0,final_t,lambdaf),time,part.x0,part.v0,final_t,lambdaf);
                    
                    //print trajectories
                    {
                        if(flag_print_POS) print_position(step,&part,time,final_t,lambdaf,shift_position_trap,shift_time);
                        
                        if(flag_print_WORK && flag_print_POS) printf("WORK_C  %d  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G\n", step,shift_time+time,dWork_smooth(part.x,time,part.x0,part.v0,final_t,lambdaf), obs_work_cycle,obs_work_engine,theor_work_engine, ((double)(cycle+1))*(-part.Dy/part.tau)*( k_potential/(2.*(k_potential+1./part.tau)*(k_potential+1./part.tau)) + part.tau*(1. - exp(-2.*final_t/(part.tau)))/8. ));
                    }

                    //calculate observables
                    /*obs_work_CHUNK += dWork(part.x,time,initial_x0,initial_v0,final_t);
                     obs_work_CHUNK_smooth += dWork_smooth(part.x,time,initial_x0,initial_v0,final_t);
                     m0_pos[step] ++;
                     m1_pos[step] += part.x;
                     m0_dwork[step] ++;
                     m1_dwork[step] += dWork(part.x,time,initial_x0,initial_v0,final_t);*/
                    
                    //update particle velocities
                    update_velocity(&part,r,time,final_t,lambdaf);
                    
                    //update particle positions
                    update_position(&part,r);

                }
                obs_work_cycle += dWork_smooth(part.x,time,part.x0,part.v0,final_t,lambdaf);
                theor_work_engine += dWork_smooth(theory_m1_x(time,part.x0,part.v0,final_t,lambdaf),time,part.x0,part.v0,final_t,lambdaf);

                m0_work_CHUNK ++;
                m1_work_CHUNK += obs_work_CHUNK;
                m0_work_CHUNK_smooth ++;
                m1_work_CHUNK_smooth += obs_work_CHUNK_smooth;
                work2_CHUNK_smooth += obs_work_CHUNK_smooth*obs_work_CHUNK_smooth;
                /*if(!((i+1)%CHUNK)){
                 obs_work = m1_work_CHUNK/m0_work_CHUNK;
                 m0_work ++;
                 m1_work += obs_work;
                 m2_work += obs_work*obs_work;
                 
                 obs_work_smooth = m1_work_CHUNK_smooth/m0_work_CHUNK_smooth;
                 m0_work_smooth ++;
                 m1_work_smooth += obs_work_smooth;
                 m2_work_smooth += obs_work_smooth*obs_work_smooth;
                 
                 obs_work2_smooth = work2_CHUNK_smooth/m0_work_CHUNK_smooth;
                 m0_work2_smooth ++;
                 m1_work2_smooth += obs_work2_smooth;
                 m2_work2_smooth += obs_work2_smooth*obs_work2_smooth;
                 
                 printf("WORK_C  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G\n", Dt, SIG, lambdaf, k_potential, final_t, u_part, a_part, m0_work_CHUNK, m0_work_CHUNK/m0_work_CHUNK, m1_work_CHUNK/m0_work_CHUNK, m0_work_CHUNK_smooth, m0_work_CHUNK_smooth/m0_work_CHUNK_smooth, m1_work_CHUNK_smooth/m0_work_CHUNK_smooth, m0_work_CHUNK_smooth, m0_work_CHUNK_smooth/m0_work_CHUNK_smooth, work2_CHUNK_smooth/m0_work_CHUNK_smooth, theory_work_m1(initial_x0,initial_v0,final_t),theory_work_m2(initial_x0,final_t));
                 
                 m0_work_CHUNK = m1_work_CHUNK = 0.;
                 m0_work_CHUNK_smooth = m1_work_CHUNK_smooth = work2_CHUNK_smooth = 0.;
                 }*/
                //shift_position_particle += part.x;
                shift_position_trap += lambdaf;
                shift_time += final_t;
                obs_work_engine += obs_work_cycle;
                if(flag_print_WORK) printf("WORK_E  %d  %.16G  %.16G  %.16G  %.16G\n", cycle+1, obs_work_cycle,obs_work_engine,theor_work_engine, ((double)(cycle+1))*(-part.Dy/part.tau)*( k_potential/(2.*(k_potential+1./part.tau)*(k_potential+1./part.tau)) + part.tau*(1. - exp(-2.*final_t/(part.tau)))/8. ));

            }
        }
        /*std_work = sqrt((m2_work/m0_work-(m1_work/m0_work)*(m1_work/m0_work))/(m0_work-1.));
        std_work_smooth = sqrt((m2_work_smooth/m0_work_smooth-(m1_work_smooth/m0_work_smooth)*(m1_work_smooth/m0_work_smooth))/(m0_work_smooth-1.));
        std_work2_smooth = sqrt((m2_work2_smooth/m0_work2_smooth-(m1_work2_smooth/m0_work2_smooth)*(m1_work2_smooth/m0_work2_smooth))/(m0_work2_smooth-1.));

        printf("WORK_F  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G\n", Dt, SIG, lambdaf, k_potential, final_t, u_part, a_part, D_part, m0_work, m0_work/m0_work, m1_work/m0_work, std_work, m0_work_smooth, m0_work_smooth/m0_work_smooth,
               m1_work_smooth/m0_work_smooth,
               std_work_smooth, m0_work2_smooth, m0_work2_smooth/m0_work2_smooth,
               m1_work2_smooth/m0_work2_smooth,
               std_work2_smooth, theory_work_m1(initial_x0,initial_v0,final_t),theory_work_m2(initial_x0,final_t));

        if(flag_print_ave_TRAJ) for(step=0;step<=n_steps;step++)
            printf("TRAJ %d  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G\n", step, step*Dt, m0_pos[step], m0_pos[step]/m0_pos[step], m1_pos[step]/m0_pos[step], theory_m1_x(step*Dt,initial_x0,initial_v0,final_t), dWork(theory_m1_x(step*Dt,initial_x0,initial_v0,final_t),step*Dt,initial_x0,initial_v0,final_t), protocol(step*Dt,initial_x0,initial_v0,final_t), protocol_smooth(step*Dt,initial_x0,initial_v0,final_t), m0_dwork[step]/m0_dwork[step], m1_dwork[step]/m0_dwork[step], dWork(m1_pos[step]/m0_pos[step],step*Dt,initial_x0,initial_v0,final_t));
    */
        printf("#Info: done.\n\n");
    }
    
    return 0;
}

void print_header (FILE *file, char*const*argv, int argc, gsl_rng *r) {
    int i;
    printf("# $Header$\n");
    printf("# Command:");
    for (i=0; i<argc; i++) printf(" %s", argv[i]);
    printf("\n");
    { time_t tm;
      tm=time(NULL);
      printf("# Info: %s", ctime(&tm));
    }
    
    { char hostname[128];
        gethostname(hostname, sizeof(hostname)-1);
        hostname[sizeof(hostname)-1]=(char)0;
        printf("# Info: Hostname: %s\n", hostname);
    }
    
    { char cwd[1024];
        cwd[0]=0;
        if(getcwd(cwd, sizeof(cwd)-1)!=NULL){
        cwd[sizeof(cwd)-1]=(char)0;
        printf("# Info: Directory: %s\n", cwd);}
    }
    
    printf("# Info: PID: %i\n", (int)getpid());

    fprintf(file,"#\n# GSL_RNG_SEED = %lu\n# GSL generator type: %s\n"
            "# GSL first value = %lu\n#\n", gsl_rng_default_seed, gsl_rng_name (r), gsl_rng_get (r));

}


void print_parameters(FILE *file, double u_part, double alpha_part, double D_part, double k, double x0_part, double v0_part, double time_f, double dt, double lambda_f){
#define print_param(p) fprintf(file,"#Info_param: %s = %g\n",#p,(double)p)
    print_param(D_part);
    print_param(u_part);
    print_param(alpha_part);
    print_param(k);
    if(!flag_stationarity){
        print_param(x0_part);
        print_param(v0_part);
    }
    
    print_param(dt);
    print_param(SIG);
    print_param(time_f);
    print_param(lambda_f);
    
    print_param(CHUNK);
    print_param(NUM_REALISATIONS);
    print_param(NUM_CYCLES);
}


void initialise_particle(struct particle *part, gsl_rng *r, double initial_x0, double initial_v0, double u_part, double D_part, double a_part){
    //initial measurement
    part->x = part->x0 = initial_x0;
    part->v = part->v0 = initial_v0;

    //parameters
    part->D = D_part;
    part->u = u_part;
    part->alpha = a_part; //tumbling rate
    part->tau = 1./a_part; //persistence time
    part->sigma = sqrt(2.*D_part*Dt); //thermal noise
    part->Dy = part->tau * u_part * u_part; //matches self-propulsion second moment
    part->sigma_y = sqrt(2.* (part->Dy) * Dt / (part->tau * part->tau) ); //AOUP noise
}

void update_velocity(struct particle *part, gsl_rng *r, double t, double time_f, double lambda_f){
    //update particle velocity given external potential and self-propulsion
    double x = part->x;
    
    //self-propulsion
#ifdef RTP
    if(gsl_ran_bernoulli(r,part->alpha*Dt/2.)) part->v *= -1.; //transmutate with prob alpha*Dt/2.
    part->w = part->v;
#elifdef AOUP
    double noise = gsl_ran_gaussian(r,part->sigma_y);
    part->v += (-(part->v)/(part->tau)*Dt + noise);
    part->w = part->v;
#else
    part->w = 0.;
#endif
    
    //external potential
    part->w -= force_potential(x,t,part->x0,part->v0,time_f,lambda_f);
}

void update_position(struct particle *part, gsl_rng *r){
    double noise = gsl_ran_gaussian(r,part->sigma); //displacement due to diffusion
    part->x += ((part->w)*Dt + noise);
}

void print_position(int step, struct particle *part, double t, double time_f, double lambda_f, double shift_lambda, double shift_time){
    if(step==0) printf("POS\n");
    printf("POS  %d  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G  %.16G\n", step,shift_time+t,shift_lambda,shift_lambda+protocol(t,part->x0,part->v0,time_f,lambda_f),shift_lambda+part->x,shift_lambda+theory_m1_x(t,part->x0,part->v0,time_f,lambda_f),part->u,part->v);
}


