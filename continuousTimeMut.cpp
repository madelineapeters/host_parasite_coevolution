//Directories
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <new>
#include <memory>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <complex>
#include <vector>
#include <sys/stat.h>
#include <random>
#include <limits>

using namespace std; //standard set of notations
//Global variables
char filename[200];


//output files
ofstream out_Pars; //global so can print to the same file in multiple functions
ofstream out_Test;
ofstream out_Data;
ofstream out_DataN;
ofstream out_erand;

//Classes (difference between C and C++)
class population; //Declare population
class metapopulation
{
    public:
            //Variables
            int kappaH,kappaP; //Maximum population size
            double phi,omega; //Migration rates for host and parasite
            double mu,nu; //Mutation rates for host and parasite
            double **betaMtx; //Vector of vectors for probability infection
            double alpha; //Mortality given infection
            double delta,gamma; //Extrinsic mortality
            population *demes;
            double *ratesVec; //event rates
            double *dVec; //vector for d terms
            int nE; //number of events
            double **eventsMtx, **eventsMtxN; //events matrix
            int nV; //number of host and parasites types in each population
            double rateTotal;
            double tmax;

            //Functions
            metapopulation();//Constructor function
            ~metapopulation();//Destructor function
            void calculateRates();
            void initializeEventsMtx();
            void initializeMeta();//Need to give output type; not passing pointer, so can't values in class
            void initializeRep();
            void simulation(char* argv1, int rep);
            void neutralEvent(int event);
        private:
};
class population
{
    public:
        //Variables
        int *state,*stateN; //state of deme (four element vector of host type 1 and 2 and parasite 1 and 2)
        metapopulation *metaPopPtr;

        //Functions
        population();//Constructor function
        ~population();//Destructor function
        void initializePop(metapopulation metaPop);//Need to give output type; not passing pointer, so can't values in class
        void resetPop(metapopulation metaPop);
        //void calcAF();

    private:
};

//(Global) Function declarations (say what functions are going to exist)
void input(metapopulation *metaPop1Ptr, char* argv1);

int chooseOutcome(vector<double> prob){
    
    // Setup the random bits
    std::random_device rd;
    std::mt19937 gen(rd());

   // Create the distribution with those weights
    std::discrete_distribution<> dist(prob.begin(), prob.end());

    // use the distribution and print the results.
    int out=dist(gen);

    return out;
}

//Main function: C++ looks for main function and runs main function
int main(int argc, char**argv) //return integer (int), first argument is integer called argc, second is array of characters called argv
{

    //Random seed
    //unsigned int seed=(unsigned int)time(NULL);
    srand(5682);//Initializing random seed
    //cout<<seed<<endl;

    //Variables
    metapopulation metaPop;

    //Setting up output files
    sprintf(filename,"./Data_%s",argv[1]);
    mkdir(filename,S_IRWXU);
    sprintf(filename,"./Data_%s/parameters_%s.csv",argv[1],argv[1]); //Makes the input file name
    out_Pars.open(filename);

    sprintf(filename,"./Data_%s/e_rand.csv",argv[1]); //Makes the input file name
    out_erand.open(filename);

    //Running program
    input(&metaPop,argv[1]); //& gives arrow, * before gives thing pointed to by arrow; (*x) is opposite of &x
    metaPop.initializeMeta(); //double check
    for(int rep=1;rep<2001;rep++){

        //Random seed
        //unsigned int seed=(unsigned int)time(NULL);
        //srand(seed);//Initializing random seed

        //metaPop.initializeRep();
        metaPop.simulation(argv[1], rep);
    }
    
    //End program
    out_Pars.close();
    out_Test.close();
    out_erand.close();

    return 0;
}

//Function definitions
void input(metapopulation *metaPopPtr,char* argv1)
{
    //Input parameter values from a file (ones that are pretty constant)
    ifstream inFile;
    sprintf(filename,"./inputMut.txt"); //Makes the input file name
    inFile.open(filename);
    string myString;
    myString="Null";

    // Look for tmax in file
    while(myString != ("(tmax):") && inFile.good()) {inFile>>myString;}
    //Once found tmax, read in next word
    inFile>>(*metaPopPtr).tmax; //Go to thing pointed to by *metaPopPtr ((*metaPopPtr) is THING POINTED TO by pointer) and take element tmax of thing pointed to
    //Read in tmax
    //cout<<"tmax: "<<(*metaPopPtr).tmax<<endl;
    out_Pars<<"tmax:, "<<(*metaPopPtr).tmax<<endl;
    
    // Look for kappaH in file
    while(myString != ("(kappaH):") && inFile.good()) {inFile>>myString;}
    //Once found kappaH, read in next word
    inFile>>(*metaPopPtr).kappaH; //Go to thing pointed to by *metaPopPtr ((*metaPopPtr) is THING POINTED TO by pointer) and take element kappaH of thing pointed to
    //Read in kappaH
    //cout<<"kappaH: "<<(*metaPopPtr).kappaH<<endl;
    out_Pars<<"kappaH:, "<<(*metaPopPtr).kappaH<<endl;
    
    // Look for kappaP in file
    while(myString != ("(kappaP):") && inFile.good()) {inFile>>myString;}
    //Once found kappaP, read in next word
    inFile>>(*metaPopPtr).kappaP; //Go to thing pointed to by *pop1Ptr and take element kappaH of thing pointed to
    //Read in kappaP
    //cout<<"kappaP: "<<(*metaPopPtr).kappaP<<endl;
    out_Pars<<"kappaP:, "<<(*metaPopPtr).kappaP<<endl;

    // Look for phi in file
    while(myString != ("(phi):") && inFile.good()) {inFile>>myString;}
    //Once found phi, read in next word
    inFile>>(*metaPopPtr).phi; //Go to thing pointed to by *pop1Ptr and take element kappaH of thing pointed to
    //Read in phi
    //cout<<"phi: "<<(*metaPopPtr).phi<<endl;
    out_Pars<<"phi:, "<<(*metaPopPtr).phi<<endl;

    // Look for omega in file
    while(myString != ("(omega):") && inFile.good()) {inFile>>myString;}
    //Once found omega, read in next word
    inFile>>(*metaPopPtr).omega; //Go to thing pointed to by *pop1Ptr and take element kappaH of thing pointed to
    //Read in omega
    //cout<<"omega: "<<(*metaPopPtr).omega<<endl;
    out_Pars<<"omega:, "<<(*metaPopPtr).omega<<endl;

    (*metaPopPtr).betaMtx=new double*[2];
    for(int r=0;r<2;r++) (*metaPopPtr).betaMtx[r]=new double[2];

    // Look for beta in file
    while(myString != ("(beta):") && inFile.good()) {inFile>>myString;}
    //Once found beta, read in next word
    //cout<<"betaMtx: "<<endl;
    out_Pars<<"betaMtx: "<<endl;
    for(int r=0;r<2;r++){
        for(int c=0;c<2;c++){
            inFile>>(*metaPopPtr).betaMtx[r][c];
            //cout<<(*metaPopPtr).betaMtx[r][c]<<",";
            out_Pars<<(*metaPopPtr).betaMtx[r][c]<<",";
        }
        //cout<<endl;
        out_Pars<<endl;
    }

    // Look for alpha in file
    while(myString != ("(alpha):") && inFile.good()) {inFile>>myString;}
    //Once found alpha, read in next word
    inFile>>(*metaPopPtr).alpha; //Go to thing pointed to by *pop1Ptr and take element kappaH of thing pointed to
    //Read in alpha
    //cout<<"alpha: "<<(*metaPopPtr).alpha<<endl;
    out_Pars<<"alpha:, "<<(*metaPopPtr).alpha<<endl;

    // Look for delta in file
    while(myString != ("(delta):") && inFile.good()) {inFile>>myString;}
    //Once found delta, read in next word
    inFile>>(*metaPopPtr).delta; //Go to thing pointed to by *pop1Ptr and take element kappaH of thing pointed to
    //Read in delta
    //cout<<"delta: "<<(*metaPopPtr).delta<<endl;
    out_Pars<<"delta:, "<<(*metaPopPtr).delta<<endl;

    // Look for gamma in file
    while(myString != ("(gamma):") && inFile.good()) {inFile>>myString;}
    //Once found gamma, read in next word
    inFile>>(*metaPopPtr).gamma; //Go to thing pointed to by *pop1Ptr and take element kappaH of thing pointed to
    //Read in gamma
    //cout<<"gamma: "<<(*metaPopPtr).gamma<<endl;
    out_Pars<<"gamma:, "<<(*metaPopPtr).gamma<<endl;

    // Look for mu in file
    while(myString != ("(mu):") && inFile.good()) {inFile>>myString;}
    //Once found mu, read in next word
    inFile>>(*metaPopPtr).mu; //Go to thing pointed to by *pop1Ptr and take element kappaH of thing pointed to
    //Read in mu
    //cout<<"mu: "<<(*metaPopPtr).mu<<endl;
    out_Pars<<"mu:, "<<(*metaPopPtr).mu<<endl;

    // Look for nu in file
    while(myString != ("(nu):") && inFile.good()) {inFile>>myString;}
    //Once found mu, read in next word
    inFile>>(*metaPopPtr).nu; //Go to thing pointed to by *pop1Ptr and take element kappaH of thing pointed to
    //Read in nu
    //cout<<"nu: "<<(*metaPopPtr).nu<<endl;
    out_Pars<<"nu:, "<<(*metaPopPtr).nu<<endl;


}
void randExp(double lambda, double* randOut)
{
    (*randOut)=-1.0/lambda*log(rand()/(double)RAND_MAX);
    while((1+(*randOut))==1)
    {
        cout<<"Epsilon failure!"<<endl; //stuck in forever loop, so redraw
        //double randMax=exp(-1.0*lambda*pow(10,-7)); //This ensures that you don't have a machine epsilon problem. This value should be near one.
        (*randOut)=-1.0/lambda*log(rand()/(double)RAND_MAX);
    }
    while(isinf((*randOut)))
    {
        cout<<"Infinity failure! "<<lambda; getchar();
        (*randOut)=-1.0/lambda*log(rand()/(double)RAND_MAX);
    }
}

//Class functions
    // Metapopulation functions
metapopulation::metapopulation()
{
    demes=NULL; betaMtx=NULL;
}
metapopulation::~metapopulation()
{

}
void metapopulation::initializeMeta()
{
    //Set nE, nV (HARD CODED)
    nE=88;
    nV=8;
    //Initialize demes
    demes= new population[2];
    demes[0].initializePop((*this));
    demes[1].initializePop((*this));
    //Initialize rates and event matrix
    ratesVec=new double[nE];
    calculateRates();
    //cout<<rateTotal<<endl;
    initializeEventsMtx();
    //Initialize demes
}
void metapopulation::calculateRates()
{
    //State order: H11, H21, P11, P21, H12, H22, P12, P22
    //state[(p-1)*4+(s-1)*2+(g-1)] H->s=1, P->s=2
    int s, ct;
    rateTotal = 0;

    //Calculate d
    double d, D;
    dVec = new double[12];
    ct=0;
    for (int u=0;u<2;u++){
        for(int i=0;i<2;i++){
            for(int j=0;j<2;j++){
                dVec[ct]=alpha*betaMtx[i][j]*(double)demes[u].state[i]*(double)demes[u].state[2+j];
                d+=dVec[ct];
                ct++;
            }
            dVec[ct]=delta*(double)demes[u].state[i];
            d+=dVec[ct];
            ct++;
        }
    }
    D=2*kappaH/d;
    ct=0;
    for (int u=0;u<2;u++){

        //Host demography
        s=1;
        for(int i=0;i<2;i++)
        {
            for(int k=0;k<2;k++)
            {
                //delta H[i,u] H[k,u]/kappaH
                ratesVec[ct]=D*delta*(double)demes[u].state[(s-1)*2+i]*(double)demes[u].state[(s-1)*2+k]/(double)kappaH;
                rateTotal+=ratesVec[ct];
                ct++;
            }
        }
        //Parasite demography
        s=2;
        for(int j=0;j<2;j++)
        {
            for(int l=0;l<2;l++)
            {
                //gamma P[j,u] P[l,u]/kappaP
                ratesVec[ct]=D*gamma*(double)demes[u].state[(s-1)*2+j]*(double)demes[u].state[(s-1)*2+l]/(double)kappaP;
                rateTotal+=ratesVec[ct];
                ct++;
            }
        }
        //Infections
        for(int i=0;i<2;i++)
        {
            for(int j=0;j<2;j++)
            {
                for(int k=0;k<2;k++)
                {
                    for(int l=0;l<2;l++)
                    {
                        //alpha beta[i,j] H[i,u] P[j,u] H[k,u]/kappaH P[l,u]/kappaP
                        ratesVec[ct]=D*alpha*betaMtx[i][j]*(double)demes[u].state[(1-1)*2+i]*(double)demes[u].state[(2-1)*2+j]*
                        ((double)demes[u].state[(1-1)*2+k]/(double)kappaH)*((double)demes[u].state[(2-1)*2+l]/(double)kappaP);
                        rateTotal+=ratesVec[ct];
                        ct++;  
                    }                
                }
            }
        }
        //Host migration
        s=1;
        for(int i=0;i<2;i++) //host that migrates from pop u
        {
            for(int k=0;k<2;k++) //host that is born in pop u
            {
                for(int m=0;m<2;m++) //host that dies in pop u'
                {
                //phi H[i,u] H[k,u]/kappaH H[m,v]/kappaH
                ratesVec[ct]=D*phi*(double)demes[u].state[(s-1)*2+i]*((double)demes[u].state[(s-1)*2+k]/(double)kappaH)*((double)demes[(u+1)%2].state[(s-1)*2+m]/(double)kappaH);
                rateTotal+=ratesVec[ct];
                ct++;
                }
            }
        }
        //Parasite migration
        s=2;
        for(int j=0;j<2;j++) //parasite that migrates from pop u
        {
            for(int l=0;l<2;l++) //parasite that is born in pop u
            {
                for(int n=0;n<2;n++) //parasites that dies in pop u'
                {
                //omega P[j,u] P[l,v]/kappaP P[n,u]/kappaP
                ratesVec[ct]=D*omega*(double)demes[u].state[(s-1)*2+j]*((double)demes[u].state[(s-1)*2+l]/(double)kappaP)*((double)demes[(u+1)%2].state[(s-1)*2+n]/(double)kappaP);
                rateTotal+=ratesVec[ct];
                ct++;
                }
            }
        }
        //Host mutation
        s=1;
        for(int i=0;i<2;i++)
        {
            //mu H[i,u]
            ratesVec[ct]=D*mu*(double)demes[u].state[(s-1)*2+i];
            rateTotal+=ratesVec[ct];
            ct++;
        }
        //Parasite mutation
        s=2;
        for(int j=0;j<2;j++)
        {
            //nu P[j,u]
            ratesVec[ct]=D*nu*(double)demes[u].state[(s-1)*2+j];
            rateTotal+=ratesVec[ct];
            ct++;
        }
    //if(u==0){cout<<"Pop 1 total rates:"<<rateTotal<<endl;}else{cout<<"Pop 1 and 2 total rates:"<<rateTotal<<endl;}
    }
}
void metapopulation::initializeEventsMtx()
{
    int s,ct;
    eventsMtx=new double*[nE];
    for (int e=0;e<nE;e++)
    {
        eventsMtx[e]=new double[nV];
        for (int v=0;v<nV;v++){
            eventsMtx[e][v]=0;
        }
    }
    ct=0;
    for (int u=0;u<2;u++){
    
        s=1;
        //Host demography (4)
        for(int i=0;i<2;i++)
        {
            for(int k=0;k<2;k++)
            {
                //i dies and k is born
                eventsMtx[ct][u*4+(s-1)*2+i]--;
                eventsMtx[ct][u*4+(s-1)*2+k]++;
                ct++;
            }
        }
        s=2;
        //Parasite demography (4)
        for(int j=0;j<2;j++)
        {
            for(int l=0;l<2;l++)
            {
                //j dies and l is born
                eventsMtx[ct][u*4+(s-1)*2+j]--;
                eventsMtx[ct][u*4+(s-1)*2+l]++;
                ct++;
            }
        }
        //Infections (16)
        for(int i=0;i<2;i++)
        {
            for(int j=0;j<2;j++)
            {
                for(int k=0;k<2;k++)
                {
                    for(int l=0;l<2;l++)
                    {

                        //i dies and k is born
                        s=1;
                        eventsMtx[ct][u*4+(s-1)*2+i]--;
                        eventsMtx[ct][u*4+(s-1)*2+k]++;
                        //j is born and l dies
                        s=2;
                        eventsMtx[ct][u*4+(s-1)*2+j]++;
                        eventsMtx[ct][u*4+(s-1)*2+l]--;
                        ct++;

                    }                
                }
            }
        }
        s=1;
        //Host migration (8)
        for(int i=0;i<2;i++)
        {
            for(int k=0;k<2;k++)
            {
                for(int m=0;m<2;m++)
                {
                    //Host H[i,u] moves from pop u to pop u'
                    eventsMtx[ct][u*4+(s-1)*2+i]--;
                    eventsMtx[ct][((u+1)%2)*4+(s-1)*2+i]++;
                    //Birth of host type H[k,u] to replace departure of H[i,u]
                    eventsMtx[ct][u*4+(s-1)*2+k]++;
                    //Death of host type H[m,u'] following arrival of host H[i,u]
                    eventsMtx[ct][((u+1)%2)*4+(s-1)*2+m]--;
                    ct++;
                }
            }
        }
        s=2;
        //Parasite migration (4)
        for(int j=0;j<2;j++)
        {
            for(int l=0;l<2;l++)
            {
                for (int m=0;m<2;m++)
                {
                    //Parasite P[j,u] moves from pop u to pop u'
                    eventsMtx[ct][u*4+(s-1)*2+j]--;
                    eventsMtx[ct][((u+1)%2)*4+(s-1)*2+j]++;
                    //Birth of parasite type P[l,u] to replace departure of P[j,u]
                    eventsMtx[ct][u*4+(s-1)*2+l]++;
                    //Death of parasite type P[m,u'] following arrival of parasite P[j,u]
                    eventsMtx[ct][((u+1)%2)*4+(s-1)*2+m]--;
                    ct++;
                }
            }
        }
        s=1;
         //Host mutation (2)
        for(int i=0;i<2;i++)
        {
            //Host H[i,u] mutates from type i to type k
            //cout<<u*4+(s-1)*2+i<<endl;
            //cout<<u*4+(s-1)*2+((i+1)%2)<<endl;
            eventsMtx[ct][u*4+(s-1)*2+i]--;
            eventsMtx[ct][u*4+(s-1)*2+((i+1)%2)]++;
            ct++;
        }
        s=2;
        //Parasite mutation (2)
        for(int j=0;j<2;j++)
        {
            //Parasite P[j,u] mutates from type j to type l
            //cout<<u*4+(s-1)*2+j<<endl;
            //cout<<u*4+(s-1)*2+((j+1)%2)<<endl;
            eventsMtx[ct][u*4+(s-1)*2+j]--;
            eventsMtx[ct][u*4+(s-1)*2+((j+1)%2)]++;
            ct++;
        }
    }

    //Printing
    out_Test.open("events.csv"); //printing out event matrix
    for (int e=0;e<nE;e++) 
    {
    //cout<<e<<",";
     for (int v=0;v<nV;v++)
     {
        //cout<<eventsMtx[e][v]<<",";
        out_Test<<eventsMtx[e][v]<<",";
     }
    //cout<<endl;
    out_Test<<endl;
    }
    out_Test.close();
    
}
void metapopulation::simulation(char* argv1, int rep)
{
    sprintf(filename,"./Data_%s/data_%s_%d.csv",argv1,argv1,rep);
    out_Data.open(filename);

    sprintf(filename,"./Data_%s/dataN_%s_%d.csv",argv1,argv1,rep);
    out_DataN.open(filename);
    
    demes[0].resetPop((*this));
    demes[1].resetPop((*this));

    double t=0; //current time
    double deltat=0; //timestep
    double randomNum;
    int e;
    int etest;
    double totRate;
    
    out_Data<<t<<",";
    for(int s=0;s<4;s++) out_Data<<demes[0].state[s]<<",";
    for(int s=0;s<3;s++) out_Data<<demes[1].state[s]<<",";
    out_Data<<demes[1].state[3]<<endl;

    out_DataN<<t<<",";
    for(int s=0;s<4;s++) out_DataN<<demes[0].state[s]<<",";
    for(int s=0;s<3;s++) out_DataN<<demes[1].state[s]<<",";
    out_DataN<<demes[1].state[3]<<endl;

    //Draw exponentially distributed number
    randExp(rateTotal, &deltat);

    //for(int c=0;c<4;c++)
    while(t<tmax)
    {
        //Pick event
        randomNum=(rand()/(double)RAND_MAX)*rateTotal;
        e=0;
        totRate=ratesVec[0];
        while(randomNum>totRate)
        {
            e++;
            totRate+=ratesVec[e];
        }

       if(rep<0){

         vector<int> eventTab(88,0);
        for (int i=0; i<2000; i++){
            
            vector<double> ratesVec2(88, 0);
            for(int z=0;z<88;z++){
                ratesVec2[z]=ratesVec[z];
            }
            etest=chooseOutcome(ratesVec2);
            eventTab[etest]=eventTab[etest]+1;
            cout<<eventTab[etest]<<endl;
        }  

        vector<double> ratesVec2(88, 0);
        for(int z=0;z<88;z++){
                ratesVec2[z]=ratesVec[z];
        }
        etest=chooseOutcome(ratesVec2);

        if(rep<0){
            for(int i=0; i<88; i++){
                out_erand<<ratesVec[i]<<",";
            }
            for(int i=0; i<87; i++){
                out_erand<<eventTab[i]<<",";
            }
            out_erand<<eventTab[87]<<endl;
        }

       } else {
           vector<double> ratesVec2(88, 0);
            for(int z=0;z<88;z++){
                ratesVec2[z]=ratesVec[z];
            }
            etest=chooseOutcome(ratesVec2);
       }
        

        //cout<<etest<<endl;

        //out_erand<<etest<<",";

        if(e>=nE){

            for(int i=0; i<88; i++){
                cout<<ratesVec[i]<<", ";
            }
            cout<<endl;
            cout<<"e"<<endl;
            cout<<e<<endl;
            cout<<"Random number"<<endl;
            cout<<randomNum<<endl;
            cout<<"totRate: "<<endl;
            cout<<totRate<<endl;
            cout<<"rateTotal: "<<endl;
            cout<<rateTotal<<endl;
            cout<<"error in event"<<endl;
            getchar();
        }
        //out_Test<<e<<endl;
        //Perform event
        for(int s=0;s<4;s++){
            demes[0].state[s]+=eventsMtx[e][s];
            demes[1].state[s]+=eventsMtx[e][4+s];
        }
        //Perform events for neutral case
        neutralEvent(e);

        //Increment t
        t+=deltat;

        //Printing
        int v1;
        v1 = rand() % 1000;
        if (v1 == 0){
        //if (t < 0.2){

            //cout<<t<<",";
             //for(int i=0; i<88; i++){
                //cout<<ratesVec[i]<<",";
            //}
            //cout <<endl;
            //Printing
            out_Data<<t<<",";
            for(int s=0;s<4;s++) out_Data<<demes[0].state[s]<<",";
            for(int s=0;s<3;s++) out_Data<<demes[1].state[s]<<",";
            out_Data<<demes[1].state[3]<<endl;

            out_DataN<<t<<",";
            for(int s=0;s<4;s++) out_DataN<<demes[0].stateN[s]<<",";
            for(int s=0;s<3;s++) out_DataN<<demes[1].stateN[s]<<",";
            out_DataN<<demes[1].stateN[3]<<endl;

        }

        //Update rate vector
        calculateRates();
        //cout<<rateTotal<<endl;

        //Draw exponentially distributed number
        randExp(rateTotal, &deltat);
        
    }
    out_Data.close();
    out_DataN.close();
}

void metapopulation::neutralEvent(int event)
{
    int* temp;
    temp=new int[8];
    for(int i=0;i<8;i++) temp[i]=0;
    int e2=event%44;
    int p=0;
    if (event > 43){p=1;}
    double randH, randP, randH2, randP2, randH3, randP3;
    randH=rand()/(double)RAND_MAX;
    randP=rand()/(double)RAND_MAX;
    randH2=rand()/(double)RAND_MAX;
    randP2=rand()/(double)RAND_MAX;
    randH3=rand()/(double)RAND_MAX;
    randP3=rand()/(double)RAND_MAX;

    //out_erand<<randH<<","<<randH2<<","<<randH3<<","<<randP<<","<<randP2<<","<<randP3<<endl;

    if(e2 < 4) //host demographic event
    {
        //Choose host to be born
        if(randH<((double)demes[p].stateN[0]/(double)kappaH)){temp[p*4+0]++;} else {temp[p*4+1]++;}
        //Choose host to die
        if(randH2<((double)demes[p].stateN[0]/(double)kappaH)){temp[p*4+0]--;} else {temp[p*4+1]--;}
    } 
    else if (e2 < 8)//parasite demographic event
    {
        //Choose parasite to be born
        if(randP<((double)demes[p].stateN[2]/(double)kappaP)){temp[p*4+2]++;} else {temp[p*4+3]++;}
        //Choose parasite to die
        if(randP2<((double)demes[p].stateN[2]/(double)kappaP)){temp[p*4+2]--;} else {temp[p*4+3]--;}   
    }
    else if (e2 < 24)//infection event
    {

        //Choose host to be born
        if(randH<((double)demes[p].stateN[0]/(double)kappaH)){temp[p*4+0]++;} else {temp[p*4+1]++;}
        //Choose host to die
        if(randH2<((double)demes[p].stateN[0]/(double)kappaH)){temp[p*4+0]--;} else {temp[p*4+1]--;}
        //Choose parasite to be born
        if(randP<((double)demes[p].stateN[2]/(double)kappaP)){temp[p*4+2]++;} else {temp[p*4+3]++;}
        //Choose parasite to die
        if(randP2<((double)demes[p].stateN[2]/(double)kappaP)){temp[p*4+2]--;} else {temp[p*4+3]--;}   

    }
    else if (e2 < 32)//host migration
    {
        //Choose host to move from focal population
        if(randH<((double)demes[p].stateN[0]/(double)kappaH)){
            temp[p*4+0]--;
            temp[((p+1)%2)*4+0]++;
        } else {
            temp[p*4+1]--;
            temp[((p+1)%2)*4+1]++;
        }

        //Choose host to be born in focal population
        if(randH2<((double)demes[p].stateN[0]/(double)kappaH)){
            temp[p*4+0]++;
        } else {
            temp[p*4+1]++;
        }

        //Choose host to die in non-focal population
        if(randH3<((double)demes[(p+1)%2].stateN[0]/(double)kappaH)){
            temp[((p+1)%2)*4+0]--;
        } else {
            temp[((p+1)%2)*4+1]--;
        }
    }
    else if (e2 < 40)//parasite migration
    {
        //Choose parasite to move from focal population
        if(randP<((double)demes[p].stateN[2]/(double)kappaP)){
            temp[p*4+2]--;
            temp[((p+1)%2)*4+2]++;
        } else {
            temp[p*4+3]--;
            temp[((p+1)%2)*4+3]++;
        }

        //Choose parasite to be born in focal population
        if(randP2<((double)demes[p].stateN[2]/(double)kappaP)){
            temp[p*4+2]++;
        } else {
            temp[p*4+3]++;
        }

        //Choose parasite to die in non-focal population
        if(randP3<((double)demes[(p+1)%2].stateN[2]/(double)kappaP)){
            temp[((p+1)%2)*4+2]--;
        } else {
            temp[((p+1)%2)*4+3]--;
        }
    }
    else if (e2 < 42)//host mutation
    {
        //Choose host to mutate in focal population
        if(randH<((double)demes[p].stateN[0]/(double)kappaH)){
            temp[p*4+0]--;
            temp[p*4+1]++;
        } else {
            temp[p*4+1]--;
            temp[p*4+0]++;
        }
    }
    else if (e2 < 44)//parasite mutation
    {
        //Choose parasite to mutate in focal population
        if(randP<((double)demes[p].stateN[2]/(double)kappaP)){
            temp[p*4+2]--;
            temp[p*4+3]++;
        } else {
            temp[p*4+3]--;
            temp[p*4+2]++;
        }
    }
    else {
        cout<<"Neutral event error"<<endl;
        getchar();
    }

    //Perform event
    for(int s=0;s<4;s++){
        demes[0].stateN[s]+=temp[s];
        demes[1].stateN[s]+=temp[4+s];
    }
}
void metapopulation::initializeRep()
{
    //Reset state back to 50:50
    for(int d=0;d<2;d++)
    {
        demes[d].state[0]=kappaH/2;
        demes[d].state[1]=kappaH - demes[d].state[0];
        demes[d].state[2]=kappaP/2;
        demes[d].state[3]=kappaP - demes[d].state[2];
    }

}


population::population() //Don't return anything, so don't need to specify void, int, etc.
{
    state=NULL;
}
population::~population() //If created something in constructor, then delete in destructor; otherwise, leave empty
{

}
void population::initializePop(metapopulation metaPop)
{
    if(metaPop.kappaH%2==1) out_Pars<<"Note kappaH is odd"<<endl;
    if(metaPop.kappaP%2==1) out_Pars<<"Note kappaP is odd"<<endl;

    state=new int[4];
    state[0]=metaPop.kappaH/2;
    state[1]=metaPop.kappaH - state[0];
    state[2]=metaPop.kappaP/2;
    state[3]=metaPop.kappaP - state[2];
    stateN=new int[4];
    stateN[0]=metaPop.kappaH/2;
    stateN[1]=metaPop.kappaH - state[0];
    stateN[2]=metaPop.kappaP/2;
    stateN[3]=metaPop.kappaP - state[2];
}
void population::resetPop(metapopulation metaPop)
{
    //Reset counts of individuals to an allele frequency of 1/2
    state[0]=metaPop.kappaH/2;
    state[1]=metaPop.kappaH - state[0];
    state[2]=metaPop.kappaP/2;
    state[3]=metaPop.kappaP - state[2];

    //Reset counts of individuals to an allele frequency of 1/2
    stateN[0]=metaPop.kappaH/2;
    stateN[1]=metaPop.kappaH - state[0];
    stateN[2]=metaPop.kappaP/2;
    stateN[3]=metaPop.kappaP - state[2];
}

//To Do:
// Iterate dynamics
// Save answer
// Iterate over replicates
