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

//Classes (difference between C and C++)
class population; //Declare population
class metapopulation
{
    public:
            //Variables
            int kappaH,kappaP; //Maximum population size
            double phi,omega; //Migration rates for host and parasite
            double **betaMtx; //Vector of vectors for probability infection
            double alpha; //Mortality given infection
            double delta,gamma; //Extrinsic mortality
            population *demes;

            //Functions
            metapopulation();//Constructor function
            ~metapopulation();//Destructor function
            void initializeMeta();//Need to give output type; not passing pointer, so can't values in class
            void simulation();
        private:
};
class population
{
    public:
        //Variables
        int nH,nP; //Realized population size
        bool *popH,*popP; //Genotypes (Could be vector<bool> pop, but wouldn't have fixed length)
        double WbarH,WbarP; //Mean fitness of host or parasite in population
        double pH,pP; //Allele frequencies
        metapopulation *metaPopPtr;

        //Functions
        population();//Constructor function
        ~population();//Destructor function
        void initializePop(metapopulation metaPop);//Need to give output type; not passing pointer, so can't values in class
        void calcAF();

    private:
};

//(Global) Function declarations (say what functions are going to exist)
void input(metapopulation *metaPop1Ptr);

//Main function: C++ looks for main function and runs main function
int main(int argc, char**argv) //return integer (int), first argument is integer called argc, second is array of characters called argv
{
    //Random seed
    unsigned int seed=(unsigned int)time(NULL);
    srand(seed);//Initializing random seed

    //Variables
    metapopulation metaPop;
    sprintf(filename,"./parameters.csv"); //Makes the input file name
    out_Pars.open(filename);

    //Setting up output files

    //Running program
    input(&metaPop); //& gives arrow, * before gives thing pointed to by arrow; (*x) is opposite of &x
    metaPop.initializeMeta();
    metaPop.demes[0].initializePop(metaPop);
    //metaPop.demes[1].initializePop(metaPop);

    //End program
    out_Pars.close();

    return 0;
}

//Function definitions
void input(metapopulation *metaPopPtr)
{
    //Input parameter values from a file (ones that are pretty constant)
    ifstream inFile;
    sprintf(filename,"./input.txt"); //Makes the input file name
    inFile.open(filename);
    string myString;
    myString="Null";
    // Look for kappaH in file
    while(myString != ("(kappaH):") && inFile.good()) {inFile>>myString;}
    //Once found kappaH, read in next word
    inFile>>(*metaPopPtr).kappaH; //Go to thing pointed to by *metaPopPtr ((*metaPopPtr) is THING POINTED TO by pointer) and take element kappaH of thing pointed to
    //Read in kappaH
    cout<<"kappaH: "<<(*metaPopPtr).kappaH<<endl;
    out_Pars<<"kappaH:, "<<(*metaPopPtr).kappaH<<endl;
    
    // Look for kappaP in file
    while(myString != ("(kappaP):") && inFile.good()) {inFile>>myString;}
    //Once found kappaP, read in next word
    inFile>>(*metaPopPtr).kappaP; //Go to thing pointed to by *pop1Ptr and take element kappaH of thing pointed to
    //Read in kappaP
    cout<<"kappaP: "<<(*metaPopPtr).kappaP<<endl;
    out_Pars<<"kappaP:, "<<(*metaPopPtr).kappaP<<endl;

    // Look for phi in file
    while(myString != ("(phi):") && inFile.good()) {inFile>>myString;}
    //Once found phi, read in next word
    inFile>>(*metaPopPtr).phi; //Go to thing pointed to by *pop1Ptr and take element kappaH of thing pointed to
    //Read in phi
    cout<<"phi: "<<(*metaPopPtr).phi<<endl;
    out_Pars<<"phi:, "<<(*metaPopPtr).phi<<endl;

    // Look for omega in file
    while(myString != ("(omega):") && inFile.good()) {inFile>>myString;}
    //Once found omega, read in next word
    inFile>>(*metaPopPtr).omega; //Go to thing pointed to by *pop1Ptr and take element kappaH of thing pointed to
    //Read in omega
    cout<<"omega: "<<(*metaPopPtr).omega<<endl;
    out_Pars<<"omega:, "<<(*metaPopPtr).omega<<endl;

    (*metaPopPtr).betaMtx=new double*[2];
    for(int r=0;r<2;r++) (*metaPopPtr).betaMtx[r]=new double[2];

    // Look for beta in file
    while(myString != ("(beta):") && inFile.good()) {inFile>>myString;}
    //Once found omega, read in next word
    cout<<"betaMtx: "<<endl;
    out_Pars<<"betaMtx: "<<endl;
    for(int r=0;r<2;r++){
        for(int c=0;c<2;c++){
            inFile>>(*metaPopPtr).betaMtx[r][c];
            cout<<(*metaPopPtr).betaMtx[r][c]<<",";
            out_Pars<<(*metaPopPtr).betaMtx[r][c]<<",";
        }
        cout<<endl;
        out_Pars<<endl;
    }
}

//Class functions
metapopulation::metapopulation()
{
    demes=NULL; betaMtx=NULL;
}
metapopulation::~metapopulation()
{

}
population::population() //Don't return anything, so don't need to specify void, int, etc.
{
    popH=NULL; popP=NULL;
}
population::~population() //If created something in constructor, then delete in destructor; otherwise, leave empty
{

}
void population::initializePop(metapopulation metaPop){

    nH=metaPop.kappaH; nP=metaPop.kappaP;
    popH=new bool[nH];
    popP=new bool[nP];
    cout<<"Hosts: "<<endl;
    for(int i=0;i<nH;i++){
        if((rand()/(double)RAND_MAX)<0.5){
            popH[i]=0;
        } else {
            popH[i]=1;
        }
        cout<<popH[i]<<",";
    }
    cout<<endl;

    cout<<"Parasites: "<<endl;
    for(int i=0;i<nP;i++){
        if((rand()/(double)RAND_MAX)<0.5) popP[i]=0;
        else popP[i]=1;
        cout<<popP[i]<<",";
    }
    cout<<endl;

    cout<<"Initial allele frequencies: "<<endl;
    calcAF();

}
void metapopulation::initializeMeta(){

    demes=new population[2];

}
void population::calcAF(){
    pH=0;pP=0;
    for(int i=0;i<nH;i++){
        pH+=popH[i];
    }
    for(int i=0;i<nP;i++){
        pP+=popP[i];
    }
    pH/=(double)nH; //type casting
    pP/=(double)nP; //type casting
    cout<<"pH: "<<pH<<" pP: "<<pP<<endl;
}

//To Do:
// Calculate fitness (expected fitness)
// Ailene will think about it... simulate gillespie algorithm?
