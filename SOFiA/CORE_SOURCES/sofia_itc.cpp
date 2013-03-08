/**********************************************************************
%  
% /// ASAR Research Group
%  
% Cologne University of Applied Sciences
% Technical University of Berlin
% Deutsche Telekom Laboratories
% WDR Westdeutscher Rundfunk
% 
% SOFiA sound field analysis
% 
% I/T/C Fast Inverse spatial Fourier Transform Core R11-1220
% 
% Copyright (C)2011 by bBrn - benjamin Bernschütz  
%                             rockzentrale 'AT' me.com
%                             +49 171 4176069 Germany  
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
% 
% p = sofia_itc(Pnm, angles, [N])
% ------------------------------------------------------------------------
% p      sound pressures (complex data)  
%        Columns : angles
%        Rows    : FFT bins   
% ------------------------------------------------------------------------
% Pnm     spatial Fourier coefficients (e.g. from SOFiA S/T/C)
%        Columns : nm coeff
%        Rows    : FFT bins   
% 
% angles  target angles [AZ1 EL1; AZ2 EL2; ... AZn ELn]
%        Columns : Angle Number 1...n
%        Rows    : AZ EL    
%          
% [N]     *** Optional: Maximum transform order 
%            If not specified the highest order available included in
%            the Pnm coefficients will be taken.
% 
% This is a pure ISFT core that does not involve extrapolation. 
% (=The pressures are refered to the original radius)
%
@ end of header
%
% CONTACT AND LICENSE INFORMATION:
% 
% /// ASAR Research Group 
%  
%     [1] Cologne University of Applied Sciences
%     [2] Technical University of Berlin 
%     [3] Deutsche Telekom Laboratories 
%     [4] WDR Westdeutscher Rundfunk 
% 
% SOFiA sound field analysis
% 
% Copyright (C)2011 bBrn - benjamin Bernschütz [1,2] et al.(§)   
% 
% Contact ------------------------------------
% Cologne University of Applied Sciences 
% Institute of Communication Systems
% Betzdorfer Street 2
% D-50679 Germany (Europe)
% 
% phone       +49 221 8275 -2496 
% cell phone  +49 171 4176069 
% mail        rockzentrale 'at' me.com 
% --------------------------------------------
% 
% This file is part of the SOFiA sound field analysis toolbox
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% 
% 
% (§) Christoph Pörschmann [1]   christoph.poerschmann 'at' fh-koeln.de
%     Sascha Spors         [2,3] sascha.spors 'at' telekom.de  
%     Stefan Weinzierl     [2]   stefan.weinzierl 'at' tu-berlin.de
% 
**********************************************************************/

#include "mex.h"
#include <iostream>
#include <complex>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#ifndef M_PI
	#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
   using namespace std;
   using namespace boost::math;
   double *grid, *PnmDataReal, *PnmDataImag;   
   int    N, Nmax, PnmDataLength;
   int    n, m, j;
   long long unsigned int ctr, ctrb, ctrc, f;
   long unsigned int FFTBlocklength;
   int    numberOfAngles, numberOfAngleInfos;
   double *AzimuthAngles,*ElevationAngles;
   double *ReturnReal,*ReturnImag; 
   complex<double> **OutputArray;
   complex<double> SHresult;
   
   #ifndef DEBUG
   mexPrintf("SOFiA I/T/C - Inverse spatial Transform Core R11-1220\n");
   #endif
   if(nrhs<2) {mexErrMsgIdAndTxt("SOFiA:ITC:notEnoughInputs",
                 "Minimum 2 Inputs required: p = sofia_itc(Pnm coefficients, angles, [N])");}
   
   //Get ARGS
         
    PnmDataReal = mxGetPr(prhs[0]);
    PnmDataImag = mxGetPi(prhs[0]);
    
    if(PnmDataImag==NULL)
        {mexErrMsgIdAndTxt("SOFiA:ITC:InputArgError",
                 "PnmData: Complex Input Data expected.");}
    
    FFTBlocklength=mxGetN(prhs[0]);            
      
    grid  = mxGetPr(prhs[1]);   
    numberOfAngles = mxGetM(prhs[1]);
    numberOfAngleInfos  = mxGetN(prhs[1]);
    
    if(numberOfAngleInfos<2)
        {mexErrMsgIdAndTxt("SOFiA:ITC:InputArgError",
         "Error: Delivered angles are not valid. Must consist of [AZ1 EL1; AZ2 EL2; ...; AZn ELn] pairs.");}
    
    PnmDataLength = mxGetM(prhs[0]);
    
    if(nrhs>2){
        N    = (int)mxGetScalar(prhs[2]); 
        Nmax = (int)sqrt((float)PnmDataLength)-1;
        if (N>Nmax)
            {mexErrMsgIdAndTxt("SOFiA:ITC:MaximumOrderExceed",
                 "Requested order too high: Maximum available order in Pnm coefficients exceeded.");}        
    }
    else {                           
        N = (int)sqrt((float)PnmDataLength)-1; 
        Nmax = N; 
    }                   
   
//     mexPrintf("Delivered FFT Taps      %ld\n",FFTBlocklength);
//     mexPrintf("Spatial Positions       %d\n",numberOfAngles);
//     mexPrintf("Transform Order         %d",N);   
    
    
    //Create Dynamic Arrays (angles)
    AzimuthAngles   = new double [numberOfAngles]; 
    ElevationAngles = new double [numberOfAngles]; 

    
    //Fill Cells GRID
    for(ctr=0; ctr<numberOfAngles; ctr++)
        {
            AzimuthAngles[ctr]   = grid[ctr];
            ElevationAngles[ctr] = grid[ctr+numberOfAngles];                           
        }        
       
    try
     {
        //Allocate Output Array Size [numberOfAngles]*[FFTBins]
        OutputArray = new complex<double>*[numberOfAngles];

        for(ctr = 0; ctr < numberOfAngles; ctr++)		
           {OutputArray[ctr] = new complex<double>[FFTBlocklength];}


        //Initializate Output Array
        for(ctr = 0; ctr < numberOfAngles; ctr++)
           {
            for(ctrb = 0; ctrb < FFTBlocklength; ctrb++)
               {
                    OutputArray[ctr][ctrb] = complex<double>(0,0);                
               }
           }        
     }
     catch(...) 
     { mexErrMsgIdAndTxt("SOFiA:ITC:OutArrayAllocation",
                 "Not able to allocate memory for the output Matrix. Maybe to large?");
     }
             
   
    //SHT TRANSFORM CORE -------------------------------------------------
        
    ctr=0;
   
     for(n=0; n<=N; n++)
     {
         for(m=-n; m<=n; m++)
         {               
             for(j=0; j<numberOfAngles; j++)
             {
                 //SHresult=complex<double>(1/(4*M_PI),0);                 
                 SHresult = spherical_harmonic(n,m,ElevationAngles[j], AzimuthAngles[j]);
                 for(f=0; f<FFTBlocklength; f++)
                 {                   
                     OutputArray[j][f] = OutputArray[j][f] + SHresult * complex<double>(PnmDataReal[ctr+f*PnmDataLength],PnmDataImag[ctr+f*PnmDataLength]);                  
                 }                 
             }
             ctr++;                
         }
     }                
     
        
    
     //Declare return ARG Matrix 
     plhs[0] = mxCreateDoubleMatrix(numberOfAngles,FFTBlocklength,mxCOMPLEX);
     ReturnReal = mxGetPr(plhs[0]);
     ReturnImag = mxGetPi(plhs[0]);
     
     //Return ARG Matrix Fill 
     ctrc=0;     
     for(ctr=0;ctr<FFTBlocklength;ctr++)
        {
        for(ctrb=0;ctrb<numberOfAngles;ctrb++)
             {
                ReturnReal[ctrc]=real(OutputArray[ctrb][ctr]);
                ReturnImag[ctrc]=imag(OutputArray[ctrb][ctr]);
                ctrc++;
              }
        }         

      delete [] AzimuthAngles;
      delete [] ElevationAngles;
      for(ctr = 0; ctr < numberOfAngles; ctr++) delete OutputArray[ctr];      
}


