/***********************************************************************
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
% S/T/C Fast Spatial Fourier Transform Core 
% 
% Copyright (C)2011 by bBrn - benjamin Bernschütz  
%                            rockzentrale 'AT' me.com
%                            +49 171 4176069 Germany  
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
%  
% Pnm = sofia_stc(N, fftData, grid)
% ------------------------------------------------------------------------
% Pnm      Spatial Fourier Coefficients 
%          Columns : nm coeff
%          Rows    : FFT bins   
% ------------------------------------------------------------------------   
% N        Maximum transform order
% 
% fftdata  Frequency domain sounfield data, e.g. from sofia_fdt()) 
%          Columns : number of spatial sampling position 
%          Rows    : FFT bins (complex sound pressure data)
% 
%          ! WARNING Cast fftdata to double if core crashes:
%            fftData = cast(fftData ,'double');
% 
% grid     Sample grid configuration
%          Columns : s=1...S spatial positions
%          Rows    : [AZ_s EL_s GridWeight_s]
%          AZ in [0...2pi] and EL [0...pi] in RAD
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
*********************************************************************/

//#define NOBOOST

#include "mex.h"
#include <iostream>
#include <complex>

#ifndef M_PI
	#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
#endif
        
#ifndef NOBOOST
        #include <boost/math/special_functions/spherical_harmonic.hpp>
#endif

//WARNING: Typecast MATLAB FFT-Output: FFTData=cast(FFTData,'double');

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
   using namespace std;
   #ifndef NOBOOST
    using namespace boost::math;
   #endif
   double *grid, *FFTdataReal, *FFTdataImag;   
   int    N, NMLocatorSize;
   int    n, m, j;
   long long unsigned int ctr, ctrb, ctrc, f;
   long unsigned int FFTBlocklength;
   int    numberOfFFTBlocks, numberOfSpatialPositionsInFFTBlock;   
   int    numberOfGridPoints, numberOfGridInfos;
   double *AzimuthAngles,*ElevationAngles, *GridWeights;
   double *ReturnReal,*ReturnImag; 
   complex<double> **OutputArray;
   complex<double> SHresult;
   
   #ifndef DEBUG
      mexPrintf("SOFiA S/T/C - Spatial Transform Core R11-1220\n");
   #endif
           
   if(nrhs!=3) {mexErrMsgIdAndTxt("SOFiA:STC:notEnoughInputs",
                 "3 Inputs required: N, FFT-Blocks, GRID-Data");}
   
   //Get ARGS
    N       = (int)mxGetScalar(prhs[0]);
        
    FFTdataReal = mxGetPr(prhs[1]);
    FFTdataImag = mxGetPi(prhs[1]);
    
    if(FFTdataImag==NULL)
        {mexErrMsgIdAndTxt("SOFiA:STC:InputArgError",
                 "FFTData: Complex Input Data expected.");}
    
    if (mxIsCell(prhs[1])==1) // Check if Input is a Cell Array
        {
            {mexErrMsgIdAndTxt("SOFiA:STC:InputArgError",
                 "FFTData: No Cell Array allowed.");}
            //numberOfFFTBlocks=mxGetN(prhs[1]);
            //mexPrintf("Multiple FFT Blocks requested: %d\n",numberOfFFTBlocks);                       
        }
    else
        {   
            numberOfFFTBlocks=1;
            numberOfSpatialPositionsInFFTBlock=mxGetM(prhs[1]);
            FFTBlocklength=mxGetN(prhs[1]);            
        }
   
    grid  = mxGetPr(prhs[2]);   
    numberOfGridPoints = mxGetM(prhs[2]);
    numberOfGridInfos  = mxGetN(prhs[2]);
    
    if(numberOfGridInfos<3)
        {mexErrMsgIdAndTxt("SOFiA:STC:InputArgError",
                 "GRID: Delivered Grid Data not valid.");}
        
    if(numberOfGridPoints!=numberOfSpatialPositionsInFFTBlock)
        {mexErrMsgIdAndTxt("SOFiA:STC:InputArgError",
                 "Number of Spatial Sampling Points in FFT/GRID is not consistent.");}

//     mexPrintf("Delivered FFT Taps           %ld\n",FFTBlocklength);
//     mexPrintf("Spatial Sampling Positions   %d\n",numberOfGridPoints);
//     mexPrintf("Transform Order              %d",N);
   
    
    
    //Create Dynamic Arrays (GRID)
    AzimuthAngles   = new double [numberOfGridPoints]; 
    ElevationAngles = new double [numberOfGridPoints]; 
    GridWeights     = new double [numberOfGridPoints]; 
    
    //Fill Cells GRID
    for(ctr=0;ctr<numberOfGridPoints;ctr++)
        {
            AzimuthAngles[ctr]   = grid[ctr];
            ElevationAngles[ctr] = grid[ctr+numberOfGridPoints]; 
            GridWeights[ctr]     = grid[ctr+2*numberOfGridPoints];                  
        }        
    
    NMLocatorSize=(N+1)*(N+1); //All n,m up to N included      
    
    try
     {
        //Allocate Output Array Size [(N+1)^2]*[FFTBins]
        OutputArray = new complex<double>*[NMLocatorSize];

        for(ctr = 0; ctr < NMLocatorSize; ctr++)		
           {OutputArray[ctr] = new complex<double>[FFTBlocklength];}


        //Initializate Output Array
        for(ctr = 0; ctr < NMLocatorSize; ctr++)
           {
            for(ctrb = 0; ctrb < FFTBlocklength; ctrb++)
               {
                    OutputArray[ctr][ctrb] = complex<double>(0,0);                
               }
           }
        
     }
     catch(...) 
     { mexErrMsgIdAndTxt("SOFiA:STC:OutArrayAllocation",
                 "Not able to allocate memory for the output Matrix. Maybe to large?");
     }
            
            
   
    //SHT TRANSFORM CORE -------------------------------------------------
        
    ctr=0;
    
     for(n=0; n<=N; n++)
     {
         for(m=-n; m<=n; m++)
         {             
             for(j=0; j<numberOfGridPoints; j++)
             {
                 SHresult=complex<double>(4*M_PI,0)*complex<double>(GridWeights[j],0);
                 
                  #ifndef NOBOOST
                          SHresult*=conj(spherical_harmonic(n,m,ElevationAngles[j],AzimuthAngles[j]));
                  #endif
                   
                 for(f=0; f<FFTBlocklength; f++)
                 {                   
                     OutputArray[ctr][f]=OutputArray[ctr][f]+SHresult*complex<double>(FFTdataReal[j+f*numberOfGridPoints],FFTdataImag[j+f*numberOfGridPoints]);                  
                 }                 
             }
             ctr++;                
         }
     }                

     
       
     //Declare return ARG Matrix 
     plhs[0] = mxCreateDoubleMatrix(NMLocatorSize,FFTBlocklength,mxCOMPLEX);
     ReturnReal = mxGetPr(plhs[0]);
     ReturnImag = mxGetPi(plhs[0]);
     
     //Return ARG Matrix Fill 
     ctrc=0;     
     for(ctr=0;ctr<FFTBlocklength;ctr++)
        {
        for(ctrb=0;ctrb<NMLocatorSize;ctrb++)
             {
                ReturnReal[ctrc]=real(OutputArray[ctrb][ctr]);
                ReturnImag[ctrc]=imag(OutputArray[ctrb][ctr]);
                ctrc++;
              }
        }         
     
      delete [] AzimuthAngles;
      delete [] ElevationAngles;
      delete [] GridWeights;
      for(ctr = 0; ctr < NMLocatorSize; ctr++) delete OutputArray[ctr];      
}

