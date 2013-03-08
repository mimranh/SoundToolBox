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
% sofia_radial.cpp (Internal C++ Header)
% 
% Copyright (C)2011 by bBrn - benjamin Bernschütz, rockzentrale 'AT' me.com  
%                         and Nils Peters, nils 'AT' icsi.berkeley.edu 
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 


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
%     Nils Peters                nils 'at' icsi.berkeley.edu
% 
**********************************************************************/

#include "sofia_radial.h"

#ifndef M_PI
	#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
#endif

// spherical functions     
std::complex<double> spbessel(int n, double kr) 
{
    std::complex<double> jn, prefactor; 
    double np;
    np=(double)n+0.5;
    prefactor=std::complex<double>(sqrt(M_PI/(2*kr)),0);
    jn = prefactor*std::complex<double>(boost::math::cyl_bessel_j(np,kr),0);  
    return jn;
}

std::complex<double> dspbessel(int n, double kr) 
{
    std::complex<double> djn;    
    djn = spbessel(n-1,kr)-spbessel(n,kr)*(std::complex<double>(n+1,0)/std::complex<double>(kr,0));   
    return djn;
}

std::complex<double> spneumann(int n, double kr) 
{
    std::complex<double> yn, prefactor;
    double np;
    np=(double)n+0.5;
    prefactor=std::complex<double>(sqrt(M_PI/(2*kr)),0);
    yn = prefactor*std::complex<double>(boost::math::cyl_neumann(np,kr),0);  
    return yn;
}   

std::complex<double> dspneumann(int n, double kr) 
{
    std::complex<double> dyn;      
    dyn = spneumann(n-1,kr)-spneumann(n,kr)*(std::complex<double>(n+1,0)/std::complex<double>(kr,0));  
    return dyn;
}   

std::complex<double> sphankel(int n, double kr) 
{
    std::complex<double> hn;  
    hn = spbessel(n,kr)-std::complex<double>(0,1)*spneumann(n, kr);
    return hn;
}   

std::complex<double> dsphankel(int n, double kr) 
{
    std::complex<double> dhn;   
    dhn = std::complex<double>(0.5,0)*(sphankel(n-1,kr)-sphankel(n+1,kr)-sphankel(n,kr)/kr);
    return dhn;
}  


//--- bn 

std::complex<double> bn_openP(int n, double krm)
{
    using namespace std;
    complex<double> bn;            
    bn  = spbessel(n,krm);
    return bn;
}

std::complex<double> bn_openPG(int n, double krm)
{
    using namespace std;
    complex<double> jn, diffjn, bn;     
    bn = complex<double>(0.5,0)*(spbessel(n,krm)-complex<double>(0,1)*dspbessel(n,krm));
    //                   ! 1/2
    return bn;
}

std::complex<double> bn_rigidP(int n, double krm, double krs)
{
    using namespace std;    
    complex<double> bn;
    bn = spbessel(n,krm)-(dspbessel(n,krs)/dsphankel(n,krs))*sphankel(n,krm);
    return bn;
}

std::complex<double> bn_rigidPG(int n, double krm, double krs) 

/* Rerence for Filter design for rigid sphere with cardioid microphones:
   P. Plessas, F. Zotter: Microphone arrays around rigid spheres for spatial recording and holography, DAGA 2010
   krm: for mic radius, krs: for sphere radius
   Implementation by Nils Peters, November 2011 */
{
    using namespace std;    
    complex<double> bn;
    bn = (spbessel(n,krm)-complex<double>(0,1)*dspbessel(n,krm) + (complex<double>(0,1)*dsphankel(n,krm)-sphankel(n,krm)) * (dspbessel(n,krs)/dsphankel(n,krs)));
    return bn;
}

std::complex<double> bn_dualOpenP(int n, double kr1, double kr2) 

/* Reference: Rafaely et al, 
   High-resolution plane-wave decomposition in an auditorium using a dual-radius scanning spherical microphone array
   JASA, 122(5), 2007

   kr1, kr2 are the kr values of the two different microphone spheres
   Implementation by Nils Peters, November 2011*/

{
    using namespace std;    
    complex<double> bn1, bn2;
    
    bn1 = bn_openP(n,kr1);
    bn2 = bn_openP(n,kr2);
    
    if (abs(bn1) >= abs(bn2))
        return bn1;
    else
        return bn2;	
}

std::complex<double> bn_npf(int n, double krm, double krs, int ac)
{
    std::complex<double> bnval;
    switch(ac)
                    {
                     case 0: bnval = bn_openP(n, krm);
                             break;

                     case 1: bnval = bn_openPG(n, krm);
                             break;

                     case 2: bnval = bn_rigidP(n, krm, krs);
                             break; 

					 case 3: bnval = bn_rigidPG(n, krm, krs);
		                     break;
                             
                     case 4: bnval = bn_dualOpenP(n, krm, krs);
		                     break;  	                    
                    }     
    
    return bnval;
}


std::complex<double> bn(int n, double krm, double krs, int ac)
{    
    std::complex<double> bnval;
    bnval = bn_npf(n, krm, krs, ac) * std::complex<double>(4*M_PI,0) * pow((std::complex<double>(0,1)),n);    
    return bnval;
}




