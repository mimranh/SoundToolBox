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
% sofia_radial.h (Internal C++ Header)
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

#ifndef SOFIA_RADIAL_H
#define SOFIA_RADIAL_H

#include <boost/math/special_functions/bessel.hpp>

std::complex<double> spbessel(int n, double kr); 
std::complex<double> dspbessel(int n, double kr); 
std::complex<double> spneumann(int n, double kr); 
std::complex<double> dspneumann(int n, double kr); 
std::complex<double> sphankel(int n, double kr); 
std::complex<double> dsphankel(int n, double kr); 
std::complex<double> bn_openP(int n, double krm);
std::complex<double> bn_openPG(int n, double krm);
std::complex<double> bn_rigidP(int n, double krm, double krs);
std::complex<double> bn_rigidPG(int n, double krm, double krs); 
std::complex<double> bn_dualOpenP(int n, double kr1, double kr2);
std::complex<double> bn(int n, double krm, double krs, int ac);
std::complex<double> bn_npf(int n, double krm, double krs, int ac);

#endif