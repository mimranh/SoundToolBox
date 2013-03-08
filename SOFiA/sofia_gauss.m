% /// ASAR Research Group
% 
%     Cologne University of Applied Sciences
%     Technical University of Berlin
%     Deutsche Telekom Laboratories
%     WDR Westdeutscher Rundfunk
% 
% SOFiA sound field analysis
% 
% SOFiA Gauss Grid R11-1220
% 
% Copyright (C)2011 by bBrn - benjamin Bernschütz 
%
% Core routines for gauss-legendre quadrature calculation:
% --------------------------------------------------------
%  Written by: Greg von Winckel - 04/13/2006
%  Contact: gregvw(at)math(dot)unm(dot)edu 
%  URL: http://www.math.unm.edu/~gregvw
% --------------------------------------------------------
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
%  
% [gridData, Npoints, Nmax] = sofia_gauss(AZnodes, ELnodes, plot)
% ------------------------------------------------------------------------     
%
% gridData           Gauss-Legendre quadrature including weigths(W):
%                    [AZ_1 EL_1 W_1;
%                     AZ_2 EL_2 W_2;
%                     ...
%                     AZ_n EL_n W_n]
%
% Npoints            Total number of nodes
% Nmax               Highest stable grid order  
%
% ------------------------------------------------------------------------
% 
% AZnodes            Number of azimutal nodes  [default = 10]
% ELnodes            Number of elevation nodes [default = 5]
% plot               Show a globe plot of the selected grid 
%                    0: Off, 1: On [default]
% 
% This function computes Gauss-Legendre quadrature nodes and weigths
% in the SOFiA/VariSphear data format.
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
%
% Core routines for gauss-legendre quadrature calculation:
% --------------------------------------------------------
%  Written by: Greg von Winckel - 04/13/2006
%  Contact: gregvw(at)math(dot)unm(dot)edu 
%  URL: http://www.math.unm.edu/~gregvw
% --------------------------------------------------------
%


function [gridData, Npoints, Nmax] = sofia_gauss(AZnodes, ELnodes, plot)

disp('SOFiA Gauss Grid R11-1220');

if nargin<3
    plot = true;
end

if nargin<2
    ELnodes = 5;
end

if nargin<1
    AZnodes = 10;
end

% Routines taken from VariSphear waveCapture
[~,T,P,W] = gauss_grid(1,ELnodes,AZnodes,1);

gridData = [P,T,W]; % Compose Grid Vector with Theta, Phi and Weights
gridData = sortrows(gridData,2);
gridData = sortrows(gridData,1);

i=1;
turnover=0;
while(1) %Sort VariSphear style
    if i>=size(gridData,1)
        break
    end
    c = find(gridData(:,1)==gridData(i,1));
    i = max(c)+1;
    if turnover == 1 
        gridData(c,:)=flipdim(gridData(c,:),1);
        turnover=0;
    else
        turnover=1;
    end
end

gridData(:,3)=gridData(:,3)/sum(gridData(:,3));

Npoints = size(gridData,1);
Nmax    = floor(sqrt(size(gridData,1)/2)-1);

if plot
    plot_grid(gridData);
end

function plot_grid(gridData)

[Xm,Ym,Zm]=sph2cart(gridData(:,1),gridData(:,2)-pi/2,1.01);

colormap Gray;

if size(Xm,1)>1500
    plot3(Xm,Ym,Zm,'marker','.','markerfacecolor','g','color','g','linestyle','none')
else
    plot3(Xm,Ym,Zm,'marker','o','markerfacecolor','g','color','g','linestyle','none')
end
axis off;
hold on;
grid off;
sphere;
axis equal;
rotate3d on;
light;
alpha(.8);
lighting phong;
camzoom(1.4);
hold off;



function [r,t,p,w]=gauss_grid(nr,nt,np,rad)

%SPHEREQUAD  Generate Gauss quadrature nodes and weights for numerically
% computing spherical volume integrals.
%
% [R,T,P,W]=SPHEREQUAD(NR,NT,NP,RAD) computes the product grid nodes in
% r, theta, and phi in spherical and the corresponding quadrature weights
% for a sphere of radius RAD>0. NR is the number of radial nodes, NT is
% the number of theta angle nodes in [0,pi], and NP is the number of phi 
% angle nodes in [0, 2*pi]. The sphere radius RAD can be set to infinity, 
% however, the functions to be integrated must decay exponentially with 
% radius to obtain a reasonable numerical approximation.
%
% Example 1: Infinite domain, theta independent
%
% f=@(R,T,P) exp(-R.^2.*(2+sin(P))); 
% [R,T,P,W]=spherequad(50,1,30,inf);
% Q=W'*f(R,T,P);
%
% Example 2: Sphere of radius 2, depends on all three
%
% f=@(R,T,P) sin(T.*R).*exp(-R.*sin(P));
% [R,T,P,W]=spherequad(24,24,24,2);
% Q=W'*f(R,T,P);
%
% Written by: Greg von Winckel - 04/13/2006
% Contact: gregvw(at)math(dot)unm(dot)edu 
% URL: http://www.math.unm.edu/~gregvw


[r,wr]=rquad(nr,2);         % radial weights and nodes (mapped Jacobi)

if rad==inf                 % infinite radius sphere
   
    wr=wr./(1-r).^4;        % singular map of sphere radius
    r=r./(1-r);
    
else                        % finite radius sphere
    
    wr=wr*rad^3;            % Scale sphere radius
    r=r*rad;
    
end

[x,wt]=rquad(nt,0); 
t=acos(2*x-1); wt=2*wt;     % theta weights and nodes (mapped Legendre)
p=2*pi*(0:np-1)'/np;        % phi nodes (Gauss-Fourier)
wp=2*pi*ones(np,1)/np;      % phi weights
[rr,tt,pp]=meshgrid(r,t,p); % Compute the product grid
r=rr(:); t=tt(:); p=pp(:);
w=reshape(reshape(wt*wr',nr*nt,1)*wp',nr*nt*np,1);


function [x,w]=rquad(N,k)

k1=k+1; k2=k+2; n=1:N;  nnk=2*n+k;
A=[k/k2 repmat(k^2,1,N)./(nnk.*(nnk+2))];
n=2:N; nnk=nnk(n);
B1=4*k1/(k2*k2*(k+3)); nk=n+k; nnk2=nnk.*nnk;
B=4*(n.*nk).^2./(nnk2.*nnk2-nnk2);
ab=[A' [(2^k1)/k1; B1; B']]; s=sqrt(ab(2:N,2));
[V,X]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
[X,I]=sort(diag(X));    
x=(X+1)/2; w=(1/2)^(k1)*ab(1,2)*V(1,I)'.^2;