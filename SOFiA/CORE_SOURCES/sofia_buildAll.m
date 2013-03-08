% /// ASAR Research Group
%  
% Cologne University of Applied Sciences
% Technical University of Berlin
% Deutsche Telekom Laboratories
% WDR Westdeutscher Rundfunk
% 
% SOFiA sound field analysis
% 
% Automatic builder
% 
% Copyright (C)2011 by Nils Peters, nils 'at' icsi.berkeley.edu    
%                 (and Benny Bernschütz, rockzentrale 'at' me.com) 
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
% void = sofia_buildAll(configuration)
% 
% configuration: [] or 'DEBUG' 
%
% This file compiles the SOFiA C++ sources. Run the file inside the
% \CORE_SOURCES folder.  
%
% Remember: You need a C++ Compiler and need to configure MATLAB to 
%           use this compiler. Further the BOOST C++ ist required.
%           Read "howtocompile.txt" for more information.


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


function sofia_buildAll(configuration)
clc

if nargin == 0
    configuration = 'RELEASE';
end

disp(' ');
disp('*** building SOFiA ***');
disp(' ');

sourcefiles = dir('*.cpp');
sources=[];
for i=1:size(sourcefiles,1)
    sources{i}  = sourcefiles(i).name;
    target_m{i} = strrep(sourcefiles(i).name,'.cpp','.m');
end

cHeaderSources = dir(['HEADER',filesep(),'*.cpp']);
cHeaders=[];

for i=1:size(cHeaderSources,1)
    cHeaders = [cHeaders, 'HEADER',filesep(),cHeaderSources(i).name,' '];
end

cHeaderPath = [pwd(),filesep(),'HEADER'];


for k = 1:length(sources);
    eval(sprintf('disp(''building %s'');',char(sources(k))));
    % mex compiling
    eval(sprintf('mex %s %s -I%s -D%s -outdir ../;',char(sources(k)),char(cHeaders), char(cHeaderPath),char(configuration)));  
        
    % header file generation: writing the header into an extra .m file 
    eval(sprintf('fir=fopen(''%s'' ,''r'',''n'',''ISO-8859-1'');',char(sources(k)))); 
    eval(sprintf('fiw=fopen(''../%s'',''w'',''n'',''ISO-8859-1'');',char(target_m(k))));
    tline = fgetl(fir); tline = fgetl(fir);
    while 1
        tline = fgetl(fir);
        if isempty(tline)
            fprintf(fiw,' \n');
            break;
        elseif char(tline(1)) == '@'
            fprintf(fiw, ' \n');        
        elseif char(tline(1)) == '%'
            fprintf(fiw, '%s\n',char(tline));
        else                
            break;
        end;
    end;            
 fclose(fir);
 fclose(fiw);
 end;

disp(' ');    
disp('***     done       ***');
disp(' ');

