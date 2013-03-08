% SOFiA example 5: Spatiotemporal resolution
% SOFiA Version  : R11-1220
% Array Dataset  : R11-1018

  clear all
  clc
  
% Read VariSphear dataset
% !!! LOCATE THE FOLDER: "EXAMPLE3_SpatioTemporal" 

  timeData = sofia_readVSAdata(); 
 
% figure(2) %Enable to get an IR overview
% clf();
% area(timeData.irOverlay');

  fftOversize = 2;
  startSample = [400 1100 1800 2500]; %Examplary values (Enable area plot)
  blockSize   = 256;

  figure(1)
  clf();
  
for ctr=1:size(startSample,2) 
      
% Here we have simply put a loop around, because this is good for
% understanding this experiment. The loop content can easily be optimized. 
   
  % Transform time domain data to frequency domain and generate kr-vector
  
  [fftData, kr, f] = sofia_fdt(timeData, fftOversize, startSample(ctr), startSample(ctr)+blockSize);
    
    % Spatial Fourier Transform
    
    Nsft = 5;
    Pnm  = sofia_stc(Nsft, fftData, timeData.quadratureGrid);

    % Radial Filters for a rigid sphere array 
    
    Nrf      = Nsft;    % radial filter order              
    maxAmp   = 10;      % Maximum modal amplification in [dB]
    ac       = 2;       % Array configuration: 2 = Rigid Sphere 
    dn       = sofia_mf(Nrf , kr, ac, maxAmp); % radial filters 

    % Make MTX  
    
    Nmtx = Nsft;
    krIndex = 90;   % Choose the kr-bin (Frequency) to display. 
    mtxData = sofia_makeMTX(Nmtx, Pnm, dn, krIndex);
  
    subplot(2,2,ctr)
    sofia_visual3D(mtxData, 0);
    view(90,0)    
        
end
  
  disp(['The plot shows responses at a frequency of ',num2str(round(10*f(krIndex))/10),'Hz']);
