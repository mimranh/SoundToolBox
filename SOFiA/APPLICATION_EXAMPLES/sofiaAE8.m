% SOFiA example 8: Superdirective impulse responses   
% SOFiA Version  : R11-1220
% Array Dataset  : R11-1018

  clear all
  clc
  
% Read VariSphear dataset
% !!! LOCATE THE FOLDER: "EXAMPLE5_RealRoom" 

  downsample = 2;
  timeData = sofia_readVSAdata(downsample);

% Transform time domain data to frequency domain and generate kr-vector
  
  [fftData, kr, f] = sofia_fdt(timeData);

% Spatial Fourier Transform
  Nsft = 4;
  Pnm  = sofia_stc(Nsft, fftData, timeData.quadratureGrid);

% Radial Filters for a rigid sphere array
  Nrf      = Nsft;   % radial filter order               
  maxAmp   = 10;     % Maximum modal amplification in [dB]
  ac       = 2;      % Array configuration: 2 = Rigid Sphere 
  
  dn                         = sofia_mf(Nrf , kr, ac, maxAmp);
  [dn, kernelSize, latency]  = sofia_rfi(dn); % Radial Filter Improvement

% Plane wave decomposition for directions given by OmegaL:
  Npdc     = Nsft;   %Plane wave decomposition order
  OmegaL   = [0 pi/2; pi pi/2];
  Y        = sofia_pdc(Npdc, OmegaL, Pnm, dn);

% Reconstruct directional impulse responses
  impulseResponses = sofia_makeIR(Y, 1/8, downsample);
  
% Remove filter latency (Be careful, this is done very roughly here for 
% demonstrtion/visualisation purpose only)

impulseResponses = impulseResponses(:, (latency*downsample):end);
  
  figure(1)
  
  tscale = linspace(0,size(impulseResponses,2)/(timeData.FS*downsample), size(impulseResponses,2));
  
  plot(tscale, impulseResponses')
  title('Directional Impulse Responses')
  xlabel('Time in s')
  ylabel('Signal Amplitude')
