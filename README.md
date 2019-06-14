# vUS-data-processing
Function description for vUS data processing based on g1 analysis

A.	vUS data processing (DV model) for in vivo data
A.1. main function
%% IQ to vUS data processing for in vivo experiment
clear all;
addpath('./Functions/');
 
SpatialMsk = questdlg('Use ULM mask for vUS spatial constrain?', ...
    'Select', ...
    'YES', 'NO', 'Cancel', 'Cancel');
%% load data 
disp(['Loading data...']);
load ('./DATA/invivoData.mat'); 
% IQ: beamformed complex quadratue data
% BB: microbubble accumuation map (resize to 25 um pixel size),
    % BB(:,:,1): up flow
    % BB(:,:,2): down flow
    % BB(:,:,3): all flow
% BBV: microbubble velocity map (resize to 25 um pixel size),
    % BBV(:,:,1): up flow
    % BBV(:,:,2): down flow
% BBVz: microbubble axial velocity map (resize to 25 um pixel size),
    % BBVz(:,:,1): up flow
    % BBVz(:,:,2): down flow
[nz,nx,nt]=size(IQ);
%% DAQ infomation and DATA processing parameter
DAQinfo.C=1540;                    % sound speed, m/s
DAQinfo.FWHM=[125 100]*1e-6;        % (X, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
DAQinfo.rFrame=5000;               % sIQ frame rate, Hz
DAQinfo.f0=16.625E6;               % Transducer center frequency, Hz
PRSSinfo.g1nT=nt;                  % g1 calculation sample number
PRSSinfo.g1nTau=100;               % maximum number of time lag
PRSSinfo.SVDrank=[25 nt];          % SVD rank [low high]
PRSSinfo.HPfC=25;                  % High pass filtering cutoff frequency, Hz
PRSSinfo.NEQ=1;                    % 0: no noise equalization; 1: apply noise equalization
%% Clutter rejection
disp('SVD Processing ...');
[sIQ, sIQHP, sIQHHP]=IQ2sIQ(IQ,DAQinfo,PRSSinfo); 
% sIQ: Singular Value Decomposition (SVD) spatiaotemporal filtered signal
% sIQHP: SVD+High Pass filtering
% sIQHHP: SVD+High Pass filtering with 70 HZ cutoff frequency
% eqNoise: obtained noise map
%% vUS data processing
switch SpatialMsk
    case 'YES'
        %% vUS data processing using the ULM mask (BB)
        clear IQ sIQ sIQHHP
        PRSSinfo.rfnScale=2; % spatial resize scale
        PRSSinfo.MskType=1; % use ULM spatial mask
        disp('vUS Processing ...(NOTE: it taks around 120 mins)');
        tic
        [F, V, Vz, SigmaVz, R]=sIQ2vUS_NP_DV(sIQHP, BB, DAQinfo,PRSSinfo);   
        toc
        %% Results plot
        [VzCmap]=Colormaps_fUS;
        xCoor=[1:nx]*0.05/PRSSinfo.rfnScale;
        zCoor=[1:nz]*0.05/PRSSinfo.rfnScale;
 
        Fig=figure;
        set(Fig,'Position',[400 400 1300 450])
        subplot(1,2,1)
        h1=imagesc(xCoor,zCoor,V(:,:,1)); % up flow
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        axis equal tight;
        hold on;
        h2=imagesc(xCoor,zCoor,V(:,:,2)); % down flow
        AlphaMsk=(abs(V(:,:,2))/3).^4;
        AlphaMsk(AlphaMsk>1)=1;
        set(h2,'AlphaData',AlphaMsk)
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        hold off
        axis equal tight;
        xlabel('x [mm]')
        ylabel('z [mm]')
        title('vUS-V [mm/s]')
        
        subplot(1,2,2)
        h1=imagesc(xCoor,zCoor,Vz(:,:,1)); % up flow
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        axis equal tight;
        hold on;
        h2=imagesc(xCoor,zCoor,Vz(:,:,2)); % down flow
        AlphaMsk=(abs(V(:,:,2))/3).^4;
        AlphaMsk(AlphaMsk>1)=1;
        set(h2,'AlphaData',AlphaMsk)
        colormap(VzCmap);
        caxis([-25 25]);
        colorbar
        hold off
        axis equal tight;
        xlabel('x [mm]')
        ylabel('z [mm]')
        title('vUS-Vz [mm/s]')
        %% save data
        save('./vUSBB.mat','-v7.3','F','V','Vz','SigmaVz','R','BB','BBV','BBVz','P');
    case 'NO'
        %% PDI processing
        disp('PDI Processing ...');
        [PDIHHP]=sIQ2PDI(sIQHHP);
        clear IQ sIQ sIQHHP
        PRSSinfo.rfnScale=1; % spatial resize scale
        PRSSinfo.MskType=0; % not use ULM spatial mask
        disp('vUS Processing ... (NOTE: it taks around 60 mins)');
        tic
        [F, V, Vz, SigmaVz, R]=sIQ2vUS_NP_DV(sIQHP,(PDIHHP).^0.5, DAQinfo, PRSSinfo); 
        toc
        %% Results plot
        [VzCmap]=Colormaps_fUS;
        Coor.x=[1:nx]*0.05/PRSSinfo.rfnScale;
        Coor.z=[1:nz]*0.05/PRSSinfo.rfnScale;
        
        Fig=figure;
        set(Fig,'Position',[400 400 1300 450])
        subplot(1,2,1)
        h1=imagesc(Coor.x,Coor.z,V(:,:,1)); % up flow
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        axis equal tight;
        hold on;
        h2=imagesc(Coor.x,Coor.z,V(:,:,2)); % down flow
        AlphaMsk=(abs(V(:,:,2))/5).^4;
        AlphaMsk(AlphaMsk>1)=1;
        set(h2,'AlphaData',AlphaMsk)
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        hold off
        axis equal tight;
        xlabel('x [mm]')
        ylabel('z [mm]')
        title('vUS-V [mm/s]')
        
        subplot(1,2,2)
        h1=imagesc(Coor.x,Coor.z,Vz(:,:,1)); % up flow
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        axis equal tight;
        hold on;
        h2=imagesc(Coor.x,Coor.z,Vz(:,:,2)); % down flow
        AlphaMsk=(abs(Vz(:,:,2))/5).^4;
        AlphaMsk(AlphaMsk>1)=1;
        set(h2,'AlphaData',AlphaMsk)
        colormap(VzCmap);
        caxis([-25 25]);
        colorbar
        hold off
        axis equal tight;
        xlabel('x [mm]')
        ylabel('z [mm]')
        title('vUS-Vz [mm/s]')
        %% plot PDI-based HSV velocity map  
        figure;
        PLOTwtV(V,(PDIHHP).^0.5,Coor,[-30 30])
        title('vUS-V [mm/s]');
        %% save data
        save('./vUS.mat','-v7.3','F','V','Vz','SigmaVz','R','PDIHHP','P');
    case 'Cancel'
end


A.2. function IQ2sIQ
% IQ to sIQ with SVD data processing, sIQ to sIQHP with high pass filtering on sIQ.
function [sIQ, sIQHP, sIQHHP, eqNoise]=IQ2sIQ(IQ,DAQinfo,PRSSinfo)
% Input:
    % IQ: complex IQ data, obtained with beamforming
    % DAQinfo.C: sound speed, m/s
    % DAQinfo.FWHM: (X, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
    % DAQinfo.rFrame: sIQ frame rate, Hz
    % DAQinfo.f0: Transducer center frequency, Hz
    % PRSSinfo.g1nt: g1 calculation sample number
    % PRSSinfo.g1nTau: maximum number of time lag
    % PRSSinfo.SVDrank: SVD rank [low high]
    % PRSSinfo.HPfC:  High pass filtering cutoff frequency, Hz
    % PRSSinfo.NEQ: do noise equalization? 0: no noise equalization; 1: apply noise equalization
% Output:
    % sIQ: Singular Value Decomposition (SVD) spatiaotemporal filtered signal
    % sIQHP: SVD+High Pass filtering
    % sIQHHP: SVD+High Pass filtering with 70 HZ cutoff frequency
    % eqNoise: obtained noise map
A.3. function sIQ2vUS_NP_DV
%% US g1 fit for in vivo data, fit negative and postive frequency signal separately
function [F, V, Vz, SigmaVz, R]=sIQ2vUS_NP_DV(sIQ, FitMsk0, DAQinfo, PRSSinfo)
% input: 
    % sIQ: bulk motion removed data
    % FitMsk0: ULM or PDI-based spatial constrain mask
    % DAQinfo: data acquistion information, including
        % DAQinfo.FWHM: (X, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
        % DAQinfo.rFrame: sIQ frame rate, Hz
        % DAQinfo.f0: Transducer center frequency, Hz
        % DAQinfo.C: Sound speed in the sample, m/s
    % PRSSinfo: data processing parameters, including 
        % PRSSinfo.g1nT: g1 calculation sample number
        % PRSSinfo.g1nTau: maximum number of time lag
        % PRSSinfo.SVDrank: SVD rank [low high]
        % PRSSinfo.HPfC:  High pass filtering cutoff frequency, Hz
        % PRSSinfo.NEQ: do noise equalization? 0: no noise equalization; 1: apply noise equalization
        % PRSSinfo.rfnScale: spatial refind scale
        % PRSSinfo.MskType: 1: use vULM data as mask; 0: PDI
 % output:
    % F: dynamic factor
    % V: total velocity, mm/s
    % Vz: axial velocity, mm/s
    % SigmaVz: axial velocity distribution, mm/s
    % R: fitting accuracy
 
B.	vUS data processing (SV model) for phantom data
B.1. Main function
%% IQ to vUS data processing for phantom data using the basic model
clear all;
addpath('./Functions/');
%% load data 
disp(['Loading data...']);
% load ('./DATA/exvivoData5a.mat');  % angled flow, preset speed 5 mm/s
load ('./DATA/exvivoData15a.mat');  % angled flow, preset speed 15 mm/s
% load ('./DATA/exvivoData9t.mat');  % transverse flow, preset speed 9 mm/s
% load ('./DATA/exvivoData25t.mat');  % transverse flow, preset speed 25 mm/s
% IQ: beamformed complex quadratue data
[nz,nx,nt]=size(IQ);
%% DAQ infomation and DATA processing parameter
DAQinfo.C=1540;                    % sound speed, m/s
DAQinfo.FWHM=[125 90]*1e-6;        % (X,Z) spatial resolution, Full Width at Half Maximum of point spread function, m
DAQinfo.rFrame=5000;               % sIQ frame rate, Hz
DAQinfo.f0=16.625E6;               % Transducer center frequency, Hz
PRSSinfo.g1nT=nt;                  % g1 calculation sample number
PRSSinfo.g1nTau=100;               % maximum number of time lag
PRSSinfo.SVDrank=[3 nt];           % SVD rank [low high]
PRSSinfo.HPfC=10;                  % High pass filtering cutoff frequency, Hz
PRSSinfo.NEQ=1;                    % 0: no noise equalization; 1: apply noise equalization
%% Clutter rejection
disp('SVD Processing ... (NOTE: it taks around 6 mins!)');
[sIQ]=IQ2sIQ(IQ,DAQinfo,PRSSinfo); 
% sIQ: SVD filtered data
%% vUS data processing
clear IQ
PRSSinfo.rfnScale=1; % spatial resize scale
disp('vUS Processing ...');
[F, V, Vz, SigmaVz, Vcz, R]=sIQ2vUS_SV(sIQ, DAQinfo,PRSSinfo);
%% Results plot
[VzCmap,VzCmapDn]=Colormaps_fUS;
Coor.x=[1:nx]*0.05/PRSSinfo.rfnScale;
Coor.z=[1:nz]*0.05/PRSSinfo.rfnScale;

Fig=figure;
set(Fig,'Position',[400 400 1700 350])
subplot(1,3,1)
h1=imagesc(Coor.x,Coor.z,abs(V)); 
colormap(VzCmapDn);
caxis([0 30]);
colorbar
axis equal tight;
xlabel('x [mm]')
ylabel('z [mm]')
title('vUS-V [mm/s]')
subplot(1,3,2)
h2=imagesc(Coor.x,Coor.z,abs(Vz)); 
colormap(VzCmapDn);
caxis([0 30]);
colorbar
axis equal tight;
xlabel('x [mm]')
ylabel('z [mm]')
title('vUS-Vz [mm/s]')
subplot(1,3,3)
h3=imagesc(Coor.x,Coor.z,abs(Vcz)); 
colormap(VzCmapDn);
caxis([0 30]);
colorbar
axis equal tight;
xlabel('x [mm]')
ylabel('z [mm]')
title('Color Doppler-Vz [mm/s]')

 B.2. function sIQ2vUS_SV
%% US g1 fit, fit all frequency signal, for single flow direction scenario
function [F, V, Vz, SigmaVz,Vcz,R]=sIQ2vUS_SV(sIQ, DAQinfo,PRSSinfo)
% input: 
    % sIQ: bulk motion removed data
    % DAQinfo: data acquistion information, including
        % DAQinfo.FWHM: (X, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
        % DAQinfo.rFrame: sIQ frame rate, Hz
        % DAQinfo.f0: Transducer center frequency, Hz
        % DAQinfo.C: Sound speed in the sample, m/s
    % PRSSinfo: data processing parameters, including 
        % PRSSinfo.g1nT: g1 calculation sample number
        % PRSSinfo.g1nTau: maximum number of time lag
        % PRSSinfo.SVDrank: SVD rank [low high]
        % PRSSinfo.HPfC:  High pass filtering cutoff frequency, Hz
        % PRSSinfo.NEQ: do noise equalization? 0: no noise equalization; 1: apply noise equalization
        % PRSSinfo.rfnScale: spatial refind scale
        % PRSSinfo.MskType: 1: use vULM data as mask; 0: PDI
 % output:
    % F: dynamic factor
    % V: total velocity, mm/s
    % Vz: axial velocity, mm/s
    % SigmaVz: axial velocity distribution, mm/s
    % Vcz: axial velocity obtained with Color Doppler, mm/s
    % R: fitting accuracy



