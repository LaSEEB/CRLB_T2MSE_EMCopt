%% File details
%
%     1. Computes RF safety metrics for Pulseq sequences 
%     a. For Pulseq sequences for deployment on Siemens scanners - 
%     computes time averaged RF power for the sequence
%     b. For Pulseq sequences for deployment on GE scanners (via TOPPE) -
%     computes the whole body SAR in W/kg
%     
% 
%    Parameters
%    ----------
%       Sample_weight : weight of the sample being imaged - double
%       ETL : Number of Echoes - double 
%       vector_flipAngle : Vector of flip Angles - double
% 
%     Returns
%     -------
%       RFpower: Time averaged RF power - whole body & head : double
%       SAR: Whole body SAR & head : double
%       T_scan_s : Total time of scan in seconds
%            
% by: TTFernandes, 2021
% adpated from: 'SAR4seq.m' in https://github.com/imr-framework/sar4seq/tree/master/sar4seq
% Copyright of the Board of Trustees of Columbia University in the City of New York

function [RFpower,SAR, T_scan_s] = SAR4seq_optFSE(Sample_weight,TE,TR,ETL,vector_flipAngle,st,nsli,res)

%% 0 - Paths 
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\Toolboxes\pulseq-master'));

%% 1 - Constants
SiemensB1fact = 1.32;     % need to explore this further - B1+ factor
GEB1fact      = 1.1725;   % need to explore this further - B1+ factor

Wbody_weight       = 103.45;   %kg - from Visible Human Male
Head_weight        = 6.024;    %kg - from Visible Human Male
Sample_head_weight = (Head_weight/Wbody_weight)*Sample_weight;

%% 2 - SAR limits
% Limits: 10 seconds must be less than twice and
%         6 minutes must be less than 4 (WB) and 3.2 (head-20)

SixMinThresh_wbg = 4;    % 6 min Glogal SAR limits Whole body (W/Kg)
TenSecThresh_wbg = 8;    % 10 second Global SAR limits Whole body (W/Kg)

SixMinThresh_hg  = 3.2;  % 6 min Glogal SAR limits - head  (W/Kg)
TenSecThresh_hg  = 6.4;  % 10 second Global SAR limits - head (W/Kg)


%% 3 - Q matrices
if(~(exist('Qmat.mat','file')))
    % ... 3.1 - Load relevant model - ONLY once per site if somebody needs to explore the Q matrix formulation ...
    dirname = uigetdir(''); %Load dir with files for EM model
    cdir = pwd;
    cd(dirname)
    load Ex.mat; model.Ex =Ex; clear Ex; 
    load Ey.mat; model.Ey =Ey; clear Ey;  
    load Ez.mat; model.Ez =Ez; clear Ez;  
    load Tissue_types.mat; model.Tissue_types =Tissue_types; clear Tissue_types;  
    load SigmabyRhox.mat;  model.SigmabyRhox =SigmabyRhox; clear SigmabyRhox;  
    load Mass_cell.mat;    model.Mass_cell =Mass_cell; clear Mass_cell; 
    cd(cdir)

    %.... 3.2 - Compute and store Q matrices once- if not done already - per model - will write relevant qmat files ...
    tic;
    Q = Q_mat_gen('Global', model,0);
    save('Qmat.mat','Q');
    toc;
else
    load('Qmat.mat','Q'); %Loads Q
end

%% 4 - Import a seq file to compute SAR for
% ... 4.1 - Parameters for pulseq ...
maxGrad        = 30;
maxSlew        = 124;
alpha_ex       = 90;                % flip angle excitation in angles
sliceThickness = st;                % slice thickness
rf_ex_phase    = pi/2;              % RF excitation Phase
rf_ref_phase   = 0;                 % RF refocusing Phase
flip_ex        = alpha_ex*pi/180;   % Flip angle excitation in rad
t_ex           = 2.5e-3;            % in ms
t_ref          = 2e-3;              % in ms TODO ver valores da philips
tend           = 0.0079*1e3;            % Spoilers gradients times in (s)

% ... 4.2 - initialize sequence ...
seq = mr.Sequence();
sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', 'MaxSlew', maxSlew, ...
              'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
              'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

%% 5 - Identify RF blocks and compute SAR - 10 seconds must be less than twice and 6 minutes must be less than 4 (WB) and 3.2 (head-20)
T_vector   = [];
SARwbg_vec = zeros(1,size(vector_flipAngle,2));
SARhg_vec  = zeros(1,size(vector_flipAngle,2));

% ... 5.1 - Impact of  rf_excitation ...
% Create rf_excitation
[rf_exc, gz]   = mr.makeSincPulse(flip_ex,'Duration',t_ex,...
    'SliceThickness',sliceThickness,...
    'apodization',0.5,'timeBwProduct',4,...
    'phaseOffset',rf_ex_phase,'system',sys);
rf_excDur = mr.calcDuration(rf_exc);
T_vector  = rf_excDur/2*1e3;
signal    = rf_exc.signal;

    % Calculate SAR - Global: Wholebody, Head, Exposed Mass
SARwbg_vec(1) = calc_SAR(Q.Qtmf,signal,Wbody_weight); %This rf could be parallel transmit as well
SARhg_vec(1)  = calc_SAR(Q.Qhmf,signal,Head_weight); %This rf could be parallel transmit as well

% ... 5.2 - Impact of rf_refocusing ...
for jj=1:size(vector_flipAngle,2)
    % Create rf_refocusing
    flip_ref       = vector_flipAngle(jj)*pi/180;  % Flip angle refocusing in rad
    [rf_refoc, gz] = mr.makeSincPulse(flip_ref, 'Duration', t_ref,...
        'SliceThickness',sliceThickness,...
        'apodization',0.5, 'timeBwProduct',4,...
        'phaseOffset',rf_ref_phase, 'use','refocusing',...
        'system',sys);
    rf_refocDur = mr.calcDuration(rf_refoc);

    % Time vector
    if jj == 1
        T_vector = [T_vector T_vector(jj)+(TE/2)];
    else
        T_vector = [T_vector T_vector(jj)+TE];
    end
    
    signal    = rf_refoc.signal;
    
    % Calculate SAR - Global: Wholebody, Head, Exposed Mass
    SARwbg_vec(1+jj) = calc_SAR(Q.Qtmf,signal,Wbody_weight); %This rf could be parallel transmit as well
    SARhg_vec(1+jj)  = calc_SAR(Q.Qhmf,signal,Head_weight); %This rf could be parallel transmit as well

end

T_vector = [T_vector T_vector(end)+tend];

% SAR vectors for all sequence, resolution and number of slices
allSARwbg_vect = repmat(SARwbg_vec,1,res*nsli);
allSARhg_vec   = repmat(SARhg_vec,1,res*nsli);

%% 6 - Obtain Time vector
RFexcDur   = rf_excDur *1e3;    % (ms)
RFrefocDur = rf_refocDur *1e3;  % (ms)

if (TE/2- RFexcDur/2 - RFrefocDur/2)<0 || (TE - RFrefocDur)<0 
    fprintf('TE is too short\n\n')
    RFpower.wbg_tavg = NaN;
    RFpower.wbg_tavg = NaN;
    SAR.SAR_10s_wbg  = NaN;
    T_scan_s         = NaN;
    return
end

% cicle for resolution & for slices
for ii=1:res*nsli-1
    T_vector = [T_vector T_vector(end) + TE + T_vector(1:ETL+1)];
end

T_vector_s = T_vector*1e-3;                                       % (s)
T_slice    = (RFexcDur + (TE/2- RFexcDur/2) + TE*(ETL-1)) * res;  % (ms) TODO
T_slice_s  = T_slice*1e-3;                                        % (s)
T_scan     = T_slice*nsli;                                        % (ms)
T_scan_s   = T_scan*1e-3;                                         % (s)
T_scan_m   = T_scan_s/60;                                         % (min)


%% 7 - Time averaged RF power - match Siemens data
% ... 7.1 - Average RFpower
RFpower.wbg_tavg = sum(allSARwbg_vect)./T_scan_s./SiemensB1fact;
RFpower.hg_tavg  = sum(allSARhg_vec)./T_scan_s./SiemensB1fact;
disp(['Time averaged RF power-Siemens is - Body: ', num2str(RFpower.wbg_tavg),'W &  Head: ', num2str(RFpower.hg_tavg), 'W']);

% ... 7.2 - instantaneous SAR
SARwbg = max(SARwbg_vec); %W/kg but Siemens reports W/lb
SARhg  = max(SARhg_vec);

SAR.wbg_predSiemens = SARwbg.* sqrt(Wbody_weight/Sample_weight)/2;   % SARwbg_pred.Siemens = SARwbg_pred.Siemens/2.20462; %Siemens reports W/lb
SAR.hg_predSiemens = SARhg.* sqrt(Head_weight/Sample_head_weight)/2;
disp(['Predicted SAR-Siemens is - Body:', num2str(SAR.wbg_predSiemens), 'W/kg &  Head: ',num2str(SAR.hg_predSiemens), 'W/kg'])

% ... 7.3 - SAR per time period 10s
timeWindV = find(T_vector_s<10);
sum_SARV  = movsum(allSARwbg_vect,timeWindV(end)-1);

SAR.SAR_10s_wbg   = max(sum_SARV.* sqrt(Wbody_weight/Sample_weight)/2);
% % figure()
% % plot(SAR.SAR_10s)
disp(['Predicted SAR-Siemens acc 10s is - Body:', num2str(SAR.SAR_10s_wbg), 'W/kg '])

% ... 7.4 - SAR per time period 6min
% % timeWindV = find(T_vector_s<60*6);
% % sum_SARV  = movsum(allSARwbg_vect,timeWindV(end)-1);

% % SAR*nsli*ETL<thres


%% 9 - Check for each instant of the time averaged SAR with appropriate time limits
% %  if(sum(SARwbg_pred.GE > TenSecThresh_wbg))
% %          error('Pulse sequence exceeding 10 second Global SAR limits, increase TR');
% %  end





