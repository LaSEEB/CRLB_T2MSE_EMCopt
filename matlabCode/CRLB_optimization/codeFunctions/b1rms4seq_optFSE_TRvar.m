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
%       b1Plus_rms: B1+rms value for the sequence designed
%       T_scan: Time of scanner in (ms)
%
% by: TTFernandes, 2021

function [b1Plus_rms,T_scan_s, TR, Trec] = b1rms4seq_optFSE_TRvar(Sample_weight,TE,TR,ETL,vector_flipAngle,st,nsli,res, display)

%% 0 - Paths 

%% 1 - Constants
SiemensB1fact = 1.32;     % need to explore this further - B1+ factor
GEB1fact      = 1.1725;   % need to explore this further - B1+ factor

Wbody_weight       = 103.45;   %kg - from Visible Human Male
Head_weight        = 6.024;    %kg - from Visible Human Male
Sample_head_weight = (Head_weight/Wbody_weight)*Sample_weight;

%% 2 - SAR limits

SystemMAX       = 5;   % units in (W) - from Visible Human Male


%% 3 - Import a seq file to compute SAR for
% ... 3.1 - Parameters for pulseq ...
maxGrad        = 30;
maxSlew        = 124;
alpha_ex       = 90;                % flip angle excitation in angles
sliceThickness = st;                % slice thickness
rf_ex_phase    = pi/2;              % RF excitation Phase
rf_ref_phase   = 0;                 % RF refocusing Phase
flip_ex        = alpha_ex*pi/180;   % Flip angle excitation in rad
t_ex           = 2.5e-3;            % in ms
t_ref          = 1.7e-3;              % in ms TODO ver valores da philips
tend           = 0.0079*1e3;        % Spoilers gradients times in (s)
gamma          = 42.54e6;           % Gyromagnetic constant (Hz/T)


% ... 4.2 - initialize sequence ...
seq = mr.Sequence();
sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', 'MaxSlew', maxSlew, ...
              'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
              'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

%% 4 - Identify RF blocks and compute SAR - 10 seconds must be less than twice and 6 minutes must be less than 4 (WB) and 3.2 (head-20)
T_vector   = [];
SARwbg_vec = zeros(1,size(vector_flipAngle,2));
SARhg_vec  = zeros(1,size(vector_flipAngle,2));

% ... 4.1 - Impact of  rf_excitation ...
% Create rf_excitation
[rf_exc, gz]   = mr.makeSincPulse(flip_ex,'Duration',t_ex,...
    'SliceThickness',sliceThickness,...
    'apodization',0.5,'timeBwProduct',4,...
    'phaseOffset',rf_ex_phase,'system',sys);
rf_excDur      = mr.calcDuration(rf_exc);
T_vector       = rf_excDur/2*1e3;
signal         = rf_exc.signal;

    % Calculate maxSignal amplitude - Global
B1plus_rf_ex   = max(signal)/gamma * 1e6;
b1plus_t_ex    = flip_ex /  (max(signal) * 2*pi) * 1e3;    % time for specific area (s): trf = flip_angle / b1_max


% ... 4.2 - Impact of rf_refocusing ...
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
    
    % Calculate maxSignal amplitude - Global
    B1plus_rf_ref     = max(signal)/gamma * 1e6;
    b1plus_t_refoc    = (vector_flipAngle(jj)*pi/180) / (max(signal) * 2*pi) * 1e3;    % time for specific area (s): trf = flip_angle / b1_max

end

T_vector = [T_vector T_vector(end)+tend];


%% 5 - Obtain Time vector
RFexcDur   = rf_excDur *1e3;    % (ms)
RFrefocDur = rf_refocDur *1e3;  % (ms)

if (TE/2- RFexcDur/2 - RFrefocDur/2)<0 || (TE - RFrefocDur)<0
    if display == 'True'
        fprintf('TE is too short\n\n')
    end
    b1Plus_rms  = NaN;
    T_scan_s      = NaN;
    return
end

% crushers
t_gs4     = 0.00198/2*1e3; % units (ms)
t_gs5     = 0.0018*1e3;    % units (ms)
t_spoiler = 5;             % units (ms)

% Get time
aux_Trec  = RFexcDur + (TE/2- RFexcDur/2) + TE * ETL;  % (ms)
aux_TRmin = aux_Trec + t_gs4 + t_gs5 + t_spoiler;      % (ms)
TRmin     = aux_TRmin * nsli;                          % (ms) 

% Respect T2 map image - recovery of longitudinal magnt big enought.
if TRmin<TR 
    T_scan  = TR*res;                               % (ms)
    Trec    = TR - aux_Trec;                        % (ms)
    TR      = TR*1e-3;                              % (s)
else
    T_scan  = TRmin * res ;                         % (ms)
    Trec    = TRmin - aux_Trec;                     % (ms)
    TR      = TRmin*1e-3;                           % (s)
end

T_scan_s   = T_scan*1e-3;                              % (s)
T_scan_m   = T_scan_s/60;                              % (min)

%% 6 - Time averaged RF power - match Siemens data
% ... 6.1 - Average RFpower
sumb1plus_rms = ( (B1plus_rf_ex )^2 * b1plus_t_ex +  ETL * (B1plus_rf_ref)^2 *b1plus_t_refoc  )  * nsli * res;
b1Plus_rms    = sqrt(sumb1plus_rms/T_scan);    % units (uT)
if display == 'True'
    disp(['Value for B1+rms is: ', num2str(b1Plus_rms),' uT']);
end

%% 7 - Check for each instant of the time averaged SAR with appropriate time limits
% %  if b1Plus_rms > 10
% %          error('B1+rms is higher than threshold - Increase TR');
% %  end





