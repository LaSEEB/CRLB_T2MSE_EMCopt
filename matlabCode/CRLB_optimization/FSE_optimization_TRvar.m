%% FSE optimization with CRLB analysis fitting EPG curves
% Dictionaries previously generated
% TTFernandes - Sept2021


% =============================
% Toolboxes:
%      - sar4seq - in: https://github.com/imr-framework/sar4seq/tree/master/sar4seq
% =============================
%% 0 - Set matlab paths and toolboxes
clear, clc

% Initial Sets
testFSEopt  = 4;

%toolboxes
% % addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\Toolboxes\sar4seq-master'))
addpath(genpath('\Toolboxes\pulseq-master'));  % add pulseq to toolboxes of matlab

%paths
file_path    = 'Define path where you added the CRLB_T2MSE_EMCopt';
file_path_dic       = [file_path '\Dictionaries'];
file_path_data      = [file_path '\Data'];
file_path_dic_save  = [file_path '\Dictionaries\FAVar_test' num2str(testFSEopt)];


addpath(genpath(file_path))


%% 1 - Tests

plotTest    = 'Fals';   % 'True' or 'Fals'
saveResults = 'True';

testSAR     = 'b1rms';  % For B1_rms test = 'b1rms' OR Loading data 'loadD' 
testDict    = 'Fals';   % Test Dictionary - 'True' OR Load Dictionary - 'Fals'
testCRLB    = 'Fals';   % Test CRLB - 'True' OR Load CRLB - 'Fals'
testMC      = 'Fals';   % Test Monte-Carlo for sanity check of CRLB model
GaussFit    = 'Fals';   % Gauss fit for the data of MC


%% 2 - Parameters to Optimize - TODO

% ... 1.1 - Parameters for model optimization ...
T2        =  40;   % T2 values for each we want to optimize FSE sequence
B1        =  1;    % B1 values for each we want to optimize FSE sequence
T1maxKnee = 1000;  % T1 max for knee (ms)

T2_range = 6:150;
B1_range = 0.6:0.02:1.2;

% ... 1.2 - Parameters of sequence ...
TR            = 5; % units (m) in Keerthivasan et al. 2019
% TR            = 7454;           % Repetition time in (ms)
nslices       = 30;             % number of slices
res           = 128;            % Resolution Nx = Ny
sliceThickn   = 2.6e-3;         % Slice Thickness (m)

% % dTE           = 8;              % Echo time (ms)
% % ETL           = 30;             % Number of Echoes
% % flipAngle     = 180;            % FLip Angle for Refocusing
% % SNR           = 40;             % SNR value

Sample_weight = 75;             % Weight of subject (kg)

% ... 1.3 - Vectors of parameters ...
vector_dTE        = 8:2:30;      % Vector of Echo SpacingTime (ms) - parecido HSa - TODO acertar ao TE possível
vector_ETL        = 6:2:30;      % Vector of # Echos - ver protocolo do HSa
vector_flipAngle  = 90:5:180;    % Vector of flip Angles (ms) - ver artigo tammir
vector_SNR        = 30;          % variance - % variar SNR de teste  [1 30 100] OR [1 5:15:150];
vector_T2         = [8 45];      % vector of different T2 values to test - [8:8:50]


fprintf('\n\n 2 - Sucessfully finished -  Parameters to Optimize \n\n')


%% 3 - Constrains - (3min)
% Control the possible combinations of parameters
tic

% ... 3.1 - Parameters generation matrix
fline          = 1;
aline          = 1;
a              = 1;
failValues     = [];
acceptedValues = [];

% ... 3.2 - Constrains
    % 3.2.1 - SAR
SAR_max     = 100; % W/kg
RFpower_max = 50;  % W
maxB1_rms   = 5;   % uT - TODO tentar provocar os constrains

    % 3.2.1 - Acquisition Time
T_acq     = 9;            % Time in min - 10/12min
SENSE     = 2;             % Acceleration factor
maxTime   = T_acq*SENSE;   % units(min)
maxTime_s = maxTime * 60;  % Time in (s)

% ... 3.3 - Test for constrains ...
if testSAR == 'b1rms'
    % -> 3.3.1 b) - test for all vector ...
    clear b1_rms T_scan
    for jj=1:length(vector_dTE)
        for ii=1:length(vector_ETL)
            for gg=1:length(vector_flipAngle)

                % - Values for RFpower, SAR and T_scan ...
                [b1_rms{jj,ii,gg},T_scan, TR_scan, Trec] = b1rms4seq_optFSE_TRvar(Sample_weight,vector_dTE(jj), TR, ...
                    vector_ETL(ii), vector_flipAngle(gg), sliceThickn, nslices, res,'Fals');

                % - check B1+rms values and max acq Time ...
                if isnan(b1_rms{jj,ii,gg})
                    Models_failed(fline,:) = [vector_dTE(jj) vector_ETL(ii) vector_flipAngle(gg) res T_scan/60 NaN Trec TR_scan];
                    fline = fline+1;
                    continue
                elseif b1_rms{jj,ii,gg} > maxB1_rms || T_scan > maxTime_s
                    Models_failed(fline,:) = [vector_dTE(jj) vector_ETL(ii) vector_flipAngle(gg) res T_scan/60 b1_rms{jj,ii,gg} Trec TR_scan];
                    fline = fline+1;
                else %Models_accepted: [TE(ms) #ETL angle TimeScan(min) b1+rms(uT) Trec(ms)]
                    Models_accepted(aline,:) = [vector_dTE(jj) vector_ETL(ii) vector_flipAngle(gg) res T_scan/60 b1_rms{jj,ii,gg} Trec TR_scan]; %[TE(ms) #ETL angle TimeScan(min)]
                    aline = aline+1;
                end
                a=a+1;
            end
                        
        end
        fprintf(['      Constrain check complete ',num2str(jj),' / ',num2str(length(vector_dTE)),'\n'])
    end
    
    % ... 3.4 - Save data ...
    if saveResults == 'True'
        cd(file_path_data)
        save(['Test',num2str(testFSEopt),'_ConstrainModels_maxB1rms',num2str(maxB1_rms),'_maxTime',num2str(maxTime_s),...
            '_mindTE',num2str(min(vector_dTE)),'_maxdTE',num2str(max(vector_dTE)), ...
            '_minETL',num2str(min(vector_ETL)),'_maxETL',num2str(max(vector_ETL)),...
            '_minFlipA',num2str(min(vector_flipAngle)),'_maxFlipA',num2str(max(vector_flipAngle)),...
            '.mat'],'Models_accepted','Models_failed','maxB1_rms','maxTime_s')
        cd(file_path)
    end
end


if testSAR == 'loadD'
    
    cd(file_path_data)
    load(['Test',num2str(testFSEopt),'_ConstrainModels_maxB1rms',num2str(maxB1_rms),'_maxTime',num2str(maxTime_s),...
          '_mindTE',num2str(min(vector_dTE)),'_maxdTE',num2str(max(vector_dTE)), ...
          '_minETL',num2str(min(vector_ETL)),'_maxETL',num2str(max(vector_ETL)),...
          '_minFlipA',num2str(min(vector_flipAngle)),'_maxFlipA',num2str(max(vector_flipAngle)),...
          '.mat'])
    cd(file_path)    
    
end   
toc

%Models_accepted: [TE(ms) #ETL angle TimeScan(min) b1+rms(uT) Trec(ms)]
fprintf('\n\n 3 - Sucessfully finished -  Constrains \n\n')

      
%% 4 - Get Dictionary with all entraces (9h14min)

tic 
% ... 4.1 - Parameters for generating signals with slr_profile.py and my_epg.py ...
T1_dic        = T1maxKnee;                  % ms
% % T2_dic        = 6:1:105;                    % ms ( 40:1:55; OR 1:1:300; )
% % B1_dic        = 0.7:0.1:1.2;                % B1 value ( 0.7:0.01:0.75; OR 0.6:0.01:1.4; )

% test 9/11
T2_dic        = vector_T2;        % ms ( 40:1:55; OR 1:1:300; )
B1_dic        = 1;                % B1 value ( 0.7:0.01:0.75; OR 0.6:0.01:1.4; )


% test 10/11
T2_dic        = 40;        % ms ( 40:1:55; OR 1:1:300; )
B1_dic        = 1;                % B1 value ( 0.7:0.01:0.75; OR 0.6:0.01:1.4; )

phase_exc     = pi/2;                       % Phase Exciatation - in (rad) - along y
FA_exc_dic    = pi/2;                       % Flip-angle Exciatation - in (rad) - along y

if testDict == 'True'    
  
    % ... test with cicles ...
    for ii=1:size(Models_accepted,1) %Models_accepted: [  TE(ms)   #ETL   angle   TimeScan(min)  ]
        numflipAngles = Models_accepted(ii,2);                                 % Number of Flipe Angles
        refoc_phase   = exp(zeros(1,Models_accepted(ii,2))./(180/pi).*1i);     % phase of refoc pulses (0 = about x axis)
        phase_refoc   = refoc_phase;                                           % Phase Refocosing -  ou zeros(1,nTE);
        FA_refoc_dic  = ones(1,Models_accepted(ii,2))*Models_accepted(ii,3);   % Flip-angle Refocosuing - in (rad) - along y
        Trec          = Models_accepted(ii,7);                                 % Recovery period (ms), from last echo to end

        % Obtain dictionary
        [All_dict_pp{ii},All_pars_pp{ii}] = dict_pars_generator_TRvar(T1_dic,T2_dic,B1_dic,...
                Models_accepted(ii,1),Models_accepted(ii,2),phase_refoc,...
                phase_exc,FA_exc_dic, FA_refoc_dic,file_path_data,...
                Trec); % 10min
        fprintf(['      Successfull models with obtained dictionary ',num2str(ii),' / ',...
            num2str(size(Models_accepted,1)),'\n'])                 
    end
    
    % ... Save Dictionary ...
    mkdir(file_path_dic_save)
    for ii=1:size(Models_accepted,1) %Models_accepted: [  TE(ms)   #ETL   angle(1)   angle(2)   TimeScan(min)  ]
        All_dict = All_dict_pp{ii};
        All_pars = All_pars_pp{ii};
        
        cd(file_path_dic_save)
        save(['Dict_TRvar_ETL',num2str(Models_accepted(ii,2)),...
            '_dTE',num2str(Models_accepted(ii,1)),...
            'ms_flipAngle',num2str(Models_accepted(ii,3)), ...
            'degr_Tseq',num2str(round(Models_accepted(ii,5))), ...
            'min_minB1',num2str(min(B1_dic)),'_maxB1',num2str(max(B1_dic)), ...
            '_minT2',num2str(min(T2_dic)),'_maxT2',num2str(max(T2_dic)),...
            '.mat'],'All_dict','All_pars')
        cd(file_path)
        fprintf(['      Successfull saved models dictionary ',num2str(ii),' / ',...
            num2str(size(Models_accepted,1)),'\n'])  
    end   
end

toc
%Models_accepted: [  TE(ms)  #ETL   angle   TimeScan(min)   b1+rms(uT)  Trec(ms)]

fprintf('\n\n 4 - Sucessfully finished -  Dictionaries \n\n')


%% 5 - CRLB (8min06s)

if testCRLB == 'True'
    
    
    if testDict == 'Fals'
        % ... Load Dictionary ...
        for ii=1:size(Models_accepted,1) %Models_accepted: [  TE(ms)   #ETL   angle(1)   angle(2)   TimeScan(min)  ]
            cd(file_path_dic_save)
            load(['Dict_TRvar_ETL',num2str(Models_accepted(ii,2)),...
                '_dTE',num2str(Models_accepted(ii,1)),...
                'ms_flipAngle',num2str(Models_accepted(ii,3)), ...
                'degr_Tseq',num2str(round(Models_accepted(ii,5))), ...
                'min_minB1',num2str(min(B1_dic)),'_maxB1',num2str(max(B1_dic)), ...
                '_minT2',num2str(min(T2_dic)),'_maxT2',num2str(max(T2_dic)),...
                '.mat'])
            cd(file_path)
            All_dict_pp{ii} = All_dict;
            All_pars_pp{ii} = All_pars;
            fprintf(['      Successfull load models dictionary ',num2str(ii),' / ',...
                num2str(size(Models_accepted,1)),'\n'])
        end
    end
    
    % ... 5.0 - Parameters ...
    T1     = T1maxKnee;
    maxETL = max(Models_accepted(:,2));

    for gg=1:size(vector_T2,2)
        for jj =1:size(vector_SNR,2) % 1.9s - Itteration by SNR values
            tic
            for ii=1:size(Models_accepted,1)
                
                % ... 5.1 - Get Signal ...
                All_dict = All_dict_pp{ii};
                All_pars = All_pars_pp{ii};
                
                S          = All_dict(:,find(All_pars(:,1) == vector_T2(gg) ...
                                               & All_pars(:,3) == B1));    % Signal S - Dictionary                                           
                auxVariab  = exp( - Models_accepted(ii,7)/T1_dic);
                aux2var    = (1 - auxVariab)/ (1 - auxVariab*S(Models_accepted(ii,2)));
                S          = aux2var.* S;                                  % Signal corrected for variable TR 
                                
                % ... 5.2 - Parameters ...
                step          = 1e-4;                                                                   % step for numerical derivative                               
                sigma         = zeros(1,size(vector_SNR,2));
                refoc_phase   = exp(zeros(1,Models_accepted(ii,2))./(180/pi).*1i);                      % phase of refoc pulses (0 = about x axis)
                phase_refoc   = refoc_phase;                                                            % Phase Refocosing -  ou zeros(1,nTE);
                FA_refoc_dic  = ones(1,Models_accepted(ii,2))*Models_accepted(ii,3);                    % Flip-angle Refocosuing - in (rad) - along y
                TRacq         = Models_accepted(ii,8);                                                  % in (s)
                sigma(jj)     = (1/vector_SNR(jj));   % Noise variance for CRLB 

                % ... 5.3 - Get CRLB num and uncertainty  ...
                [CRLB_num{ii,jj,gg},uCRLB(ii,jj,gg)] = CRLB_epg_optFSE_TRvar(S, phase_exc, refoc_phase, ...
                    B1, T1,  vector_T2(gg), Models_accepted(ii,1), Models_accepted(ii,2), ...
                    TRacq, Models_accepted(ii,7), FA_refoc_dic, sigma(jj), ...
                    step, maxETL, plotTest, file_path_data);
                fprintf(['      Test ',num2str(ii),' / ',num2str(size(Models_accepted,1)),'\n'])              
% %                 fprintf(['      uCRLB ',num2str(uCRLB(ii,jj,gg)),'\n'])

            end
            toc
        end
        fprintf(['\n\n      Successfull CRLB model obtained for different T2 values ',num2str(gg),...
            ' / ',num2str(size(vector_T2,1)),'\n'])
    end
    
    % save CRLB data
    cd(file_path_data)
    save(['CRLBvalues_test',num2str(testFSEopt),...
        '_TRvar_ModelsAccepted',num2str(size(Models_accepted,1)),...
        '_minSNR',num2str(min(vector_SNR)),'_maxSNR',num2str(max(vector_SNR)), ...
        '_minT2-',num2str(min(vector_T2)),'_maxT2-',num2str(max(vector_T2)),...
        '.mat'],'CRLB_num','uCRLB')
    cd(file_path)
    
elseif testCRLB == 'Fals'
    % load CRLB data    
    cd(file_path_data)    
    load(['CRLBvalues_test',num2str(testFSEopt),...
        '_TRvar_ModelsAccepted',num2str(size(Models_accepted,1)),...
        '_minSNR',num2str(min(vector_SNR)),'_maxSNR',num2str(max(vector_SNR)), ...
        '_minT2-',num2str(min(vector_T2)),'_maxT2-',num2str(max(vector_T2)),...
        '.mat'])
    cd(file_path)
end


fprintf('\n\n 5 - Sucessfully finished - CRLB uncertainty\n\n')


%% 6 - Get best model fit
clear aux_results_dTE aux_results_ETL aux_results_flipA a b c uni_dTE_x uni_dTE_y ...
    uni_ETL_x uni_ETL_y uni_flipA_x uni_flipA_y results_dTE results_ETL results_flipA ...
    Results
    
% ... 6.1 - Parameters ...
SNR_val   = 30;     % Possible SNR values = 1 30 100
T2_val    = 8;      % units (ms) - Possible T2 values = 8    16    24    32    40    48
testdTE   = 8;     % units (ms)
testETL   = 12;     % echoes
testFlipA = 150;    % units (º)

SNRindx    = find(vector_SNR == SNR_val);
aux_T2indx = find(vector_T2 == T2_val);
T2indx     = (size(Models_accepted,2)+aux_T2indx);
T2thrs     = 10*T2_val;


% ... 6.2 - Get  Results ...
%Models_accepted: [  TE(ms)  #ETL   angle  res  TimeScan(min)   b1+rms(uT)  Trec(ms)]
Results    = Models_accepted; 
% % for gg=1:size(vector_T2,2)
% %     Results  = [Results uCRLB(:,SNRindx,gg)];
% % end

for gg=1:size(vector_T2,2)
    for ii=1:size(Models_accepted,1)
        new_uCRLB(ii,SNRindx,gg) = uCRLB(ii,SNRindx,gg)*sqrt(Models_accepted(ii,5)*60);
    end
end

for gg=1:size(vector_T2,2)
    Results  = [Results new_uCRLB(:,SNRindx,gg)];
end

intervPlot = [min(min(log(Results(:,T2indx:end)))) max(max(log(Results(:,T2indx:end))))]; % for caxis

% ... 6.3 - Get matrix for plot ...
a=1;b=1;c=1;
for ii=1:size(Models_accepted,1)    
    % Test for fixed Echo Time (dTE)
    if Results(ii,1) == testdTE  % Test for fixed dTE
        aux_results_dTE(a,:) = Results(ii,[2 3 T2indx]);
        a = a+1;
    end
    % Test for fixed Echo train Length (ETL)
    if Results(ii,2) == testETL
        aux_results_ETL(b,:) = Results(ii,[1 3 T2indx]);
        b = b+1;
    end
    % Test for fixed flip Angle
    if Results(ii,3) == testFlipA
        aux_results_flipA(c,:) = Results(ii,[1 2 T2indx]);
        c = c+1;
    end    
end

% ... 6.4 - Get matrix for plot ...
uni_dTE_x   = unique(aux_results_dTE(:,1));
uni_dTE_y   = unique(aux_results_dTE(:,2));
uni_ETL_x   = unique(aux_results_ETL(:,1));
uni_ETL_y   = unique(aux_results_ETL(:,2));
uni_flipA_x = unique(aux_results_flipA(:,1));
uni_flipA_y = unique(aux_results_flipA(:,2));

results_dTE   = NaN(length(uni_dTE_x),length(uni_dTE_y));
results_ETL   = NaN(length(uni_ETL_x),length(uni_ETL_y));
results_flipA = NaN(length(uni_flipA_x),length(uni_flipA_y));

for i=1:length(aux_results_dTE)
    for g=1:length(uni_dTE_x)
        for j=1:length(uni_dTE_y)
            if uni_dTE_x(g)==aux_results_dTE(i,1) && uni_dTE_y(j)==aux_results_dTE(i,2)
                results_dTE(g,j) = aux_results_dTE(i,3);
            end
        end
    end
end
for i=1:length(aux_results_ETL)
    for g=1:length(uni_ETL_x)
        for j=1:length(uni_ETL_y)
            if uni_ETL_x(g)==aux_results_ETL(i,1) && uni_ETL_y(j)==aux_results_ETL(i,2)
                results_ETL(g,j) = aux_results_ETL(i,3);
            end
        end
    end
end
for i=1:length(aux_results_flipA)
    for g=1:length(uni_flipA_x)
        for j=1:length(uni_flipA_y)
            if uni_flipA_x(g)==aux_results_flipA(i,1) && uni_flipA_y(j)==aux_results_flipA(i,2)
                results_flipA(g,j) = aux_results_flipA(i,3);
            end
        end
    end
end
clear  aux_results_dTE aux_results_ETL aux_results_flipA


% ... 6.5 - Plots ...
figure()
% plot for fixed dTE
subplot(221)
imagesc(log(results_dTE))
hold on
grid on
title(['log CRLB | T2 = ',num2str(vector_T2(aux_T2indx)),' ms | fixed dTE = ',num2str(testdTE),' ms'],'fontsize',20)
xlabel(['flipAngle (º)'],'fontsize',20)
ylabel(['# Echoes'],'fontsize',20)
set(gca,'FontSize',15)
xticks(1:1:size(uni_dTE_y,1)*2)
yticks(1:1:size(uni_dTE_x,1)*2)
xticklabels(uni_dTE_y)
yticklabels(uni_dTE_x)
caxis(intervPlot)

% plot for fixeed ETL
subplot(222)
imagesc(log(results_ETL))
hold on
grid on
title(['log CRLB | T2 = ',num2str(vector_T2(aux_T2indx)),' ms | fixed ETL = ',num2str(testETL)],'fontsize',20)
xlabel(['flipAngle (º)'],'fontsize',20)
ylabel(['dTE (ms)'],'fontsize',20)
set(gca,'FontSize',15)
xticks(1:1:size(uni_ETL_y,1)*2)
yticks(1:1:size(uni_ETL_x,1)*2)
xticklabels(uni_ETL_y)
yticklabels(uni_ETL_x)
caxis(intervPlot)

% plot for fixeed flip Angle
subplot(223)
imagesc(log(results_flipA))
hold on
grid on
title(['log CRLB | T2 = ',num2str(vector_T2(aux_T2indx)),'ms | fixed fA = ',num2str(testFlipA),'º'],'fontsize',20)
xlabel(['# Echoes'],'fontsize',20)
ylabel(['dTE (ms)'],'fontsize',20)
set(gca,'FontSize',15)
xticks(1:1:size(uni_flipA_y,1)*2)
yticks(1:1:size(uni_flipA_x,1)*2)
xticklabels(uni_flipA_y)
yticklabels(uni_flipA_x)
caxis(intervPlot)


% ... 6.7 - Values of minimized

%SNR = 30 for different T2 values
fprintf(['\n            ... Optimized Results with CRLB variance ... \n\n'])
for gg=1:size(vector_T2,2)
    [minValue(gg), indxResult] = (min( (Results(:,T2indx-aux_T2indx+gg)) ) ) ;
    fprintf(['       CRLB min ' ,num2str(log(minValue(gg))),' variance model for: ', ...
        'T2 = ', num2str(vector_T2(gg)),...        
        ' ms, dTE = ', num2str(Results(indxResult,1)),...
        ' ms, ETL = ', num2str(Results(indxResult,2)),...
        ', flipAngle = ', num2str(Results(indxResult,3)),'º\n'])
end



