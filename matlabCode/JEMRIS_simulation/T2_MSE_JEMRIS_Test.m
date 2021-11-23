% Tests T2 MSE estimation JEMRIS Phantom
% TTFernandes 10_11_2021

%% 0 - Set matlab paths and toolboxes
clear all
clc

PC = 1; % PC=1 - My Computer OR PC=2 - Seia

addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\T2_knee')); % add to directory

%% 0.5 - Tests
plotTest    = 'Fals';
saveResults = 'True';

maskTest  = 'ready';     % 'getIt' - Get mask            OR 'ready' - load mask
dictTest  = 'getIt';     % 'getIt' - Get Dictionary      OR 'ready' - load Dictionary
methodDic = 'JUST_DIC';  % 'DICT_SLR' consider the slice profile & 'JUST_DIC' just build the dictionary
tempMatch = 'getIt';     % 'getIt' - Get ind_parameters  OR 'ready' - load ind_parameters
testSNR   = 'TrueL';     % 'TrueH' OR 'TrueL' OR 'False'

%% 1 - load Data MSE

% ... 1.1 - get file ...
myCD = ('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\MSE_Jemris'); % Define your Directory
dir_dic = [myCD, '\Dictionaries'];
mkdir(dir_dic)

% ... 1.2.0 - Load data ... example: D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\MSE_Jemris\test4\Test4_JEMRISPhantom_Image__Echoes6_NxNy128_flipAngle150_TE8.mat
cd(myCD);
str = {'select file'};
s   = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',str);
if s==1
    [file,file_folder]=uigetfile;
else
    folder=uigetdir;
end
cd(file_folder);
load(file);

% ... 1.2.2 - Defining parameters ...
nechoes      = Echos;                % Number of echos
ESP          = TE;                   % Inter-echo time space in (s)
MSE_TE       = ESP:ESP:ESP*nechoes;  % vector w/ Echoes
nslice       = 1;                    % Number of slices
Ny           = Nx;                   % Ny
sliceTest    = 1;                    % Slice for Testing
slice_inf    = 1;                    % Slice inf
slice_sup    = 1;                    % Slice sup
gamma        = 42.576e6;             % Gyro-magnetic ratio Hz/T
T2optimiz    = optimSeqT2;

if testSNR == 'TrueL'
    imgMSE_final = abs(im_l);    
elseif testSNR == 'TrueH'
    imgMSE_final = abs(im_h);    
elseif testSNR == 'False'
    imgMSE_final = abs(im);
    nreps = 1;
end
    
% ... 1.3 - Parameters ...
slice_inf = min(sliceTest);
slice_sup = max(sliceTest);

nTE       = nechoes; % Number of Echos of TEs
TE_first  = TE;      % First TE - in (s)
ESP       = TE;      % Spacing Time between Echoes - in ()
TE_vector = MSE_TE;  % TE vector - in (s)

nl = Nx; % Image size

FA_exc    = pi/2;   % Flip-angle Exciatation - in (rad) - along y
phase_exc = pi/2;   % Phase Exciatation - in (rad) - along y

FA_refoc    = flipA*pi/180;                             % Flip-angle Refocosing - in (rad) - along x
phase_refoc = exp(zeros(1,nTE)./(180/pi).*1i); % Phase Refocosing -  ou zeros(1,nTE);

% ... 1.4 - time definition ...
t           = TE_vector';   % in ()

fprintf('\n\n 1 - Sucessfully finished - Load\n\n')


%% 2 - inputs
% ... 2.1 - Create and display images montage

montage_imgMSE_disp = permute(imgMSE_final(:,:,slice_inf:slice_sup,1),[1 2 4 3]); %1st echo

% ... 2.2 - segment JEMRIS phantom for 5 selected slices ...
if maskTest == 'getIt'    
    figure; imshow(imgMSE_final(:,:,1,1),[]); % display image 1st echo random slice with good cartilage display
    roiMSE       = roipoly;                   % select roi with cartilage region only
    phantom_mask_m  = roiMSE;
    phantom_mask    = roiMSE(:)';
    
    cd(file_folder)
    save(['z_phantom_mask.mat'],'phantom_mask')

elseif maskTest == 'ready'
    cd(file_folder)
    load('z_phantom_mask.mat')
end
fprintf('\n\n 2 - Sucessfully finished - Get Inputs & Mask\n\n')


%% 3 - Build dictionary - bi-exponencial T2 estimation
tic 

% ... 3.1 - Parameters
% ESP - echo Spacing
T1 = 1000;                  %ms
T2 = [0.1 1:6 6.1:0.1:10 10:1:43 43.1:0.1:47 47:1:80];                %ms
B1 = 0.6:0.01:1.4;
% ESP - Spacing Time between Echoes - in (ms)

% ... 3.2 - Get Dictionary ...
if dictTest == 'getIt'
    if methodDic == 'DICT_SLR'
        [dict_phantom_shortTE,pars] = dict_pars_generator(T1,T2,B1,ESP, ...
            nTE,phase_refoc,phase_exc,FA_exc, FA_refoc*180/pi, methodDic); % 10min
    elseif methodDic == 'JUST_DIC'
        [dict_phantom_shortTE,pars] = dict_pars_generator(T1,T2,B1,ESP, ...
            nTE,phase_refoc,phase_exc,FA_exc, FA_refoc, methodDic); % 10min
    end
    
    col_T2 = pars(:,1);     % indexes of T2
    col_B1 = pars(:,3);     % indexes of B1
    col_B1 = col_B1*100;
    
    % normalize dictionary
    for ii=1:size(dict_phantom_shortTE,2)
        Dict_phantom_shortTE_norm(:,ii) = dict_phantom_shortTE(:,ii)./norm(dict_phantom_shortTE(:,ii));
    end
    figure(); plot(abs(dict_phantom_shortTE))
    % Save Dictionary
    cd(dir_dic)
    save(['Dictionary_MSE_test',num2str(test),'_JEMRIS_',methodDic,'.mat'],'Dict_phantom_shortTE_norm','col_T2','col_B1')
    cd(myCD)
elseif dictTest == 'ready'
    cd(dir_dic)
    load(['Dictionary_MSE_test',num2str(test),'_JEMRIS_',methodDic,'.mat'])    
end

toc

fprintf('\n\n 3 - Sucessfully finished - Build dictionary - bi-exponencial T2 estimation\n\n')



%% 4 - T2 Dictionary Matching (selected 5 slices only)
% tic

% ... 4.1 - Parameters ...
nslices     = 1;
X           = zeros(nechoes,nl*nl,nslices);
T2_dict_map = zeros(nl,nl,nslices);
n_slice     = 1;
cd_results  = [file_folder, '\Results'];

if tempMatch == 'getIt'
    for rep=1:nreps
        % ... 4.2 - Reshape data ...
        X_MSE = squeeze(imgMSE_final(:,:,:,rep));
        X_MSE = permute(X_MSE,[3 1 2]);
        X_MSE = reshape(X_MSE,[nechoes nl*nl]);
        X (:,:,1) = abs(X_MSE);
        
        % ... 4.3 - template_matching of the dictionary ...
        tic
        for i = 1:size(X,2)
            if phantom_mask(1,i) == 0
                [ind_param(1,i,rep)] = 1;
            else
                [ind_param(1,i,rep)] = template_match(Dict_phantom_shortTE_norm,X(:,i));
            end
            fprintf(['      Successfull template match ',num2str(i),' / ',...
                num2str(size(X,2)),'\n'])
        end
        toc
        
        % ... 4.4 - Save Results ...
        mkdir(cd_results); cd(cd_results)
        save(['z_ind_param_Matlab_SNR',testSNR,'_',num2str(SNRval),'.mat'],'ind_param')
        cd(file_folder)
    end
elseif tempMatch == 'ready'
    mkdir(cd_results)
    cd(cd_results)
    load(['z_ind_param_Matlab_SNR',testSNR,'_',num2str(SNRval),'.mat'])
end


% ... 4.5 - estimate of the parameters ...
figure()
for rep=1:nreps
    T2_dict             = col_T2(ind_param(:,:,rep));
    T2_dict_map (:,:,rep) = reshape((T2_dict),[nl,nl]); % [Nx, Ny, nSlices,
    
    
    if plotTest == 'True'
        figure(22); plot(T2_dict);
        figure(23); plot(T2_dict_map(:,:,rep));
    end
    
    % .. get maxValue ..
    maxVal = max(max(max(max((double(abs(T2_dict_map)))))))+2;
    
    % .. Plots ..
    subplot(1,nreps,rep);imshow(T2_dict_map(:,:,rep),[]);colormap hot;
end
% toc
fprintf('\n\n 4 - Sucessfully finished - Build dictionary - bi-exponencial T2 estimation\n\n')

if plotTest == 'True'
    figure();imshow(T2_dict_map(:,:,1),[]);colormap hot;
    caxis([0 55])
    colorbar
    title(['T2 map estimation for T2 = ',num2str(optimSeqT2)])
end
%% 5 - ROI analysis
% ... 5.1 - Mask for inner vial ...
if maskTest == 'getIt'    
    figure; imshow(imgMSE_final(:,:,1,1),[]);       % display image 1st echo random slice with good cartilage display
    roiMSE_inner         = roipoly;               % select roi with cartilage region only
    phantom_mask_m_inner = roiMSE_inner;
    phantom_mask_m_out   = phantom_mask_m - phantom_mask_m_inner;    
    phantom_mask_inner   = phantom_mask_m_inner(:)';
    phantom_mask_out     = phantom_mask_m_out(:)';    
    % save
    cd(file_folder)
    save(['z_phantom_mask.mat'],'phantom_mask','phantom_mask_inner','phantom_mask_out', ...
                               'phantom_mask_m','phantom_mask_m_inner','phantom_mask_m_out')    
else % load
    cd(file_folder)
    load('z_phantom_mask.mat')    
end

% ... 5.2 - Get T2 Distribution for Vials ...
se = strel('disk',2);
aux_phantom_mask_m_inner = imerode(phantom_mask_m_inner,se);
aux_phantom_mask_m_out   = imerode(phantom_mask_m_out,se);

for rep=1:nreps
    T2_inner(:,:,rep) = aux_phantom_mask_m_inner.*T2_dict_map(:,:,rep);
    T2_out(:,:,rep)   = aux_phantom_mask_m_out.*T2_dict_map(:,:,rep);
end

T2_inner = mean(T2_inner,3);
T2_out   = mean(T2_out,3);

if plotTest == 'True'   
    figure();subplot(121);imshow(aux_phantom_mask_m_inner,[]);subplot(122);imshow(aux_phantom_mask_m_out,[]);colormap hot;

    figure();subplot(121);imshow(T2_inner(:,:,1),[]);subplot(122);imshow(T2_out(:,:,1),[]);colormap hot;
end

T2_inner_dist = nonzeros(T2_inner(:)');
T2_out_dist   = nonzeros(T2_out(:)');

T2m_in     = mean(T2_inner_dist);
T2std_in   = std(T2_inner_dist);
T2medi_in  = median(sort(T2_inner_dist));
T2iqr_in   = iqr(sort(T2_inner_dist));

T2m_out    = mean(T2_out_dist);
T2std_out  = std(T2_out_dist);
T2medi_out = median(sort(T2_out_dist));
T2iqr_out  = iqr(sort(T2_out_dist));

% ... 5.3 - Get boxplot ...
boxP = [T2_inner_dist; T2_out_dist];
grp = [T2_low*ones(size(T2_inner_dist,1),1); T2_high*ones(size(T2_out_dist,1),1)];
figure();boxplot(boxP,grp);
xlabel('Different Vials','fontsize',20)
ylabel('T2 (ms)','fontsize',20)
ylim([0 50])
a = gca; % get the current axis;
a.Position(3) = 0.6;
annotation('textbox', [0.4, 0.7, 0.1, 0.1], 'String', "med = " + T2medi_out + ' | iqr = ' + T2iqr_out)
hold on
a.Position(3) = 0.9;
annotation('textbox', [0.4, 0.1, 0.1, 0.1], 'String', "med = " + round(T2medi_in,2) + ' | iqr = ' + T2iqr_in)
title(['Optimized sequence for T2 = ',num2str(T2optimiz),'ms'])

fprintf(['\n\n Value for mu T2_8  = ', num2str(T2m_in),'  &  Std = ',num2str(T2std_in),' \n'])
fprintf([' Value for mu T2_45 = ', num2str(T2m_out),'  &  Std = ',num2str(T2std_out),' \n'])

fprintf('\n\n 5 - Sucessfully finished - ROI Analysis & Plots \n\n')
