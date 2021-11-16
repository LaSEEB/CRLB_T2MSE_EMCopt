%% '.xml' JEMRIS - MSE
% Test MSE in JEMRIS

%% 0 - Get paths
clear all
clc
% % close all

PC = 1;
if PC == 1
    addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\Toolboxes\Matlab'));
    cd('C:\Users\filia\Documents\Simulation')
    jempath
end

%% 1 - Run JEMRIS_seq
test = 6;

dir_test = ['D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\MSE_Jemris\test',num2str(test)];
mkdir(dir_test)
cd(dir_test)

% JEMRIS_seq

%% 1.5 - MSE parameters
testMSE  = test;
testPlot = 'Fals';
maskTest = 'Fals';

if test == 4
    Echos   = 6;
    TE      = 8;    % in ms
    Ny      = 128;
    Nx      = Ny;
    flipA   = 150;
    T2_low  = 8;
    T2_high = 45;
    SNRval  = 30;
    nreps   = 3;
    optimSeqT2   = 8;
    echoTest = 3;
    
elseif test == 5
    Echos   = 6;
    TE      = 12;    % in ms
    Ny      = 128;
    Nx      = Ny;
    flipA   = 165;
    T2_low  = 8;
    T2_high = 45;
    SNRval  = 30;%27.5770;
    nreps   = 2;    
    optimSeqT2   = 45;
    echoTest = 2;

elseif test == 6
    Echos   = 6;
    TE      = 8;    % in ms
    Ny      = 128;
    Nx      = Ny;
    flipA   = 150;
    T2_low  = 8;
    T2_high = 45;
    SNRval  = 30;
    nreps   = 3;
    optimSeqT2   = 8;
    echoTest = 3;    
end
%% 2 - Run JEMRIS_sim - Run in Pangeia
% % JEMRIS_sim

%  change directory: cd C:\Users\filia\Documents\PhD\Projetos\qMRI\Reconstruction\MSE\test2
%  jemris simu.xml   ou    mpirun -np 15 pjemris simu.xml
%                        - para correr em paralelo. Neste caso utilizando processadores
%% 3 - Load Data
name_simu = ['test_',num2str(test),'mse.mat'];
read_h5signals(name_simu,'signals.h5')
load(name_simu)
kdata_aux = complex(M(:,1),M(:,2));   % data obtained with simulation

%% 4 - Re-organize data
kdata = zeros(Nx,Ny,Echos);
aux = reshape(kdata_aux,Echos*Nx,Ny);
for i=1:Ny
    for ech=1:Echos
        kdata(i,:,ech) = aux((ech-1)*Nx+1:ech*Nx,i);
    end
end
clear aux
% % kdata = zeros(Nx,Ny,Echos);
% % for ech=1:Echos
% %     for i=2:Ny
% %         itt = Echos*Ny*(i-1)+1+Ny*(ech-1);
% %         kdata(:,i,ech) = kdata_aux(itt:itt+Ny-1); %Dim [Nx*Ny,Echos]=(32*32,4)
% %     end
% % end


%% 5 - Recon with FFT - with added noise
% ... 5.1 - Recon with FFT...
for ech=1:Echos
    im(:,:,ech) = ifftshift(ifft2(fftshift(kdata(:,:,ech))));
end

% ... 5.2 - plot das figuras ...
if testPlot == 'True'
    figure()
    axisBound = [min(min(min(abs(im)))) max(max(max(abs(im))))];
    for ech=1:Echos
        subplot(round(Echos/2),round(Echos/2)+1,ech)
        %     figure;imshow(abs(kdata_recon),[])  % plot do espaco_k
        imshow(abs(im(:,:,ech)),[]);
        hold on
        title(['Echo ',num2str(ech)])
        caxis(axisBound)
    end
end
%% 6 - Add Noise - Get Sigma from Image Space - Apply in K-space
% ... 6.1 - Get Mask ...
if maskTest == 'True'
    % all mask
    figure; imshow(abs(im(:,:,1)),[]); % display image 1st echo random slice with good cartilage display
    roiMSE          = roipoly;                   % select roi with cartilage region only
    phantom_mask_m  = roiMSE;
    phantom_mask    = roiMSE(:)';
    
    % Inner mask
    figure; imshow(abs(im(:,:,1)),[]);       % display image 1st echo random slice with good cartilage display
    roiMSE_inner         = roipoly;               % select roi with cartilage region only
    phantom_mask_m_inner = roiMSE_inner;
    phantom_mask_m_out   = phantom_mask_m - phantom_mask_m_inner;
    phantom_mask_inner   = phantom_mask_m_inner(:)';
    phantom_mask_out     = phantom_mask_m_out(:)';
    
    % save
    cd(file_folder)
    save(['z_phantom_mask.mat'],'phantom_mask','phantom_mask_inner','phantom_mask_out', ...
        'phantom_mask_m','phantom_mask_m_inner','phantom_mask_m_out')
elseif maskTest == 'Fals'
    load('z_phantom_mask.mat')
end
se = strel('disk',2);
aux_phantom_mask_m_inner = imerode(phantom_mask_m_inner,se);
aux_phantom_mask_m_out   = imerode(phantom_mask_m_out,se);

% ... 6.2 - Obtain signma from SNR to T2_high ...
img_echo1_T2h = nonzeros(abs(im(:,:,echoTest)).*aux_phantom_mask_m_out);   % use echo 1
sigma_h       = (1/SNRval)*(mean(abs(img_echo1_T2h)));          % sigma from SNR
%     figure();imshow(img_echo1_T2h,[])

% ... 6.3 - Obtain signma from SNR to T2_low ...
img_echo1_T2l = nonzeros(abs(im(:,:,echoTest)).*aux_phantom_mask_m_inner);  % use echo 1
sigma_l       = (1/SNRval)*(mean(abs(img_echo1_T2l)));           % sigma from SNR
%     figure();imshow(img_echo1_T2l,[])

for rep=1:nreps
    % ... 6.4 - Signal with noise added to k-space ...    
    kdata_h(:,:,:,rep) = kdata + sigma_h.*(randn(size(kdata))+1i*randn(size(kdata)))./sqrt(2)*Ny;
    kdata_l(:,:,:,rep) = kdata + sigma_l.*(randn(size(kdata))+1i*randn(size(kdata)))./sqrt(2)*Ny;
    
    % ... 6.5 - Recon images with noise ...
    for ech=1:Echos
        im_h(:,:,ech,rep) = ifftshift(ifft2(fftshift(kdata_h(:,:,ech))));
    end
    for ech=1:Echos
        im_l(:,:,ech,rep) = ifftshift(ifft2(fftshift(kdata_l(:,:,ech))));
    end
end

% ... 6.6 - Measuring SNReff ...


aux_im_h_out = nonzeros(abs(im_h(:,:,1,1)).*aux_phantom_mask_m_out);     
aux_im_h_in  = nonzeros(abs(im_h(:,:,1,1)).*aux_phantom_mask_m_inner);   
aux_im_l_out = nonzeros(abs(im_l(:,:,1,1)).*aux_phantom_mask_m_out);     
aux_im_l_in  = nonzeros(abs(im_l(:,:,1,1)).*aux_phantom_mask_m_inner);   

SNReff.h_out = mean(aux_im_h_out)/std(aux_im_h_out);
SNReff.h_in  = mean(aux_im_h_in)/std(aux_im_h_in);
SNReff.l_out = mean(aux_im_l_out)/std(aux_im_l_out);
SNReff.l_in  = mean(aux_im_l_in)/std(aux_im_l_in);

if testPlot == 'True'
    figure();
    plot(squeeze(abs(im_h(30,30,:,1))),'ro--')
    hold on
    plot(squeeze(abs(im_l(30,30,:,1))),'g+--')
    hold on
    plot(squeeze(abs(im(30,30,:))),'bx--')
    hold on
    legend('im_h','im_l','im')
    
    % % figure;
    % % subplot(221);imshow(abs(im_l(:,:,1,1)-im_h(:,:,1,1)),[]); % display image 1st echo random slice with good cartilage display
    % % subplot(222); imshow(abs(im(:,:,1)-im_l(:,:,1,1)),[]); % display image 1st echo random slice with good cartilage displa
    % % subplot(223); imshow(abs(im(:,:,1)-im_h(:,:,1,1)),[]); % display image 1st echo random slice with good cartilage display
    
    scaletest = [min(min(abs(im(:,:,echoTest)))) max(max(abs(im(:,:,echoTest))))];
    figure;
    subplot(221);imshow(abs(im(:,:,echoTest)),[]); % display image 1st echo random slice with good cartilage display
    hold on; caxis([scaletest])
    subplot(222); imshow(abs(im_l(:,:,echoTest)),[]); % display image 1st echo random slice with good cartilage displa
    hold on; caxis([scaletest])
    subplot(223); imshow(abs(im(:,:,echoTest)-im_h(:,:,echoTest)),[]); % display image 1st echo random slice with good cartilage display
    hold on; caxis([scaletest])
end
%% 6 - Save Images with and w/ NOISE

save(['Test',num2str(test),'_JEMRISPhantom_Image_',...
    '_Echoes',num2str(Echos),...
    '_NxNy',num2str(Nx), ...
    '_flipAngle',num2str(flipA), ...
    '_TE',num2str(TE)],...
    'im','im_h','im_l','Echos','Nx','flipA',...
    'TE','test','T2_low','T2_high','SNRval','SNReff','nreps','optimSeqT2')
SNReff
