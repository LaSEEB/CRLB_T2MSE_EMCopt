function [uMC_fit, uMC_std] = MC_epg_optFSE(nreps, sig, Dict_norm, All_pars, sigma, T2_dic, T2, CRLB_num, SNR,GaussFit, plotTest)
% Obtain the variance expected for a specific model/set of parameters
%
% Functions used:
%
% Inputs:
%   nreps:  	number of reps for MC
%   sig:        Signal
%   Dict_norm:  All dictionary normalized
%   All_pars:   Parameters from dictionary
%   sigma:      sigma measured for SNR value
%   T2_dic:     T2 values
%   T2:         T2 true value
%   CRLB_num:   value for CRLB for comparison
%   SNR:        SNR value
%   GaussFit:   'True' or 'Fals'
%   plotTest:   'True' or 'Fals'
%
% Ouputs:
%   uMC:       MC uncertainty
%

%%
% ... 5.3 - Get noisy measurment & Template Match ...
for rep=1:nreps
    sig_noisy      = sig+complex(randn(size(sig,1),1)*sigma,randn(size(sig,1),1)*sigma);
    sig_noisy      = sig_noisy./norm(sig_noisy);
    
    [tempMatch, ~] = templatematch(Dict_norm,sig_noisy);
    T2_est(rep)    = All_pars(tempMatch,1);
    %             B1_est(rep) = All_pars(tempMatch,3);
    
    sinais(:,rep)=sig_noisy;    % store noisy signals if needed
end

% ... 5.4 - histograms & Normal distribution ...
%     std(B1_est)/mean(B1_est)
[counts,x] = hist(T2_est,T2_dic);             % Histogram of MC
y          = normpdf(T2_dic,T2,CRLB_num(2));  % For CRLB variance

% ... 5.5 - Check standard deviation values of T2 distributions ...
prob          = counts./trapz(x,counts);                                % probability as area under the curve
[mu,sg,curve] = mygaussfit(x,prob,[mean(T2_est) std(T2_est)],'noplot'); % Gaussian fit
uMC_fit       = (sg^2)/(T2^2);                                          % Uncertainty of MC
uMC_std       = (std(T2_est)^2)/(T2^2);                                 % Uncertainty of MC


% ... 5.5 - Figures ...
if plotTest == 'True'
    T2plot = T2_dic(1:T2*2);
    figure, hist(T2_est,T2plot)
    hold on
    title(['T2 estimate histogram - SNR=',num2str(SNR)])
    plot(T2plot,y(T2plot)/max(y(T2plot))*max(counts),'r')
end


% ... 5.7 - See standard deviation for each echo ...
% -----------------------------------------------------------
if GaussFit == 'True'
    figure, bar(x,prob,'hist')
    hold on, plot(x,curve,'r','LineWidth',2), title(['Mu: ' num2str(mu) ' | ' 'Std: ' num2str(sg)])
    xlim([mu-4*sg mu+4*sg])
    
    for i=1:ETL
        i
        sinais = abs(sinais);
        [counts,x] = hist(sinais(i,:));
        prob = counts./trapz(x,counts);
        mu_0 = mean(sinais(i,:));
        sg_0 = std(sinais(i,:));
        [mu,sg,curve] = mygaussfit(x,prob,[mu_0 sg_0],'noplot');
        
        if plotTest == 'True'
            figure(10)
            subplot(ETL/2,2,i)
            bar(x,prob,'hist')
            hold on, plot(x,curve,'r','LineWidth',2), title(['Mu: ' num2str(mu) ' | ' 'Std: ' num2str(sg)])
            xlim([mu-4*sg mu+4*sg]), title(['Mu: ' num2str(mu) ' | ' 'Std: ' num2str(sg)])
        end
    end
    if plotTest == 'True'
        T2plot = T2_dic(1:T2*2);
        [mu1, sigma1] = normfit(T2_est);
        pd = fitdist(T2_est,'normal');
        
        y = pdf(pd,T2plot);
        hold on, plot(T2plot,y,'g--','LineWidth',2)
    end
end




end

