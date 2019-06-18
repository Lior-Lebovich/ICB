%% Fig. 3: The Poisson network model

gamma_sigma = 1;
gamma_k = 0.1330;
nu_hat = 1.2556; 
baselineTheta = 0.65;
figure;


%% Fig. 3a: stimulus-dependent distribution of firing rates

dL_normpdf = 5;
x_normPdf = -40:0.01:40;
normPdf_dL_0 = normpdf(x_normPdf, log(nu_hat)+0*gamma_k, gamma_sigma);
normPdf_dL_P = normpdf(x_normPdf, log(nu_hat)+dL_normpdf*gamma_k, gamma_sigma);
normPdf_dL_M = normpdf(x_normPdf, log(nu_hat)-dL_normpdf*gamma_k, gamma_sigma);
subplot(2,4,5);
plot(x_normPdf, normPdf_dL_0, 'b-', ...
    x_normPdf, normPdf_dL_P, 'm-', ...
    x_normPdf, normPdf_dL_M, 'y-', 'lineWidth', 1 )
ylim([0,0.42]); xlim([log(nu_hat)-4,log(nu_hat)+4]); axis square; box off;
xlabel('log(FR)'); ylabel('PDF')


%% Fig. 3b: Example trial

nNeurons = 100000; % number of neurons in one population
threshold = 1*baselineTheta*sqrt(2*nNeurons);
Z_U = randn(nNeurons,1); % neuron-specific and trial-independent quenched noise
Z_D = randn(nNeurons,1); 
dL_L = 0;
nu_U = nu_hat * exp(gamma_sigma*repmat(Z_U,[1,length(dL_L)])); % because dL = 0.
nu_D = nu_hat * exp(gamma_sigma*repmat(Z_D,[1,length(dL_L)]));
nu_U_sumNeurons = sum(nu_U);
nu_D_sumNeurons = sum(nu_D);
activity_by_dt = []; activity_U_in_dt = []; activity_D_in_dt = [];
dt = 0.0001;
activity_by_dt(1) = 0;
i = 1; activity_U_in_dt(1) = 0; activity_D_in_dt(1) = 0; activity = 0;
while (activity<threshold)&&(activity>-threshold)
    i = i+1;
    activity_U_in_dt(i) = poissrnd(dt*nu_U_sumNeurons);
    activity_D_in_dt(i) = poissrnd(dt*nu_D_sumNeurons);
    activity = activity + activity_U_in_dt(i) - activity_D_in_dt(i);
    activity_by_dt(i) = activity;
end
DT = length(activity_by_dt)*dt;
subplot(2,4,[2:4, 6:8]);
plot(0:dt:(DT-dt), activity_by_dt, 'k-',...
    [0, DT+0.05],[threshold, threshold],'k--',...
    [0, DT+0.05],[-threshold, -threshold],'k--',...
    [0, DT+0.05],[0, 0],'k-'); xlim([0, 0.2])
yticks( [-threshold, 0, threshold] );
yticklabels( {'-\theta', '0', '\theta'} );
xlim( [0, ceil( (DT-dt)*100 ) / 100] );


%% psychometric curve for one simulated network:

%{
Z_U = randn(nNeurons,1);
Z_D = randn(nNeurons,1); % use Z_D=Z_U to model symmetrical populations;
dL_L = linspace(-0.1,0.1,41); 
nTrials = 1000;
nu_U = nu_hat * exp(gamma_k*repmat(dL_L,[length(Z_U),1])) .* ...
    exp(gamma_sigma*repmat(Z_U,[1,length(dL_L)]));
nu_D = nu_hat * exp(gamma_k*repmat(dL_L,[length(Z_D),1])) .* ...
    exp(gamma_sigma*repmat(Z_D,[1,length(dL_L)]));
nu_U_sumNeurons = sum(nu_U);
nu_D_sumNeurons = sum(nu_D);
A_byTrialType = nu_U_sumNeurons - flip(nu_D_sumNeurons); % for DDM first passage eqiations
C2_byTrialType = nu_U_sumNeurons + flip(nu_D_sumNeurons); % for DDM first passage eqiations
Prob_U_by_dL_L = zeros(1,length(dL_L));
meanDT = zeros(1,length(dL_L));
for l = 1:length(dL_L)
    l
    DT = zeros(1,nTrials);
    trial_decision = zeros(nTrials,1);
    for trial = 1:nTrials;
        activity_by_dt = []; activity_U_in_dt = []; activity_D_in_dt = []; 
        dt = 0.001;
        activity_by_dt(1) = 0;
        i = 1; activity_U_in_dt(1) = 0; activity_D_in_dt(1) = 0; activity = 0;
        while (activity<threshold)&&(activity>-threshold)
            i = i+1;
            activity_U_in_dt(i) = poissrnd(dt*nu_U_sumNeurons(l));
            activity_D_in_dt(i) = poissrnd(dt*nu_D_sumNeurons(length(dL_L)+1-l));
            activity = activity + activity_U_in_dt(i) - activity_D_in_dt(i);
            activity_by_dt(i) = activity;
        end
        DT(trial) = length(activity_by_dt)*dt;
        if activity>0
            trial_decision(trial) = 1;
        elseif activity<0
            trial_decision(trial) = -1;
        else
            trial_decision(trial) = 0;
        end
    end
    meanDT(l) = mean(DT);
    Prob_U_by_dL_L(l) = sum(trial_decision>0)/nTrials;
end
figure;
st_err = sqrt((Prob_U_by_dL_L.*(1-Prob_U_by_dL_L))./ nTrials) ;
errorbar(dL_L,Prob_U_by_dL_L,st_err,...
'MarkerEdgeColor','black',...
'MarkerSize',3,'Marker','s','LineStyle','none','lineWidth',0.6,...
'color','black','MarkerFaceColor','none');
hold on;
% plot DDM first passage using simulated drift (A) and diffusion (c^2):
plot(dL_L,1./ (1+exp(-2*threshold*A_byTrialType./C2_byTrialType)),'k-','lineWidth',1);
xlim([dL_L(1),dL_L(end)])
xlabel('\DeltaL/L'); ylabel('p_{Up}')
%}
