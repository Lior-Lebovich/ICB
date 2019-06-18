%% Fig. 4: ICBs in the Poisson network model

% model parameters:
gamma_sigma = 1;
gamma_k = 0.1330;
nu_hat = 1.2556; 
nNeurons = 100000; % number of neurons in a population 

% theta parameters (columns of Fig. 4):
baselineTheta = 0.65;
theta_34 = baselineTheta * 0.75; 
theta_3 = baselineTheta  * 1; 
theta_27 = baselineTheta * 1.75; 

% A_magnitude (below) denotes nu_U/nu_D (i.e., the ratio of sum of 
% rates in the two populations). To see how this ratio is obtained and
% results in a psychometric curve - visit the bottom comment in 'fig3.m'.
A_mid = 0.99989477;
A_low = 0.9947;
A_high = 1.0091;

dL_L = linspace(-0.21,0.21,1001); % relative offset
p = [10.^([-10 -9 -8 -7]) 0.0001:0.0001:(1-0.0001) 1-(10.^([-7 -8 -9 -10]))]; % becuase the PDF is defined for eps..1-eps

figure;


%% Figure 4a - Psychometric curves of three networks:

psycho_edge = 0.135;
p_mat = nan(3,3);
j = 0;
for theta = [theta_34 theta_3 theta_27]
    j = j+1;
    subplot(3,3,j);
    plot(dL_L, (1+exp(-2*theta*sqrt(2*nNeurons)*(1-2./(1+A_low*exp(2*gamma_k*dL_L))))).^-1, 'b-',...
    dL_L, (1+exp(-2*theta*sqrt(2*nNeurons)*(1-2./(1+A_mid*exp(2*gamma_k*dL_L))))).^-1, 'k-',...
    dL_L, (1+exp(-2*theta*sqrt(2*nNeurons)*(1-2./(1+A_high*exp(2*gamma_k*dL_L))))).^-1, 'r-',...
    'lineWidth',1.5);
    xlim([-psycho_edge*1.05,psycho_edge*1.05]); ylim([-0.05,1]); axis square; box off;
    title(['\theta = ' num2str(theta/baselineTheta) '\theta*'])
    xlabel('\DeltaL/L'); ylabel('p_{Up}');
    % for marks on ICB PDF:
    p_mat(j,:) = (1+exp(-2*theta*sqrt(2*nNeurons)*(1-2./(1+[A_low, A_mid, A_high])))).^-1;
end


%% Fig. 4b -  ICB distributions:

j = 0;
col = {'b', 'k', 'r'};
for theta = [theta_34 theta_3 theta_27]
    j = j+1;
    subplot(3,3,j+3);
    threshold = theta*sqrt(2*nNeurons);
    pr_P = (  1  /  ( theta * sqrt( 8 * pi * ( exp(gamma_sigma^2) - 1 ) ) )  ) * ...
        ( 1 ./ ( p .* (1-p) ) ) .* ...
        exp(  - ( ( log(p) - log(1-p) ).^2 )  /  ( 8*(theta^2)*( exp(gamma_sigma^2) - 1 ) )  );
    bias = (2*p) - 1;
    pr_Bias = 0.5 * pr_P; 
    area(bias, pr_Bias, 'FaceColor',[.5,.5,.5],'EdgeColor','none' )
    axis square; box off; ylim([0,max([2.5, max(pr_Bias)])])
    xlabel('ICB'); ylabel('PDF'); hold on;
    % mark ICBs of 3 psychometric curves:
    for k = 1:3
        text( -1 + 2 * p_mat(j,k), ...
            0.5 * pr_P( -1 + min( find( p > p_mat(j,k) ) ) ), ...
            '\downarrow', 'color', col{k}, 'FontSize',15, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontWeight', 'bold'); hold on;
    end
end


%% Figure 4c - Distribution of reaction times:

j = 0;
for theta = [theta_34 theta_3 theta_27]
    j = j+1;
    subplot(3,3,j+6);
    i = i+1;
    kkk = 1:2:10000;
    tEnd = 10^1;
    tt = 0.001:0.001:tEnd;
    otherPart = (...
        pi * (1/(2*(theta^2))) * nu_hat * exp(0.5*((gamma_sigma)^2)) *...
        (1 ./ sqrt( 1 + (tt * nu_hat * exp(0.5*((gamma_sigma)^2)) * (exp((gamma_sigma)^2) - 1)) ) ) .*...
        exp(0.5 * (theta^2) ./ (...
            ( tt * nu_hat * exp(0.5*((gamma_sigma)^2)) ) + ...
            (1/(exp((gamma_sigma)^2)-1)) ...
            )) ...
        );
    inf_partA_Sin = kkk .* (2*(mod(kkk,4)==1)-1);
    inf_partB_inExp = (kkk.^2) * (pi^2) * (1/8) * nu_hat * (1/(theta^2)) * ...
            exp(0.5*((gamma_sigma)^2));
    inf_sum = zeros(size(tt));
    for gg = 1:length(tt)
        inf_sum(gg) = inf_partA_Sin * exp( -tt(gg) * inf_partB_inExp  )';
    end
    pdf_t = otherPart .* inf_sum;
    area(tt, pdf_t, 'FaceColor',[.5,.5,.5],'EdgeColor','none' );
    xlim([-0.1,3]); ylim([0, ceil(max(pdf_t))])
    axis square; box off; 
    xlabel('RT [s]'); ylabel('PDF');
end
