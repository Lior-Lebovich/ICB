%% supp Fig. S4: Poisson model parameters:


% the lambda vector for which we plot all 4 subplots:
lambda_vect = logspace( -1, 0, 1000 );
% define 3 particular values for theta_tilda and lambda (see equation 5):
midThresh = 0.65;
bigThresh = midThresh * (7/4);
lowThresh = midThresh * (3/4);
lambda_byThresh = @(thresh)( 1 ./ ( thresh * sqrt( 8 * ( exp(1) - 1 ) ) ) );
lambda_midThresh = lambda_byThresh(midThresh); 
lambda_bigThresh = lambda_byThresh(bigThresh); 
lambda_lowThresh = lambda_byThresh(lowThresh); 
col_mid = [.5 .5 .5];
col_big = [0 .5 0];
col_low = [.93 .69 .13];

figure; 


% computing std(p):
E_p = 0.5;
% E[p^2] (2nd moment): 
hist_p_Poisson_p2 =  @(p,Lambda) (  (Lambda/sqrt(pi)) * ( 1 ./ ( p .* (1-p) ) ) ...
    .* exp( -1 * ( Lambda ^ 2 ) * ( ( log(p) - log(1-p) ).^2 ) ) .* (p.^2)  );
E_p2 = nan(size(lambda_vect));
for i = 1:length(lambda_vect)
    lambda = lambda_vect(i);
    E_p2(i) = integral( @(p)hist_p_Poisson_p2(p,lambda), ...
        1e-10,1-(1e-10), 'AbsTol', 1e-13 );
end
variance_by_lambda = E_p2 - ( E_p ^ 2 );
% computing std(ICB):
std_midThresh = sqrt( integral( @(p)hist_p_Poisson_p2(p,lambda_midThresh), ...
    1e-10,1-(1e-10), 'AbsTol', 1e-13 ) - ( E_p ^ 2 ) );
std_bigThresh = sqrt( integral( @(p)hist_p_Poisson_p2(p,lambda_bigThresh), ...
    1e-10,1-(1e-10), 'AbsTol', 1e-13 ) - ( E_p ^ 2 ) );
std_lowThresh = sqrt( integral( @(p)hist_p_Poisson_p2(p,lambda_lowThresh), ...
    1e-10,1-(1e-10), 'AbsTol', 1e-13 ) - ( E_p ^ 2 ) );
% plot std(ICB PDF) as a function of lambda (note that V[ICB]=4*V[2p-1]):
subplot(2,2,4); 
loglog( lambda_vect, 2*sqrt(variance_by_lambda), 'k-', 'lineWidth', 1 ); hold on;
plot( lambda_midThresh, 2*std_midThresh, 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_mid ); hold on;
plot( lambda_lowThresh, 2*std_lowThresh, 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_low ); hold on;
plot( lambda_bigThresh, 2*std_bigThresh, 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_big ); hold on;
xlabel('\lambda'); ylabel('std(ICB PDF)'); box off; 
axis square; grid on;
ggg = gca;
set(gca, 'FontName', 'Times'); box off;
ggg.XMinorTick = 'on';
ggg.YMinorTick = 'on';
ylim([0.3,1]); yticks([.3,.4,.6,.8,1]);


% plot lamda as a function of theta_tilda:
subplot(2,2,2);
min_theta = ( max(lambda_vect) * sqrt(8 * ( exp(1) - 1 ) ) ) ^ -1;
max_theta = ( min(lambda_vect) * sqrt(8 * ( exp(1) - 1 ) ) ) ^ -1;
theta_vect = linspace(min_theta, max_theta, 10000);
loglog( theta_vect, lambda_byThresh( theta_vect ), ...
    'k-', 'lineWidth', 1 ); hold on;
plot( midThresh, lambda_midThresh, 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_mid ); hold on;
plot( lowThresh, lambda_lowThresh, 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_low ); hold on;
plot( bigThresh, lambda_bigThresh, 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_big ); hold on;
ylabel('\lambda'); xlabel('$$\tilde{\theta}$$', 'Interpreter', 'LaTeX'); box off;
axis square; grid on;
xlim( [0.2,3] ); 
ggg = gca;
set(gca, 'FontName', 'Times'); box off;
ggg.XMinorTick = 'on';
ggg.YMinorTick = 'on';
xticks([.2, 0.5, 1, 3]);


% plot lamda as a function of gamma*sigma:
subplot(2,2,1);
gamSig_byLambdaAndThresh = @( lambda, thresh ) ( ...
    sqrt( log( 1 + ( 1 / ( thresh*sqrt(8)*lambda )^2 ) ) ) );
lambda_byGamSig = @( GamSig, thresh )( 1 ./ ( thresh * sqrt( 8 * ( exp( GamSig.^2 ) - 1 ) ) ) );
min_gammaSig = gamSig_byLambdaAndThresh( max(lambda_vect), lowThresh );
max_gammaSig = gamSig_byLambdaAndThresh( min(lambda_vect), lowThresh );
gammaSig_vect = linspace( min_gammaSig, max_gammaSig, 10000 );
loglog( gammaSig_vect, lambda_byGamSig( gammaSig_vect, lowThresh ), ...
    'Color', col_low, 'lineWidth', 1 ); hold on;
loglog( 1, lambda_byGamSig( 1, lowThresh ), 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_low ); hold on;
min_gammaSig = gamSig_byLambdaAndThresh( max(lambda_vect), midThresh );
max_gammaSig = gamSig_byLambdaAndThresh( min(lambda_vect), midThresh );
gammaSig_vect = linspace( min_gammaSig, max_gammaSig, 10000 );
loglog( gammaSig_vect, lambda_byGamSig( gammaSig_vect, midThresh ), ...
    'Color', col_mid, 'lineWidth', 1 ); hold on;
loglog( 1, lambda_byGamSig( 1, midThresh ), 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_mid ); hold on;
min_gammaSig = gamSig_byLambdaAndThresh( max(lambda_vect), bigThresh );
max_gammaSig = gamSig_byLambdaAndThresh( min(lambda_vect), bigThresh );
gammaSig_vect = linspace( min_gammaSig, max_gammaSig, 10000 );
loglog( gammaSig_vect, lambda_byGamSig( gammaSig_vect, bigThresh ), ...
    'Color', col_big, 'lineWidth', 1 ); hold on;
loglog( 1, lambda_byGamSig( 1, bigThresh ), 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_big  ); hold on;
ylabel('\lambda'); xlabel('\gamma\sigma'); box off;
axis square; grid on;
xlim([0.3, 2]); xticks([.3, .5, 1, 2]);
ggg = gca;
set(gca, 'FontName', 'Times'); box off;
ggg.XMinorTick = 'on';
ggg.YMinorTick = 'on';


% plot lamda as a function of nu_hat (independent):
subplot(2,2,3);
loglog( [0.1,5], ones(2,1) * lambda_byThresh( lowThresh ), ...
    'Color', col_low, 'lineWidth', 1 ); hold on;
loglog( 1.2556, lambda_byThresh( lowThresh ), 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_low ); hold on;
loglog( [0.1,5], ones(2,1) * lambda_byThresh( midThresh ), ...
    'Color', col_mid, 'lineWidth', 1 ); hold on;
loglog( 1.2556, lambda_byThresh( midThresh ), 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_mid ); hold on;
loglog( [0.1,5], ones(2,1) * lambda_byThresh( bigThresh ), ...
    'Color', col_big, 'lineWidth', 1 ); hold on;
loglog( 1.2556, lambda_byThresh( bigThresh ), 'marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', col_big ); hold on;
ylabel('\lambda'); xlabel('$\bar{\nu}$','Interpreter','Latex'); box off;
ylim([10^-1,10]); xlim([0.3,2]); axis square; grid on;
ggg = gca;
set(gca, 'FontName', 'Times'); box off;
yticks( logspace(-1,0,2) ); ylim([10^-1,1]);
ggg.XMinorTick = 'on';
ggg.YMinorTick = 'on';
xticks([0.3, 0.5, 1, 2]);
