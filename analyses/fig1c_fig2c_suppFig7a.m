%% Model comparison using DIC: fig 1c, fig 2c, supp fig 7a: 


%% load DIC data:

% Rows of DIC_mat: 1. Bisection task, 2. g=2.7, 3. g=3, 4. g=3.4,
% 5-14. motor pairs (1-10). 
% Columns of DIC_mat: 1. unbiased drift and IC, 2. biased drift, unbiased
% IC, 3. biased IC, unbiased drift, 4. biased drift and IC
% DIC_mat_A and DIC_mat_B are the same DIC_mat, for the non-motor datasets.

addpath('analyses','data');

load('DIC_data_informative_40K');

nDatasets = 14; % bisection, g=2.7, g=3, g=3.4, 10 motors (each set was analyzed separately) 
dataName = {'Vertical bisection', 'g = 2.7', 'g = 3', 'g = 3.4', 'Motor task'};
forModelsCols = [0 0 0; 0.0431 0.5020 0.2510; 0.4980 0.1843 0.5451];
modelsNameForDIC = {'b. drift','c. IC','a. drift + IC'};
nModelsInclude = 4;
plotDataOrder = [1, 5, 4, 3, 2];


%% Run over datasets and DDMs, compute SEM of DIC and plot:

figure;

for dataNum = 1:5
    
    tit = dataName{dataNum};
    subplot( 5, 1, plotDataOrder(dataNum) );
    
    % If motor task then compute for each DDM the average and SEM of 
    % deltaDIC over all pairs. Else, compute deltaDIC based on three 
    % repetitions of the fitting procedures.
    % deltaDIC = DIC of a DDM variant - baseline DIC of a DDM with unbiased
    % drift and IC.
    if dataNum == 5 
        DIC_vect = mean( DIC_mat(5:end,2:(nModelsInclude)) - DIC_mat(5:end,1) );
        DIC_vect_sem = std( DIC_mat(5:end,2:(nModelsInclude)) - DIC_mat(5:end,1) ) ...
            ./ sqrt( length( dataNum:nDatasets ) );
    else
        forDIC_vect = [DIC_mat(dataNum,2:(nModelsInclude)); ...
            DIC_mat_A(dataNum,2:(nModelsInclude)); ...
            DIC_mat_B(dataNum,2:(nModelsInclude))];
        forDIC_vect1 = [DIC_mat(dataNum,1); DIC_mat_A(dataNum,1); ...
            DIC_mat_B(dataNum,1)];
        DIC_vect = mean( forDIC_vect - forDIC_vect1, 'omitnan' );
        DIC_vect_sem = std( forDIC_vect - forDIC_vect1, 'omitnan' ) ./ sqrt( size(forDIC_vect1,1) );
    end
    
    % plot:
    bar_plot = barh( categorical(modelsNameForDIC), DIC_vect);
    bar_plot.FaceColor = 'flat';
    bar_plot.CData(1,:) = forModelsCols(1,:);
    bar_plot.CData(2,:) = forModelsCols(2,:);
    bar_plot.CData(3,:) = forModelsCols(3,:);
    xlabel('\Delta DIC'); box off;
    title(tit);
    hold on;
    errorbar( DIC_vect, categorical(modelsNameForDIC), ...
        DIC_vect_sem, 'horizontal', 'lineStyle', 'none', ...
        'lineWidth', 1, 'color', [.5 .5 .5])
    g = gca;
    g.XDir = 'reverse';    
    
end 
