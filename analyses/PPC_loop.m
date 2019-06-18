
%% Reading and analyzing PPC data for all dataset (5) and all models (4):

addpath('analyses'); addpath('data');

PPC_link = 'hddm related code and data\';
dataName_vect = {'BL', 'g27', 'g3', 'g34', ...
    'motor1', 'motor2', 'motor3', 'motor4', 'motor5', ...
    'motor6', 'motor7', 'motor8', 'motor9', 'motor10'}; 
modelName_vect = {'none', 'drift', 'ic', 'both'}; 

% running over data (5):
for dataNum = 1:14
    dataName = dataName_vect(dataNum);
    % #samples in PPC per participant and trial:
    if ( strcmp(dataName, 'g27') || strcmp(dataName, 'g3') || ...
            strcmp(dataName, 'g34') )
        reps = 20; 
    else 
        reps = 100;
    end
    % #participants in data:
    if strcmp(dataName, 'BL') 
        nSubjects = 100;
        nTrials = 20;
    elseif ( strcmp(dataName, 'g27') || strcmp(dataName, 'g3') || ...
            strcmp(dataName, 'g34') )
        nSubjects = 200;
        nTrials = 500;
    else
        nSubjects = 20;
        nTrials = 20;
    end
    
    % running over models (4):
    for modelNum = 1:4
        modelName = modelName_vect(modelNum);
        dataLink = PPC_link + string(dataName) + ...
            '\informative priors 40K\' + string(modelName) + '\ppc_data.csv';
        PPC_data =  csvread( dataLink ,1 ,1 );
        
        % computing observed and expected probabilities, overall and in 5 
        % RT quantiles:
        pE_pO_qE1to5_qO1to5_A = PPC_in_model(PPC_data,reps,nSubjects);
        PPC_analysis_struct.( dataName_vect{dataNum} ).( modelName_vect{modelNum} ) ...
            = pE_pO_qE1to5_qO1to5_A;

        % RT matrices for observed and PPC-simulated data:
        [RT_observed_A,RT_expected_A,RT_observed_real_A,RT_expected_real_A] = ...
            DT_data(PPC_data,reps,nSubjects,nTrials); 
        PPC_RT.( dataName_vect{dataNum} ).( modelName_vect{modelNum} ) =  RT_expected_A;
        PPC_RT_real.( dataName_vect{dataNum} ).( modelName_vect{modelNum} ) =  RT_expected_real_A;
        if modelNum == 4
            PPC_RT.( dataName_vect{dataNum} ).observed = RT_observed_A;
            PPC_RT_real.( dataName_vect{dataNum} ).observed = RT_observed_real_A;
            PPC_RT_both.( dataName_vect{dataNum} ).RT_observed = RT_observed_A;
            PPC_RT_both.( dataName_vect{dataNum} ).RT_expected = RT_expected_A;
        end
    end
end

save('data\PPC_loop_data.mat');
