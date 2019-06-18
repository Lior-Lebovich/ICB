function [RT_observed,RT_expected,RT_observed_real,RT_expected_real] = ...
    DT_data(PPC_data,reps,nSubjects,nTrials)

    RT_observed = nan(nSubjects,nTrials);
    RT_expected = nan(nSubjects,nTrials*reps);
    RT_observed_real = nan(nSubjects,nTrials);
    RT_expected_real = nan(nSubjects,nTrials*reps);
    
    for ii = 1:nSubjects
        % Abs so no negative RTs (coded as such for lower boundary responses)
        subj_data = abs(  PPC_data( PPC_data(:, 6 ) == (ii-1), : )  );
        RT_expected(ii,:) = [subj_data( :, 3 )' nan(1,nTrials*reps-size( subj_data, 1))];
        subj_nTrials = size( subj_data, 1) / reps;
        RT_observed(ii,:) = [subj_data( 1:subj_nTrials, 7 )' nan(1,nTrials-subj_nTrials)];
        
        subj_data_real = PPC_data( PPC_data(:, 6 ) == (ii-1), : );
        RT_expected_real(ii,:) = [subj_data_real( :, 3 )' nan(1,nTrials*reps-size( subj_data_real, 1))];
        RT_observed_real(ii,:) = [subj_data_real( 1:subj_nTrials, 7 )' nan(1,nTrials-subj_nTrials)];
    end


end

