function [pE_pO_qE1to5_qO1to5, mRTbyQuantile_O, mRTbyQuantile_E ] = ...
    PPC_in_model(PPC_data,reps,nSubjects)

    pE_pO_qE1to5_qO1to5 = nan(nSubjects,13);
    mRTbyQuantile_O = nan(nSubjects,5); %%%
    mRTbyQuantile_E = nan(nSubjects,5); %%%
    for ii = 1:nSubjects
        pE_pO_qE1to5_qO1to5( ii, 13 ) = ii - 1;
        % Abs so no negative RTs (coded as such for lower boundary responses)
        subj_data = abs(  PPC_data( PPC_data(:, 6 ) == (ii-1), : )  );
        % Calculating p_expected (given PPC) and p_observed (given data)
        pE_pO_qE1to5_qO1to5( ii, 1:2 ) = mean( subj_data(:,[4,8]) );
        % Calculating quantile RT for expected (PPC base):
        samples_in_q = size( subj_data, 1) / 5;
        RT_p_expected = sortrows( subj_data( :, [3,4] ), [1] );
        for jj = 1:5
            pE_pO_qE1to5_qO1to5( ii, 2+jj ) = ...
                mean( RT_p_expected( (1+(jj-1)*samples_in_q):(jj*samples_in_q), 2 ) );
            if nargout == 3 %%%
                mRTbyQuantile_E(ii, jj) = ...
                    mean( RT_p_expected( (1+(jj-1)*samples_in_q):(jj*samples_in_q), 1 ) );
            end %%%
        end
        % Calculating quantile RT for observed (data base):
        % Note the the 5th quantile is more noisy since it contains between 0-5
        % trials. For the case in which it has zerp trials - we will later 
        % omit the participant from the relevant analysis (conditional plot).
        subj_nTrials = size( subj_data, 1) / reps;
        trials_in_q = ceil( subj_nTrials / 5 );
        RT_p_observed = sortrows( subj_data( 1:subj_nTrials, [7,8]), [1] );
        for jj = 1:5
            pE_pO_qE1to5_qO1to5( ii, 7+jj ) = ...
                mean( RT_p_observed( (1+(jj-1)*trials_in_q):...
                min( [(jj*trials_in_q),subj_nTrials])  , 2 ) );
            if nargout == 3 %%%
                mRTbyQuantile_O(ii, jj) = ...
                    mean( RT_p_observed( (1+(jj-1)*trials_in_q):...
                    min( [(jj*trials_in_q),subj_nTrials])  , 1 ) );
            end %%%
        end
    end

end
