%% p_prev: Comparing impossible trial to the previous trial of the same type:

% Here we compute the upper and lower 90% bounds of p_prev.
% Given pUp[before impossible] and pUp[impossible], we compute the pdf.

numImpTrials = 20;
possiblePprev = 0:(1/numImpTrials):1;
CI90_okay_Up_rPrev_cImp = nan( numImpTrials + 1 );
CI90_okay_Down_rPrev_cImp = nan( numImpTrials + 1 );

for numUpPos = 0:numImpTrials % running over pUp[before impossible]
    % For the PDF of pPrev, order instantiations of the same 
    % pUp[before impossible] can be reduces to a general order. We will 
    % then compute the PDF by considering possible orders of the
    % impossible decisions. 
    genPrevVect = [ones(1,numUpPos), -1*ones(1,numImpTrials-numUpPos)]; 
    for numUpImp = 0:numUpPos %:numImpTrials % running over pUp[impossible]
        count_p_prev = zeros( 1, numImpTrials+1 );
        pdf_p_prev = zeros( 1, numImpTrials+1 );
        % For p_up = 0 or 1, the PDF of p_prev has only one non-zero value:
        if abs(10-numUpPos)==10 || abs(10-numUpImp)==10
            genImpVect = [ones(1,numUpImp), ...
                -1*ones(1,numImpTrials-numUpImp)];
            p_prev = sum( (genPrevVect .* genImpVect) == 1 ) / numImpTrials;
            pdf_p_prev( numImpTrials*p_prev + 1 ) = 1;
            PprevPDF.(['impUp' num2str(numUpImp) '_prevUp' num2str(numUpPos)]) = ...
                pdf_p_prev;
            PprevPDF.(['impUp' num2str(numUpPos) '_prevUp' num2str(numUpImp)]) = ...
                pdf_p_prev;
            CI90_okay_Up.(['impUp' num2str(numUpImp) '_prevUp' num2str(numUpPos)]) = ...
                p_prev;
            CI90_okay_Up.(['impUp' num2str(numUpPos) '_prevUp' num2str(numUpImp)]) = ...
                p_prev;
            CI90_okay_Down.(['impUp' num2str(numUpImp) '_prevUp' num2str(numUpPos)]) = ...
                p_prev;
            CI90_okay_Down.(['impUp' num2str(numUpPos) '_prevUp' num2str(numUpImp)]) = ...
                p_prev;
            CI_okay_UpPos.(['impUp' num2str(numUpPos) '_prevUp' num2str(numUpImp)]) = ...
                p_prev;
            CI_okay_UpPos.(['impUp' num2str(numUpImp) '_prevUp' num2str(numUpPos)]) = ...
                p_prev;
            CI_okay_DownPos.(['impUp' num2str(numUpPos) '_prevUp' num2str(numUpImp)]) = ...
                p_prev;
            CI_okay_DownPos.(['impUp' num2str(numUpImp) '_prevUp' num2str(numUpPos)]) = ...
                p_prev;
        % For other p_up, we compute:
        else
            %deno1st_vect = max([numUpPos-numUpImp,0]):min([numImpTrials-numUpImp,numUpPos]);
            %min_up = min([numUpPos,numUpImp]); 
            orders = numUpImp:-1:max( [0, numUpImp-numImpTrials+numUpPos] );
            for i = 1:length(orders) % running over possible orders (from high to low..)
                genImpVect = [ones(1,orders(i)), ...
                    -1*ones(1,numImpTrials-numUpImp), ...
                    ones(1,numUpImp-orders(i)) ]; % for computing p_prev only
                p_prev = sum( (genPrevVect .* genImpVect) == 1 ) / numImpTrials;
                count_p_prev( numImpTrials*p_prev + 1 ) = count_p_prev( numImpTrials*p_prev + 1 ) + ...
                    nchoosek( numUpPos, orders(i) ) * ...
                    nchoosek( numImpTrials-numUpPos, numUpImp-orders(i) );
            end
            pdf_p_prev = count_p_prev / nchoosek( numImpTrials, numUpImp );
            % verifying that the pdf of p_prev sums to 1:
            if round( sum( 1e5 * pdf_p_prev ) ) / 1e5 ~= 1
                disp([numUpPos numUpImp])
                disp('Oh no!!');
                disp(pdf_p_prev);
                disp( sum(pdf_p_prev) );
            end
            PprevPDF.(['impUp' num2str(numUpImp) '_prevUp' num2str(numUpPos)]) = ...
                pdf_p_prev;
            PprevPDF.(['impUp' num2str(numUpPos) '_prevUp' num2str(numUpImp)]) = ...
                pdf_p_prev;
            % caculating the lowest, nonsignificant pPrev 
            pPrev_down = max( [0, find(cumsum(pdf_p_prev) <= 0.05)] ) / numImpTrials;
            CI90_okay_Down_rPrev_cImp(numUpPos+1,numUpImp+1) = pPrev_down;
            CI90_okay_Down_rPrev_cImp(numUpImp+1,numUpPos+1) = pPrev_down;
            CI90_okay_Down.(['impUp' num2str(numUpImp) '_prevUp' num2str(numUpPos)]) = ...
                pPrev_down;
            CI90_okay_Down.(['impUp' num2str(numUpPos) '_prevUp' num2str(numUpImp)]) = ...
                pPrev_down;
            CI_okay_DownPos.(['impUp' num2str(numUpPos) '_prevUp' num2str(numUpImp)]) = ...
                ( min( find(pdf_p_prev > 0) ) - 1 ) / numImpTrials;
            CI_okay_DownPos.(['impUp' num2str(numUpImp) '_prevUp' num2str(numUpPos)]) = ...
                ( min( find(pdf_p_prev > 0) ) - 1 ) / numImpTrials;
            
            % caculating the highest, nonsignificant pPrev
            pPrev_up = ( numImpTrials - max( [0, find(cumsum(flip(pdf_p_prev)) <= 0.05)] ) ) / numImpTrials;
            CI90_okay_Up_rPrev_cImp(numUpPos+1,numUpImp+1) = pPrev_up;
            CI90_okay_Up_rPrev_cImp(numUpImp+1,numUpPos+1) = pPrev_up;
            CI90_okay_Up.(['impUp' num2str(numUpImp) '_prevUp' num2str(numUpPos)]) = ...
                pPrev_up;
            CI90_okay_Up.(['impUp' num2str(numUpPos) '_prevUp' num2str(numUpImp)]) = ...
                pPrev_up;
            CI_okay_UpPos.(['impUp' num2str(numUpPos) '_prevUp' num2str(numUpImp)]) = ...
                ( max( find(pdf_p_prev > 0) ) - 1 ) / numImpTrials;
            CI_okay_UpPos.(['impUp' num2str(numUpImp) '_prevUp' num2str(numUpPos)]) = ...
                ( max( find(pdf_p_prev > 0) ) - 1 ) / numImpTrials;
        end
    end
end

clearvars -except PprevPDF possiblePprev CI90_okay_Up CI90_okay_Down...
    CI90_okay_Up_rPrev_cImp CI90_okay_Down_rPrev_cImp CI_okay_UpPos CI_okay_DownPos

save( 'data/prob4pPrevNotBiggerThanP2_data_SHORT.mat' )

clear all