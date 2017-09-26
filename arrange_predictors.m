function [ predictors ] = arrange_predictors( correlation_matrix_xcorr, correlation_matrix_MI, correlation_matrix_JE, correlation_matrix_TE, correlation_matrix_GTE, correlation_matrix_SC, correlation_matrix_PP, correlation_matrix_PP_samebin )
%Arrange predictors from all the correlation matrices

    n_n=size(correlation_matrix_xcorr,1); % Anzahl der Zellen (vermutlich immer 30)
    predictors=zeros(n_n^2-n_n,1);
    c=1;
    for s=1:n_n
        for t=1:n_n
            if s~=t % Wenn die Werte auf der Diagonalen liegen werden sie ignoriert und nicht mit in die predictors geschrieben
                predictors(c,1)=correlation_matrix_xcorr(s,t);
                predictors(c,2)=correlation_matrix_MI(s,t);
                predictors(c,3)=correlation_matrix_JE(s,t);
                predictors(c,4)=correlation_matrix_TE(s,t);
                predictors(c,5)=correlation_matrix_GTE(s,t);
                predictors(c,6)=correlation_matrix_SC(s,t);
                predictors(c,7)=correlation_matrix_PP(s,t);
                predictors(c,8)=correlation_matrix_PP_samebin(s,t);
                c=c+1;
            end
        end
    end



end

