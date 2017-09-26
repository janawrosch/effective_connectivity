%scipt descriction
%This function calulates the transfere Entropy value for a cell pair.
%Input parameter:
% data: binary array, rows = cells, columns = time frames
% k : time frames in the past which are taken into account for estimating
% the future, mostly 2
% bin_size : reduces the data, to reduce the error of time resolution:
% mostly 3

% Output parameter:
% results : Matrix containg TE values for each connections.
%           Rows= source cell, columns= target cell
% results_list: contains connections with TE values (vector of results)
% results_vec_bin_true : contains TE values in order of results_list and
%                        and '1'&'0' for true or false connection(only for simulated data)
% temporary parameters:
% counting_list: list of pattern appearing in data

% Edit by Vicky von Einem, Sep 16



function [results]= auswTE16(data)

%% Taking the time
dispstat('', 'init')
dispstat('TE calculation started', 'timestamp','keepthis');

bin_size =3;
k=2;

% combining bins of the number of bin_size
num_col = size(data,2); % Number of coloumns
num_rows = size(data,1);
results = zeros(num_rows,num_rows);
result_list=[0,0,0];
columns_tail = mod(num_col,bin_size); %Rest bigger than bin_size
if  columns_tail ~= 0 % if rest is bigger than 0, than I have to shorten data
    data_end = num_col-columns_tail;
else
    data_end = num_col;
end

%%
% reduce data by binding bins of one bin window together
words_array=zeros(num_rows,data_end/bin_size);
res=0;
%%%% bin_size reduction jumping, not overlapping
for i=0:(data_end/bin_size)-1
words_array(:,i+1)=sum(data(:,1+i*bin_size:i*bin_size+bin_size),2)>=1;
end

num_col=size(words_array,2);
%%
m=1;
% h = waitbar(0,'Please wait...');
steps = 1000;
for z= 1:num_rows % going through rows
%     waitbar(z / num_rows)
    for p = 1 :num_rows-z
        % counting_list=[xyx1,counts(xyx1),P(xyx1),
        %                xy,counts(xy),P(x1|xy),
        %                x1x,counts(x1x),P(x1|x),
        %                x,counts(x)]
        % P(x1|xy)=P(xyx1)/P(xy)= col3/(col5./num_col)
        %
        counting_list=[0,0,0,0,0,0,0,0,0,0,0,0];
        counting_list_r=[0,0,0,0,0,0,0,0,0,0,0,0];
        for s = 1:num_col-k
            l=1;
            words_k_x=0;
            words_k_y=0;
            words_k_x1=0;
            
            %reverse
            wordsr_k_x=0;
            wordsr_k_y=0;
            wordsr_k_x1=0;
            words_k_x1=  words_array(z,k+s);
            wordsr_k_x1 = words_array(z+p,k+s);
            l=1;
            for o=k-1:-1:0
                words_k_x = words_k_x + words_array(z,s+o)*l;
                words_k_y = words_k_y + words_array(z+p,s+o)*l;
                l=l*10;
            end
            wordsr_k_y = words_k_x;
            wordsr_k_x = words_k_y;
            
            xyx1 = words_k_x*10^(k+1) ...
                +words_k_y*10^(1)+words_k_x1;
            xy = words_k_x*10^k+ words_k_y;
            x1x = words_k_x1*10^k+ words_k_x;
            
            %reverse
            rxyx1 = wordsr_k_x*10^(k+1)...
                +wordsr_k_y*10^(1)+wordsr_k_x1;
            rxy = wordsr_k_x*10^k+ wordsr_k_y;
            rx1x = wordsr_k_x1*10^k+ wordsr_k_x;
            
            
            %% Generating counting list for one pair of cells ( Forwards ) %%%%%
            [counting_list]=generating_countinglist( counting_list,xyx1,xy,x1x,words_k_x);
            
            %% Generating Couintinng_list_r for one pari of cells ( Backwards) %%%
            [counting_list_r]=generating_countinglist(counting_list_r,rxyx1,rxy,rx1x,wordsr_k_x);
        end % The ending of running though columns
        
        counting_list = sortrows(counting_list); % sorts rows, to make counting of appearnce quicker
        counting_list_r = sortrows(counting_list_r);
        
        %%   %%%%%%%%%%%%%%%%%%%%Forwards counting of conditions %%%%%%
        %%counting X(n)in col 11
        counting_list= counting_sorted_list_xy(counting_list,10,11);
        %%% counting X(n)Y(n)in col 5
        counting_list= counting_sorted_list_xy(counting_list,4,5);
        %%% counting X(n+1)X(n) in col 8
        counting_list= sortrows(counting_list,7);
        counting_list= counting_sorted_list_xy(counting_list,7,8);
        %%
        
        counting_list(:,12)=(counting_list(:, 2)./(num_col-k)).*log2((counting_list(:,2)./...
            counting_list(:,5))./(counting_list(:,8)./counting_list(:, 11)));
        summer=0;
        summe= sum(counting_list(:,12));
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%% Backwards counting of conditions%%%%%%%%%%%%%
        %%%counting X(n)in col 11
        counting_list_r= counting_sorted_list_xy(counting_list_r,10,11);
        %%%%counting X(n)Y(n)in col 5
        counting_list_r= counting_sorted_list_xy(counting_list_r,4,5);
        %%%counting X(n+1)X(n) in col 8
        counting_list_r= sortrows(counting_list_r,7);
        counting_list_r= counting_sorted_list_xy(counting_list_r,7,8);
        
        
        %% Caltculating TE
        counting_list_r(:,12)=(counting_list_r(:, 2)./(num_col-k)).*log2((counting_list_r(:,2)./...
            counting_list_r(:,5))./(counting_list_r(:,8)./counting_list_r(:, 11)));
        
        rsumme=0;
        rsumme= sum(counting_list_r(:,12));
        
        results(z+p,z)=summe;    %TE value for connection y-->x, 2-->1
        results(z,z+p)= rsumme;  %TE value for connection x-->y, 1-->2 (reverse)
    
  
    
end
dispstat(sprintf('Currently calculating source cell %i of %i', z ,size(data,1)), 'timestamp')
end
dispstat('TE calculation finished', 'timestamp', 'keepthis')
% close(h)
end % end function

