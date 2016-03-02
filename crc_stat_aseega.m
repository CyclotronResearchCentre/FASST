function [comp_res] = crc_stat_aseega(CM,flags)

% Three levels of comparison,
%   matches on 2 conditions : wake and sleep
%   matches on 2 conditions : REM  - NREM
%   matches on 5 conditions : wake - S1 - S2 - S3&S4 - REM
% NB: MT and unscorable windows are not included in statistic
% 
% Po = Observed agreement : number of matches over total number of windows scored
% kappa of Cohen = (P0 - Pe) / (1 - Pe) with Pe an expected probablity
%
% With FASST chosen as gold standard, and for each state:
%   TPR = True Positive Ration (sensitivity), TPR = matches in state S over
%   the total number of windows scored as S with FASST
%
%   PPV = Positive Predictive Ratio, PPV = matches in state S over the
%   total number of windows scored as S with ASEEGA

%%

% if more than one CM, they are summed
CM = sum(CM,3);
N  = sum(sum(CM(1:6,1:6))); % Number of elements in the  sleep stages
N_n_rem = sum(sum(CM(3:6,3:6))); % Number of elements in REM and NREM (S2 to S5 => 3,4,5,6)

% Comparison with 2 conditions / wake(1) and sleep(2 to 6)
cond_c2 = [CM(1,1) sum(CM(1,2:6)); sum(CM(2:6,1)) sum(sum(CM(2:6,2:6)))];
Po(1) = sum(diag(cond_c2))/N;
Pe(1) = 1/N^2*sum(sum(cond_c2).*sum(cond_c2,2)'); 
k(1) = (Po(1)-Pe(1))/(1-Pe(1));

% Comparison with 2 conditions / NREM(3 to 5)-REM(6)
cond_c3 = [sum(sum(CM(3:5,3:5))) sum(CM(3:5,6));...
     sum(CM(6,3:5)) CM(6,6)];
Po(2) = sum(diag(cond_c3))/N_n_rem;
Pe(2) = 1/N_n_rem^2*sum(sum(cond_c3).*sum(cond_c3,2)'); 
k(2) = (Po(2) - Pe(2))/(1-Pe(2));
TPR = diag(cond_c3)./sum(cond_c3,2);
PPV = diag(cond_c3)./sum(cond_c3)';

% Comparison with 5 conditions / wake(1)-S1(2)-S2(3)-S3&S4(4to5)-REM(6)
cond_c5 = [CM(1,1) CM(1,2) CM(1,3) sum(CM(1,4:5)) CM(1,6);...
     CM(2,1) CM(2,2) CM(2,3) sum(CM(2,4:5)) CM(2,6);...
     CM(3,1) CM(3,2) CM(3,3) sum(CM(3,4:5)) CM(3,6); ...
     sum(CM(4:5,1)) sum(CM(4:5,2)) sum(CM(4:5,3)) sum(sum(CM(4:5,4:5))) sum(CM(4:5,6)); ... 
     CM(6,1) CM(6,2) CM(6,3) sum(CM(6,4:5)) CM(6,6)];
Po(3) = sum(diag(cond_c5))/N;
Pe(3) = 1/N^2*sum(sum(cond_c5).*sum(cond_c5,2)'); 
k(3) = (Po(3) - Pe(3))/(1-Pe(3));

TPR_tot = diag(cond_c5)./sum(cond_c5,2);
PPV_tot = diag(cond_c5)./sum(cond_c5)';

if flags.addpie
    label = {'W_a','S1_a','S2_a','S3,4_a','REM_a'};
    label2 = {'NREM_a', 'REM_a'};
    [r c v]=find(cond_c5);
    explode = eye(5);
    figure;
    subplot(2,3,1)
    pie(cond_c5(1,c(r==1)), explode(1,c(r==1)),label(c(r==1))),title(['W : ', num2str(sum(cond_c5(1,:))), ' scoring windows in FASST'])
    subplot(2,3,2)
    pie(cond_c5(2,c(r==2)),explode(2,c(r==2)),label(c(r==2))),title(['S1 : ', num2str(sum(cond_c5(2,:))), ' scoring windows in FASST'])
    subplot(2,3,3)
    pie(cond_c5(3,c(r==3)),explode(3,c(r==3)),label(c(r==3))),title(['S2 : ', num2str(sum(cond_c5(3,:))), ' scoring windows in FASST'])
    subplot(2,3,4)
    pie(cond_c5(4,c(r==4)),explode(4,c(r==4)),label(c(r==4))),title(['S3 & S4 : ', num2str(sum(cond_c5(4,:))), ' scoring windows in FASST'])
    subplot(2,3,5)
    pie(cond_c5(5,c(r==5)),explode(5,c(r==5)),label(c(r==5))),title(['REM : ', num2str(sum(cond_c5(5,:))), ' scoring windows in FASST'])

    [r c v]=find(cond_c3);
    figure;
    subplot(121),pie(cond_c3(1,c(r==1)),label2(c(r==1))),title(['NREM (S2-3-4) : ', num2str(sum(cond_c3(1,:))), ' scoring winsows in FASST'])
    subplot(122),pie(cond_c3(2,c(r==2)),label2(c(r==2))),title(['REM : ', num2str(sum(cond_c3(2,:))), ' scoring winsows in FASST'])
end

% artifact pages
MTf = CM(7,:);
MTa = CM(:,7);
MTT = sum([MTf MTa']);
Po(4) = CM(7,7)/MTT;

%add in structure comp_res
comp_res.CM = cond_c5;
comp_res.Po = Po;
comp_res.kappa = k;
comp_res.tpr = TPR';
comp_res.ppv = PPV';
comp_res.tpr_tot = TPR_tot';
comp_res.ppv_tot = PPV_tot';
