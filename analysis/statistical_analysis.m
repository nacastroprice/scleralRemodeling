% anova2 for area analysis

rep  = 1;
%         exp1    exp2    exp3
% treated 
% untreated
format long
% aream = [0.00198642   0.0015094;6.88708e-06  5.15011e-06;1.17804e-05  5.21541e-06 ];
aream = [0.00198642   0.0015094;6.88708e-06  5.15011e-06;1.17804e-05  5.21541e-06 ];

% 0.00198642   0.0015094
% 6.88708e-06  5.15011e-06
% 1.17804e-05  5.21541e-06  

arease = [5.6858e-05    0.0001718;5.0993e-09    7.2701e-09;9.7054e-09    7.4887e-09];
%SE
% 5.6858e-05    0.0001718
% 5.0993e-09    7.2701e-09
% 9.7054e-09    7.4887e-09

[~,~,stats] = anova2(aream(2:end,:),rep);




y = aream(2:end,:);
group = {'T','U'};
p = anova1(y,group);

%% This is the analysis with vector of % change in area for 2043 and 2025

load("Z:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\analysis\dataanalysisnorm_2043.mat")
load("Z:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\analysis\dataanalysisnorm_2025.mat")

diff = length(dataanalysisnorm_2043) - length(dataanalysisnorm_2025);
dataanalysisnorm_2043 = dataanalysisnorm_2043(1:end-diff,:);

dat = dataanalysisnorm_2043(1:100:end,:);
dat(:,3:4) = dataanalysisnorm_2025(1:100:end,:);
rep = 1;

group = {'T','U'};
p = anova1(dat,group);
[i,p,stats] = anova2(dat,rep);


%% For aoctool

x1 = [1:726];
x2 = x1;
x = [x1,x2]';
group1 = ones(1,length(x1));
group2 = group1 * 2;
groups = [group1,group2]';
y = [dat(:,1);dat(:,2)];

[h,atab,ctab,stats] = aoctool(x,y,groups);

multcompare(stats,0.05,'on','','s')











