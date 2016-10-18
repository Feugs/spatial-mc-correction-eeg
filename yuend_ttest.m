function [h_yuen, p, CI, t_yuen, diff, se, tcrit, df] = yuend_ttest(cond_1_data, cond_2_data, percent, alpha)
%
% Computes t_yuen (Yuen's T statistic) to compare the trimmed means of two
% dependent groups. Each group needs the same number of data points.
%  
% 
% Inputs:
%
%   cond_1_data             vector of observations in group/condition 1
%   cond_2_data             vector of observations in group/condition 2
%   percent                 percent trimming, must be between 0 and 50 
%                           Default = 20
%   alpha                   nominal alpha level. Default = 0.05
%
% Outputs:
%
%   t_yuen.tstat            Yuen T statistic. t_yuen is distributed approximately
%                           as Student's t with estimated degrees of freedom, df.
%   diff                    difference between trimmed means of cond_1_data 
%                           and cond_2_data.
%   se                      standard error
%   CI                      confidence interval around the difference
%   p                       p-value
%   tcrit                   1 - alpha / 2 quantile of the Student's t distribution 
%                           with adjusted degrees of freedom                         
%   df                      degrees of freedom
%
% See Wilcox (2012), Introduction to Robust Estimation and Hypothesis
% Testing (3rd Edition), page 195-198 for a description of the Yuen
% procedure for dependent groups.
%
% _________________________________________________________________________
% 
% Copyright (c) 2016 Daniel Feuerriegel
%
%
% This file is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
% Modified from the limo_yuend_ttest function in the LIMO Toolbox:
% Pernet, C.R., Chauveau, N., Gaspar, C. & Rousselet, G.A. 
% LIMO EEG: a toolbox for hierarchical LInear Modeling of EletroEncephaloGraphic data.
% Computational Intelligence and Neuroscience, Volume 2011 (2011), Article ID 831409, 
% 11 pages, doi:10.1155/2011/831409
% Original version copyright (C) LIMO Team 2010 under a GNU Lesser General Public License (LGPL)
%
% original version: GAR, University of Glasgow, Dec 2007
% 3D, standalone version: GAR, University of Glasgow, June 2010
% GAR fixed bug for covariance / se
%


% Defaults (if not specified in input arguments)
if nargin < 4 
    alpha = .05;
end
if nargin < 3
    percent = 20;
end

if isempty(cond_1_data) || isempty(cond_2_data) 
    error('yuend_ttest:InvalidInput', 'data vectors cannot have length = 0');
end

if (percent >= 100) || (percent < 0)
    error('yuend_ttest:InvalidPercent', 'PERCENT must be between 0 and 50.');
end

if percent >= 50
    error('yuend_ttest:InvalidPercent', 'PERCENT cannot be 50 or higher, use a method for medians instead.');
end

% number of trials
n_cond_1 = length(cond_1_data);
n_cond_2 = length(cond_2_data);

if n_cond_1 ~= n_cond_2
    error('yuend_ttest:InvalidInput', 'input data vectors must be the same size.');
else
    n = n_cond_1;
end

g = floor((percent / 100) * n); % number of items to winsorize and trim
n_trimmed = n - 2 .* g; % effective sample size after trimming

% winsorise cond_1_data (the first dataset)
cond_1_sorted = sort(cond_1_data);
loval = cond_1_sorted(g + 1);
hival = cond_1_sorted(n - g);
winsorized_cond_1 = cond_1_data;
winsorized_cond_1(winsorized_cond_1 <= loval) = loval;
winsorized_cond_1(winsorized_cond_1 >= hival) = hival;

winsorized_var_cond_1 = var(winsorized_cond_1, 0); % Calculate the winsorized variance

% winsorise cond_2_data (the second dataset)
cond_2_sorted = sort(cond_2_data);
loval = cond_2_sorted(g + 1);
hival = cond_2_sorted(n - g);
winsorized_cond_2 = cond_2_data;
winsorized_cond_2(winsorized_cond_2 <= loval) = loval;
winsorized_cond_2(winsorized_cond_2 >= hival) = hival;

winsorized_var_cond_2 = var(winsorized_cond_2, 0); % Calculate the winsorized variance

% yuen's estimate of standard errors for cond_1_data and cond_2_data
d_cond_1 = (n - 1) .* winsorized_var_cond_1;
d_cond_2 = (n - 1) .* winsorized_var_cond_2;

% covariance of winsorized samples
tmp = cov(winsorized_cond_1, winsorized_cond_2);
winsorized_covariance = tmp(1, 2);

winsorized_cov_cond_1_2 = (n - 1) .* winsorized_covariance;


% trimmed means
mean_cond_1 = mean(cond_1_sorted(g + 1 : n - g));
mean_cond_2 = mean(cond_2_sorted(g + 1 : n - g));

diff = mean_cond_1 - mean_cond_2; % Calculate difference in trimmed means

df = n_trimmed - 1; % Calculate degrees of freedom
se = sqrt( (d_cond_1 + d_cond_2 - 2 .* winsorized_cov_cond_1_2) ./ (n_trimmed .* (n_trimmed - 1)) ); % Calculate standard error

t_yuen = diff./se; % Calculate yuen's t

p = 2 * (1 - tcdf(abs(t_yuen), df)); % 2-tailed probability

if p < alpha;
    h_yuen = 1;
else
    h_yuen = 0;
end

tcrit = tinv(1 - alpha ./ 2, df); % 1-alpha./2 quantile of Student's distribution with df degrees of freedom
 
CI(1) = diff - tcrit .* se; 
CI(2)= diff + tcrit .* se;

end