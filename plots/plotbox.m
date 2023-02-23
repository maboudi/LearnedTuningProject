function plotbox(group, dv, groupLabels, varargin)

% dv: dependent variable

if nargin < 4
   clr = 'k';
else
    clr = varargin{1};
end

hold on

boxplot(dv, group)
set(gca, 'XTickLabels', groupLabels', 'box', 'off')

set(findobj(gcf,'LineStyle','--'),'LineStyle','-')

a = get(get(gca,'children'),'children');   % Get the handles of all the objects

t = get(a{1},'Tag');   % List the names of all the objects 

for kk = 1 : length(t)

    obj = a{1}(kk); 

    set(obj, 'Color', clr); 
end

