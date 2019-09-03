% Linear Ranking
function [new_pop, winners]=linear_ranking(pop, fitness)
% Ranking Selection
[p, n]=size(pop);
new_pop=zeros(p, n);
format long
fitness=fitnessRank(fitness, 1.4, p);
fitness=1.0./fitness; % new fitness in order to select the shortest paths
fitness=fitness.*10000; % "normalize" fitness value
winners=zeros(p/2, n);
fitness_cumperc=cumsum(fitness./sum(fitness));
selected=zeros(1,p/2);
i=1;
while length(nonzeros(selected)) < p/2
    idx=find(fitness_cumperc>=rand);
    if ismember(idx(1), selected)
        continue
    end
    winners(i, :)=pop(idx(1), :);
    selected(i)=idx(1);
    i=i+1;
end
new_pop(1:p/4, :)=winners(1:p/4, :);
new_pop(p/2+1:3*p/4, :)=winners(p/4+1:p/2, :);

function y=fitnessRank(pos, sp, nind)
    y=2-sp+2*(sp-1)*((pos-1)/(nind-1));
