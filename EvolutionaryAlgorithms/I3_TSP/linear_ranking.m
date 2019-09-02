% Linear Ranking
function [new_pop, winners_r1]=linear_ranking(pop, fitness)
% Ranking Selection
[p, n]=size(pop);
new_pop=zeros(p, n);
format long
fitness=fitnessRank(fitness, 1.4, p);
fitness=1.0./fitness; % new fitness in order to select the shortest paths
fitness=fitness.*10000; % "normalize" fitness value
winners_r1=zeros(p/2, n);
fitness_cumperc=cumsum(fitness./sum(fitness));
roulette=rand(1,p/2);
for i=1:length(roulette)
    idx=find(fitness_cumperc>=roulette(i));
    idx=idx(1);
    winners_r1(i, :)=pop(idx, :);
end
new_pop(1:p/4, :)=winners_r1(1:p/4, :);
new_pop(p/2+1:3*p/4, :)=winners_r1(p/4+1:p/2, :);

function y=fitnessRank(pos, sp, nind)
    y=2-sp+2*(sp-1)*((pos-1)/(nind-1));
