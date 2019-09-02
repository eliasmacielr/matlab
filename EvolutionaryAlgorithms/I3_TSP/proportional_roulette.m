% -- Proportional Roulette
function [new_pop, winners_r1]=proportional_roulette(pop, fitness)
% -- Roulette Selection - Round One
[p, n]=size(pop);
new_pop=zeros(p, n);
format long
fitness=1.0./fitness; % new fitness in order to select the shortest paths
fitness=fitness.*100000; % "normalize" fitness value
winners_r1=zeros(p/2, n);
tmp_fitness=zeros(1, p/2);
fitness_cumperc=cumsum(fitness./sum(fitness));
roulette=rand(1,p/2);
for i=1:length(roulette)
    idx=find(fitness_cumperc>=roulette(i));
    idx=idx(1);
    winners_r1(i, :)=pop(idx, :);
    tmp_fitness(i)=fitness(idx);
end
% -- Roulette Selection - Round Two
winners=zeros(p/4, n);
fitness_cumperc=cumsum(tmp_fitness./sum(tmp_fitness));
roulette=rand(1,p/4);
for i=1:length(roulette)
    idx=find(fitness_cumperc>=roulette(i));
    idx=idx(1);
    winners(i, :)=winners_r1(idx, :);
end
new_pop(1:p/4, :)=winners;
new_pop(p/2+1:3*p/4, :)=winners;
