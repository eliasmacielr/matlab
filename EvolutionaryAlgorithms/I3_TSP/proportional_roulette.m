% -- Proportional Roulette
function [new_pop, winners]=proportional_roulette(pop, fitness)
% -- Roulette Selection - Round One
[p, n]=size(pop);
new_pop=zeros(p, n);
format long
fitness_roulette=n./fitness; % new fitness in order to select the shortest paths
fitness_roulette=fitness_roulette.*n; % "normalize" fitness_roulette value
winners=zeros(p/2, n);
fitness_roulette_cumperc=cumsum(fitness_roulette./sum(fitness_roulette));
roulette=rand(1,p/2);
for i=1:length(roulette)
    idx=find(fitness_roulette_cumperc>=roulette(i));
    winners(i, :)=pop(idx(1), :);
end
new_pop(1:p/4, :)=winners(1:p/4, :);
new_pop(p/2+1:3*p/4, :)=winners(p/4+1:p/2, :);
