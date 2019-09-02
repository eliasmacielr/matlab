% Binary Tournament
function [new_pop, winners_r1]=binary_tournament(pop, fitness)
% Tournament Selection - Round One
[p, n]=size(pop);
new_pop=zeros(p, n);
ts_r1=randperm(p);
winners_r1=zeros(p/2, n);
tmp_fitness=zeros(1, p/2);
for i=2:2:p
    if fitness(ts_r1(i-1)) > fitness(ts_r1(i))
        winners_r1(i/2, :)=pop(ts_r1(i), :);
        tmp_fitness(i/2)=fitness(ts_r1(i));
    else
        winners_r1(i/2, :)=pop(ts_r1(i-1), :);
        tmp_fitness(i/2)=fitness(ts_r1(i-1));
    end
end
% Tournament Selection - Round Two
ts_r2=randperm(p/2);
winners=zeros(p/4, n);
for i=2:2:p/2
    if tmp_fitness(ts_r2(i-1)) > tmp_fitness(ts_r2(i))
        winners(i/2, :)=winners_r1(ts_r2(i), :);
    else
        winners(i/2, :)=winners_r1(ts_r2(i-1), :);
    end
end
new_pop(1:p/4, :)=winners;
new_pop(p/2+1:3*p/4, :)=winners;
