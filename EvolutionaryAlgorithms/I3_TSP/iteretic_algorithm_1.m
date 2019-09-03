% -- The genetic algorithm search function is implemented using the MATLAB
% -- Code shown below
function new_pop=iteretic_algorithm_1(pop, fitness, mutate_rate)
% -- Selection
% [new_pop, winners]=binary_tournament(pop, fitness);
[new_pop, winners]=proportional_roulette(pop, fitness);
% [new_pop, winners]=linear_ranking(pop, fitness);
% -- Crossover
new_pop=order_1(new_pop, winners);
% new_pop=pmx(new_pop, winners);
% new_pop=mox(new_pop, winners);
% -- Mutate
% new_pop=twors(new_pop, mutate_rate);
new_pop=rsm(new_pop, mutate_rate);
% new_pop=psm(new_pop, mutate_rate);
