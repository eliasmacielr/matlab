function result=fitness(population)
% Returns the fitness for each row in a population
result=sum(population, 2);

function result=distribution(population)
% Takes the population data and returns the population
% data with each genotype paired with its fraction of
% the total fitness of the population
genotypes=noduplicates(population);
total_fitness=sum(fitness(genotypes));
result=[(fitness(genotypes)/total_fitness), genotypes];
