function result = select(popWithDistrib)
% Select some genotypes from the population,
% and possibly mutates them.
selector = rand;
total_prob = 0;
% Default to last in case of rounding error
genotype = popWithDistrib(end,2:end);
for i = 1:size(popWithDistrib,1)
    total_prob = total_prob + popWithDistrib(i,1);
    if total_prob > selector
        genotype = popWithDistrib(i,2:end);
        break;
    end
end
result = mutate(genotype);
