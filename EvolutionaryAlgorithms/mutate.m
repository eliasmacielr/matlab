function result=mutate(genotype)
% Possibly mutates a genotype
result=abs(genotype - rand(size(genotype,1), size(genotype,2))<0.03);
