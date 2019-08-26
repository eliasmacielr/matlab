function result=reproduce(population)
% Returns next generation of a population
fitted_pop=distribution(population);
for i=1:(size(population,1)/2)
    x=select(fitted_pop);
    y=select(fitted_pop);
    [x,y]=crossover(x,y);
    result(2*i-1,:)=x;
    result(2*i,:)=y;
end
