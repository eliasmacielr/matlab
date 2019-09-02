% RSM
function new_pop=rsm(new_pop, mutate_rate)
[p, n]=size(new_pop);
mutate=randperm(p/2);
num_mutate=round(mutate_rate*p/2);
for i=1:num_mutate
    ndx=ceil(n*sort(rand(1, 2)));
    while ndx(1)==ndx(2)
        ndx=ceil(n*sort(rand(1, 2)));
    end
    new_pop(p/2+mutate(i), ndx(1):ndx(2))= ...
        fliplr(new_pop(p/2+mutate(i), ndx(1):ndx(2)));
end
