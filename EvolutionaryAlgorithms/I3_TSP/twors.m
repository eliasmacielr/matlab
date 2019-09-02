% TWORS
function new_pop=twors(new_pop, mutate_rate)
[p, n]=size(new_pop);
mutate=randperm(p/2);
num_mutate=round(mutate_rate*p/2);
for i=1:num_mutate
    ndx=ceil(n*sort(rand(1, 2)));
    while ndx(1)==ndx(2)
        ndx=ceil(n*sort(rand(1, 2)));
    end
    if mod(ndx(2)-ndx(1)+1, 2) ~= 0
        ndx(2)=ndx(2)-1;
    end
    subarray=new_pop(p/2+mutate(i), ndx(1):ndx(2));
    idxs=randperm(length(subarray));
    for j=1:((ndx(2)-ndx(1)+1)/2)
        subarray([idxs(j) idxs(j*2)])=subarray([idxs(j*2) idxs(j)]);
    end
    new_pop(p/2+mutate(i), ndx(1):ndx(2))=subarray;
end
