% ORDER 1
function new_pop=order_1(new_pop, winners)
[p, n]=size(new_pop);
crossover=randperm(p/2);
children=zeros(p/4, n);
for i=2:2:p/2
    parent1=winners(crossover(i-1), :);
    parent2=winners(crossover(i), :);
    child=parent2;
    ndx=ceil(n*sort(rand(1, 2)));
    while ndx(1)==ndx(2)
        ndx=ceil(n*sort(rand(1, 2)));
    end
    tmp=parent1(ndx(1):ndx(2));
    for j=1:length(tmp)
        child(child==tmp(j))=0;
    end
    child=[child(1:n) tmp];
    child=nonzeros(child)';
    children(i/2, :)=child;
end
new_pop(p/4+1:p/2, :)=children;
new_pop(3*p/4+1:p, :)=children;
