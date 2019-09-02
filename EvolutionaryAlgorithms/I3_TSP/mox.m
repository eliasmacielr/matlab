% MOX
function new_pop=mox(new_pop, winners_r1)
[p, n]=size(new_pop);
crossover=randperm(p/2);
children=zeros(p/4, n);
children2=zeros(p/4, n);
for i=2:2:p/2
    parent1=1:n;
    parent2=winners_r1(crossover(i), :);
    ndx=ceil(n*sort(rand(1)));
    child1=[parent1(1:ndx-1) intersect(parent1(ndx:n), parent2, 'stable')];
    child2=[parent2(1:ndx-1) intersect(parent2(ndx:n), parent1, 'stable')];
    children(i/2, :)=child1;
    children2(i/2, :)=child2;
end
new_pop(p/4+1:p/2, :)=children;
new_pop(3*p/4+1:p, :)=children2;
