% -- The genetic algorithm search function is implemented using the MATLAB
% -- Code shown below
function new_pop=iteretic_algorithm(pop, fitness, mutate_rate)
[p, n]=size(pop);
% -- Tournament Selection - Round One
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
% -- Tournament Selection - Round Two
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
% -- Crossover
crossover=randperm(p/2);
children=zeros(p/4, n);
% When using MOX
children2=zeros(p/4, n);
for i=2:2:p/2
    % Original code
%     parent1=winners_r1(crossover(i-1), :);
%     parent2=winners_r1(crossover(i), :);
%     child=parent2;
%     ndx=ceil(n*sort(rand(1, 2)));
%     while ndx(1)==ndx(2)
%         ndx=ceil(n*sort(rand(1, 2)));
%     end
%     tmp=parent1(ndx(1):ndx(2));
%     for j=1:length(tmp)
%         child(child==tmp(j))=0;
%     end
%     %child=[child(1:ndx(1)) tmp child(ndx(1)+1:n)];
%     % ORDER 1
%     child=[child(1:n) tmp];
%     child=nonzeros(child)';
%     children(i/2, :)=child;
    %
    % PMX
%     parent1=winners_r1(crossover(i-1), :);
%     parent2=1:n;
%     child=zeros(1,n);
%     ndx=ceil(n/2*sort(rand(1, 2)));
%     ndx(2)=ndx(1)+n/2;
%     child(ndx(1):ndx(2))=parent1(ndx(1):ndx(2));
%     tmp=setdiff(parent2(ndx(1):ndx(2)),parent1(ndx(1):ndx(2)));
%     % try to fit tmp right to the left of the copied part from parent1
%     if length(tmp) > ndx(1)-1 % tmp does not fit at the beginning
%         child(1:ndx(1)-1)=tmp(1:ndx(1)-1);
%         tmp_2=tmp(ndx(1):end);
%         child(ndx(2)+1:ndx(2)+length(tmp_2))=tmp_2;
%         child(ndx(2)+length(tmp_2)+1:n)= ...
%             parent2(ndx(2)+length(tmp_2)+1:n);
%     else % tmp fits to the left of the copied part
%         idx=ndx(1)-length(tmp);
%         child(idx:ndx(1)-1)=tmp;
%         child(1:idx-1)=parent2(1:idx-1);
%         child(ndx(2)+1:n)=parent2(ndx(2)+1:n);
%     end
%     children(i/2, :)=child;
    %
    % MOX
    parent1=1:n;
    parent2=winners_r1(crossover(i), :);
    ndx=ceil(n*sort(rand(1)));
    child1=[parent1(1:ndx-1) intersect(parent1(ndx:n), parent2, 'stable')];
    child2=[parent2(1:ndx-1) intersect(parent2(ndx:n), parent1, 'stable')];
    children(i/2, :)=child1;
    children2(i/2, :)=child2;
    %
end
% new_pop(p/4+1:p/2, :)=children;
% new_pop(3*p/4+1:p, :)=children;
% When using MOX (comment the two lines above)
new_pop(p/4+1:p/2, :)=children;
new_pop(3*p/4+1:p, :)=children2;
% -- Mutate
mutate=randperm(p/2);
num_mutate=round(mutate_rate*p/2);
for i=1:num_mutate
    ndx=ceil(n*sort(rand(1, 2)));
    while ndx(1)==ndx(2)
        ndx=ceil(n*sort(rand(1, 2)));
    end
    % TWORS
%     if mod(ndx(2)-ndx(1)+1, 2) ~= 0
%         ndx(2)=ndx(2)-1;
%     end
%     subarray=new_pop(p/2+mutate(i), ndx(1):ndx(2));
%     idxs=randperm(length(subarray));
%     for j=1:((ndx(2)-ndx(1)+1)/2)
%         subarray([idxs(j) idxs(j*2)])=subarray([idxs(j*2) idxs(j)]);
%     end
%     new_pop(p/2+mutate(i), ndx(1):ndx(2))=subarray;
    % RSM
    new_pop(p/2+mutate(i), ndx(1):ndx(2))= ...
        fliplr(new_pop(p/2+mutate(i), ndx(1):ndx(2)));
    % PSM
%     subarray=new_pop(p/2+mutate(i), ndx(1):ndx(2));
%     new_pop(p/2+mutate(i), ndx(1):ndx(2))= ...
%         subarray(randperm(length(subarray)));
end
