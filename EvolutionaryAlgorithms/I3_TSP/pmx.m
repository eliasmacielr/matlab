% PMX
function new_pop=pmx(new_pop, winners_r1)
[p, n]=size(new_pop);
crossover=randperm(p/2);
children=zeros(p/4, n);
for i=2:2:p/2
    parent1=winners_r1(crossover(i-1), :);
    parent2=1:n;
    child=zeros(1,n);
    ndx=ceil(n/2*sort(rand(1, 2)));
    ndx(2)=ndx(1)+n/2;
    child(ndx(1):ndx(2))=parent1(ndx(1):ndx(2));
    tmp=setdiff(parent2(ndx(1):ndx(2)),parent1(ndx(1):ndx(2)));
    % try to fit tmp right to the left of the copied part from parent1
    if length(tmp) > ndx(1)-1 % tmp does not fit at the beginning
        child(1:ndx(1)-1)=tmp(1:ndx(1)-1);
        tmp_2=tmp(ndx(1):end);
        child(ndx(2)+1:ndx(2)+length(tmp_2))=tmp_2;
        child(ndx(2)+length(tmp_2)+1:n)= ...
            parent2(ndx(2)+length(tmp_2)+1:n);
    else % tmp fits to the left of the copied part
        idx=ndx(1)-length(tmp);
        child(idx:ndx(1)-1)=tmp;
        child(1:idx-1)=parent2(1:idx-1);
        child(ndx(2)+1:n)=parent2(ndx(2)+1:n);
    end
    children(i/2, :)=child;
end
new_pop(p/4+1:p/2, :)=children;
new_pop(3*p/4+1:p, :)=children;
