function [x,y]=crossover(x,y)
% Possibly takes some information from one genotype and
% swaps it with information from another genotype
if rand < 0.6
    gene length=size(x,2);
    % site is between 2 and gene length
    site=ceil(rand * (gene_length-1))+1;
    tmp=x(site:gene_length);
    x(site:gene_length)=y(site:gene_length);
    y(site:gene_length)=tmp;
end
