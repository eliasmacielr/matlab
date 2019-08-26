fplot(@(x)x+10.*sin(5.*x)+7.*cos(4.*x))

% Let's create a random starting population of size 10
initPop=initializega(10,[0 9],'gaDemo1Eval');
% To have a look at this population the initPop is plotted as shown in 
% Figure 3.14
hold on
plot (initPop(:,1),initPop(:,2),'g+')

% Now the GA is run for one generation
[x, endPop]=ga([0 9],'gaDemo1Eval',[],initPop,[1e-6 1 1],'maxGenTerm',1,...
    'normGeomSelect',[0.08],['arithXover'],[2 0],'nonUnifMutation',...
    [2 1 3]);

% The resulting population is plotted and the plot is shown in Figure 3.15
plot (endPop(:,1),endPop(:,2),'ro')

% Next the GA is run for 25 generations
[x endPop]=ga([0 9],'gaDemo1Eval',[],initPop,[1e-6 1 1],'maxGenTerm',...
    25,'normGeomSelect',[0.08],['arithXover'],[2],'nonUnifMutation',...
    [2 25 3]);

disp(x)

% The resulting population is plotted and the plot is shown in Figure 3.16.
plot(endPop(:,1),endPop(:,2),'y*')
