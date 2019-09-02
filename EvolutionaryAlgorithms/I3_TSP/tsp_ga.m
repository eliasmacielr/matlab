function varargout=tsp_ga(varargin)
% the function tsp ga finds a (near) optimal solution
% to the Traveling
% Salesman Problem by setting up a Genetic Algorithm
% (GA) to search for the
% shortest path (least distance needed to travel to
% each city exactly once)

% -- Declare the parameters used in TSP
narginchk(0, 9);
num_cities=50;
cities=10*rand(num_cities, 2);
pop_size=10*num_cities;
num_iter=1000;
mutate_rate=0.1;
show_progress=1;
show_results=0;

% -- Declare the Process Inputs
cities_flag=0;
option_flag=0;
for var=varargin
    if option_flag
        if ~isfloat(var{1}), error(['Invalid value for option' ...
                upper(option)]);
        end
        switch option
            case 'popsize', pop_size=4*ceil(real(var{1}(1))/4);
                option_flag=0;
            case 'mrate', mutate_rate=min(abs(real(var{1}(1))), 1);
                option_flag=0;
            case 'numiter', num_iter=round(real(var{1}(1)));
                option_flag=0;
            otherwise, error(['Invalid option ' upper(option)])
        end
    elseif ischar(var{1})
        switch lower(var{1})
            case '-noplot', show_progress=0;
            case '-results', show_results=1;
            otherwise, option=lower(var{1}); option_flag=1;
        end
    elseif isfloat(var{1})
        if cities_flag
            error('CITIES or NUM CITIES may be specified, but not both');
        end
        if length(var{1})==1
            num_cities=round(real(var{1}));
            if num_cities < 2
                error('NUM CITIES must be an integer greater than 1');
            end
            cities=10*rand(num_cities, 2); cities_flag=1;
        else
            cities=real(var{1});
            [num_cities, nc]=size(cities); cities_flag=1;
            if or(num_cities < 2, nc ~=2)
                error('CITIES must be an Nx2 matrix of floats, with N > 1')
            end
        end
    else
        error('Invalid input argument.')
    end
end

% -- Construction of the Distance Matrix by using Distance measurement
% -- formula
dist_matx=zeros(num_cities);
for ii=2:num_cities
    for jj=1:ii-1
        dist_matx(ii, jj)=sqrt(sum((cities(ii, :)-cities(jj,:)).^2));
        dist_matx(jj, ii)=dist_matx(ii, jj);
    end
end

% -- The cities and the distance matrix are plotted. The plot is shown in
% -- Figure 3.20
if show_progress
    figure(1)
    subplot(2, 2, 1)
    plot(cities(:,1), cities(:,2), 'b.')
    if num_cities < 75
        for c=1:num_cities
            text(cities(c, 1), cities(c, 2), [' ' num2str(c)], 'Color', ...
            'k', 'FontWeight', 'b')
        end
    end
    title([num2str(num_cities) ' Cities'])
    subplot(2, 2, 2)
    imagesc(dist_matx)
    title('Distance Matrix')
    colormap(flipud(gray))
end

% -- Initalize Population in a random manner
pop=zeros(pop_size, num_cities);
pop(1, :)=(1:num_cities);
for k=2:pop_size
    pop(k, :)=randperm(num_cities);
end

% -- Calculation of the best route
if num_cities < 25
    display_rate=1;
else
    display_rate=10;
end
fitness=zeros(1, pop_size);
best_fitness=zeros(1, num_iter);
for iter=1:num_iter
    for p=1:pop_size
        d=dist_matx(pop(p, 1), pop(p, num_cities));
        for city=2:num_cities
            d=d+dist_matx(pop(p, city-1), pop(p, city));
        end
        fitness(p)=d;
    end
    [best_fitness(iter), index]=min(fitness);
    best_route=pop(index, :);
    if and(show_progress, ~mod(iter, display_rate))
        figure(1)
        subplot(2, 2, 3)
        route=cities([best_route best_route(1)], :);
        plot(route(:, 1), route(:, 2)', 'b.-')
        title(['Best GA Route (dist=' num2str(best_fitness(iter)) ')'])
        subplot(2, 2, 4)
        plot(best_fitness(1:iter), 'r', 'LineWidth', 2)
        axis([1 max(2, iter) 0 max(best_fitness)*1.1])
    end
    % Genetic Algorithm Search
    pop=iteretic_algorithm(pop, fitness, mutate_rate);
end
% Plotting the best fitness. The plot is shown in
% Figure 3.20
if show_progress
    figure(1)
    subplot(2, 2, 3)
    route=cities([best_route best_route(1)], :);
    plot(route(:, 1), route(:, 2)', 'b.-')
    title(['Best GA Route (dist=' num2str(best_fitness(iter)) ')'])
    subplot(2, 2, 4)
    plot(best_fitness(1:iter), 'r', 'LineWidth', 2)
    title('Best Fitness')
    xlabel('Generation')
    ylabel('Distance')
    axis([1 max(2, iter) 0 max(best_fitness)*1.1])
end
if show_results
    figure(2)
    imagesc(dist_matx)
    title('Distance Matrix')
    colormap(flipud(gray))
    figure(3)
    plot(best_fitness(1:iter), 'r', 'LineWidth', 2)
    title('Best Fitness')
    xlabel('Generation')
    ylabel('Distance')
    axis([1 max(2, iter) 0 max(best_fitness)*1.1])
    figure(4)
    route=cities([best route best route(1)], :);
    plot(route(:, 1), route(:, 2)', 'b.-')
    for c=1:num_cities
        text(cities(c, 1), cities(c, 2), [' ' num2str(c)], 'Color', ...
            'k', 'FontWeight', 'b')
    end
    title(['Best GA Route (dist=' num2str(best_fitness(iter)) ')'])
end

[~, indx]=min(best_route);
best_ga_route=[best_route(indx:num_cities) best_route(1:indx-1)];
if best_ga_route(2) > best_ga_route(num_cities)
    best_ga_route(2:num_cities)=fliplr(best_ga_route(2:num_cities));
end
varargout{1}=cities(best_ga_route, :);
varargout{2}=best_ga_route;
varargout{3}=best_fitness(iter);
