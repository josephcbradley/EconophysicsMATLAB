function fig = plot_graph(g, groups, group_colours, group_names)
%% Description 
%% Inputs

%% Ouputs 

%% Setup 


%% Calculation
    
%check that num nodes is equal to group list size
if numnodes(g) ~= length(groups) categorical
    error("The length of `groups` must be the same as the number of nodes in the graph.")
end

%check that group_colors is mx3
if size(group_colours, 2) ~= 3
    error("`group_colours` must have three columns")
end

%check that there are enough group colors 
if max(groups) > size(group_colours, 1)
    error("Too many groups and not enough colors!")
end


N = length(groups);
%allocate colors for graph
node_colours = zeros(N, 3);
for i = 1:N
    s = groups(i);
    node_colours(i, :) = group_colours(s, :);
end

fig = figure;
plot(g, 'Layout', 'force', 'NodeColor', node_colours);
hold on 
%how many groups are present? 
present_groups = unique(groups);
n_present_groups = length(present_groups);
%initialise an object for every present group
h=gobjects(n_present_groups,1);
%use s to track h groups
s = 1;
n_groups = size(group_colours, 1);
%for every present group, add a an invisible scatter plot so that we can
%make a legend
for i = 1:n_groups %for all possible groups
    if ~isempty(find(present_groups == i, 1)) %if the group is present
        %find its colour
        group_colour = group_colours(i, :);
        %add invisible scatter
        h(s)=scatter(nan,nan, 1, group_colour, "filled");
        %increment to keep track of h objects
        s = s + 1;
    end
end
legend(h, group_names(present_groups), 'Location', 'southoutside', 'NumColumns', 2, 'interpreter', 'latex')

end