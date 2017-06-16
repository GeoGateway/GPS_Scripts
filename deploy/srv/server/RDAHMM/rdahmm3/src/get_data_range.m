function [min_data, range] = get_data_range(data);

[rows, cols] = size(data);

% Flip data matrix if it looks like it is oriented the wrong way
if rows >= cols
  data = data';
end

[rows, cols] = size(data);

if nargin <= 5
  min_data = min(data,[],2);
  max_data = max(data,[],2);
  range = max_data - min_data;
end
