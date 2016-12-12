function range_norm_data = range_norm(data, min_data, range);

[rows, cols] = size(data);

range_norm_data = (data - min_data * ones(1, cols)) ./ (range * ones(1,cols));
