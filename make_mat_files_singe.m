expressionPath = '/home/ubuntu/CPTAC_single_cell_splits_subset_pca_t/';

% Use dir function with wildcard to get file information
files = dir([expressionPath '/*']);

% Loop through each file and process it
for i = 1:length(files)
    % Construct the full filename with path
    if ~files(i).isdir
	filename = fullfile(expressionPath, files(i).name);
    	fid = fopen(filename, 'r');
    	%data_line = fgetl(fid); % Read a line
	%fclose(fid);
	%fid = fopen(filename, 'r');
	data = textscan(fid, [formatSpec repmat('%q', 1, inf)], 'Delimiter', '\t');
	%data2 = cellfun(@(x) x(2:end), data(2:end), 'UniformOutput', false);
	disp(data)
	size(data)
	%disp(data2)
	%size(data2)
	data3 = cat(2, data2{:});
    	fclose(fid);

   	% Read the pseudotime file.
   	%pFilename = join(["pseudotime_", split(filename, "_pca.tsv")(1), ".csv"]);
   	%pFilename
   end  
end
