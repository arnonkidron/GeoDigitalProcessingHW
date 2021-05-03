function urls = askUserForMultipleOFFfiles()
    % ask the user to select the files 
    [filenames, paths] = uigetfile('*.off', ...
        'Select files to process', ...
        'MultiSelect', 'on');
    urls = fullfile(paths, filenames);
    
    % if the user has selected one file, put it in a cell, for consistency
    if(~iscell(urls))
        urls = {urls};
    end

end
