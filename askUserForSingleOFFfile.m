function url = askUserForSingleOFFfile()
    [filename, path] = uigetfile('*.off');
    url = fullfile(path, filename);
end

