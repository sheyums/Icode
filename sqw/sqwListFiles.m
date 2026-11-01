function files = sqwListFiles(files)
%SQWLISTFILES  Count files of a group, in the order every fly index refers to (sorted, as Python's sorted glob).
%
%   files = sqwListFiles(files)
%
%   files  a folder (searched recursively for *.txt) or a list of files.
%          Returns a sorted column string array. Fly k of a group is files(k);
%          sqw_ref's 0-based fly index i is files(i + 1).

files = string(files);
if isscalar(files) && isfolder(files)
    d = dir(fullfile(files, '**', '*.txt'));
    files = string(fullfile({d.folder}, {d.name}));
end
files = sort(files(:));
if isempty(files), error('sqw:noFiles', 'No count files given.'); end
end
