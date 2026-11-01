function x = readCounts(filePath)
%READCOUNTS  1-s movement counts of one fly from a text file (sqw_ref.read_counts).
%
%   x = readCounts(filePath)
%
%   filePath  text file holding one non-negative integer per line, one line
%             per second, starting at ZT0
%   x         (nSeconds x 1) double, the counts
%
%   The whole file must parse as numbers: a header line, a non-numeric token
%   or a negative / fractional value is an error, not something to skip. A
%   file that is silently read short shifts every later minute's clock time.
%
%   See also SLEEPBOUTS, SQWSYMBOLS.

arguments
    filePath (1,1) string {mustBeFile}
end

fid = fopen(filePath, 'r');
if fid < 0
    error('readCounts:open', 'Cannot open %s', filePath);
end
closer = onCleanup(@() fclose(fid));
[x, n] = fscanf(fid, '%f');
if ~feof(fid)
    pos = ftell(fid);
    error('readCounts:parse', ['%s: stopped after %d numbers at byte %d; ' ...
        'the file holds something other than one number per line.'], filePath, n, pos);
end
x = x(:);
if isempty(x)
    error('readCounts:empty', '%s holds no counts.', filePath);
end
if any(~isfinite(x)) || any(x < 0) || any(x ~= round(x))
    error('readCounts:values', '%s: counts must be non-negative integers.', filePath);
end
end
