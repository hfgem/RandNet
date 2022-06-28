function out = loaddatastruct(animaldir, animalprefix, datatype, days)

% out = loaddatastruct(animaldir, animalprefix, datatype)
% out = loaddatastruct(animaldir, animalprefix, datatype, days)
%
% Load the components of a data cell array of combines then into one
% variable.  Datatype is a string with the base name of the files.  If DAYS
% is omitted, all files are loaded.  Otherwise, only the files for the
% specified days will be included.

if (nargin < 4)
    days = [];
end
out = [];
datafiles = dir([animaldir, animalprefix, datatype,'*']);

% Ryan -- code to warn user if no files found
if isempty(datafiles)
	errormessage = ['Looked in %s for animal %s''s ' ...
				'datatype ''''%s'''' and failed to find it' '\n\n'];
	warning(errormessage,animaldir,animalprefix,datatype);
end

for i = 1:length(datafiles)
    if isempty(days)
        load([animaldir,datafiles(i).name]);
        eval(['out = datavaradd(out,',datatype,');']);
    else
        s = datafiles(i).name;
        fileday = str2num(s(strfind(s,datatype)+length(datatype):strfind(s,'.')-1));  %get the experiment day from the filename
        if (isempty(fileday))|(ismember(fileday,days))
            load([animaldir,datafiles(i).name]);
            if strcmp(datatype,'ripplesep1'), % file ripplesep1 is already loaded
                datatype = 'ripples';
            end
            eval(['out = datavaradd(out,',datatype,');']);
        end
    end
end

%--------------------------------------
function out = datavaradd(origvar, addcell)

out = origvar;
for i = 1:length(addcell)
    if (~isempty(addcell{i}))
        out{i} = addcell{i};
    end
end
