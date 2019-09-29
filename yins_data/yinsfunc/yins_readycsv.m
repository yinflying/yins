function ret = yins_readycsv(fname, type)
% read ycsv format to struct or matrix. 
% args:
%       fname: ycsv file name
%       type: (optional), output data type, 'mat' or 'struct', default 'struct'

fprintf('opening file: %s ...\n', fname);
fid = fopen(fname);
if(fid == -1)
    disp('file open error');
    ret = -1;
    return;
end
header = fgetl(fid);
header_count = 0;
while ischar(header)
    header_count = header_count + 1;
    if(~isempty(header) && strcmp(header(1), '>'))
        property = strsplit(strtrim(header(2:end)), {',','\s+'}, ...
            'DelimiterType', 'RegularExpression');
        break;
    end
    header = fgetl(fid);
end
fclose(fid);

% reading csv file
data = csvread(fname, header_count, 0);

if(exist('type', 'var') && strcmp(type, 'mat'))
    ret = data;
    return;
end

% csvdata to struct
for i = 1:length(property)
    eval(['ret.' property{i} '= data(:,' num2str(i) ');']);
end

if(isfield(ret, 'status'))
    ret.dec_status = zeros(length(ret.status), 16);
    for i = 1:length(ret.status)
        s = sscanf(dec2bin(ret.status(i),16),'%1i');
        ret.dec_status(i, :) = (flipud(s))';
    end
end
    