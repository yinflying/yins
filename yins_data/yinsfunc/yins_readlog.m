function D = yins_readlog(fname)
fprintf('opening file: %s ... [', fname);
fid = fopen(fname);
if(fid == -1)
    disp('file open error');
    D = -1;
    return;
end
j = 0;
tline = fgetl(fid);
D = struct;
name_list = {};data_list = {};
while ischar(tline)
    data = strsplit(tline, {'\s+'}, 'DelimiterType', 'RegularExpression');
    if(strcmp(data{1}(1:5), 'TRACE'))
        for i = 2:2:length(data)
            pos = find_namepos(data{i}, name_list);
            if(pos > 0)
                data_list{pos}(end+1,1) = str2double(data{i+1});
            else
                name_list{end+1} = data{i};
                data_list{end+1} = str2double(data{i+1});
            end
        end
    end
    tline = fgetl(fid);
    j = j + 1;
    if(j == 1)
        fprintf('%7i]\n',j);
    end
    if(mod(j,1000) == 0)
        fprintf('\b\b\b\b\b\b\b\b\b');
        str = sprintf('%7i]',j);
        fprintf('%s\n',str);
    end
end
fprintf('\b\b\b\b\b\b\b\b');
str = sprintf('%7i]',j);
fprintf('%s\n',str);
fclose(fid);
for i = 1:length(name_list)
    eval([ 'D.' name_list{i} ' = data_list{' num2str(i) '};']);
end

function pos = find_namepos(name,  name_list)
for i = 1:length(name_list)
    if(strcmp(name, name_list{i}) == true)
        pos = i;return;
    end
end
pos = -1;return;

