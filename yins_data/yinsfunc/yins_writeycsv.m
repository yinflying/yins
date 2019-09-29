function yins_writeycsv(data, fname)
fprintf('write struct data to ycsv format file: %s\n', fname);
fid = fopen(fname','w');
S = fieldnames(data);
fprintf(fid, '>');
for i = 1:length(S)-1
    if(strcmp(S{i}, 'week')) fprintf(fid, '%4s,', S{i});
    else fprintf(fid, '%12s,', S{i});
    end
end
fprintf(fid, '%12s', S{end});
fprintf(fid, '\n');

cell_data = struct2cell(data);
data_num = 0;
field_num = length(cell_data);
for i = 1:field_num
    if (data_num < length(cell_data{i}))
        data_num = length(cell_data{i});
    end
end

for i = 1:data_num
    for j = 1:field_num
        if(strcmp(S{j}, 'week')) fprintf(fid, '%5.0f', cell_data{j}(i));
        elseif(strcmp(S{j}, 'sec'))  fprintf(fid,  '%12.3f', cell_data{j}(i));
        elseif(strcmp(S{j}, 'status'))  fprintf(fid,  '%12.0f', cell_data{j}(i));
        elseif(strcmp(S{j}, 'ext_status'))  fprintf(fid, '%12.0f', cell_data{j}(i));
        else fprintf(fid, '%12.7f', cell_data{j}(i));
        end
        if(j ~= field_num) fprintf(fid, ',');end
    end
    fprintf(fid, '\n');
end

fclose(fid);
end 

