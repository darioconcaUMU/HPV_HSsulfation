old = cd;
selPath = uigetdir(cd, 'Select folder');
if selPath == 0
    warning("No folder was selected. Analisis is cancelled");
    return
end
cd(selPath)
extract_IR;
cd(old);

function extract_IR

files = dir;
fileNameCumul = {files(:).name};

direct = [files(:).isdir];

for i=3:length(fileNameCumul)
    if direct(i) == 1 && fileNameCumul{i}(1)~='$' && fileNameCumul{i}(1)~='.'
        cd(fileNameCumul{i})
        display(fileNameCumul{i});
        extract_IR
    elseif ~isempty(regexp(fileNameCumul{i},'\w*cumulative.mat','ONCE'))
        
        clear spot_record metadata tracks track_props duration cum_duration...
        arrival cum_arrival particle_frame frame_props surf_coverage...
        gof_analysis fit_association fit_dissociation
    
        L = 1;
        duration = [];
        currentfile = fileNameCumul{i};
        load(currentfile);
        table_xls{L,1} = currentfile;
        L = L+1;
        
        if isempty(cum_duration)
            table_xls{L,1} = 'Time last dissociation';
            L = L+1;
            table_xls{L,1} = 'Initial number';
            L = L+1;
            table_xls{L,1} = 'Value at 2500s';
            L = L+1;
            table_xls{L,1} = 'Value End';
            L = L+1;
            table_xls{L,1} = 'Irr. Frac. 2500';
            L = L+1;
            table_xls{L,1} = 'Irr. Frac. End';
        else 
            value1 = cum_duration(1);
            temp = find(duration>2500,1);
            if isempty(temp)
                value2500 = cum_duration(end);
            else
                value2500 = cum_duration(temp-1);
            end
            valueEnd = cum_duration(end);
            table_xls{L,1} = 'Time last dissociation';
            table_xls{L,2} = duration(end);
            L = L+1;
            table_xls{L,1} = 'Initial number';
            table_xls{L,2} = value1;
            L = L+1;
            table_xls{L,1} = 'Value at 2500s';
            table_xls{L,2} = value2500;
            L = L+1;
            table_xls{L,1} = 'Value End';
            table_xls{L,2} = valueEnd;
            L = L+1;
            table_xls{L,1} = 'Irr. Frac. 2500';
            table_xls{L,2} = value2500/value1*100;
            L = L+1;
            table_xls{L,1} = 'Irr. Frac. End';
            table_xls{L,2} = valueEnd/value1*100;
        end        
        writecell(table_xls,'analysis irreversible fraction.xlsx');
    end
end
cd('..');
return

end