old = cd;
selPath = uigetdir(cd, 'Select folder');
if selPath == 0
    warning("No folder was selected. Analisis is cancelled");
    return
end
cd(selPath)
extract_fit;
cd(old);

function extract_fit

files = dir;
fileNameFolder = {files(:).name};

direct = [files(:).isdir];

for i=3:length(fileNameFolder)
    
    if direct(i) == 1 && fileNameFolder{i}(1)~='$' && fileNameFolder{i}(1)~='.'
        cd(fileNameFolder{i})
        display(fileNameFolder{i});
        extract_fit
    elseif ~isempty(regexp(fileNameFolder{i},'\w*cumulative.mat','ONCE'))
        
        clear spot_record metadata tracks track_props duration cum_duration...
            arrival cum_arrival particle_frame frame_props surf_coverage...
            gof_analysis fit_association fit_dissociation
        
        L = 1;
        currentfile = fileNameFolder{i};
        load(currentfile);
        
        % Check if some analysis is not performed due to the lack of
        % datapoints
        if isempty(cum_arrival)
            fit_ass = [];
        else
            fit_ass = fit_association;
        end
        if isempty(cum_duration)
            fit_dis = [];
        else
            fit_dis = fit_dissociation;
        end
        if exist('fit_dissociation_2exp','var')
            fit_2exp = fit_dissociation_2exp;
            gof_2exp = gof_2edx;
        else
            fit_2exp = [];
            gof_2exp = [];
        end
        
        surf_coverage.slope = surf_cov_slope;
        
        table_xls = [];
        
        [table_xls, L] = build_table(table_xls, fit_ass, fit_dis, fit_2exp,...
            gof_analysis, gof_2exp, surf_coverage, currentfile, L);
        L = L+2;
        
        if exist('fraction','var')
            
            boundary = ["min", string(intensity_fraction), "max"];
            
            if isempty(fraction{1}) % compatibility with an earlier problem
                fraction(1) = [];
            end
            
            for j = 1:length(fraction)
                
                if isempty(fraction{j}.cum_arrival)
                    fit_ass = [];
                else
                    fit_ass = fraction{j}.fit_association;
                end
                if isempty(fraction{j}.cum_duration)
                    fit_dis = [];
                else
                    fit_dis = fraction{j}.fit_dissociation;
                end
                if isfield(fraction{j}, 'fit_dissociation_2exp')
                    fit_2exp = fraction{j}.fit_dissociation_2exp;
                    gof_2exp = fraction{j}.gof_2edx;
                else
                    fit_2exp = [];
                    gof_2exp = [];
                end
                
                title_sec = strcat("Intensity fraction between ", boundary(j) , " and ", boundary(j+1));
                fraction{j}.surf_coverage.slope = fraction{j}.surf_cov_slope;
                
                [table_xls, L] = build_table(table_xls, fit_ass, fit_dis, fit_2exp,...
                    gof_analysis, gof_2exp, fraction{j}.surf_coverage, title_sec, L);
                
                L = L+2;
            end
            
        end
        
        writecell(table_xls,'analysis cumulative.xlsx');
    end
end
cd('..');
return

end

function [table_xls, L] = build_table(table_xls, fit_ass, fit_dis, fit_2exp,...
    gof_analysis, gof_2exp, surf_coverage, currentfile, L)

table_xls{L,1} = currentfile;
L = L+1;
table_xls{L,1} = 'k_on';
if isempty(fit_ass)
    table_xls{L,2} = [];
    L = L+1;
else
    table_xls{L,2} = fit_ass.p1;
    table_xls{L,3} = gof_analysis.association.arr_rate_95conf(1);
    table_xls{L,4} = gof_analysis.association.arr_rate_95conf(2);
    table_xls{L,5} = 'R_sq';
    table_xls{L,6} = gof_analysis.association.rsquare;
end
L = L+1;
table_xls{L,1} = 'k_off';
if isempty(fit_dis)
    table_xls{L,2} = [];
    L = L+1;
    table_xls{L,1} = 'Irreversible fraction';
    table_xls{L,2} = [];
else
    table_xls{L,2} = fit_dis.b;
    table_xls{L,3} = gof_analysis.dissociation.koff_95conf(1);
    table_xls{L,4} = gof_analysis.dissociation.koff_95conf(2);
    table_xls{L,5} = 'R_sq';
    table_xls{L,6} = gof_analysis.dissociation.rsquare;
    L = L+1;
    table_xls{L,1} = 'Irreversible fraction';
    irr_fraction = fit_dis.c/(fit_dis.a+fit_dis.c)*100;
    table_xls{L,2} = irr_fraction;
    table_xls{L,3} = gof_analysis.dissociation.irr_frac_95conf(1);
    table_xls{L,4} = gof_analysis.dissociation.irr_frac_95conf(2);
    if ~isempty(fit_2exp)
        L = L+1;
        table_xls{L,1} = 'k_off_1';
        table_xls{L,2} = fit_2exp.b;
        table_xls{L,3} = fit_2exp.b - gof_2exp.dissociation.koff_95conf(1);
        table_xls{L,4} = fit_2exp.b + gof_2exp.dissociation.koff_95conf(1);
        rev_value = fit_2exp.a+fit_2exp.c+fit_2exp.e;
        fraction_1 = fit_2exp.a/rev_value*100;
        table_xls{L,5} = 'Frac_1';
        table_xls{L,6} = fraction_1;
        perc_error = gof_2exp.dissociation.frac1_95conf(1)*fraction_1/fit_2exp.a;
        table_xls{L,7} = fraction_1 - perc_error;
        table_xls{L,8} = fraction_1 + perc_error;
        table_xls{L,9} = 'R_sq';
        table_xls{L,10} = gof_2exp.dissociation.rsquare;
        L = L+1;
        table_xls{L,1} = 'k_off_2';
        table_xls{L,2} = fit_2exp.d;
        table_xls{L,3} = fit_2exp.d - gof_2exp.dissociation.koff_95conf(2);
        table_xls{L,4} = fit_2exp.d + gof_2exp.dissociation.koff_95conf(2);
        fraction_2 = fit_2exp.c/rev_value*100;
        table_xls{L,5} = 'Frac_2';
        table_xls{L,6} = fraction_2;
        perc_error = gof_2exp.dissociation.frac1_95conf(2)*fraction_2/fit_2exp.c;
        table_xls{L,7} = fraction_2 - perc_error;
        table_xls{L,8} = fraction_2 + perc_error;
        table_xls{L,9} = 'R_sq';
        table_xls{L,10} = gof_2exp.dissociation.rsquare;
        L = L+1;
        table_xls{L,1} = 'Irreversible fraction';
        irr_fraction_2exp = fit_2exp.e/rev_value*100;
        table_xls{L,2} = irr_fraction_2exp;
        table_xls{L,3} = gof_2exp.dissociation.irr_frac_95conf(1);
        table_xls{L,4} = gof_2exp.dissociation.irr_frac_95conf(2);
    end
end
L = L+1;
table_xls{L,1} = 'Surface coverage';
if ~isempty(surf_coverage)
    table_xls{L,2} = surf_coverage.mean;
    table_xls{L,3} = surf_coverage.slope;
end

end