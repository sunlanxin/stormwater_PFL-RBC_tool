function TVGM_to_swmm(input,inp_Precip)
% ------------------------------------------------------------------------------------------------------------
% This function couples TVGM-Urban and SWMM model in a loose way.
% The TVGM-Urban flows is writen as a sequence of timeseries to replace the contents of the [TIMESERIES] module in the SWMM.inp file,
% and modifies START_DATE, START_TIME, REPORT_START_DATE, REPORT_START_TIME, END_DATE, END_TIME.
% The [INFLOWS] module of the SWMM.inp file use this timeseries during the control simulation 
% ------------------------------------------------------------------------------------------------------------

add = inp_Precip; 
START_DATE = inp_Precip(1,2);
START_TIME = inp_Precip(1,3);
REPORT_START_DATE = inp_Precip(1,2);
REPORT_START_TIME = inp_Precip(1,3);
END_DATE = inp_Precip(end,2);
END_TIME = inp_Precip(end,3);
oldinp = fopen(input,'rt+');% 需要读取的文件
timeseries = fopen('timeseries.txt','wt+');

Title = ["[TIMESERIES]", ";;Name	Date	Time	Value", ";;-------------- ------------------ ------------------"];
fprintf(timeseries,'%s\n',Title);    

% save the original information in [TIMESERIES]
%i = 0;
%k1 = 0;
%while ~feof(oldinp)
%    tline= fgetl(oldinp);
%    if tline == -1
%        break;
%    end
%    %disp(tline)
%    i=i+1;
%    if contains(tline,'[TIMESERIES]')
%        k1 = 1;
%        
%    elseif contains(tline,'[REPORT]')
%        k1 = 0;
%    end
%    
%    if k1 && tline~=""
%        fprintf(timeseries,'%s\n',tline);
%    end
%end
%fprintf(timeseries,'%s\n',';');

for i = 1 : length(add)
    tline = strjoin(add(i,:),'\t');
    fprintf(timeseries,'%s\n',tline);
end
fprintf(timeseries,'%s\n',' ');

frewind(oldinp)
frewind(timeseries)
i = 0;
k1 = 0;
while ~feof(oldinp)
    tline = fgetl(oldinp);
    i = i+1;
    if tline == -1
        break
    elseif contains(tline,'[TIMESERIES]')
        k1 = 1;
        while ~feof(timeseries)
            newData = fgetl(timeseries);
            if newData == -1
                break
            else
                newline{i} = newData;
            end
            i=i+1;
        end
    elseif contains(tline,'[REPORT]')
        k1 = 0;
    end
    
    if k1
        i = i-1;
    else
        if contains(tline,'START_DATE') && ~contains(tline,'REPORT')
            newline{i} = convertStringsToChars(strcat("START_DATE","             ",START_DATE));
        elseif contains(tline,'START_TIME') && ~contains(tline,'REPORT')
            newline{i} = convertStringsToChars(strcat("START_TIME","             ",START_TIME));
        elseif contains(tline,'REPORT_START_DATE')
            newline{i} = convertStringsToChars(strcat("REPORT_START_DATE","             ",REPORT_START_DATE));
        elseif contains(tline,'REPORT_START_TIME')
            newline{i} = convertStringsToChars(strcat("REPORT_START_TIME","             ",REPORT_START_TIME));
        elseif contains(tline,'END_DATE')
            newline{i} = convertStringsToChars(strcat("END_DATE","             ",END_DATE));
        elseif contains(tline,'END_TIME')
            newline{i} = convertStringsToChars(strcat("END_TIME","             ",END_TIME));
        else
            newline{i} = tline;
        end
    end
end
fclose(oldinp);
delete(input);
newinp = fopen(input,'wt+');
for k = 1:i-1
    fprintf(newinp,'%s\n',newline{k});
end
fclose(newinp);
fclose(timeseries);
delete("timeseries.txt");
end