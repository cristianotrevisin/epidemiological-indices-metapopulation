% Read from table (downloaded from 
% https://www.ecdc.europa.eu/en/publications-data/data-virus-variants-covid-19-eueea)
% public archive
variants_raw = readtable("ECDC.csv");

% Select Italy
variants_raw = variants_raw(strcmp(variants_raw.country,'Italy')==1,:);

% Remove columns
variants_raw = variants_raw(:,[3 4 9 12]);

t0_2020 = datetime(2020, 1, 5, 'Format', 'dd-MMMM-yyyy');  % First day of the year
t0_2021 = datetime(2021, 1, 10, 'Format', 'dd-MMMM-yyyy');  % First day of the year
t0_2022 = datetime(2022, 1, 9, 'Format', 'dd-MMMM-yyyy');  % First day of the year

% Add true date
for i = 1:height(variants_raw);
    raw_val = variants_raw{i,1}; raw_val = raw_val{1};
    year = str2num(raw_val(1:4));
    week = str2num(raw_val(6:7));
    if year == 2020
        variants_raw.date(i) = t0_2020+(week-1)*7;
    elseif year == 2021
        variants_raw.date(i) = t0_2021+(week-1)*7;
    elseif year == 2022
        variants_raw.date(i) = t0_2022+(week-1)*7;
    end
end

variants_raw = variants_raw(variants_raw.date < datetime(2022,10,25),:);
variants_raw = variants_raw(strcmp(variants_raw.source,'TESSy')==1,:);

variants_raw = variants_raw(:,[5 3 4]);

variants_raw.percent_variant = str2double(variants_raw.percent_variant);

% remove unwanted variants
% retaining WT, alpha, beta, gamma, delta, omicrons BA.1-5
% WT: Other
% Alpha: B.1.1.7
% Beta: B.1.351
% Gamma: P.1
% Delta: B.1.617.2
% Omicron BA.1, BA.2 BA.4, BA.5
% Other
% 'UNK' removed (part of "Other")
variants_raw(strcmp(variants_raw.variant,'UNK'),:) = [];
% 'B.1.1.529' removed (too low)
variants_raw(strcmp(variants_raw.variant,'B.1.1.529'),:) = [];
% 'B.1.427/B.1.429' removed (too low)
variants_raw(strcmp(variants_raw.variant,'B.1.427/B.1.429'),:) = [];
% 'B.1.525' removed (too low)
variants_raw(strcmp(variants_raw.variant,'B.1.525'),:) = [];
% 'B.1.616' removed (too low)
variants_raw(strcmp(variants_raw.variant,'B.1.616'),:) = [];
% 'B.1.620' removed (too low)
variants_raw(strcmp(variants_raw.variant,'B.1.620'),:) = [];
% 'B.1.621' removed (too low)
variants_raw(strcmp(variants_raw.variant,'B.1.621'),:) = [];
% 'BA.2.75' removed (too low)
variants_raw(strcmp(variants_raw.variant,'BA.2.75'),:) = [];
% 'BA.2.86' removed (too low)
variants_raw(strcmp(variants_raw.variant,'BA.2.86'),:) = [];
% 'BQ.1'  integrated into BA.5 (later);
% 'C.37' removed (too low)
variants_raw(strcmp(variants_raw.variant,'C.37'),:) = [];
% 'P.3' removed (too low)
variants_raw(strcmp(variants_raw.variant,'P.3'),:) = [];
% 'XBB' (et similia) removed 
variants_raw(strcmp(variants_raw.variant,'XBB'),:) = [];
variants_raw(strcmp(variants_raw.variant,'XBB.1.5-like'),:) = [];
variants_raw(strcmp(variants_raw.variant,'XBB.1.5-like+F456L'),:) = [];

% Create prevalence file
number_of_week = unique(variants_raw.date);

variants = array2table(zeros(length(number_of_week),1+length(unique(variants_raw.variant))));
variants.Properties.VariableNames = {'Date','B.1.1.7', ...
           'B.1.351','B.1.617.2', 'BA.1','BA.2','BA.4', 'BA.5', 'BQ.1', 'Other', 'P.1'};
variants.Date = number_of_week;
variants_list = unique(variants_raw.variant);
% Assign
for i = 1:length(number_of_week)
    temp = variants_raw(variants_raw.date==number_of_week(i),:);
    for v = 1:length(variants_list);
        try
        variants(i,variants_list(v)) = temp(strcmp(temp.variant,variants_list(v))==1,3);
        catch; 
        end
    end
end

% Inglobate BQ.1 in BA.5 due to similar properties
variants.("BA.5") = variants.("BA.5") +variants.("BQ.1");
variants = variants(:,[1 10 2 3 11 4:8]);

% Check that sum equals 1 at all times
variants(:,2:end) = array2table(table2array(variants(:,2:end)) * 1 ./ sum(table2array(variants(:,2:end)),2));

% Interpolate to get daily values
% Change date of first row to 24 feb 2020, first day of data. This is
% neutral anyway because, at that time, Wild Type is 100% (so data aren't
% impacted)
variants(1,1) = array2table(datetime(2020,02,24));
variants = table2timetable(variants);
variants = retime(variants,'daily','linear');

% Correct and simplify
% 1. Alpha variant is seemingly detected in early 2020. This is probably
% caused by testing errors, and it is brought to 0 in our dataset.
variants(1:40,"Other") = array2table(ones(40,1)); variants(1:40, 'B.1.1.7') = array2table(zeros(40,1));
% 2. Beta variant was never very high in Italy. It is hence removed, and
variants = variants(:,[1 2 4:end]);
% 3. On August 22, 2021, the "Other category" reaches 0, meaning that the
% wild type has gone extinct and other "Other" cases are attributable to
% variants. After this date, this category is removed, and cases are
% distributed among other variants
variants(546:959,"Other") = array2table(zeros(959-546+1,1));
% 4. In August 2022, the Delta variant resurges a bit, but not in a
% statistically significative manner. For the sake of simplification, this
% resurgence is cancelled
variants(911:959, "B.1.617.2") = array2table(zeros(959-911+1,1));
% 5. In summer 2021, BA.2 was detected. This is brought to 0
variants(1:666, "BA.2") = array2table(zeros(666,1));
% 6. Alpha spur in December 2021. This is also cancelled
variants(647:959, "B.1.1.7") = array2table(zeros(959-647+1,1));
% 7. Smaller Gamma resurgence in Oct 2021
variants(596:643, "P.1") = array2table(zeros(643-596+1,1));
% Cleaning complete. Check that sum equals 1 at all times
variants(:,1:end) = array2table(table2array(variants) ./ sum(table2array(variants),2));

% 14 day smoothing
variants(:,1:end) = array2table(smoothdata(table2array(variants), 1, 'movmean', [13 0]));

%figure;
%bar(variants.Date,table2array(variants)','stacked')
save("variants_italy.mat",'variants')
