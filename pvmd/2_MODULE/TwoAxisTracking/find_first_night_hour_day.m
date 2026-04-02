% Script to find the first hour of night after sunset every day in a given
% location % return the indices of the hour out of 8760 but for the days
% choosen!! min one day

% works only for full year data
% Change location of need it 

function [first_night_hours]=find_first_night_hour_day(location,init_day,end_day,init_month,end_month)

start_day =init_day;
%start_day =1;

%end_day = end_day;
%end_day = 2;

start_month=init_month;
%start_month=7;

end_month=end_month;

%end_month=7;


%IN Delft for now  % ADD ALL to path before 

%cd('C:\Users\orest\Desktop\SET\Thesis\Toolbox Code\toolbox-m-code-develop_updated license\data\Weather\Locations')

load(location)


DHI=weather_data(:,8);


first_night_hours =zeros(365,1); % Array to store the first night hour of each day

for day = 1:365 % Assuming non-leap year with 365 days
    day_start_index = (day - 1) * 24 + 1; % Starting index of the day
    day_dhi = DHI(day_start_index : day_start_index + 23); % DHI values for the day
    
    sunset_hour = find(day_dhi > 0, 1, 'last'); % Find the last hour of daylight (sunset hour)
    
    if ~isempty(sunset_hour) && sunset_hour < 24
        night_hour = sunset_hour+1; % Find the first hour after sunset with DHI = 0
        
        if ~isempty(night_hour)
            first_night_hours(day) = day_start_index-1  + night_hour; % Store the hour index after sunset % I put -1 here it works!!
        end
    end
end
% % Display the first night hour for each day
% for i = 1:length(night_hours)
%     fprintf('Day %d: First night hour = %d\n', ceil(night_hours(i)/24), mod(night_hours(i), 24));
% end


% Get the starting and ending dates as serial date numbers
start_date = datenum(2023, start_month, start_day);
end_date = datenum(2023, end_month, end_day);

% Get the day numbers within the range
day_numbers = (start_date:end_date) - datenum(2023, 1, 1) + 1; 

%disp(day_numbers);
 
% return only the first night hours of the days you input
 first_night_hours=first_night_hours(day_numbers);

%end



