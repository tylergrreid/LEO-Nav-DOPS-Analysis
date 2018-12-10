function UTC_date = GPS2UTC(GPS_week,t_sec)
%% DESCRIPTION
%
%  This function takes a gps time based date and converts it into its
%  equivalent utc date by finding the appropriate number of leap seconds to
%  subtract from the gps date.
%
%% INPUT
%
%  GPS_week = GPS week number
%
%  t_sec    = second of GPS week (0<t<24*7*3600 = 604800)
%
%% OUPUT
%
%  UTC_date = UTC date vector in the form [yyyy mm dd HH MM SS.FFF]
%
%% IMPLEMENTATION

% define the 1st epoch, Jan 6 1980 00:00:00 
epoch1 = datenum([1980 01 06 00 00 00]);

% total number of seconds from 1st epoch
total_secs = t_sec + GPS_week*7*24*3600;

% add this number of seconds to the first epoch date
GPS_date = epoch1 + total_secs/3600/24;

% determine the number of leap seconds which correspond to the given date
leap_dates = [...
    'Jan 6 1980'
    'Jul 1 1981'
    'Jul 1 1982'
    'Jul 1 1983'
    'Jul 1 1985'
    'Jan 1 1988'
    'Jan 1 1990'
    'Jan 1 1991'
    'Jul 1 1992'
    'Jul 1 1993'
    'Jul 1 1994'
    'Jan 1 1996'
    'Jul 1 1997'
    'Jan 1 1999'
    'Jan 1 2006'
    'Jan 1 2009'
    'Jul 1 2012'
    'Jul 1 2015'];

% compute leap dates to serial date form
serial_leap = datenum(leap_dates);

% find number of leap seconds
leapsec = 0;
for i = 2:length(serial_leap)
    if GPS_date >= serial_leap(i)
        leapsec = leapsec + 1;
    end
end

% compute the UTC date
temp = addtodate(GPS_date,-leapsec,'second');

% convert to date vector
UTC_date = datevec(temp);

