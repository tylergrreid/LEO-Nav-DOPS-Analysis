function alm_param=read_yuma(filename)
%*************************************************************************
%*     Copyright c 2013 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Todd Walter at:      *
%*     twalter@stanford.edu                                              *
%*************************************************************************
%
%FUNCTION ALM_PARAM=READ_YUMA(FILENAME) reads in a yuma almanac file and stores
%  the results in a matlab matrix, alm_param.
%   FILENAME contains the name of the file with yuma almanac parameters
%   ALM_PARAM is a matrix whose rows correspond to satellites and with columns:
%    1    2    3    4     5    6       7       8          9     10  11 
%   PRN ECCEN TOA INCLIN RORA SQRT_A R_ACEN ARG_PERIG MEAN_ANOM AF0 AF1
%    -    -   sec  rad    r/s m^(1/2) rad     rad        rad     s  s/s
%
% see also ALM2XYZ

%based on code by Jon Nichols
%last modification March 12, 2013 Todd Walter

alm_param=zeros(64,11);

if nargin < 1
  error('you must specify a filename or week number')
end
%if ~ischar(filename) 
%  filename=['YUMA' num2str(filename) '.TXT'];
%end
i=1;
while i<=size(filename,2)
  if iscell(filename)
    fid=fopen(filename{i});
    if fid==-1
      error(['no such file ' filename{i}]);
    end       
  else
    fid=fopen(filename);
    i = size(filename,2);
    if fid==-1
      error(['no such file ' filename]);
    end    
  end

  %read in file
  while fgets(fid)~=-1

    %get prn number
    str = fgets(fid);
    prn = str2double(str(28:end));
    alm_param(prn,1) = prn;
    fgets(fid);

    %get eccentricity
    str = fgets(fid);
    alm_param(prn,2) = str2double(str(28:end));

    %get time of applicability
    str = fgets(fid);
    alm_param(prn,3) = str2double(str(28:end));

    %get inclination angle
    str = fgets(fid);
    alm_param(prn,4) = str2double(str(28:end));

    %get rate of right ascention
    str = fgets(fid);
    alm_param(prn,5) = str2double(str(28:end));

    %get square root of semi-major axis
    str = fgets(fid); 
    alm_param(prn,6) = str2double(str(28:end));

    %get right ascention
    str = fgets(fid); 
    alm_param(prn,7) = str2double(str(28:end));

    %get argument of perigee
    str = fgets(fid); 
    alm_param(prn,8) = str2double(str(28:end));

    %get mean anomaly
    str = fgets(fid); 
    alm_param(prn,9) = str2double(str(28:end));

    %get Af0
    str = fgets(fid); 
    alm_param(prn,10) = str2double(str(28:end));

    %get Af1
    str = fgets(fid); 
    alm_param(prn,11) = str2double(str(28:end));

    %get week number
    str = fgets(fid); 
    %alm_param(prn,3) = alm_param(prn,3) + 604800.0*str2double(str(28:end));
    fgets(fid);

  end;

  fclose(fid);
  i = i+1;
end

%save only nonzero rows
i = find(alm_param(:,1));
alm_param = alm_param(i,:);
