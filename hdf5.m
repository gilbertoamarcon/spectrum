clear all;
close all;
clc;

INPUT_FILE='/home/gil/sea/all.h5';
INPUT_DSET='/glider/bob/';

OUTPUT_FILE='/home/gil/sea/relevance.h5';
OUTPUT_DSET='/glider/relevance/';

BINS = 0:1:50;
DST = 0.25;
DSAMPLE = 100;

% Loading data
datetime	= h5read(INPUT_FILE,strcat(INPUT_DSET,'datetime'));
longitude	= h5read(INPUT_FILE,strcat(INPUT_DSET,'longitude'));
latitude	= h5read(INPUT_FILE,strcat(INPUT_DSET,'latitude'));
depth		= h5read(INPUT_FILE,strcat(INPUT_DSET,'depth'));
temperature	= h5read(INPUT_FILE,strcat(INPUT_DSET,'temperature'));

% Cropping
beg			= 35000;
lim			= 50000;
datetime	= datetime(beg:lim);
longitude	= longitude(beg:lim);
latitude	= latitude(beg:lim);
depth		= depth(beg:lim);
temperature	= temperature(beg:lim);


% Ad-hoc temperature filering
for i=1:length(temperature)
	if temperature(i) < 5
		temperature(i) = NaN;
	end
end









% Remove NaNs and Inf
disp('filter')
cleanup_mat = cleanup([datetime longitude latitude depth temperature]);
datetime	= cleanup_mat(:,1);
longitude	= cleanup_mat(:,2);
latitude	= cleanup_mat(:,3);
depth		= cleanup_mat(:,4);
temperature	= cleanup_mat(:,5);

% Depth-separated data
temperatures = depth_separate_temperature(depth,temperature,BINS);

% Compute array of distances
disp('compute_distance_array')
distance = compute_distance_array(longitude, latitude);

% Resample variables to a fixed sampling rate, relative to distance
disp('distance_resample')
[resampled_distance,vars]	= distance_resample(distance, DST, [datetime longitude latitude depth temperature temperatures]);
resampled_datetime			= vars(:,1);
resampled_longitude			= vars(:,2);
resampled_latitude			= vars(:,3);
resampled_depth				= vars(:,4);
resampled_temperature		= vars(:,5);
resampled_temperatures		= vars(:,6:end);



% Gaussian filter
smooth_temps = gaussian_filter(resampled_temperatures,DST,5000,3);

% Difference filter
sobel = difference_filter(smooth_temps);


disp('plot')
hold on;
axis([min(resampled_distance) max(resampled_distance) 5 50])
scatter(resampled_distance,resampled_depth,50,resampled_temperature,'filled');
for i=1:length(BINS)-1
	% plot(resampled_distance,resampled_temperatures(:,i)+BINS(i),'-k')
	% plot(resampled_distance,smooth_temps(:,i)+BINS(i),'-b')
	plot(resampled_distance,1*sobel(:,i)+BINS(i),'-r')
end
legend('Temperature (Celsius)','Relevance Factor');
xlabel('Horizontal Odometry (meters)');
ylabel('Depth (meters)');
c = colorbar;
c.Label.String = 'Temperature (Celsius)';



odat = [];
olon = [];
olat = [];
odep = [];
orel = [];
for i=1:length(BINS)-1
	odat = [olon; resampled_datetime];
	olon = [olon; resampled_longitude];
	olat = [olat; resampled_latitude];
	odep = [odep; ones(size(resampled_datetime))*BINS(i)];
	orel = [orel; sobel(:,i)];
end

odat = odat(1:DSAMPLE:end);
olon = olon(1:DSAMPLE:end);
olat = olat(1:DSAMPLE:end);
odep = odep(1:DSAMPLE:end);
orel = orel(1:DSAMPLE:end);

% scatter3(olon,olat,odep,[],orel,'.');

dset_len = length(odat);

delete(OUTPUT_FILE);
h5create(OUTPUT_FILE,strcat(OUTPUT_DSET,'datetime'),[1 dset_len]);
h5create(OUTPUT_FILE,strcat(OUTPUT_DSET,'longitude'),[1 dset_len]);
h5create(OUTPUT_FILE,strcat(OUTPUT_DSET,'latitude'),[1 dset_len]);
h5create(OUTPUT_FILE,strcat(OUTPUT_DSET,'depth'),[1 dset_len]);
h5create(OUTPUT_FILE,strcat(OUTPUT_DSET,'relevance'),[1 dset_len]);

h5write(OUTPUT_FILE,strcat(OUTPUT_DSET,'datetime'),odat');
h5write(OUTPUT_FILE,strcat(OUTPUT_DSET,'longitude'),olon');
h5write(OUTPUT_FILE,strcat(OUTPUT_DSET,'latitude'),olat');
h5write(OUTPUT_FILE,strcat(OUTPUT_DSET,'depth'),odep');
h5write(OUTPUT_FILE,strcat(OUTPUT_DSET,'relevance'),orel');




% Remove NaNs and Inf
function mat = cleanup(mat)
	mat(any(isnan(mat), 2), :) = [];
	mat(any(isinf(mat), 2), :) = [];
end



function out = depth_separate_temperature(depth,in,bins)
	len = length(in);
	num_bins = length(bins)-1;
	out = zeros(len,num_bins);
	for i=1:len
		for j=1:num_bins
			if depth(i) >= bins(j) && depth(i) < bins(j+1)
				out(i,j) = in(i);
			else
				out(i,j) = NaN;
			end
		end
	end
end


% Compute distance (in meters) betweein lat-lon horizontal positions
function distance = lat_lon_distance(lat1, lon1, lat2, lon2)
	R = 6378.137;
	dLat = lat2 * pi / 180 - lat1 * pi / 180;
	dLon = lon2 * pi / 180 - lon1 * pi / 180;
	a = sin(dLat/2) * sin(dLat/2) + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dLon/2) * sin(dLon/2);
	c = 2 * atan2(sqrt(a), sqrt(1-a));
	d = R * c;
	distance = d * 1000; 
end

% Compute array of distances
function distance = compute_distance_array(longitude, latitude)
	len			= length(longitude);
	distance	= zeros(len,1);
	for i=2:len
		distance(i) = distance(i-1) + lat_lon_distance(latitude(i-1),longitude(i-1),latitude(i),longitude(i));
	end
end

% Resample variable to a fixed sampling rate, relative to distance
function [new_distance,new_vars] = distance_resample(distance, avg_dst, vars)

	num_vars = size(vars,2);

	old_len = length(distance);
	new_len = floor(max(distance)/avg_dst);

	new_distance	= linspace(0,max(distance),new_len);
	new_vars		= zeros(new_len,num_vars);

	% For each variable
	for v=1:num_vars

		ja = 1;
		jb = 1;
		while isnan(vars(ja,v)) || isinf(vars(ja,v))
			ja = ja+1;
			jb = jb+1;
		end
		for i=2:new_len
			dst = new_distance(i);
			while jb < old_len && (distance(jb) < dst || isnan(vars(jb,v)) || isinf(vars(jb,v)))
				if ~isnan(vars(jb,v)) && ~isinf(vars(jb,v))
					ja = jb;
				end
				jb = jb+1;
			end
			if ja == jb || isnan(vars(jb,v)) || isinf(vars(jb,v))
				new_vars(i,v) = vars(ja,v);
			else
				alpha = (dst-distance(ja))/(distance(jb)-distance(ja));
				new_vars(i,v) = (1-alpha)*vars(ja,v)+alpha*vars(jb,v);
			end
		end

	end
	new_vars(1,:) = new_vars(2,:);

end

% Remove NaNs and Inf
function smooth_temps = gaussian_filter(inputs,fs,window_width,stdev)
	wsize = window_width/fs;
	gaussian = gausswin(wsize,stdev);
	gaussian = gaussian/sum(gaussian);
	inputs = padarray(inputs,wsize,'replicate');
	smooth_temps = filter(gaussian,1,inputs);
	smooth_temps = smooth_temps(wsize+1+0.5*wsize:end-wsize+0.5*wsize,:);
end

function out = difference_filter(in)
	len = length(in);
	out = zeros(size(in));
	for i=2:len
		out(i,:) = abs(in(i,:)-in(i-1,:));
	end
	out = out/max(max(out));
end

