clear all;
close all;
clc;

DSET='/glider/bob/';
FILE='/home/gil/sea/ship-pat/bob.h5';
BINS = 0:1:50;
DST = 0.50;

% Loading data
datetime	= h5read(FILE,strcat(DSET,'datetime'));
longitude	= h5read(FILE,strcat(DSET,'longitude'));
latitude	= h5read(FILE,strcat(DSET,'latitude'));
depth		= h5read(FILE,strcat(DSET,'depth'));
temperature	= h5read(FILE,strcat(DSET,'temperature'));

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
[resampled_distance,vars]	= distance_resample(distance, DST, [depth temperature temperatures]);
resampled_depth				= vars(:,1);
resampled_temperature		= vars(:,2);
resampled_temperatures		= vars(:,3:end);




smooth_temps = gaussian_filter(resampled_temperatures,DST,5000,3);


% sobel filtering
disp('sobel')
sobel = sobel_filter(smooth_temps(:,20));



disp('plot')
hold on;
scatter(resampled_distance,resampled_depth,[],resampled_temperature,'.');
plot(resampled_distance,resampled_temperatures(:,20),'-k')
plot(resampled_distance,smooth_temps(:,20),'-b')
plot(resampled_distance,1000*sobel,'-r')








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

function out = sobel_filter(in)
	len = length(in);
	out = zeros(len,1);
	for i=2:len
		out(i) = abs(in(i)-in(i-1));
	end
	out = out/norm(out);
end

