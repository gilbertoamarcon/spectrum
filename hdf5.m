clear all 
clc

DSET='/glider/bob/';
FILE='/home/gil/sea/ship-pat/bob.h5';

longitude	= h5read(FILE,strcat(DSET,'longitude'));
latitude	= h5read(FILE,strcat(DSET,'latitude'));
depth		= h5read(FILE,strcat(DSET,'depth'));

temperature	= h5read(FILE,strcat(DSET,'temperature'));
datetime	= h5read(FILE,strcat(DSET,'datetime'));

[datetime, longitude, latitude, depth, temperature] = filter(datetime, longitude, latitude, depth, temperature);

[distance, avg_dst] = compute_distance_array(longitude, latitude, depth);

disp(avg_dst)

[temperature,distance] = distance_resample(distance, avg_dst, temperature);

% plot(distance,temperature);
scatter(distance,temperature,'.');

% spectrogram(temperature);

function [datetime, longitude, latitude, depth, temperature] = filter(datetime, longitude, latitude, depth, temperature)
	mat = [datetime, longitude latitude depth temperature];
	mat(any(isnan(mat), 2), :) = [];
	mat(any(isinf(mat), 2), :) = [];
	datetime	= mat(:,1);
	longitude	= mat(:,2);
	latitude	= mat(:,3);
	depth		= mat(:,4);
	temperature	= mat(:,5);
end

function [new_temperature,new_distance] = distance_resample(distance, avg_dst, temperature)

	old_len = length(distance);
	new_len = floor(max(distance)/avg_dst);

	new_temperature	= zeros(new_len,1);
	new_distance	= linspace(0,max(distance),new_len);

	ja = 1;
	jb = 1;
	for i=2:new_len
		dst = new_distance(i);
		while distance(jb) < dst;
			ja = jb;
			jb = jb+1;
		end
		alpha = (dst-distance(ja))/(distance(jb)-distance(ja));
		new_temperature(i) = (1-alpha)*temperature(ja)+alpha*temperature(jb);
	end

	
end

function [distance, avg_dst] = compute_distance_array(longitude, latitude, depth)
	len = length(longitude);
	distance = zeros(len,1);
	ind_dsts = zeros(len,1);
	for i=2:len
		hor_dist = lat_lon_distance(latitude(i-1),longitude(i-1),latitude(i),longitude(i));
		pos_0 = [0 depth(i-1)];
		pos_1 = [hor_dist depth(i)];
		ind_dsts(i) = norm(pos_1 - pos_0);
		distance(i) = distance(i-1) + ind_dsts(i);
	end
	avg_dst = mean(ind_dsts(i));
end

function distance = lat_lon_distance(lat1, lon1, lat2, lon2)
	R = 6378.137;
	dLat = lat2 * pi / 180 - lat1 * pi / 180;
	dLon = lon2 * pi / 180 - lon1 * pi / 180;
	a = sin(dLat/2) * sin(dLat/2) + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dLon/2) * sin(dLon/2);
	c = 2 * atan2(sqrt(a), sqrt(1-a));
	d = R * c;
	distance = d * 1000; 
end