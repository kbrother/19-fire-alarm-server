% ~/MATLAB/R2022a/bin/matlab -nodisplay -nodesktop -r "main SFAS91";
% ps -A | grep MATLAB | awk '{print $1}' | xargs kill -9 $1
% mosquitto_sub -h 143.248.221.190 -u kaist_fire -P 'KaistFire!123' -t /PRE_S/SFAS91 > sfas91_ours_0613.txt
function main(topic)
addpath(genpath('lib'));
addpath(genpath('sofia'));
addpath(genpath('sofia/mylib'));
addpath(genpath('MQTT'))
warning('off')
disp(topic);
%javaaddpath('MQTT/jar/org.eclipse.paho.client.mqttv3-1.1.0.jar')
%javaaddpath('MQTT/mqttasync.jar')

init_vars();
brokerAddress = 'tcp://143.248.221.190';
port = 1883;
connection(topic, brokerAddress, port);
end

% function for subscribe
function connection(topic, brokerAddress, port)
	global mqClient
	mqClient = mqttclient(brokerAddress, Port=port, Username="kaist_fire", Password="KaistFire!123");
	mysub = subscribe(mqClient, strcat('/CFD/', topic), Callback=@cb);
	disp(mysub);
end

% function that initialize topic list
function init_vars
	global alpha 
	global alpha_sum	
	global cnt

	cnt = 0;
	alpha = 0.6;
	alpha_sum = 0;			
end

function initialize_tensor(msg, num_channel, window_len)
	% Create tensors
	global X_last
	global init_tensor	
	global init_tensor_mean
	global init_tensor_std
	global weighted_sum
	global weighted_sq_sum
	global U
	global sigma
	global Yt
	global W
	global L 
	global B
	global F
	global real_num_id	

	fprintf('number of channel in this topic: %d\n', num_channel);	
	init_tensor = cell(num_channel, 1);
	init_tensor_mean = cell(num_channel, 1);
	init_tensor_std = cell(num_channel, 1);
	window_sum = cell(num_channel, 1);
	window_sq_sum = cell(num_channel, 1);
	U = cell(num_channel, 2);		
	sigma = cell(num_channel, 1);
	Yt = cell(num_channel, 1);
	W = cell(num_channel, 1);
	L = cell(num_channel, 1);
	B = cell(num_channel, 1);
	F = cell(num_channel, 1);
	X_last = cell(num_channel, 1);
	real_num_id = cell(num_channel, 1);	
			
	for i=1:num_channel										
		num_id = max(msg(i).id);
		real_num_id{i} = 0;
		% fprintf("num id: %d\n", num_id);
		for j=1:num_id			
			if msg(i).err(j) == 0
				real_num_id{i} = real_num_id{i} + 1;
			end
		end

		if real_num_id{i} > 0
			init_tensor{i} = tenzeros([real_num_id{i} 7 window_len]);
			init_tensor_mean{i} = tenzeros([real_num_id{i} 7 window_len]);
			init_tensor_std{i} = tenzeros([real_num_id{i} 7 window_len]);	
			weighted_sum{i} = zeros(real_num_id{i}, 7);
			weighted_sq_sum{i} = zeros(real_num_id{i}, 7);					
		end	
	end 		
end

function fill_tensor(msg, num_channel, time_stamp, window_len)
	global X_last
	global init_tensor			
	global init_tensor_mean
	global init_tensor_std		
	global alpha
	global alpha_sum
	global weighted_sum
	global weighted_sq_sum
	global real_num_id

	alpha_sum = alpha*alpha_sum + 1;
	for i=1:num_channel				
		num_id = max(msg(i).id);
		if real_num_id{i} == 0
			continue
		end						

		curr_matrix = X_last{i};				
		weighted_sum{i} = alpha * weighted_sum{i} + curr_matrix;					
		weighted_sq_sum{i} = alpha * weighted_sq_sum{i} + curr_matrix .* curr_matrix;	

		init_tensor{i}(:, :, time_stamp) = X_last{i};					
		curr_mean = weighted_sum{i} / alpha_sum;
		curr_std = sqrt(weighted_sq_sum{i} / alpha_sum - curr_mean .* curr_mean);				
		init_tensor_mean{i}(:, :, time_stamp) = curr_mean;
		init_tensor_std{i}(:, :, time_stamp) = curr_std;									
	end	
end

function parse_traffic(msg, num_channel)	
	global X_last;
	global init_tensor				
	global alpha	
	global weighted_sum
	global weighted_sq_sum
	global real_num_id

	for i=1:num_channel	
		if real_num_id{i} > 0				
			idx_cnt = 1;			
			num_id = max(msg(i).id);										
			curr_matrix = zeros(real_num_id{i}, 7);									
			for j=1:num_id
				if msg(i).err(j) == 0
					curr_matrix(idx_cnt, 1) = msg(i).temp(j);
					curr_matrix(idx_cnt, 2) = msg(i).hum(j);
					curr_matrix(idx_cnt, 3) = msg(i).pm1(j);
					curr_matrix(idx_cnt, 4) = msg(i).pm2(j);
					curr_matrix(idx_cnt, 5) = msg(i).pm10(j);
					curr_matrix(idx_cnt, 6) = msg(i).co2(j);
					curr_matrix(idx_cnt, 7) = msg(i).co(j);	
					idx_cnt = idx_cnt + 1;
				end
			end	
			X_last{i} = curr_matrix;																			
		end				
	end	
end

function normalize_init_tensor(msg, num_channel)
	global init_tensor
	global init_tensor_mean
	global init_tensor_std
	global real_num_id

	for i=1:num_channel				
		if real_num_id{i} > 0		
			curr_mask = find(init_tensor_std{i} < 1e-3);
			if ~isempty(curr_mask)
				init_tensor_std{i}(curr_mask) = 1e-3;			
			end
			init_tensor{i} = init_tensor{i} - init_tensor_mean{i};		
			init_tensor{i} = init_tensor{i} ./ init_tensor_std{i};		
		end		
	end	
end

function normalize_matrix(msg, num_channel)
	global tensor_window
	global Yt
	global weighted_sum
	global weighted_sq_sum
	global real_num_id
	global alpha_sum
	global X_last

	for i=1:num_channel		
		if real_num_id{i} > 0 				
			avg_mat = weighted_sum{i} / alpha_sum; 					
			std_mat = sqrt(weighted_sq_sum{i} / alpha_sum - avg_mat .* avg_mat);					
			std_mat(std_mat < 1e-3) = 1e-3;
			
			Yt{i} = X_last{i} - avg_mat;
			Yt{i} = Yt{i} ./ std_mat;					
		end		
	end
end

function denorm_Xhat = denormalize_matrix(Xhat, channel_id, outlier)
	global tensor_window
	global Yt
	global weighted_sum
	global weighted_sq_sum
	global real_num_id
	global alpha_sum
	global alpha
	global X_last;

	assert(real_num_id{channel_id} > 0);
	avg_mat = weighted_sum{channel_id} / alpha_sum;
	std_mat = sqrt(weighted_sq_sum{channel_id} / alpha_sum - avg_mat .* avg_mat);
	std_mat(std_mat < 1e-3) = 1e-3;

	denorm_O = outlier .* std_mat;
	denorm_O(abs(denorm_O) < 500) = 0;
	denorm_Xhat = Xhat .* std_mat + avg_mat;	

	weighted_sum{channel_id} = alpha * weighted_sum{channel_id} + (X_last{channel_id} - denorm_O);
	weighted_sq_sum{channel_id} = alpha * weighted_sq_sum{channel_id} + (X_last{channel_id}-denorm_O) .* (X_last{channel_id}-denorm_O);	
end

function cb(topic, data)
	global init_tensor	
	global U
	global O
	global L
	global B
	global F 
	global Yt
	global sigma
	global W
	global real_num_id;	
	global sfas53_cnt;
	global alpha_sum
	global alpha
	global cnt
	global mqClient
	global weighted_sq_sum
	global weighted_sum

	lambda1 = 0.001;
	lambda3 = 10;
	mu = 0.01;
	window_len = 100;	
	rank = 20;
	period = 43200;

	% temp hum pm1 pm2 pm10 co2 co
	%fireMQTT = mqtt('tcp://143.248.221.190', 'Username', 'kaist_fire', 'Port', 1883, 'Password', "KaistFire!123");    		
	jsondata = jsondecode(data);
    topic = split(topic, "/");
	topic = topic{3};

	% Update count
	if cnt <= window_len
		cnt = cnt + 1;
	end
		
	% Initialize tensor
	num_channel = size(jsondata.cfd, 1);
    if cnt == 1		
		initialize_tensor(jsondata.cfd, num_channel, window_len);	   
		disp("init finish");
    end

	% Fill init tensor	
	parse_traffic(jsondata.cfd, num_channel);		
	if cnt <= window_len							
		fill_tensor(jsondata.cfd, num_channel, cnt, window_len);		
		disp("fill finish");
	end
	
	% Intialize sofia
	if cnt == window_len
		% Normalize tensor
		normalize_init_tensor(jsondata.cfd, num_channel);				
		for i=1:num_channel		
			% Initialization
			if real_num_id{i} > 0
				omega_temp = logical(double(tenones([real_num_id{i} 7 window_len])));		
				%disp(init_tensor{i});						
				[U_init, X_hat, O, ~] = sofia_init(init_tensor{i}, omega_temp, rank, lambda1, lambda3);							
				U{i, 1} = U_init{1};											
				U{i, 2} = U_init{2};				
				W_init = U_init{3};						
				for n=1:2
					weights = sqrt(sum(U{i, n}.^2, 1));
					U{i, n} = U{i, n} ./ weights;
					W_init = W_init .* weights;
				end;
				W{i} = W_init(end, :);				

				% HW fitting		
				%disp(W_init);
				[~,L{i},B{i},F{i}] = hw_add_add_fit(W_init);								

				% Initialize error scale tensor
				Ysz = size(init_tensor{i});
				sigma{i} = 0.1*ones([Ysz(1), Ysz(2)]);										
			end
		end		
		disp("sofia init finish");
	end	

	% Dynamic updates	
	if cnt > window_len				
		normalize_matrix(jsondata.cfd, num_channel);						
		output_json.pre = cell(num_channel, 1);
		for i=1:num_channel		
			num_id = max(jsondata.cfd(i).id);
			curr_output.ch = i;	

			if real_num_id{i} > 0				
				Omega_temp = ones(real_num_id{i}, 7);
				U_N = hw_add_add_forecast(L{i}, B{i}, 1);						
				Ythat = U{i, 1} * diag(U_N) * U{i, 2}';								

				Rt = Yt{i} - Ythat;								
				cRt = huber(Rt./sigma{i}) .*sigma{i};								
				sigma{i} = sigma_update(sigma{i}, Rt, Omega_temp, 0.01);
										
				cRt = Omega_temp .* cRt;			
				G = cell(3, 1);
				
				G{1} = cRt * U{i, 2} * diag(U_N);
				G{2} = cRt' * U{i, 1} * diag(U_N);													
				G{3} = (khatrirao(U{i, 1},U{i, 2})' * reshape(cRt,[],1))' + lambda1 * (W{i} - U_N);
						
				for n=1:2
					G{n} = G{n} * min(1, mu * sqrt(rank)/norm(G{n}, 'fro')); % What?
					U{i, n} = U{i, n} + mu * G{n};
				end					

				G{3} = G{3} * min(1, mu * sqrt(rank)/norm(G{3}, 'fro')); % What?
				U_N = U_N + mu * G{3};					

				for n=1:2
					weights = sqrt(sum(U{i, n}.^2, 1));
					U{i, n} = U{i, n} ./ weights;
					U_N = U_N .* weights;
				end
							
				[L{i}, B{i}] = hw_add_add_update(U_N, L{i}, B{i}, F{i});				
				%disp("hw add add update");
				W{i} = U_N;

				X_hat = double(full(ktensor({U{i, 1}, U{i, 2}, U_N})));
				outlier = Yt{i} - (Ythat + cRt);			
				denorm_X_hat = denormalize_matrix(X_hat, i, outlier);								
				denorm_X_hat(denorm_X_hat < 0) = 1e-3;
				idx_cnt = 1;

				if num_id == 1
					curr_output.id = {1};																								
					curr_output.temp = {denorm_X_hat(1, 1)}; 
					curr_output.hum = {denorm_X_hat(1, 2)}; 
					curr_output.pm1 = {denorm_X_hat(1, 3)}; 
					curr_output.pm2 = {denorm_X_hat(1, 4)}; 
					curr_output.pm10 = {denorm_X_hat(1, 5)}; 
					curr_output.co2 = {denorm_X_hat(1, 6)}; 
					curr_output.co = {denorm_X_hat(1, 7)};
					curr_output.err = {jsondata.cfd(i).err}; 																		
				else	
					curr_output.temp = zeros(num_id, 1);
					curr_output.hum = zeros(num_id, 1);
					curr_output.pm1 = zeros(num_id, 1);
					curr_output.pm2 = zeros(num_id, 1);
					curr_output.pm10 = zeros(num_id, 1);
					curr_output.co2 = zeros(num_id, 1);
					curr_output.co = zeros(num_id, 1);
					
					curr_output.id = [1:num_id];	
					for j=1:num_id					
						if jsondata.cfd(i).err(j) == 0						
							curr_output.temp(j) = denorm_X_hat(idx_cnt, 1); 
							curr_output.hum(j) = denorm_X_hat(idx_cnt, 2); 
							curr_output.pm1(j) = denorm_X_hat(idx_cnt, 3); 
							curr_output.pm2(j) = denorm_X_hat(idx_cnt, 4); 
							curr_output.pm10(j) = denorm_X_hat(idx_cnt, 5); 
							curr_output.co2(j) = denorm_X_hat(idx_cnt, 6); 
							curr_output.co(j) = denorm_X_hat(idx_cnt, 7); 							
							idx_cnt = idx_cnt + 1;
						else
							curr_output.temp(j) = jsondata.cfd(i).temp(j);
							curr_output.hum(j) = jsondata.cfd(i).hum(j);
							curr_output.pm1(j) = jsondata.cfd(i).pm1(j); 
							curr_output.pm2(j) = jsondata.cfd(i).pm2(j); 
							curr_output.pm10(j) = jsondata.cfd(i).pm10(j); 
							curr_output.co2(j) = jsondata.cfd(i).co2(j); 
							curr_output.co(j) = jsondata.cfd(i).co(j);	
						end											
					end	
					curr_output.err = jsondata.cfd(i).err;	
				end
				curr_output.sw_v = jsondata.cfd(i).sw_v;
				curr_output.tm = jsondata.cfd(i).tm;
				output_json.pre{i} = curr_output;						
			elseif num_id == 1
				curr_output.id = {1};	
				curr_output.temp = {jsondata.cfd(i).temp(1)};
				curr_output.hum = {jsondata.cfd(i).hum(1)};
				curr_output.pm1 = {jsondata.cfd(i).pm1(1)}; 
				curr_output.pm2 = {jsondata.cfd(i).pm2(1)}; 
				curr_output.pm10 = {jsondata.cfd(i).pm10(1)}; 
				curr_output.co2 = {jsondata.cfd(i).co2(1)}; 
				curr_output.co = {jsondata.cfd(i).co(1)};	
				curr_output.err = {jsondata.cfd(i).err};
				curr_output.sw_v = jsondata.cfd(i).sw_v;
				curr_output.tm = jsondata.cfd(i).tm;
				output_json.pre{i} = curr_output;
			else
				output_json.pre{i} = jsondata.cfd(i);
			end	
		end				
		output_json = jsonencode(output_json);					
		write(mqClient, strcat('/PRE_S/', topic), output_json);						
		alpha_sum = alpha*alpha_sum + 1;
	end				
end


function new = sigma_update(old, Rt, Omegat, phi)
	rho = biweight(Rt./old); 
	old_2 = old.^2;
	
	new = phi * rho .* old_2 + (1-phi) * old_2;
	new = sqrt(new);
	new = Omegat.*new + (1-Omegat).*old;
end
	
