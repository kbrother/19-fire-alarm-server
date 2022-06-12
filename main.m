%  ~/MATLAB/R2021a/bin/matlab -nodisplay -nodesktop -r "run('main.m');"
% ps -A | grep MATLAB | awk '{print $1}' | xargs kill -9 $1
% mosquitto_sub -h 143.248.221.190 -u kaist_fire -P 'KaistFire!123' -t /PRE_S/SFAS81 > sfas81.txt
addpath(genpath('lib'));
addpath(genpath('sofia'));
addpath(genpath('sofia/mylib'));
addpath(genpath('MQTT'))
javaaddpath('MQTT/jar/org.eclipse.paho.client.mqttv3-1.1.0.jar')
javaaddpath('MQTT/mqttasync.jar')

fireMQTT = mqtt('tcp://143.248.221.190', 'Username', 'kaist_fire', 'Port', 1883, 'Password', "KaistFire!123");
declare_tensor();
fireSub = subscribe(fireMQTT, '/CFD/#', 'Callback', @cb);

% function that initialize topic list
function declare_tensor
	global init_tensor
	global init_tensor_mean
	global init_tensor_std
	global weighted_sum
	global weighted_sq_sum
	global X_last
	global U
	global O
	global L
	global B
	global F
	global sigma
	global output_json
	global Yt
	global W
	global real_num_id
	global alpha 
	global alpha_sum
	global sfas53_cnt
	global cnt

	sfas53_cnt = 0;
	alpha = 0.99;
	%alpha_sum = 0;	
	alpha_sum = struct;
	init_tensor = struct;
	init_tensor_mean = struct;
	init_tensor_std = struct;
	tensor_window = struct;
	weighted_sum = struct;
	weighted_sq_sum = struct;
	init_tensor_mean = struct;
	U = struct;
	O = struct;
	L = struct;
	B = struct;
	F = struct;
	ouput_json = struct;
	sigma = struct;
	Yt = struct;
	W = struct;
	real_num_id = struct;	
	cnt = struct;
end

function initialize_tensor(topic, msg, num_channel, window_len)
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
	if ~isfield(init_tensor, topic)
		init_tensor.(topic) = cell(num_channel, 1);
		init_tensor_mean.(topic) = cell(num_channel, 1);
		init_tensor_std.(topic) = cell(num_channel, 1);
		window_sum.(topic) = cell(num_channel, 1);
		window_sq_sum.(topic) = cell(num_channel, 1);
		U.(topic) = cell(num_channel, 2);		
		sigma.(topic) = cell(num_channel, 1);
		Yt.(topic) = cell(num_channel, 1);
		W.(topic) = cell(num_channel, 1);
		L.(topic) = cell(num_channel, 1);
		B.(topic) = cell(num_channel, 1);
		F.(topic) = cell(num_channel, 1);
		X_last.(topic) = cell(num_channel, 1);
		real_num_id.(topic) = cell(num_channel, 1);
	end 					
			
	for i=1:num_channel										
		num_id = max(msg(i).id);
		real_num_id.(topic){i} = 0;
		% fprintf("num id: %d\n", num_id);
		for j=1:num_id			
			if msg(i).err(j) == 0
				real_num_id.(topic){i} = real_num_id.(topic){i} + 1;
			end
		end

		if real_num_id.(topic){i} > 0
			init_tensor.(topic){i} = tenzeros([real_num_id.(topic){i} 7 window_len]);
			init_tensor_mean.(topic){i} = tenzeros([real_num_id.(topic){i} 7 window_len]);
			init_tensor_std.(topic){i} = tenzeros([real_num_id.(topic){i} 7 window_len]);	
			weighted_sum.(topic){i} = zeros(real_num_id.(topic){i}, 7);
			weighted_sq_sum.(topic){i} = zeros(real_num_id.(topic){i}, 7);					
		end	
	end 		
end

function fill_tensor(topic, msg, num_channel, time_stamp)
	global X_last
	global init_tensor			
	global init_tensor_mean
	global init_tensor_std		
	global alpha
	global alpha_sum
	global weighted_sum
	global weighted_sq_sum
	global real_num_id
	
	for i=1:num_channel				
		num_id = max(msg(i).id);
		if real_num_id.(topic){i} == 0
			continue
		end		
		idx_cnt = 1;	
		for j=1:num_id				
			if msg(i).err(j) ~= 0
				continue
			end			
			init_tensor.(topic){i}(idx_cnt, 1, time_stamp) = msg(i).temp(j);
			init_tensor.(topic){i}(idx_cnt, 2, time_stamp) = msg(i).hum(j);
			init_tensor.(topic){i}(idx_cnt, 3, time_stamp) = msg(i).pm1(j);
			init_tensor.(topic){i}(idx_cnt, 4, time_stamp) = msg(i).pm2(j);
			init_tensor.(topic){i}(idx_cnt, 5, time_stamp) = msg(i).pm10(j);
			init_tensor.(topic){i}(idx_cnt, 6, time_stamp) = msg(i).co2(j);
			init_tensor.(topic){i}(idx_cnt, 7, time_stamp) = msg(i).co(j);			
			%fprintf("fill tensor end\n"); 					
			idx_cnt = idx_cnt + 1;
		end			

		curr_matrix = double(init_tensor.(topic){i}(:, :, time_stamp));
		X_last.(topic){i} = curr_matrix;		
		weighted_sum.(topic){i} = alpha * weighted_sum.(topic){i} + curr_matrix; 		
		weighted_sq_sum.(topic){i} = alpha * weighted_sq_sum.(topic){i} + curr_matrix .* curr_matrix;				
		curr_mean = weighted_sum.(topic){i} / alpha_sum.(topic);
		curr_std = sqrt(weighted_sq_sum.(topic){i} / alpha_sum.(topic) - curr_mean .* curr_mean);				
		init_tensor_mean.(topic){i}(:, :, time_stamp) = curr_mean;
		init_tensor_std.(topic){i}(:, :, time_stamp) = curr_std;					
	end	
end

function update_window(topic, msg, num_channel)	
	global X_last;
	global init_tensor			
	global init_tensor_mean
	global init_tensor_std		
	global alpha	
	global weighted_sum
	global weighted_sq_sum
	global real_num_id

	for i=1:num_channel			
		if real_num_id.(topic){i} > 0				
			idx_cnt = 1;
			num_id = max(msg(i).id);
			curr_matrix = zeros(real_num_id.(topic){i}, 7);			
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
			
			X_last.(topic){i} = curr_matrix;			
			weighted_sum.(topic){i} = alpha * weighted_sum.(topic){i} + curr_matrix;					
			weighted_sq_sum.(topic){i} = alpha * weighted_sq_sum.(topic){i} + curr_matrix .* curr_matrix;					
		end		
	end	
end

function normalize_init_tensor(msg, topic, num_channel)
	global init_tensor
	global init_tensor_mean
	global init_tensor_std
	global real_num_id

	for i=1:num_channel				
		if real_num_id.(topic){i} > 0		
			curr_mask = find(init_tensor_std.(topic){i} < 1e-3);
			if ~isempty(curr_mask)
				init_tensor_std.(topic){i}(curr_mask) = 1e-3;			
			end
			init_tensor.(topic){i} = init_tensor.(topic){i} - init_tensor_mean.(topic){i};		
			init_tensor.(topic){i} = init_tensor.(topic){i} ./ init_tensor_std.(topic){i};		
		end		
	end	
end

function normalize_matrix(msg, topic, num_channel)
	global tensor_window
	global Yt
	global weighted_sum
	global weighted_sq_sum
	global real_num_id
	global alpha_sum
	global X_last

	for i=1:num_channel		
		if real_num_id.(topic){i} > 0 				
			avg_mat = weighted_sum.(topic){i} / alpha_sum.(topic); 					
			std_mat = sqrt(weighted_sq_sum.(topic){i} / alpha_sum.(topic) - avg_mat .* avg_mat);					
			std_mat(std_mat < 1e-3) = 1e-3;
			
			Yt.(topic){i} = X_last.(topic){i} - avg_mat;
			Yt.(topic){i} = Yt.(topic){i} ./ std_mat;					
		end		
	end
end

function denorm_Xhat = denormalize_matrix(Xhat, topic, channel_id)
	global tensor_window
	global Yt
	global weighted_sum
	global weighted_sq_sum
	global real_num_id
	global alpha_sum

	assert(real_num_id.(topic){channel_id} > 0);
	avg_mat = weighted_sum.(topic){channel_id} / alpha_sum.(topic);
	std_mat = sqrt(weighted_sq_sum.(topic){channel_id} / alpha_sum.(topic) - avg_mat .* avg_mat);
	std_mat(std_mat < 1e-3) = 1e-3;
	denorm_Xhat = Xhat .* std_mat + avg_mat;	
end

function cb(topic, msg)
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

	lambda1 = 0.001;
	lambda3 = 10;
	mu = 0.01;
	window_len = 100;
	num_gd = 1;
	rank = 20;

	% temp hum pm1 pm2 pm10 co2 co
	fireMQTT = mqtt('tcp://143.248.221.190', 'Username', 'kaist_fire', 'Port', 1883, 'Password', "KaistFire!123");    
	jsondata = jsondecode(msg);
    topic = split(topic, "/");
	topic = topic{3};
	if topic == "SFAS53"
		sfas53_cnt = sfas53_cnt + 1;
		if rem(sfas53_cnt, 2) == 1
			return
		end
	end
	
	if ~isfield(cnt, topic)
		cnt.(topic) = 1;
	else
		cnt.(topic) = cnt.(topic) + 1;
	end

	if cnt.(topic) == 1
		alpha_sum.(topic) = 0;			
	end	
	alpha_sum.(topic) = alpha*alpha_sum.(topic) + 1;
	
	num_channel = size(jsondata.cfd, 1);
    if cnt.(topic) == 1		
		initialize_tensor(topic, jsondata.cfd, num_channel, window_len);	   
		disp("init finish");
    end

	if cnt.(topic) <= window_len
		fill_tensor(topic, jsondata.cfd, num_channel, cnt.(topic));
		disp("fill finish");
	end
	
	if cnt.(topic) > window_len
		update_window(topic, jsondata.cfd, num_channel);				
		normalize_matrix(jsondata.cfd, topic, num_channel);		
		
		% Dynamic updates			
		output_json.pre = cell(num_channel, 1);
		for i=1:num_channel		
			num_id = max(jsondata.cfd(i).id);
			curr_output.ch = i;	

			if real_num_id.(topic){i} > 0				
				Omega_temp = ones(real_num_id.(topic){i}, 7);
				U_N = hw_add_add_forecast(L.(topic){i}, B.(topic){i}, 1);						
				Ythat = U.(topic){i, 1} * diag(U_N) * U.(topic){i, 2}';								

				Rt = Yt.(topic){i} - Ythat;								
				cRt = huber(Rt./sigma.(topic){i}) .*sigma.(topic){i};								
				sigma.(topic){i} = sigma_update(sigma.(topic){i}, Rt, Omega_temp, 0.01);
										
				cRt = Omega_temp .* cRt;			
				G = cell(3, 1);

				for j=1:num_gd
					G{1} = cRt * U.(topic){i, 2} * diag(U_N);
					G{2} = cRt' * U.(topic){i, 1} * diag(U_N);													
					G{3} = (khatrirao(U.(topic){i, 1},U.(topic){i, 2})' * reshape(cRt,[],1))' + lambda1 * (W.(topic){i} - U_N);
							
					for n=1:2
			            G{n} = G{n} * min(1, mu * sqrt(rank)/norm(G{n}, 'fro')); % What?
            			U.(topic){i, n} = U.(topic){i, n} + mu * G{n};
					end					

					G{3} = G{3} * min(1, mu * sqrt(rank)/norm(G{3}, 'fro')); % What?
            		U_N = U_N + mu * G{3};					
				end

				for n=1:2
					weights = sqrt(sum(U.(topic){i, n}.^2, 1));
					U.(topic){i, n} = U.(topic){i, n} ./ weights;
					U_N = U_N .* weights;
				end
							
				[L.(topic){i}, B.(topic){i}] = hw_add_add_update(U_N, L.(topic){i}, B.(topic){i}, F.(topic){i});				
				%disp("hw add add update");
				W.(topic){i} = U_N;

				X_hat = double(full(ktensor({U.(topic){i, 1}, U.(topic){i, 2}, U_N})));
				denorm_X_hat = denormalize_matrix(X_hat, topic, i);								
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
		publish(fireMQTT, strcat('/PRE_S/', topic), output_json);	
		publish(fireMQTT, strcat('/PRE_T/', topic), msg);
	end	
	
	if cnt.(topic) == window_len
		% Normalize tensor
		normalize_init_tensor(jsondata.cfd, topic, num_channel);				
		for i=1:num_channel		
			% Initialization
			if real_num_id.(topic){i} > 0
				omega_temp = logical(double(tenones([real_num_id.(topic){i} 7 window_len])));
				%disp(init_tensor.(topic){i});				
				[U_init, X_hat.(topic), O.(topic), ~] = sofia_init(init_tensor.(topic){i}, omega_temp, rank, lambda1, lambda3);							
				U.(topic){i, 1} = U_init{1};											
				U.(topic){i, 2} = U_init{2};				
				W_init = U_init{3};					
				for n=1:2
					weights = sqrt(sum(U.(topic){i, n}.^2, 1));
					U.(topic){i, n} = U.(topic){i, n} ./ weights;
					W_init = W_init .* weights;
				end
				W.(topic){i} = W_init(window_len, :);

				% HW fitting				
				[~,L.(topic){i},B.(topic){i},F.(topic){i}] = hw_add_add_fit(W_init);								

				% Initialize error scale tensor
				Ysz = size(init_tensor.(topic){i});
				sigma.(topic){i} = 0.1*ones([Ysz(1), Ysz(2)]);							
			end
		end		
	end	
end


function new = sigma_update(old, Rt, Omegat, phi)
	rho = biweight(Rt./old); 
	old_2 = old.^2;
	
	new = phi * rho .* old_2 + (1-phi) * old_2;
	new = sqrt(new);
	new = Omegat.*new + (1-Omegat).*old;
end
	
