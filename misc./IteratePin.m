%Iterate on a single pin
function IteratePin()

    
    %Local variables
    file_name = 'remap2pin';
	file_nums = [2 31 39 42 47 50]; 
	
	%Fit HMM to each file
	for i = 1:length(file_nums),
        
        disp(sprintf('Loading Pin %d...',file_nums(i)));
        
		%Generate Filename
		if(file_nums(i) < 10),
            file = sprintf('%s0%d_spks.mat',file_name,file_nums(i));
        else
            file = sprintf('%s%d_spks.mat',file_name,file_nums(i));
        end;
        
        %Load file data
        load(file);
        
        %Train HMM
        for j = 1:50,
         
            disp(sprintf('Run %d...',j));   
            [P,Q,pi_i,LL,observ,gamma] = HMM_train_Mult_2(spk_trains(2:size(spk_trains,1),250:2000,:),time_vec(250:2000),[1 64],5,1);
            file = sprintf('Pin%d_%d.mat',file_nums(i),-LL);
            save(file,'P','Q','pi_i','LL','gamma');
            
        end;
        
             
	end;

return;