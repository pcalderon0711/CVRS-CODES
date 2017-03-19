function writeTextFile(paramsetID,P0,params,param_names,resnorm,exitflag,output,N,T,ts,W,code)
%WRITETEXTFILE Summary of this function goes here
%   Detailed explanation goes here

cd('txts');
mkdir(code);
cd(code);

% create name for the file with the counter in text format (e.g.,
% results#.txt)
mynewfilename = strcat(code,'.txt');

% create file for writing  
myfile = fopen(mynewfilename,'w');

% print the test subject name (patient name)
fprintf(myfile,'%14s %12s\n','code : ',code);
fprintf(myfile,'\n');

% print the set of parameters to be identified with the corresponding
% initial values
fprintf(myfile,'%31s\n','parameter set to be identified,:');
for i=1:length(paramsetID)
    fprintf(myfile,'%6s %12.4f\n',paramsetID{i},P0(i));    
end
    
fprintf(myfile,'\n');

% print the complete set of parameters with the nominal values
fprintf(myfile,'%23s\n','full set of parameters:');
for i=1:length(params)
    fprintf(myfile,'%12s %12.5f\n',param_names{i},params(i));
end

fprintf(myfile,'\n');

% print the number of data points
fprintf(myfile,'%17s %4.2f\n','discretization N: ',N);
fprintf(myfile,'\n');

% print the simulation time
fprintf(myfile,'%17s [%4.2f %4.2f]\n','simulation time: ',0,T);
fprintf(myfile,'\n');

% print the rest-exer transition time
fprintf(myfile,'%17s %4.2f\n','transition time: ',ts);
fprintf(myfile,'\n');

% print the workload
fprintf(myfile,'%17s %4.2f\n','workload: ',W);
fprintf(myfile,'\n');

% % print the simulation time
% fprintf(myfile,'%17s [%4.2f %4.2f]\n','simulation time: ',simtime);
% fprintf(myfile,'\n');
% 
% % print the time used for parameter identification
% fprintf(myfile,'%21s [%4.2f %4.2f]\n','identification time: ',IDtime);
% fprintf(myfile,'\n');

% residuum (averaged)
fprintf(myfile,'%10s %4.8f\n','resnorm: ',resnorm);
fprintf(myfile,'\n');

% % residuum (averaged)
% fprintf(myfile,'%10s %4.8f\n','residuum: ',residual);
% fprintf(myfile,'\n');

% exitflag with which fminsearch exited
fprintf(myfile,'%10s %d\n','exitflag: ',exitflag);
fprintf(myfile,'\n');

% 
fprintf(myfile,'%17s %d\n','# of iterations: ',output.iterations);
fprintf(myfile,'\n');

fclose(myfile); % closing myfile

cd ../..
end

