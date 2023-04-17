

riekf = RIEKF; 



a = [0.737820573000000;
-0.498077751000000;
-9.21242504700000]; 


g = [-0.130113259000000;
0.543296754000000;
-0.187950999000000];

dt = 0.01; 

riekf.prediction(a, g, 0.01)

% obj = riekf
% 
% blkdiag(obj.cov_g, obj.cov_a, obj.cov_gb, obj.cov_ab);


riekf.nonholonomic()
riekf.correction(a)