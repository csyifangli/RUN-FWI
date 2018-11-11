%% Setting up dimensions
dims.dy =     10; % [m]
dims.dx =     10; % [m]
dims.dt = 1.0e-3; % [s]
 
dims.ny = 201; % Cells in y-direction
dims.nx = 301; % Cells in x-direction
dims.nt = 801; % Amount of time steps
 
%% Model dimensions
dims.modely = 100:150;
dims.modelx = 100:200;
dims.my = length(dims.modely);
dims.mx = length(dims.modelx);
 
%% Source locations
sx = min(dims.modelx):max(dims.modelx);
sy = min(dims.modely)*ones(1,length(sx));
dims.srcPos = sy + dims.ny*sx;
 
%% Receiver locations
rx = min(dims.modelx):max(dims.modelx);
ry = min(dims.modely)*ones(1,length(rx));
dims.recPos = ry+dims.ny*rx;
 
%% Creating background model
bg = zeros(dims.ny,dims.nx,'single');
bg(:) = 2.0e3;         % [m/s] - Background
bg(115:end,:) = 2.3e3; % [m/s] - Layer
 
%% Making initial fields
% using extended boundaries to avoid reflection from boundaries
u = zeros(dims.ny , dims.nx ,  'single');  
uNext = zeros(dims.ny , dims.nx , 'single');
uPrev = zeros(dims.ny , dims.nx , 'single');
 
 
 
%% Begin iteration
model = bg;     % Starting model
dims.ds = 10;   % Grid point distance between sources
maxIter = 30;   % Maximum number of iterations per frequency
freqs = [4,6,8,10,12];  % Frequencies to use in inversion
 


errVec = zeros(1,maxIter*length(freqs));
STlength = zeros(1,maxIter*length(freqs));
  
 

 
it = 1; tic;
for f = freqs
    %% Generating ricker source signature wavelet 
    source = rickerWave(f,dims);
    %% Load true recording
    load (['trueRec_',num2str(f),'Hz.mat']);   
    stepLength = 512;
    for i = 1:maxIter
        fprintf('Iteration %i\t' , i);
        %% Calculate gradient ## IMPLEMENT ##
        [grad, err] = calculateGradient(u,uNext,uPrev,dims,model,source,trueRec)
            
        %% Taper gradient
        G = taperGradient(grad);
        figure(1)
        imagesc(G(dims.modely,dims.modelx));
        title('Gradient');
        axis('image');
        colorbar();
        drawnow();
 
        %% Calculate step length ## IMPLEMENT ##
        [stepLength,newErr] = calculateStepLength(u,uNext,uPrev,dims,model,source,trueRec,G,err,stepLength);
        fprintf("Error:%6.4f\n",newErr);
 
        %% Update model
        model = model + stepLength*G;
        

        
        figure(2)
        imagesc(model(dims.modely,dims.modelx),[1800,2700]);
        title('Model');
        axis('image');
        colorbar(); 
        drawnow();
 
        errVec(it) = newErr; 
        STlength(it) = stepLength;
        figure(3);
        subplot(2,1,1);semilogy(errVec,'-','LineWidth',2);
        title('Error')
        xlabel('Number Of Iteration')
        ylabel('Error')
        subplot(2,1,2);semilogy(STlength,'-','LineWidth',2);
        title('Step Length')
        xlabel('Number Of Iteration')
        ylabel('Step Length')
        drawnow();
        it = it + 1;
        toc
        
    end
end
 
