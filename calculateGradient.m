function [gradient, err] = calculateGradient(u,uNext,uPrev,dims,model,source,trueRec) 
%% Depending on how you solve the problem you most probably need more input variables
    recording = zeros(dims.nt,length(dims.recPos),'single');
    gradient  = zeros(dims.ny,dims.nx,'single');
    forwardField = zeros(dims.my,dims.mx,dims.nt,'single'); 
    adjointField = zeros(dims.my,dims.mx,dims.nt,'single');
    err = 0;
    for s = 1:dims.ds:length(dims.srcPos)
        %% Run forward simulation on background model
        u(:) = 0.0;
        uNext(:) = 0.0;
        uPrev(:) = 0.0;
        for t = 1:dims.nt
            % Solve wave equation
            uNext = solveWaveEqn(u,uNext,uPrev,dims,model,source,dims.srcPos(s),t);
            uPrev=u;
            u=uNext;
            % Record traces
            recording(t,:) = u(dims.recPos);
            % Save fprward field for use in correlation
            forwardField(:,:,t) = u(dims.modely,dims.modelx);
        end

        %% Calculate difference and error
        chi = recording-trueRec(:,:,s); % Difference between Synthetic model and true model
        chi = flipud(chi);  %reverse in time and use this residuals as new Source
        err = err+sum(abs(chi(:)));
        
        %% Run adjoint simulation
        u(:)  = 0.0;
        uNext(:) = 0.0;
        uPrev(:) = 0.0;
        
        for t = 1:dims.nt
           uNext = solveWaveEqn(u,uNext,uPrev,dims,model,chi,dims.recPos,t); %using wave Eq. for  back propagation
           uPrev=u;
           u=uNext;
            adjointField(:,:,dims.nt-t+1) = u(dims.modely,dims.modelx);
        end
        %% Correlate
        for t = 2:dims.nt-1
            
           % making sure that both the forward and adjoint fields are aligned in time
           
            derivAdj = (adjointField(:,:,t+1) - adjointField(:,:,t))/dims.dt;
            
            derivForw = (forwardField(:,:,t) - forwardField(:,:,t-1))/dims.dt;
            
    
             gradient(dims.modely,dims.modelx) = gradient(dims.modely,dims.modelx) + derivAdj.*derivForw;  
            
            
        end
    end
end
 

