function [stepLength,newErr] = calculateStepLength(u,uNext,uPrev,dims,model,source,trueRec,G,oldErr,stepLength);
%% Depending on how you solve the problem you most probably need more input variables                      
    recording = zeros(dims.nt,length(dims.recPos),'single');
     stepLength = 512;
    newErr = inf;
    while (newErr > oldErr)
        newErr = 0;
        stepLength = stepLength/2;
        newmodel = model + stepLength*G;
        for s = 1:dims.ds:length(dims.srcPos)
                                  u(:) = 0.0;
                                  uNext(:) = 0.0;
                                  uPrev(:) = 0.0;
            for t = 1:dims.nt
               %  Solve wave equation using test model update
 
            uNext = solveWaveEqn(u,uNext,uPrev,dims,newmodel,source,dims.srcPos(s),t);
            uPrev=u;
            u=uNext;
                %  Record traces
                recording(t,:) = u(dims.recPos);
            end
            %% Calculate new error and check against old
 chi = recording-trueRec(:,:,s);
 newErr=newErr+sum(abs(chi(:)));
 
        end
        
    end
end

