function [ dF_F ] = relative_fluorescence( meanstack )

dispstat('', 'init')
n_n=size(meanstack,1);
n_frames=size(meanstack,2);
t0_s=0.2;
t1_s=0.75;
t2_s=3;
f=27.33; % Hz recording

t0=round(t0_s*f);
t1=round(t1_s*f);
t2=round(t2_s*f);




td_baseline=zeros(size(meanstack));
w_integral=zeros(n_frames,1);
dF_F=td_baseline;
F_bar=td_baseline;
rel_f=td_baseline;

for cell=1:n_n
    for frame=t1:n_frames-t1
        F_bar(cell, frame)=1/t1*trapz(meanstack(cell,frame-round(t1/2):frame+round(t1/2)));
    end
    for frame=t2+2:size(F_bar,2)
        td_baseline(cell,frame)=min(F_bar(cell,frame-t2-1:frame-1));
    end
    for frame=1:n_frames
        rel_f(cell,frame)=(meanstack(cell,frame)-td_baseline(cell,frame))/td_baseline(cell,frame);
    end
end
rel_f(:,1:t2+t1)=zeros(n_n,t2+t1);
rel_f(:,end-t1+1:end)=zeros(n_n,t1);


n_tau=0:1:n_frames;
w=exp(-(abs(n_tau))./t0);
w_integral=cumtrapz(w);
for cell=1:n_n
    for frame=1:n_frames
        S(1:frame)=rel_f(cell, frame:-1:1).*w(1:frame);
        dF_F(cell,frame)=trapz(S(1:frame))/w_integral(frame+1);
    end
    
    %     figure
    % s1=subplot(3,1,1);
    % plot(meanstack(cell,:))
    % title('raw fluorescence')
    % s2=subplot(3,1,2);
    % plot(rel_f(cell,:))
    % title('rel_f')
    % s3=subplot(3,1,3);
    % plot(dF_F(cell,:))
    % title('dF_F')
    % linkaxes([s1,s2,s3], 'x')
    %
    
    
    dispstat(sprintf('Prepared fluorescence trace for cell %i of %i', cell, n_n))
end



dispstat('Fluorescence traces prepared', 'keepthis', 'timestamp')

