function [G,m_synth,tt,VEL] = tomo_old(X,ngrid)
% randomly generates a velocity grid and raypaths
% and outputs the G matrix (data/model relationship kernel), 
% column vector of randomly-generated model parameters,
% and the travel times for each raypath (no error added)
%
% X is the number o raypaths
% ngrid is the number of grids to a side

% random velocity matrix (.75+/-.1 units/sec)
VEL=(rand(ngrid,ngrid)-.5)./5+.75;

% making into slowness
SLOW=1./VEL;

% square for now (not necessary though)
[k o]=size(SLOW);
if k~=o
    make the grid square please
end

% plotting the velocity grid
figure
pcolor(0:ngrid,0:ngrid,[VEL zeros(ngrid,1); zeros(1,ngrid) 0])
axis ij
caxis([0.65 0.8501])
a=colorbar;
ylabel(a,'Velocity (units/sec)','fontsize',14,'fontweight','bold')
title('Velocity Grid and Raypaths','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold')

% randomizing position on a side
Y=k*rand(X,2);

% randomly choosing a side for start/end
for i=1:X
    milk=round(4*rand(1)+.5);
    juice=round(4*rand(1)+.5);
    
    % start/end shouldn't be on same side
    while milk==juice
        juice=round(4*rand(1)+.5);
    end
    
    % start position
    if milk==1
        A(i,:)=[0 , Y(i,1)];
    elseif milk==2
        A(i,:)=[k , Y(i,1)];
    elseif milk==3
        A(i,:)=[Y(i,1), 0];
    else
        A(i,:)=[Y(i,1), k];
    end
    
    % end position
    if juice==1
        B(i,:)=[0 , Y(i,2)];
    elseif juice==2
        B(i,:)=[k , Y(i,2)];
    elseif juice==3
        B(i,:)=[Y(i,2), 0];
    else
        B(i,:)=[Y(i,2), k];
    end
end

% plotting the raypaths
Y=[A B];
hold on
for i=1:X
    plot(Y(i,1:2:3),Y(i,2:2:4),'w','linewidth',4)
end

% Finding dimensions and checking for problems
m=size(A,1);
n=size(B,1);
if m~=n | m==0
    you have too many/few start or end points
end
if A(:,1)>k | A(:,2)>o | A(:,1)<0 | A(:,2)<0
    your start points are outside the range
end
if B(:,1)>k | B(:,2)>o | B(:,1)<0 | B(:,2)<0
    your end points are outside the range
end
for i=1:m
    if A(i,1)==B(i,1) & A(i,1)-round(A(i,1))==0
        illegal horizontal raypath
    elseif A(i,2)==B(i,2) & A(i,2)-round(A(i,2))==0
        illegal vertical raypath
    elseif A(i,1)==B(i,1) & A(i,2)==B(i,2)
        not a ray
    end
end


% Converting the slowness matrix to a column vector
m_synth=zeros(k*o,1);
for i=1:k
    for j=1:o
        m_synth((i-1)*o+j)=SLOW(i,j);
    end
end


% Creating the length matrix 
G=zeros(m,k*o);
check=zeros(1,1);
for i=1:m
    % cleaning up for the next raypath
    l=zeros(1,1);
    w=zeros(1,1);
    L=zeros(k,o);
    
    % making start/end point positions scalars
    a1=A(i,1);
    b1=A(i,2);
    c1=B(i,1);
    d1=B(i,2);
    
    % finding the slopes (x/y and y/x)
    if a1==c1
        s1=k;
        r1=0;
    elseif b1==d1
        s1=0;
        r1=o;
    else
        s1=(d1-b1)/(c1-a1);
        r1=1/s1;
    end
    
    % defining the starting grid and the step directions
    if c1>a1
        u1=ceil(a1+1e-10);
        t1=ceil(c1+1-1e-10);
        y1=1;
    elseif c1==a1
        u1=ceil(a1);
        t1=floor(c1);
        y1=0;
    else
        u1=ceil(a1);
        t1=floor(c1);
        y1=-1;
    end
    if d1>b1
        v1=ceil(b1+1e-10);
        w1=ceil(d1+1-1e-10);
        z1=1;
    elseif d1==b1
        v1=ceil(b1);
        w1=floor(d1);
        z1=0;
    else
        v1=ceil(b1);
        w1=floor(d1);
        z1=-1;
    end
    
    %finding how far from the grid edge the start point is
    if z1==1
        p1=abs(b1-floor(b1));
    elseif z1==-1
        p1=abs(b1-ceil(b1));
    end
    if y1==1
        q1=abs(a1-floor(a1));
    elseif y1==-1
        q1=abs(a1-ceil(a1));
    end
    
    %finding each grid containing the ray and assigning the
    %corresponding length
    j1=1;
    while u1~=t1 & v1~=w1 & (sum(l)-abs(b1-d1)<-1e-10 | s1==0) & (sum(w)-abs(a1-c1)<-1e-10 | r1==0)
        
        %finding the horizontal length in the grid
        l(j1)=abs(s1);
        if l(j1)>1
            l(j1)=1;
        end
        if j1==1
            if l(j1)+p1>1
                l(j1)=l(j1)-(l(j1)+p1-1);
            end
        end
        if j1>1
            if sum(l(1:j1))+p1-ceil(sum(l(1:j1-1))+p1+1e-10)>0
                l(j1)=l(j1)-sum(l(1:j1))-p1+ceil(sum(l(1:j1-1))+p1+1e-10);
            end
        end
        if z1*(d1-(z1*sum(l)+b1))<0
            l(j1)=l(j1)+z1*(d1-(z1*sum(l)+b1));
        end
        
        %finding w(j) from l(j)
        w(j1)=abs(r1)*l(j1);
        if l(j1)==0
            w(j1)=1;
        end
        if j1==1
            if w(j1)+q1>1
                w(j1)=w(j1)-(w(j1)+q1-1);
                l(j1)=abs(s1)*w(j1);
            end
        end
        if j1>1
            if sum(w(1:j1))+q1-ceil(sum(w(1:j1-1))+q1+1e-10)>0
                w(j1)=w(j1)-sum(w(1:j1))-q1+ceil(sum(w(1:j1-1))+q1+1e-10);
                l(j1)=abs(s1)*w(j1);
            end
        end
        if y1*(c1-(y1*sum(w)+a1))<0
            w(j1)=w(j1)+y1*(c1-(y1*sum(w)+a1));
            l(j1)=abs(s1)*w(j1);
        end
        
        %now setting the length of the grid
        L(u1,v1)=sqrt(l(j1)^2+w(j1)^2);
        
        %moving to the next grid
        if abs((z1*sum(l)+b1)-round(z1*sum(l)+b1))<1e-10
            v1=v1+z1;
        end
        if abs((y1*sum(w)+a1)-round(y1*sum(w)+a1))<1e-10
            u1=u1+y1;
        end
        
        j1=j1+1;
    end
    
    %vectorizing the L matrix into a row in the G matrix
    for i1=1:k
        for j1=1:o
            G(i,(i1-1)*o+j1)=L(i1,j1);
        end
    end
    
    %checking that the length matrix was correct
    check(i)=abs(sum(G(i,:))-sqrt((a1-c1)^2+(b1-d1)^2));
    if check(i)>1e-10
        check
        uhoh
    end
end

%calculating the travel times for each ray
tt=G*m_synth;

end
